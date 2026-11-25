"""
Estimated read counts for the Shotgun flow.

YaMAS runs MetaPhlAn with its default output (relative abundance, %).
This module adds a second MetaPhlAn pass over the mapping file that the
first run already saved (fastq/<sample>.bowtie2.bz2), this time with
`-t rel_ab_w_read_stats`. The slow alignment step is not repeated. It merges
MetaPhlAn's per-clade read estimates into a counts table that sits next to
the existing % table.

Nothing that exists today changes: the % profiles, the % tables and the
HUMAnN input stay exactly as they are. Only new files are added:

  qza/<sample>_readstats.txt           MetaPhlAn read-stats profile, per sample
  export/<id>_final_counts.txt         merged counts, same layout as <id>_final.txt
  export/<id>_final_counts_table.csv   same layout as <id>_final_table.csv
  export/<id>_reads_summary.csv        reads processed / mapped, per sample

The counts are MetaPhlAn's *estimated* number of reads per clade (column
`estimated_number_of_reads_from_the_clade`), not exact per-taxon read counts.
"""
import os
import re
import csv

import pandas as pd

from .utilities import run_cmd

READSTATS_SUFFIX = "_readstats.txt"   # must NOT end with "_profile.txt" (those get merged into the % table)
COUNT_COLUMN = "estimated_number_of_reads_from_the_clade"


def add_read_counts(dir_path, dataset_id, threads, target_index, db_param):
    """Write the counts tables next to the % table.

    Never raises: if anything goes wrong it prints a warning, and the rest of
    the Shotgun flow (and its existing outputs) continues exactly as before.
    """
    try:
        _add_read_counts(dir_path, dataset_id, threads, target_index, db_param)
    except Exception as e:
        print(f"[Counts] Warning: read counts were not created ({e}). "
              f"The % tables are not affected.")


def _add_read_counts(dir_path, dataset_id, threads, target_index, db_param):
    fastq_dir = os.path.join(dir_path, "fastq")
    qza_dir = os.path.join(dir_path, "qza")
    export_dir = os.path.join(dir_path, "export")

    counts_txt = os.path.join(export_dir, f"{dataset_id}_final_counts.txt")
    counts_csv = os.path.join(export_dir, f"{dataset_id}_final_counts_table.csv")
    summary_csv = os.path.join(export_dir, f"{dataset_id}_reads_summary.csv")
    percent_txt = os.path.join(export_dir, f"{dataset_id}_final.txt")

    # Remove counts tables from an earlier run, so a failed run never leaves stale counts behind.
    for old in (counts_txt, counts_csv, summary_csv):
        if os.path.exists(old):
            os.remove(old)

    # Same samples as the % table: one per qza/*_profile.txt
    profiles = sorted(f for f in os.listdir(qza_dir) if f.endswith("_profile.txt"))
    if not profiles:
        print("[Counts] No MetaPhlAn profiles found - skipping read counts.")
        return

    print(f"[Counts] Estimating read counts for {len(profiles)} sample(s)...")
    samples = []    # (column name, database line, header numbers, {clade: count})
    skipped = []
    for profile in profiles:
        prefix = profile[:-len("_profile.txt")]
        # Column name exactly as MetaPhlAn's merge names it in the % table
        sample = os.path.splitext(profile)[0].replace('_profile', '')
        map_file = os.path.join(fastq_dir, f"{prefix}.bowtie2.bz2")
        stats_file = os.path.join(qza_dir, f"{prefix}{READSTATS_SUFFIX}")

        if not os.path.exists(map_file):
            print(f"[Counts] Warning: no mapping file for {sample} ({map_file}) - skipped.")
            skipped.append(sample)
            continue

        if os.path.exists(stats_file):
            os.remove(stats_file)
        cmd = (
            f"metaphlan {map_file} --input_type bowtie2out -t rel_ab_w_read_stats "
            f"--nproc {threads} -o {stats_file} "
            f"--index {target_index}{db_param}"
        )
        run_cmd([cmd], suppress_warnings=["Failed to launch x86-64-v3"])

        if not os.path.exists(stats_file):
            print(f"[Counts] Warning: MetaPhlAn did not create read counts for {sample} - skipped.")
            skipped.append(sample)
            continue

        try:
            version, numbers, relab, counts = _read_stats_file(stats_file)
        except (ValueError, IndexError) as e:
            print(f"[Counts] Warning: could not read {stats_file} ({e}) - skipped.")
            skipped.append(sample)
            continue

        # Safety check: the counts must come from the same profiling as the % table.
        try:
            same = _same_relab(_read_relab(os.path.join(qza_dir, profile)), relab)
        except (ValueError, IndexError):
            same = False
        if not same:
            print(f"[Counts] Warning: the read counts of {sample} do not match its % profile "
                  f"(the mapping file may come from a different run).")

        samples.append((sample, version, numbers, counts))

    if skipped:
        print(f"[Counts] Note: no read counts for {len(skipped)} sample(s): {', '.join(skipped)}")
    if not samples:
        print("[Counts] No read counts were created.")
        return

    versions = {version for _, version, _, _ in samples}
    if len(versions) > 1:
        print("[Counts] Warning: the samples were profiled with different MetaPhlAn databases - "
              "counts table not created (MetaPhlAn applies the same rule to the % table).")
        return

    no_estimate = sum(v is None for _, _, _, counts in samples for v in counts.values())

    # Merge the same way MetaPhlAn merges the % profiles: a clade missing from a sample = 0
    table = pd.DataFrame({sample: pd.Series(counts, dtype="float64")
                          for sample, _, _, counts in samples})
    table = _align_to_percent_table(table, percent_txt)
    table = table.fillna(0).astype("int64")
    table.index.name = "clade_name"
    if no_estimate:
        print(f"[Counts] Note: {no_estimate} value(s) had no read estimate from MetaPhlAn ('-') "
              f"and were written as 0.")

    with open(counts_txt, "w") as out:
        out.write(versions.pop() + "\n")
        table.to_csv(out, sep="\t")
    _txt_to_csv(counts_txt, counts_csv)

    numbers_by_sample = {sample: numbers for sample, _, numbers, _ in samples}
    with open(summary_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["sample", "reads_processed", "estimated_reads_mapped_to_known_clades"])
        for sample in table.columns:
            numbers = numbers_by_sample[sample]
            writer.writerow([sample, numbers.get("reads_processed", ""), numbers.get("mapped", "")])

    print(f"[Counts] Read counts for {len(samples)} sample(s) written to {counts_txt}")


def _read_metaphlan_file(path):
    """Split a MetaPhlAn output file into (comment lines, column names, data rows)."""
    comments, rows = [], []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                comments.append(line)
            else:
                rows.append(line.split("\t"))
    if not comments:
        raise ValueError("no header lines")
    columns = comments[-1].lstrip("#").split("\t")
    return comments, columns, rows


def _clade(name):
    # With no detected taxa, MetaPhlAn writes "UNCLASSIFIED" in the % profile
    # but "unclassified" in the read-stats profile.
    return "UNCLASSIFIED" if name.upper() == "UNCLASSIFIED" else name


def _read_stats_file(path):
    """Read one `-t rel_ab_w_read_stats` profile.

    Returns (database line, {"reads_processed", "mapped"}, {clade: %}, {clade: count or None}).
    """
    comments, columns, rows = _read_metaphlan_file(path)
    if "relative_abundance" not in columns or COUNT_COLUMN not in columns:
        raise ValueError("not a read-stats profile")
    i_rel = columns.index("relative_abundance")
    i_count = columns.index(COUNT_COLUMN)

    numbers = {}
    for c in comments:
        m = re.match(r"#(\d+) reads processed", c)
        if m:
            numbers["reads_processed"] = int(m.group(1))
        m = re.match(r"#estimated_reads_mapped_to_known_clades:(\d+)", c)
        if m:
            numbers["mapped"] = int(m.group(1))

    relab, counts = {}, {}
    for row in rows:
        clade = _clade(row[0])
        relab[clade] = float(row[i_rel])
        value = row[i_count] if len(row) > i_count else "-"
        counts[clade] = int(value) if value.isdigit() else None   # "-" = MetaPhlAn gave no estimate
    return comments[0], numbers, relab, counts


def _read_relab(profile_path):
    """{clade: %} from an existing default (rel_ab) profile."""
    _, columns, rows = _read_metaphlan_file(profile_path)
    i_rel = columns.index("relative_abundance")
    return {_clade(row[0]): float(row[i_rel]) for row in rows}


def _same_relab(a, b):
    return set(a) == set(b) and all(abs(a[c] - b[c]) < 1e-4 for c in a)


def _align_to_percent_table(table, percent_txt):
    """Put rows and columns in the same order as the % table, so the two tables line up."""
    if not os.path.exists(percent_txt):
        return table
    with open(percent_txt) as f:
        lines = f.read().splitlines()
    if len(lines) < 2:
        return table
    pct_columns = lines[1].split("\t")[1:]
    pct_rows = [line.split("\t")[0] for line in lines[2:] if line]
    known_rows, known_cols = set(pct_rows), set(pct_columns)
    rows = [r for r in pct_rows if r in table.index] + [r for r in table.index if r not in known_rows]
    cols = [c for c in pct_columns if c in table.columns] + [c for c in table.columns if c not in known_cols]
    return table.reindex(index=rows, columns=cols)


def _txt_to_csv(input_file, output_file):
    """Same conversion as metaphlan_txt_csv() in create_visualization.py, for the counts file.
    Kept separate so that the existing function stays untouched."""
    with open(input_file, 'r') as txt:
        lines = txt.readlines()
    if not lines:
        return

    headers = lines[0].strip().split('\t')
    data = [line.strip().split('\t') for line in lines[1:]]
    headers = [h.replace('|', ',') for h in headers]
    data = [[e.replace('|', ',') for e in row] for row in data]
    transposed = list(map(list, zip(*data)))

    with open(output_file, 'w', newline='') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(headers)
        writer.writerows(transposed)
