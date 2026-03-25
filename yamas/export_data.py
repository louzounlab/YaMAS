from __future__ import annotations

import csv
import datetime
import glob
import os
import pickle
import pandas as pd
from collections import defaultdict

from .utilities import ReadsData, run_cmd

nodes_names = []


# ---------------------------------------------------------------------------
# Helper: reference database path
# ---------------------------------------------------------------------------

def _reference_db_path(classifier_path: str):
    """Return the classifier/reference FASTA path.

    For the new standalone pipeline the user supplies a FASTA reference
    database (e.g. the SILVA or GreenGenes FASTA formatted for vsearch
    --sintax).  The argument name is kept as 'classifier' for CLI
    backwards-compatibility.
    """
    if not os.path.exists(classifier_path):
        raise FileNotFoundError(
            f"Reference database not found at: {classifier_path}\n"
            "Please supply a vsearch-compatible FASTA reference "
            "(e.g. silva_138_99_sintax.fasta or gg_13_8_99_sintax.fasta)."
        )
    return classifier_path


# ---------------------------------------------------------------------------
# Step 1 – Quality trimming with cutadapt  (replaces qiime dada2 trimming)
# ---------------------------------------------------------------------------

def trim_reads(reads_data: ReadsData, left, right, threads: int = 12):
    """Trim/truncate FASTQ reads with cutadapt.

    *left*  = bases to trim from the 5' end
    *right* = length to truncate to

    Can be:
    - Single-end: integer (e.g., "13" or 13)
    - Paired-end: comma-separated (e.g., "13,13" or "150,200")
    """
    def parse_param(param_str, paired):
        """Parse trim/trunc parameter, handling both single and comma-separated formats."""
        param_str = str(param_str).strip()
        if ',' in param_str:
            parts = [int(p.strip()) for p in param_str.split(',')]
            if len(parts) != 2:
                raise ValueError(f"Expected 1 or 2 values, got {len(parts)}")
            return parts[0], parts[1]
        else:
            val = int(param_str)
            if paired:
                raise ValueError(
                    f"Paired-end data requires comma-separated values "
                    f"(e.g., 'trim' as '13,13' not '13'). "
                    f"To process as single-end, use --export-as-single flag."
                )
            return val, None

    paired = reads_data.fwd and reads_data.rev
    fastq_dir = os.path.join(reads_data.dir_path, "fastq")
    trimmed_dir = os.path.join(reads_data.dir_path, "qza", "trimmed")
    os.makedirs(trimmed_dir, exist_ok=True)

    fastq_files = sorted(
        f for f in os.listdir(fastq_dir)
        if f.endswith(".fastq") or f.endswith(".fq")
    )

    if paired:
        trim_f, trim_r = parse_param(left, paired=True)
        trunc_f, trunc_r = parse_param(right, paired=True)
        
        r1_files = [f for f in fastq_files if "_1.fastq" in f or "_1.fq" in f]
        for r1 in r1_files:
            r2 = r1.replace("_1.fastq", "_2.fastq").replace("_1.fq", "_2.fq")
            if r2 not in fastq_files:
                continue
            out1 = os.path.join(trimmed_dir, r1)
            out2 = os.path.join(trimmed_dir, r2)
            cmd = (
                f"cutadapt"
                f" -u {trim_f} -U {trim_r}"
                f" -l {trunc_f} -L {trunc_r}"
                f" -j {threads}"
                f" -o {out1} -p {out2}"
                f" {os.path.join(fastq_dir, r1)} {os.path.join(fastq_dir, r2)}"
            )
            run_cmd([cmd])
    else:
        trim_val, _ = parse_param(left, paired=False)
        trunc_val, _ = parse_param(right, paired=False)
        
        # NEW: In single-end mode, use ONLY forward reads to prevent garbage ASVs
        # from reverse-complement sequences inflating the feature table.
        # When PE fails and we Fall back to SE, reverse reads are discarded.
        r1_files = [f for f in fastq_files if "_1.fastq" in f or "_1.fq" in f]
        
        # Fallback: if data was genuinely single-end (no _1/_2 suffixes), process all
        if not r1_files:
            r1_files = fastq_files
        
        for fq in r1_files:
            out = os.path.join(trimmed_dir, fq)
            cmd = (
                f"cutadapt"
                f" -u {trim_val}"
                f" -l {trunc_val}"
                f" -j {threads}"
                f" -o {out}"
                f" {os.path.join(fastq_dir, fq)}"
            )
            run_cmd([cmd])

    return trimmed_dir


# ---------------------------------------------------------------------------
# Step 2 – Dereplication + denoising + chimera removal with vsearch
#           (replaces qiime dada2 denoise + qiime vsearch cluster)
# ---------------------------------------------------------------------------

def dereplicate_and_denoise(reads_data: ReadsData, trimmed_dir: str, threads: int = 12):
    """Merge PE (if needed), dereplicate, denoise, cluster at 99%, remove chimeras."""
    paired = reads_data.fwd and reads_data.rev
    qza_dir = os.path.join(reads_data.dir_path, "qza")
    merged_dir = os.path.join(qza_dir, "merged")
    os.makedirs(merged_dir, exist_ok=True)

    trimmed_files = sorted(
        f for f in os.listdir(trimmed_dir)
        if f.endswith(".fastq") or f.endswith(".fq")
    )

    # --- 2a. Merge paired ends (if paired) ---
    if paired:
        r1_files = [f for f in trimmed_files if "_1.fastq" in f or "_1.fq" in f]
        for r1 in r1_files:
            r2 = r1.replace("_1.fastq", "_2.fastq").replace("_1.fq", "_2.fq")
            sample = r1.replace("_1.fastq", "").replace("_1.fq", "")
            merged_out = os.path.join(merged_dir, f"{sample}_merged.fastq")
            cmd = (
                f"vsearch --fastq_mergepairs {os.path.join(trimmed_dir, r1)}"
                f" --reverse {os.path.join(trimmed_dir, r2)}"
                f" --fastqout {merged_out}"
                f" --threads {threads}"
            )
            run_cmd([cmd])
        input_dir = merged_dir
    else:
        input_dir = trimmed_dir

    # --- 2b. Concatenate all samples into one FASTA for global dereplication ---
    all_fasta = os.path.join(qza_dir, "all_reads.fasta")
    input_files = sorted(glob.glob(os.path.join(input_dir, "*.fastq")) +
                         glob.glob(os.path.join(input_dir, "*.fq")))

    # Convert FASTQ -> FASTA with quality filtering, relabelling by sample
    for fq in input_files:
        basename = os.path.splitext(os.path.basename(fq))[0]
        sample_label = basename.replace("_merged", "")
        tmp_fasta = fq + ".fasta"
        run_cmd([
            f"vsearch --fastq_filter {fq}"
            f" --fastq_maxee 1.0"
            f" --fastaout {tmp_fasta}"
            f" --relabel {sample_label}."
            f" --fasta_width 0"
        ])

    # Cat all per-sample FASTAs
    tmp_fastas = " ".join(f + ".fasta" for f in input_files)
    run_cmd([f"cat {tmp_fastas} > {all_fasta}"])

    # --- 2c. Dereplicate ---
    derep_fasta = os.path.join(qza_dir, "derep.fasta")
    run_cmd([
        f"vsearch --derep_fulllength {all_fasta}"
        f" --output {derep_fasta}"
        f" --sizeout --minuniquesize 2"
        f" --fasta_width 0"
    ])

    # --- 2d. Denoise (unoise3 algorithm – equivalent to DADA2 denoising) ---
    denoised_fasta = os.path.join(qza_dir, "denoised.fasta")
    run_cmd([
        f"vsearch --cluster_unoise {derep_fasta}"
        f" --centroids {denoised_fasta}"
        f" --sizein --sizeout"
        f" --minsize 2"
    ])

    # --- 2e. Chimera removal (de novo) ---
    nonchim_fasta = os.path.join(qza_dir, "nonchim.fasta")
    run_cmd([
        f"vsearch --uchime3_denovo {denoised_fasta}"
        f" --nonchimeras {nonchim_fasta}"
        f" --sizein --sizeout"
        f" --fasta_width 0"
    ])

    # --- 2f. Cluster at 99% identity ---
    rep_seqs = os.path.join(qza_dir, "rep-seqs-dn-99.fasta")
    run_cmd([
        f"vsearch --cluster_size {nonchim_fasta}"
        f" --id 0.99"
        f" --centroids {rep_seqs}"
        f" --sizein --sizeout"
        f" --fasta_width 0"
    ])

    # --- 2g. Build OTU table by mapping all reads back to centroids ---
    otu_map = os.path.join(qza_dir, "otu_map.txt")
    run_cmd([
        f"vsearch --usearch_global {all_fasta}"
        f" --db {rep_seqs}"
        f" --id 0.99"
        f" --otutabout {otu_map}"
        f" --threads {threads}"
    ])

    return rep_seqs, otu_map


# ---------------------------------------------------------------------------
# Step 3 – Taxonomy classification (replaces qiime feature-classifier)
# ---------------------------------------------------------------------------

def assign_taxonomy(reads_data: ReadsData, data_type, classifier_path: str, threads: int = 12):
    """Run vsearch --sintax for taxonomy assignment."""
    ref_db = _reference_db_path(classifier_path)
    qza_dir = os.path.join(reads_data.dir_path, "qza")
    rep_seqs = os.path.join(qza_dir, "rep-seqs-dn-99.fasta")
    tax_out = os.path.join(qza_dir, "taxonomy_raw.tsv")

    run_cmd([
        f"vsearch --sintax {rep_seqs}"
        f" --db {ref_db}"
        f" --tabbedout {tax_out}"
        f" --sintax_cutoff 0.8"
        f" --threads {threads}"
    ])
    return tax_out


# ---------------------------------------------------------------------------
# Step 4 – Filter taxonomy (replaces qiime taxa filter + feature-table filter)
# ---------------------------------------------------------------------------

def _load_otu_table(otu_map_path: str) -> pd.DataFrame:
    """Load the vsearch OTU table (TSV with '#OTU ID' as first col).

    pandas comment='#' would silently skip the header row and promote the
    first data row to column names, so we read without comment handling and
    strip the leading '#' from the first column name manually.
    """
    df = pd.read_csv(otu_map_path, sep='\t')
    df.rename(columns={df.columns[0]: df.columns[0].lstrip('#').strip()}, inplace=True)
    return df

def _load_taxonomy(tax_raw_path: str) -> pd.DataFrame:
    """Parse the vsearch --sintax tabbedout into a DF with columns
    ['Feature ID', 'Taxon', 'Confidence'].

    sintax tabbedout columns:
      0: query ID (may have ;size=N)
      1: full taxonomy with per-rank scores, e.g. k:Bacteria(0.89),p:...
      2: strand
      3: taxonomy above cutoff (no scores), e.g. k:Bacteria,p:Bacteroidetes
    Confidence is the score of the deepest rank that passed the sintax_cutoff,
    obtained by looking up the last rank of col3 in the scored col1 string.
    """
    import re
    rows = []
    with open(tax_raw_path, 'r') as fh:
        for line in fh:
            parts = line.strip().split('\t')
            feat_id = parts[0].split(';')[0]
            full_taxon = parts[1] if len(parts) > 1 else ""
            # col 3 is the taxonomy above sintax_cutoff (clean, no scores)
            taxon = parts[3] if len(parts) >= 4 and parts[3] else full_taxon
            # Build a score lookup from col1: {rank_label: score}
            score_map = {}
            for m in re.finditer(r'([a-z]:[^(,]+)\(([0-9.]+)\)', full_taxon):
                score_map[m.group(1).strip()] = m.group(2)
            # Confidence = score of the deepest rank in the accepted taxonomy (col3)
            accepted_ranks = [r.strip() for r in taxon.split(',') if r.strip()]
            confidence = score_map.get(accepted_ranks[-1], "") if accepted_ranks else ""
            rows.append({"Feature ID": feat_id, "Taxon": taxon, "Confidence": confidence})
    return pd.DataFrame(rows)

def clean_taxonomy(reads_data: ReadsData, otu_map: str, tax_raw: str, min_samples: int = 3, min_frequency: int = 10):
    """Filter OTUs: remove mitochondria/chloroplast, apply prevalence/abundance filters."""
    exports_dir = os.path.join(reads_data.dir_path, "export")
    os.makedirs(exports_dir, exist_ok=True)

    otu_df = _load_otu_table(otu_map)
    tax_df = _load_taxonomy(tax_raw)

    # Identify OTU ID column
    id_col = otu_df.columns[0]  # '#OTU ID' or 'OTU ID'

    # Remove mitochondria / chloroplast by taxonomic annotation
    exclude_ids = set()
    for _, row in tax_df.iterrows():
        taxon_str = str(row["Taxon"]).lower()
        if "mitochondria" in taxon_str or "chloroplast" in taxon_str:
            exclude_ids.add(row["Feature ID"])
    otu_df = otu_df[~otu_df[id_col].isin(exclude_ids)]

    # Filter by min prevalence (present in >= min_samples samples).
    # Cap min_samples to the number of samples in this dataset so that
    # single-sample (or few-sample) runs are not filtered to zero OTUs.
    sample_cols = [c for c in otu_df.columns if c != id_col]
    effective_min_samples = min(min_samples, len(sample_cols))
    prevalence = (otu_df[sample_cols] > 0).sum(axis=1)
    otu_df = otu_df[prevalence >= effective_min_samples]

    # Filter by min total frequency
    total_freq = otu_df[sample_cols].sum(axis=1)
    otu_df = otu_df[total_freq >= min_frequency]

    return otu_df, tax_df


# ---------------------------------------------------------------------------
# Step 5 – Export OTU table
# ---------------------------------------------------------------------------

def export_otu(reads_data: ReadsData, otu_df: pd.DataFrame):
    exports_dir = os.path.join(reads_data.dir_path, "export")
    os.makedirs(exports_dir, exist_ok=True)
    otu_tsv = os.path.join(exports_dir, "otu.tsv")
    otu_df.to_csv(otu_tsv, sep='\t', index=False)
    return otu_tsv


# ---------------------------------------------------------------------------
# Step 6 – Export taxonomy table
# ---------------------------------------------------------------------------

def export_taxonomy(reads_data: ReadsData, tax_df: pd.DataFrame, otu_df: pd.DataFrame):
    exports_dir = os.path.join(reads_data.dir_path, "export")
    os.makedirs(exports_dir, exist_ok=True)
    id_col = otu_df.columns[0]
    kept_ids = set(otu_df[id_col].tolist())
    tax_filtered = tax_df[tax_df["Feature ID"].isin(kept_ids)]
    tax_tsv = os.path.join(exports_dir, "taxonomy.tsv")
    tax_filtered.to_csv(tax_tsv, sep='\t', index=False)
    return tax_tsv


# ---------------------------------------------------------------------------
# Step 7 – Phylogeny: mafft + fasttree  (replaces qiime alignment mafft,
#           qiime phylogeny fasttree, qiime phylogeny midpoint-root)
# ---------------------------------------------------------------------------

def export_phylogeny(reads_data: ReadsData, otu_df: pd.DataFrame, threads: int = 12):
    """Align rep-seqs with mafft, build tree with FastTree, midpoint-root with BioPython."""
    qza_dir = os.path.join(reads_data.dir_path, "qza")
    exports_dir = os.path.join(reads_data.dir_path, "export")
    os.makedirs(exports_dir, exist_ok=True)

    rep_seqs_full = os.path.join(qza_dir, "rep-seqs-dn-99.fasta")
    id_col = otu_df.columns[0]
    kept_ids = set(otu_df[id_col].tolist())

    # Filter rep-seqs to only include OTUs that survived filtering.
    # Strip ;size=... annotations from headers: semicolons are Newick tree
    # terminators, so leaving them in causes BioPython to misparse the tree.
    rep_seqs_filtered = os.path.join(qza_dir, "rep-seqs-filtered.fasta")
    with open(rep_seqs_full, 'r') as fin, open(rep_seqs_filtered, 'w') as fout:
        write = False
        seq_id = None
        for line in fin:
            if line.startswith('>'):
                seq_id = line[1:].strip().split(';')[0]
                write = seq_id in kept_ids
                if write:
                    fout.write(f'>{seq_id}\n')
            elif write:
                fout.write(line)

    # mafft alignment
    aligned = os.path.join(qza_dir, "aligned-rep-seqs.fasta")
    run_cmd([f"mafft --thread {threads} --auto {rep_seqs_filtered} > {aligned}"])

    # FastTree
    unrooted_nwk = os.path.join(exports_dir, "fasttree-tree.nwk")
    run_cmd([f"FastTree -gtr -nt {aligned} > {unrooted_nwk}"])

    # Midpoint-root with BioPython
    from Bio import Phylo as BioPhylo
    tree = BioPhylo.read(unrooted_nwk, "newick")
    tree.root_at_midpoint()
    rooted_nwk = os.path.join(exports_dir, "tree.nwk")
    BioPhylo.write(tree, rooted_nwk, "newick")

    return rooted_nwk


# ---------------------------------------------------------------------------
# Step 8 – Convert TSV exports to CSV
# ---------------------------------------------------------------------------

def convert_to_csv(reads_data: ReadsData):
    exports_dir = os.path.join(reads_data.dir_path, "export")

    otu_tsv = os.path.join(exports_dir, "otu.tsv")
    otu_csv = os.path.join(exports_dir, "otu.csv")
    with open(otu_tsv, 'r', newline='') as tsv, open(otu_csv, 'w', newline='') as csv_out:
        writer = csv.writer(csv_out)
        reader = csv.reader(tsv, delimiter='\t')
        for row in reader:
            writer.writerow(row)

    tax_tsv = os.path.join(exports_dir, "taxonomy.tsv")
    tax_csv = os.path.join(exports_dir, "taxonomy.csv")
    with open(tax_tsv, 'r', newline='') as tsv, open(tax_csv, 'w', newline='') as csv_out:
        writer = csv.writer(csv_out)
        reader = csv.reader(tsv, delimiter='\t')
        for row in reader:
            writer.writerow(row)


# ---------------------------------------------------------------------------
# Step 9 – Pad OTU table for tree tips
# ---------------------------------------------------------------------------

def export_otu_padding_for_tree(reads_data: ReadsData):
    from Bio import Phylo as BioPhylo

    tree_file = os.path.join(reads_data.dir_path, 'export', 'tree.nwk')
    otu_path = os.path.join(reads_data.dir_path, 'export', 'otu.csv')
    otu_padding_path = os.path.join(reads_data.dir_path, 'export', 'otu_padding.csv')
    otu = pd.read_csv(otu_path)
    id_col = otu.columns[0]
    asv_list = otu[id_col].tolist()
    tree = BioPhylo.read(tree_file, "newick")

    tree_tips = []
    _collect_terminal_names(tree.root, tree_tips)

    in_tree_not_asv = [i for i in tree_tips if i not in asv_list]
    print(f"Adding {len(in_tree_not_asv)} ASVs from tree to CSV...")

    with open(otu_path, 'r', newline='') as orig, open(otu_padding_path, 'w', newline='') as pad:
        reader = csv.reader(orig)
        writer = csv.writer(pad)
        for row in reader:
            writer.writerow(row)

    with open(otu_padding_path, 'a', newline='') as pad:
        writer = csv.writer(pad)
        for i in in_tree_not_asv:
            new_row = [i] + [0] * len(otu.columns[1:])
            writer.writerow(new_row)


def _collect_terminal_names(clade, result: list):
    if clade.is_terminal():
        if clade.name:
            result.append(clade.name)
    else:
        for sub_clade in clade.clades:
            _collect_terminal_names(sub_clade, result)


# ---------------------------------------------------------------------------
# Main entry point  (signature unchanged for CLI compatibility)
# ---------------------------------------------------------------------------

def export(output_dir: str, data_type, trim, trunc, classifier_file_path: str, threads: int = 12, force_single_end: bool = False):
    threads = int(threads)
    print(f"\n{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')}")
    print(f"### Exporting {data_type} ###")

    reads_data: ReadsData = pickle.load(open(os.path.join(output_dir, "reads_data.pkl"), "rb"))
    _reference_db_path(classifier_file_path)

    # NEW: Override paired-end detection if user passes --export-as-single
    if force_single_end and reads_data.rev:
        reads_data = ReadsData(reads_data.dir_path, fwd=True, rev=False)
        print("Forcing single-end processing (--export-as-single flag set)")

    os.makedirs(os.path.join(reads_data.dir_path, "export"), exist_ok=True)

    print("Trimming reads with cutadapt...")
    trimmed_dir = trim_reads(reads_data, left=trim, right=trunc, threads=threads)

    print("Dereplicating, denoising & clustering (vsearch)...")
    rep_seqs, otu_map = dereplicate_and_denoise(reads_data, trimmed_dir, threads=threads)

    print("Assigning taxonomy (vsearch --sintax)...")
    tax_raw = assign_taxonomy(reads_data, data_type, classifier_file_path, threads=threads)

    print("Filtering taxonomy & OTU table...")
    otu_df, tax_df = clean_taxonomy(reads_data, otu_map, tax_raw)

    print("Exporting OTU table...")
    export_otu(reads_data, otu_df)

    print("Exporting taxonomy...")
    export_taxonomy(reads_data, tax_df, otu_df)

    print("Building phylogeny (mafft + FastTree)...")
    export_phylogeny(reads_data, otu_df, threads=threads)

    print("Converting to CSV...")
    convert_to_csv(reads_data)

    print("Padding OTU table for tree...")
    export_otu_padding_for_tree(reads_data)

    print("Export finished.")