import os
import shutil
from pathlib import Path
from typing import Union, List, Optional
from .utilities import run_cmd


def run_humann_pipeline(dir_path: Union[str, Path], dataset_id: str, threads: int = 8):
    """
    Runs HUMAnN once on all samples combined:
      1. Concatenate all FASTQ files (all samples, both paired reads) into one file.
      2. Merge all per-sample MetaPhlAn profiles into one combined profile.
      3. Run HUMAnN once on the combined input.
      4. Rename outputs to:
           {dataset_id}_merged_genefamilies.tsv
           {dataset_id}_merged_pathabundance.tsv
           {dataset_id}_merged_pathcoverage.tsv
           {dataset_id}_humann.log
      5. Clean up temp files.
    """
    base_dir = Path(dir_path)
    fastq_dir = base_dir / "fastq"
    qza_dir = base_dir / "qza"
    humann_dir = base_dir / "humann_results"
    humann_dir.mkdir(parents=True, exist_ok=True)

    print(f"[HUMAnN] Starting single-run pipeline for dataset: {dataset_id}")

    # ── 1. Collect all FASTQ files ────────────────────────────────────────────
    all_fastq = sorted(
        list(fastq_dir.glob("*.fastq")) + list(fastq_dir.glob("*.fq"))
    )
    if not all_fastq:
        print(f"[HUMAnN] No FASTQ files found in {fastq_dir}. Aborting.")
        return

    print(f"[HUMAnN] Found {len(all_fastq)} FASTQ file(s): {[f.name for f in all_fastq]}")

    # ── 2. Concatenate all FASTQs into one combined file ──────────────────────
    combined_fastq = humann_dir / f"{dataset_id}_combined_input.fastq"
    fastq_list = " ".join(str(f) for f in all_fastq)
    print(f"[HUMAnN] Concatenating all FASTQs → {combined_fastq.name} ...")
    run_cmd([f"cat {fastq_list} > {combined_fastq}"])
    print(f"[HUMAnN] Concatenation complete.")

    # ── 3. Collect and merge all MetaPhlAn profiles ───────────────────────────
    profile_files = sorted(qza_dir.glob("*_profile.txt"))
    if not profile_files:
        print(f"[HUMAnN] No MetaPhlAn profiles found in {qza_dir}. Aborting.")
        combined_fastq.unlink(missing_ok=True)
        return

    print(f"[HUMAnN] Found {len(profile_files)} profile(s): {[p.name for p in profile_files]}")

    combined_profile = humann_dir / f"{dataset_id}_combined_profile.txt"

    if len(profile_files) == 1:
        shutil.copy(profile_files[0], combined_profile)
        print(f"[HUMAnN] Single profile — copied as combined profile.")
    else:
        profile_list = " ".join(str(p) for p in profile_files)
        merge_cmd = f"merge_metaphlan_tables.py {profile_list} > {combined_profile}"
        print(f"[HUMAnN] Merging profiles → {combined_profile.name} ...")
        run_cmd([merge_cmd])
        print(f"[HUMAnN] Profile merge complete.")

    if not combined_profile.exists() or combined_profile.stat().st_size == 0:
        print(f"[HUMAnN] Combined profile is missing or empty. Aborting.")
        combined_fastq.unlink(missing_ok=True)
        return

    # ── 4. Run HUMAnN once on the combined input ──────────────────────────────
    log_file = humann_dir / f"{dataset_id}_humann.log"

    cmd_parts = [
        "humann",
        f"--input {combined_fastq}",
        f"--output {humann_dir}",
        f"--taxonomic-profile {combined_profile}",
        f"--threads {threads}",
        "--input-format fastq",
        "--remove-temp-output",
    ]

    full_cmd = " ".join(cmd_parts) + f" > {log_file} 2>&1"
    print(f"[HUMAnN] Running: {full_cmd}")
    run_cmd([full_cmd])
    print(f"[HUMAnN] HUMAnN run complete.")

    # ── 5. Rename HUMAnN outputs to dataset-level names ──────────────────────
    # HUMAnN names outputs after the input stem, e.g.:
    #   {dataset_id}_combined_input_genefamilies.tsv
    #   {dataset_id}_combined_input_pathabundance.tsv
    #   {dataset_id}_combined_input_pathcoverage.tsv
    input_stem = combined_fastq.stem  # e.g. "SRR041654_combined_input"

    rename_map = {
        f"{input_stem}_genefamilies.tsv": f"{dataset_id}_merged_genefamilies.tsv",
        f"{input_stem}_pathabundance.tsv": f"{dataset_id}_merged_pathabundance.tsv",
        f"{input_stem}_pathcoverage.tsv":  f"{dataset_id}_merged_pathcoverage.tsv",
    }

    for src_name, dst_name in rename_map.items():
        src = humann_dir / src_name
        dst = humann_dir / dst_name
        if src.exists():
            src.rename(dst)
            print(f"[HUMAnN] Renamed: {src_name} → {dst_name}")
        else:
            print(f"[HUMAnN] Warning: expected output not found: {src_name}")

    # ── 6. Clean up temp combined files ──────────────────────────────────────
    for tmp in [combined_fastq, combined_profile]:
        if tmp.exists():
            tmp.unlink()
            print(f"[HUMAnN] Removed temp file: {tmp.name}")

    print(f"[HUMAnN] Pipeline complete. Final files in {humann_dir}:")
    for f in sorted(humann_dir.iterdir()):
        if f.is_file():
            print(f"  {f.name}")