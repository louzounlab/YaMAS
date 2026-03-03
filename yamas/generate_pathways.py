import os
from pathlib import Path
from typing import Union
from collections import defaultdict
from .utilities import run_cmd

def run_humann_pipeline(dir_path: Union[str, Path], dataset_id: str, threads: int = 8):
    """
    Runs HUMAnN on each sample individually (concatenating Paired-end if needed), 
    then joins the tables so the final output is exactly 4 files: 3 merged tables and 1 log file.
    """
    base_dir = Path(dir_path)
    fastq_dir = base_dir / "fastq"
    qza_dir = base_dir / "qza"
    humann_dir = base_dir / "humann_results"
    humann_dir.mkdir(parents=True, exist_ok=True)

    print(f"[HUMAnN] Starting pipeline for dataset: {dataset_id}")

    all_fastq = sorted(list(fastq_dir.glob("*.fastq")) + list(fastq_dir.glob("*.fq")))
    if not all_fastq:
        print(f"[HUMAnN] No FASTQ files found in {fastq_dir}. Aborting.")
        return

    log_file = humann_dir / f"{dataset_id}_humann.log"
    if log_file.exists():
        log_file.unlink()

    # 1. Group files by sample to handle Paired-end data correctly
    samples = defaultdict(list)
    for fq in all_fastq:
        name = fq.name
        # Strip common paired-end suffixes to get the base sample name
        for suffix in ["_1_1.fastq", "_1_2.fastq", "_1.fastq", "_2.fastq", ".fastq", ".fq"]:
            if name.endswith(suffix):
                sample_name = name.replace(suffix, "")
                break
        samples[sample_name].append(fq)

    print(f"[HUMAnN] Found {len(samples)} sample(s) from {len(all_fastq)} FASTQ file(s). Running HUMAnN per sample...")

    # 2. Run HUMAnN on each sample
    for sample_name, fq_list in samples.items():
        fq_list.sort() # Ensure _1 comes before _2
        
        # Determine the correct MetaPhlAn profile path
        # It might have a "_1" suffix depending on how MetaPhlAn named it during paired execution
        profile_path = qza_dir / f"{sample_name}_profile.txt"
        if not profile_path.exists():
            profile_path = qza_dir / f"{sample_name}_1_profile.txt"

        input_fastq = fq_list[0]
        temp_concat = False
        
        # If Paired-end (2 files), concatenate them into a temp file for HUMAnN
        if len(fq_list) > 1:
            input_fastq = humann_dir / f"{sample_name}_temp_concat.fastq"
            fq_str = " ".join(str(f) for f in fq_list)
            print(f"  -> Concatenating paired reads for {sample_name}...")
            run_cmd([f"cat {fq_str} > {input_fastq}"], strict=False)
            temp_concat = True

        cmd_parts = [
            "humann",
            f"--input {input_fastq}",
            f"--output {humann_dir}",
            f"--threads {threads}",
            "--input-format fastq",
            "--remove-temp-output"
        ]
        
        # Skip taxonomic prescreen if profile exists
        if profile_path.exists():
            cmd_parts.append(f"--taxonomic-profile {profile_path}")
        else:
            print(f"  -> Warning: Profile not found for {sample_name}. HUMAnN will run MetaPhlAn internally.")

        # Execute and APPEND (>>) to the master log
        full_cmd = " ".join(cmd_parts) + f" >> {log_file} 2>&1"
        print(f"  -> Processing {sample_name}...")
        run_cmd([full_cmd], strict=False)
        
        # Clean up temp concat file to save disk space
        if temp_concat and input_fastq.exists():
            input_fastq.unlink()

    # 3. Join the individual tables into 3 merged matrices
    print("\n[HUMAnN] Joining tables across all samples...")
    tables_to_join = ["genefamilies", "pathabundance", "pathcoverage"]
    
    for table in tables_to_join:
        out_merged = humann_dir / f"{dataset_id}_merged_{table}.tsv"
        join_cmd = (
            f"humann_join_tables --input {humann_dir} "
            f"--output {out_merged} --file_name {table}"
        )
        run_cmd([join_cmd], strict=False)

    # 4. Cleanup: Remove individual sample TSVs, keeping only the 3 merged tables + log
    print("[HUMAnN] Cleaning up intermediate sample files...")
    for table in tables_to_join:
        for f in humann_dir.glob(f"*_{table}.tsv"):
            if not f.name.startswith(f"{dataset_id}_merged_"):
                f.unlink()

    # 5. Fix column headers in the final merged files
    print("\n[HUMAnN] Formatting column headers in final tables...")
    for table in tables_to_join:
        merged_file = humann_dir / f"{dataset_id}_merged_{table}.tsv"
        if merged_file.exists():
            with open(merged_file, 'r') as f:
                lines = f.readlines()
            
            if lines:
                headers = lines[0].strip('\n').split('\t')
                clean_headers = []
                for h in headers:
                    if h.startswith("#"):
                        clean_headers.append(h)  # Keep "# Pathway" or "# Gene Family"
                    else:
                        # מחיקת החלקים המיותרים בלבד:
                        clean_h = h.replace("_temp_concat", "").replace("-RPKs", "")
                        clean_headers.append(clean_h)
                
                lines[0] = '\t'.join(clean_headers) + '\n'
                
                with open(merged_file, 'w') as f:
                    f.writelines(lines)

    print(f"\n[HUMAnN] Pipeline complete! Final output files in {humann_dir}:")
    for f in sorted(humann_dir.iterdir()):
        if f.is_file():
            print(f"  - {f.name}")