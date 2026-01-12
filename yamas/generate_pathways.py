import os
import glob
from pathlib import Path
from typing import Union, List, Optional
from .utilities import run_cmd

def run_humann_pipeline(dir_path: Union[str, Path], dataset_id: str, threads: int = 8):
    """
    High-level orchestrator that finds samples and runs HUMAnN on them.
    """
    base_dir = Path(dir_path)
    fastq_dir = base_dir / "fastq"
    qza_dir = base_dir / "qza" # Where MetaPhlAn profiles are stored
    humann_dir = base_dir / "humann_results"
    humann_dir.mkdir(parents=True, exist_ok=True)

    print(f"[HUMAnN] Starting pipeline for dataset: {dataset_id}")

    # 1. Identify Samples from FASTQ files
    # Group by prefix to handle paired reads (e.g., sample_1.fastq, sample_2.fastq)
    fastq_files = sorted(list(fastq_dir.glob("*.fastq")) + list(fastq_dir.glob("*.fq")))
    samples = {}

    for f in fastq_files:
        if "_1.fastq" in f.name or "_1.fq" in f.name:
            sample_name = f.name.replace("_1.fastq", "").replace("_1.fq", "")
            if sample_name not in samples: samples[sample_name] = []
            samples[sample_name].append(f)
        elif "_2.fastq" in f.name or "_2.fq" in f.name:
            sample_name = f.name.replace("_2.fastq", "").replace("_2.fq", "")
            if sample_name not in samples: samples[sample_name] = []
            samples[sample_name].append(f)
        else:
            # Single end or unknown format
            sample_name = f.stem
            samples[sample_name] = [f]

    print(f"[HUMAnN] Found {len(samples)} samples to process.")

    # 2. Process each sample
    for sample_name, files in samples.items():
        print(f"\n[HUMAnN] Processing sample: {sample_name}")
        
        # A. SMART PROFILE FINDER
        # Check for various naming patterns (e.g. Sample_profile.txt vs Sample_1_profile.txt)
        candidates = [
            qza_dir / f"{sample_name}_profile.txt",
            qza_dir / f"{sample_name}_1_profile.txt",
            qza_dir / f"{files[0].stem}_profile.txt"
        ]
        
        profile_path = None
        for cand in candidates:
            if cand.exists():
                profile_path = cand
                break
        
        if not profile_path:
            print(f"[HUMAnN] Warning: No taxonomic profile found for {sample_name}. Checked: {[p.name for p in candidates]}")
            continue

        # B. Prepare Input (Concatenate if paired)
        input_for_humann = files[0]
        temp_cat_file = None

        files.sort() # Ensure _1 comes before _2
        if len(files) == 2:
            # Concatenate paired reads for HUMAnN
            temp_cat_file = humann_dir / f"{sample_name}_merged.fastq"
            print(f"[HUMAnN] Merging paired reads to {temp_cat_file}...")
            # Simple concatenation: cat file1 file2 > merged
            run_cmd([f"cat {files[0]} {files[1]} > {temp_cat_file}"])
            input_for_humann = temp_cat_file
        
        # C. Run HUMAnN
        try:
            _run_single_humann(
                input_file=input_for_humann,
                output_dir=humann_dir,
                meta_profile=profile_path,
                threads=threads
            )
        except Exception as e:
            print(f"[HUMAnN] Error processing {sample_name}: {e}")
        finally:
            # Cleanup merged file
            if temp_cat_file and temp_cat_file.exists():
                os.remove(temp_cat_file)

def _run_single_humann(
    input_file: Path,
    output_dir: Path,
    meta_profile: Path,
    threads: int
):
    """
    Low-level wrapper to execute the HUMAnN shell command.
    """
    # Define output logs
    log_file = output_dir / f"{input_file.stem}_humann.log"
    
    # Construct command
    # --input-format fastq is explicitly safer
    # --taxonomic-profile is critical to bypass MetaPhlAn re-run
    cmd = [
        "humann",
        f"--input {input_file}",
        f"--output {output_dir}",
        f"--taxonomic-profile {meta_profile}",
        f"--threads {threads}",
        "--input-format fastq",
        "--remove-temp-output" # Clean up intermediate bowtie/diamond files
    ]

    full_cmd = " ".join(cmd) + f" > {log_file} 2>&1"
    
    print(f"[HUMAnN] Executing: {full_cmd}")
    run_cmd([full_cmd])
    print(f"[HUMAnN] Completed {input_file.name}")