import os
import shutil
from pathlib import Path
from typing import Union
from .utilities import run_cmd

def _validate_host_db(db_path: Union[str, Path]) -> Path:
    db_dir = Path(db_path)
    # Check root
    bt2_files = list(db_dir.glob("*.bt2")) + list(db_dir.glob("*.bt2l"))
    if bt2_files: return db_dir

    # Check subdirs
    for sub in ["human_genome", "bowtie2", "human_genome/bowtie2"]:
        check_path = db_dir / sub
        if check_path.exists():
            bt2_files = list(check_path.glob("*.bt2")) + list(check_path.glob("*.bt2l"))
            if bt2_files: return check_path

    raise FileNotFoundError(f"[generate_dehost] Host DB not found in {db_dir}")

def run_dehost_pipeline(dir_path: Union[str, Path], threads: int = 8):
    base_dir = Path(dir_path)
    fastq_dir = base_dir / "fastq"
    knead_out_base = base_dir / "knead_out"
    knead_out_base.mkdir(parents=True, exist_ok=True)
    
    host_db_path = os.environ.get("YAMAS_HOST_DB")
    if not host_db_path:
        print("Warning: YAMAS_HOST_DB not set. Skipping.")
        return

    try:
        valid_db_path = _validate_host_db(host_db_path)
    except FileNotFoundError as e:
        print(e)
        raise

    fastq_files = sorted(list(fastq_dir.glob("*.fastq")) + list(fastq_dir.glob("*.fq")))
    samples = {}
    for f in fastq_files:
        if "_1.fastq" in f.name:
            s = f.name.replace("_1.fastq", "")
            samples.setdefault(s, []).append(f)
        elif "_2.fastq" in f.name:
            s = f.name.replace("_2.fastq", "")
            samples.setdefault(s, []).append(f)
        else:
            samples.setdefault(f.stem, []).append(f)

    print(f"[KneadData] Found {len(samples)} samples.")

    for sample_name, files in samples.items():
        sample_out_dir = knead_out_base / sample_name
        sample_out_dir.mkdir(parents=True, exist_ok=True)

        # FIXED COMMAND: Bypass trimming to avoid Java memory crashes completely
        cmd = [
            "kneaddata",
            f"--threads {threads}",
            f"--reference-db {valid_db_path}",
            f"--output {sample_out_dir}",
            "--bypass-trim", # Skip Trimmomatic (Solves memory crash + argument error)
            "--bypass-trf"   # Skip TRF
        ]

        files.sort()
        if len(files) == 2:
            cmd.append(f"--input1 {files[0]}")
            cmd.append(f"--input2 {files[1]}")
        elif len(files) == 1:
            cmd.append(f"--input {files[0]}")
        
        full_cmd = " ".join(cmd) + f" > {sample_out_dir}/kneaddata.log 2>&1"
        
        print(f"[KneadData] Running safe mode: {full_cmd}")
        run_cmd([full_cmd])
        
        clean_fastq_dir = base_dir / "fastq_clean"
        clean_fastq_dir.mkdir(parents=True, exist_ok=True)
        
        # Robust copy: Look for the specific paired output files
        for f in sample_out_dir.glob("*.fastq"):
            # Only copy the final clean outputs (exclude intermediates)
            if "paired_" in f.name and "contam" not in f.name:
                shutil.copy2(str(f), str(clean_fastq_dir / f.name))
            elif len(files) == 1 and "kneaddata.fastq" in f.name:
                shutil.copy2(str(f), str(clean_fastq_dir / f.name))