from __future__ import annotations

import json
import os
import shutil
from pathlib import Path
from typing import Iterable, List, Dict, Optional, Union
from .utilities import run_cmd  # << like HUMAnN

# ------------------------- internal helpers -------------------------
def _ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)

def _is_fastq(p: Path) -> bool:
    n = p.name.lower()
    return n.endswith(".fastq") or n.endswith(".fastq.gz") or n.endswith(".fq") or n.endswith(".fq.gz")

def _detect_pairs(fastq_dir: Path) -> Dict[str, List[Path]]:
    """Group FASTQ files by sample key, pairing _R1/_R2 or _1/_2 when possible."""
    files = sorted([p for p in fastq_dir.iterdir() if p.is_file() and _is_fastq(p)])
    groups: Dict[str, List[Path]] = {}
    for f in files:
        stem = f.name.replace(".gz", "")
        key = stem
        for tok in ("_R1.fastq", "_R2.fastq", "_1.fastq", "_2.fastq",
                    "_R1.fq", "_R2.fq", "_1.fq", "_2.fq",
                    ".fastq", ".fq"):
            key = key.replace(tok, "")
        groups.setdefault(key, []).append(f)
    return groups

def _load_db_from_config() -> Optional[str]:
    p = Path("~/.yamas/config.json").expanduser()
    if p.exists():
        try:
            return json.loads(p.read_text()).get("clean_db")
        except Exception:
            return None
    return None

def _resolve_host_db(host_db: Optional[Union[str, Path]]) -> Path:
    candidate = host_db or os.environ.get("YAMAS_HOST_DB") or _load_db_from_config()
    if not candidate:
        raise ValueError(
            "[generate_dehost] Host DB path is not set. "
            "Set env YAMAS_HOST_DB, or create ~/.yamas/config.json with "
            '{"clean_db": "/path/to/kneaddata_database/output"}, '
            "or pass host_db=... to run_dehost_pipeline()."
        )
    return Path(candidate)

def _validate_host_db(db_dir: Path) -> None:
    """Require a directory containing Bowtie2 index files (*.bt2/ *.bt2l)."""
    if not db_dir.exists() or not db_dir.is_dir():
        raise FileNotFoundError(
            f"[generate_dehost] Host DB directory not found: {db_dir}\n"
            "Expected a directory created by 'kneaddata_database --download human_genome bowtie2 <ROOT>'."
        )
    idx = list(db_dir.glob("*.bt2")) + list(db_dir.glob("*.bt2l")) + list(db_dir.glob("*.bt2*"))
    if not idx:
        raise FileNotFoundError(
            f"[generate_dehost] No Bowtie2 index files (*.bt2/ *.bt2l) under: {db_dir}"
        )


# ------------------------- public API -------------------------
def run_dehost_pipeline(
    base_dir: Union[str, Path],
    threads: int = 8,
    host_db: Optional[Union[str, Path]] = None,
    run_fastqc: bool = True,
    bypass_trf: bool = False,
    trimmomatic_adapters: Optional[Union[str, Path]] = None,
    extra_kneaddata_args: Optional[Iterable[str]] = None,
) -> List[Path]:
    """
    Run KneadData for all samples in <base_dir>/fastq, using a single-shell command via run_cmd.
    - Per-sample outputs: <base_dir>/knead_out/<sample>/
    - Cleaned copies:     <base_dir>/fastq_clean/
    Returns: sorted list of cleaned FASTQ paths.
    """
    base = Path(base_dir)
    fastq_dir = base / "fastq"
    out_root  = base / "knead_out"
    clean_dir = base / "fastq_clean"

    _ensure_dir(out_root)
    _ensure_dir(clean_dir)

    if not fastq_dir.exists():
        raise FileNotFoundError(f"[generate_dehost] FASTQ folder not found: {fastq_dir}")

    db_dir = _resolve_host_db(host_db)
    _validate_host_db(db_dir)

    pairs = _detect_pairs(fastq_dir)
    cleaned: List[Path] = []

    for key, files in pairs.items():
        if not files:
            continue

        out_dir = out_root / key
        _ensure_dir(out_dir)

        # Build the KneadData command 
        cmd_parts: List[str] = [
            "kneaddata",
            f"--threads {threads}",
            f"--reference-db {db_dir}",
            f"--output {out_dir}",
        ]
        if run_fastqc:
            cmd_parts += ["--run-fastqc-start", "--run-fastqc-end"]
        if bypass_trf:
            cmd_parts += ["--bypass-trf"]
        if trimmomatic_adapters:
            cmd_parts += [
                "--trimmomatic trimmomatic",
                f"--trimmomatic-options ILLUMINACLIP:{trimmomatic_adapters}:2:30:10"
            ]
        if extra_kneaddata_args:
            cmd_parts += list(extra_kneaddata_args)

        files_sorted = sorted(files)
        if len(files_sorted) >= 2:
            r1, r2 = files_sorted[0], files_sorted[1]
            cmd_parts += [f"--input1 {r1}", f"--input2 {r2}"]
        else:
            cmd_parts += [f"--unpaired {files_sorted[0]}"]

        log_file = out_dir / "kneaddata.log"
        full_cmd = " ".join(cmd_parts) + f" > {log_file} 2>&1"
        print(f"[KneadData] Running: {full_cmd}")
        run_cmd([full_cmd])

        # Collect cleaned outputs and copy to fastq_clean/
        for pattern in ("*_kneaddata*.fastq", "*_kneaddata*.fastq.gz",
                        "*_kneaddata*.fq",   "*_kneaddata*.fq.gz"):
            for p in out_dir.glob(pattern):
                dest = clean_dir / p.name
                if dest.exists():
                    dest.unlink()
                shutil.copy2(p, dest)
                cleaned.append(dest)

    return sorted(cleaned)


def combine_for_humann(base_dir: Union[str, Path]) -> Path:
    """
    Concatenate all cleaned FASTQs into one file for HUMAnN compatibility.
    Output: <base_dir>/fastq_clean/combined.clean.fastq
    """
    base = Path(base_dir)
    clean_dir = base / "fastq_clean"
    _ensure_dir(clean_dir)

    combined = clean_dir / "combined.clean.fastq"
    if combined.exists():
        combined.unlink()

    with open(combined, "wb") as w:
        for ext in (".fastq", ".fastq.gz", ".fq", ".fq.gz"):
            for p in sorted(clean_dir.glob(f"*{ext}")):
                if p.name == combined.name:
                    continue
                w.write(p.read_bytes())
    return combined