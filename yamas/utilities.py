import os
import subprocess
from dataclasses import dataclass

@dataclass(frozen=True)
class ReadsData:
    dir_path: str
    fwd: bool = True
    rev: bool = False

def run_cmd(command: list, strict: bool = False, suppress_warnings: list = None):
    """Executes a shell command.
    
    Args:
        suppress_warnings: List of substrings. Any stderr line containing one of
                           these strings will be silently dropped.
    """
    import sys
    cmd_str = " ".join(command)
    # print(f"Running: {cmd_str}") # Optional: Uncomment for debugging
    if suppress_warnings:
        result = subprocess.run(cmd_str, shell=True, stderr=subprocess.PIPE, text=True)
        for line in result.stderr.splitlines():
            if not any(pattern in line for pattern in suppress_warnings):
                print(line, file=sys.stderr)
        exit_code = result.returncode
    else:
        exit_code = os.system(cmd_str)
    if exit_code != 0:
        msg = f"Command '{cmd_str}' returned non-zero exit status {exit_code}."
        if strict:
            raise RuntimeError(msg)
        print(f"Warning: {msg}")

def qiime2_version():
    """
    Legacy function kept for reference. No longer needed since QIIME2
    dependency has been removed. Returns a placeholder string.
    """
    return "qiime2-removed"

def download_classifier_url():
    """Return a URL hint for downloading a vsearch-compatible reference database."""
    return ("https://www.arb-silva.de/download/archive/qiime/ or "
            "ftp://greengenes.microbio.me/greengenes_release/ "
            "(download a FASTA formatted for vsearch --sintax)")

def check_conda_qiime2():
    """No-op. QIIME2 is no longer required — kept for backward compatibility."""
    pass