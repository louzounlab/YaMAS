import os
import subprocess
from dataclasses import dataclass

@dataclass(frozen=True)
class ReadsData:
    dir_path: str
    fwd: bool = True
    rev: bool = False

def run_cmd(command: list):
    """Executes a shell command."""
    cmd_str = " ".join(command)
    # print(f"Running: {cmd_str}") # Optional: Uncomment for debugging
    exit_code = os.system(cmd_str)
    if exit_code != 0:
        print(f"Warning: Command '{cmd_str}' returned non-zero exit status {exit_code}.")

def qiime2_version():
    """
    Attempts to find QIIME2 version safely.
    Returns a default fallback if not found to prevent crashes in Python 3.13.
    """
    try:
        # Use subprocess.getoutput to actually capture the string.
        # os.system (used in legacy code) only returns exit codes, not output.
        o = subprocess.getoutput("conda env list")
        for env_1 in o.split("\n"):
            for env_2 in env_1.split("/"):
                if "qiime2" in env_2 and " " not in env_2:
                    return env_2
    except Exception:
        pass
    return "qiime2-unknown"

def download_classifier_url():
    version = qiime2_version()
    # Handle cases where version parsing fails
    ver_str = version.split('-')[1] if '-' in version else "2023.2"
    return f"https://data.qiime2.org/{ver_str}/common/gg-13-8-99-nb-classifier.qza"

def check_conda_qiime2():
    """
    Refactored to WARN rather than CRASH.
    Allows Shotgun pipeline to run in non-QIIME environments (Python 3.13).
    """
    conda_prefix = os.environ.get("CONDA_PREFIX", "")
    if "qiime2" not in os.path.split(conda_prefix)[-1]:
        print("Warning: 'qiime2' environment not detected. "
              "Shotgun analysis will proceed, but 16S features requiring QIIME2 will fail.");