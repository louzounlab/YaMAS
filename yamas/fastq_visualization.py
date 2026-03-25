import os

import pandas as pd
import csv
import gzip
import os.path
import pickle
import datetime
from tqdm import tqdm
from .utilities import run_cmd, ReadsData
import json
import shutil
import os


def check_metadata(metadata_path):
    metadata = pd.read_csv(metadata_path, sep='\t')
    if 'barcode' in metadata.columns:
        return "yes"
    return "no"


def _read_barcodes_from_metadata(metadata_path):
    """Read barcode -> sample-id mapping from metadata TSV."""
    df = pd.read_csv(metadata_path, sep='\t')
    if 'barcode' not in df.columns:
        raise ValueError("Metadata must have a 'barcode' column.")
    id_col = df.columns[0]  # First column is typically sample-id
    return dict(zip(df['barcode'], df[id_col]))


def _reverse_complement(seq):
    comp = str.maketrans('ACGTacgt', 'TGCAtgca')
    return seq.translate(comp)[::-1]


def demux_emp_single(dir_path, fastq_path, barcode_path, metadata_path):
    """Demultiplex EMP single-end reads by streaming barcodes and sequences in lockstep.

    Both the barcode file and the sequence file are read one record at a time
    (4 lines each), so memory usage is O(1) regardless of file size.
    Matches each barcode (forward and reverse-complement) against the
    metadata mapping to assign reads to per-sample FASTQ files.

    This replaces qiime demux emp-single.
    """
    fastq_dir = os.path.join(dir_path, "fastq")
    os.makedirs(fastq_dir, exist_ok=True)

    barcode_map = _read_barcodes_from_metadata(metadata_path)

    # Build forward + reverse-complement barcode lookup
    full_map = {}
    for barcode, sample_id in barcode_map.items():
        full_map[barcode] = sample_id
        full_map[_reverse_complement(barcode)] = sample_id

    bc_opener = gzip.open if barcode_path.endswith('.gz') else open
    fq_opener = gzip.open if fastq_path.endswith('.gz') else open

    sample_files = {}
    matched = 0
    unmatched = 0

    try:
        with bc_opener(barcode_path, 'rt') as bfh, fq_opener(fastq_path, 'rt') as fq:
            while True:
                # Read one FASTQ record from barcode file
                bc_header = bfh.readline()
                if not bc_header:
                    break
                bc_seq = bfh.readline().strip()
                bfh.readline()  # +
                bfh.readline()  # qual

                # Read one FASTQ record from sequence file
                header = fq.readline()
                if not header:
                    break
                seq = fq.readline()
                plus = fq.readline()
                qual = fq.readline()

                sample_id = full_map.get(bc_seq)

                if sample_id:
                    if sample_id not in sample_files:
                        out_path = os.path.join(fastq_dir, f"{sample_id}_1.fastq")
                        sample_files[sample_id] = open(out_path, 'w')
                    sample_files[sample_id].write(header + seq + plus + qual)
                    matched += 1
                else:
                    unmatched += 1
    finally:
        for fh in sample_files.values():
            fh.close()

    print(f"Demultiplexed {matched} reads into {len(sample_files)} samples. {unmatched} unmatched.")
    return fastq_dir


def run_qc_visualization(dir_path, threads=8):
    """Run FastQC + MultiQC on demultiplexed FASTQs. Replaces qiime demux summarize."""
    fastq_dir = os.path.join(dir_path, "fastq")
    vis_dir = os.path.join(dir_path, "vis")
    fastqc_out = os.path.join(vis_dir, "fastqc")
    os.makedirs(fastqc_out, exist_ok=True)

    fastq_files = [
        os.path.join(fastq_dir, f)
        for f in sorted(os.listdir(fastq_dir))
        if f.endswith(".fastq") or f.endswith(".fq")
    ]
    if not fastq_files:
        print("Warning: No FASTQ files found for QC.")
        return None

    file_list = " ".join(fastq_files)
    run_cmd([f"fastqc -t {threads} -o {fastqc_out} {file_list}"])
    run_cmd([f"multiqc {fastqc_out} -o {vis_dir} -f -n multiqc_report"])

    report_path = os.path.join(vis_dir, "multiqc_report.html")
    return report_path


def detect_reads_data(dir_path):
    """Detect if data is single or paired-end by inspecting the fastq/ directory."""
    fastq_dir = os.path.join(dir_path, "fastq")
    fastq_files = sorted([f for f in os.listdir(fastq_dir) if f.endswith(".fastq") or f.endswith(".fq")])
    has_r2 = any("_2.fastq" in f or "_2.fq" in f for f in fastq_files)
    return ReadsData(dir_path, fwd=True, rev=has_r2)


def fastq_visualization(fastq_path, barcode_path, metadata_path, data_type, verbose_print):
    verbose_print("\n")
    verbose_print(datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S'))
    dir_path = os.path.commonpath([os.path.dirname(fastq_path), os.path.dirname(metadata_path), os.path.dirname(barcode_path)])

    verbose_print("\n")
    verbose_print('Checking metadata...', end=" ")
    if check_metadata(metadata_path) == "no":
        verbose_print("The 'barcode' column does not exist in metadata.tsv. check and try again.")
        return
    verbose_print('Done.')
    verbose_print("\n")
    os.makedirs(os.path.join(dir_path, "qza"), exist_ok=True)
    os.makedirs(os.path.join(dir_path, "vis"), exist_ok=True)

    verbose_print("Find ALL NEW data in the directory you created:", dir_path)

    if data_type == '16S' or data_type == '18S':

        verbose_print("\n")
        verbose_print(
            f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start 'Demultiplex reads with cutadapt' (1/2)")
        demux_emp_single(dir_path, fastq_path, barcode_path, metadata_path)
        verbose_print(
            f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish 'Demultiplex reads with cutadapt' (1/2)")
        verbose_print("\n")

        verbose_print(
            f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start 'QC visualization (FastQC + MultiQC)' (2/2)")
        vis_file_path = run_qc_visualization(dir_path)
        verbose_print(
            f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish 'QC visualization (FastQC + MultiQC)' (2/2)")

        reads_data = detect_reads_data(dir_path)
        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish creating visualization\n")

        pickle.dump(reads_data, open(os.path.join(reads_data.dir_path, "reads_data.pkl"), "wb"))

        print(f"QC report is located in {vis_file_path}\n"
              f"It is highly recommended to start the trimming from 13")
        
        if reads_data.fwd and reads_data.rev:
            print(f"Note: The data has both forward and reverse reads.\n"
                  f"Therefore, you must give the parameters 'trim' and 'trunc' of export() "
                  f"as comma-separated values (e.g. 13,13 150,200).")
        else:
            print(f"Note: The data has only a forward read.\n"
                  f"Therefore, you must give the parameters 'trim' and 'trunc' of export() "
                  f"as single integer values.")

        return reads_data.dir_path

    else:
        print("YaMAS doesnt support downloading shotgun, yet.")
