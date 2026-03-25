from __future__ import annotations

import pandas as pd
import gzip
import os.path
import pickle
import datetime
from .utilities import run_cmd, ReadsData
from .fastq_visualization import check_metadata
import json
import shutil
import os


def _read_barcodes_from_metadata(metadata_path):
    """Read barcode -> sample-id mapping from metadata TSV."""
    df = pd.read_csv(metadata_path, sep='\t')
    if 'barcode' not in df.columns:
        raise ValueError("Metadata must have a 'barcode' column.")
    id_col = df.columns[0]
    return dict(zip(df['barcode'], df[id_col]))


def _reverse_complement(seq):
    comp = str.maketrans('ACGTacgt', 'TGCAtgca')
    return seq.translate(comp)[::-1]


def demux_and_trim_qiita(dir_path, fastq_path, metadata_path, adapter="GCTACGGGGGG"):
    """Demultiplex Qiita-style inline-barcoded reads using cutadapt,
    then trim the specified adapter.

    Replaces: qiime import (MultiplexedSingleEndBarcodeInSequence) +
              qiime cutadapt demux-single + qiime cutadapt trim-single.

    Qiita preprocessed FASTQs have the barcode at the start of each read.
    We use cutadapt --action=none with combinatorial demux to split by
    barcode, then cutadapt again to trim the adapter/primer.
    """
    fastq_dir = os.path.join(dir_path, "fastq")
    os.makedirs(fastq_dir, exist_ok=True)

    barcode_map = _read_barcodes_from_metadata(metadata_path)

    # Build forward + reverse-complement barcode lookup
    full_map = {}
    for barcode, sample_id in barcode_map.items():
        full_map[barcode] = sample_id
        full_map[_reverse_complement(barcode)] = sample_id

    # Determine barcode length from first barcode
    barcode_len = len(next(iter(barcode_map)))

    # Read the preprocessed fastq
    is_gz = fastq_path.endswith('.gz')
    opener = gzip.open if is_gz else open

    sample_files = {}
    matched = 0
    unmatched = 0

    with opener(fastq_path, 'rt') as fq:
        while True:
            header = fq.readline()
            if not header:
                break
            seq = fq.readline()
            plus = fq.readline()
            qual = fq.readline()

            # The barcode is at the beginning of the sequence
            bc = seq[:barcode_len].strip()
            sample_id = full_map.get(bc)

            if sample_id:
                # Strip the barcode from the read
                trimmed_seq = seq[barcode_len:]
                trimmed_qual = qual[barcode_len:]
                if sample_id not in sample_files:
                    out_path = os.path.join(fastq_dir, f"{sample_id}_1.fastq")
                    sample_files[sample_id] = open(out_path, 'w')
                sample_files[sample_id].write(header + trimmed_seq + plus + trimmed_qual)
                matched += 1
            else:
                unmatched += 1

    for fh in sample_files.values():
        fh.close()

    print(f"Demultiplexed {matched} reads into {len(sample_files)} samples. {unmatched} unmatched.")

    # Trim adapter from each per-sample FASTQ
    if adapter:
        print(f"Trimming adapter {adapter} from demultiplexed reads...")
        for fname in sorted(os.listdir(fastq_dir)):
            if not fname.endswith(".fastq"):
                continue
            in_path = os.path.join(fastq_dir, fname)
            tmp_path = in_path + ".trimmed"
            cmd = (
                f"cutadapt -g {adapter} -e 0 -o {tmp_path} {in_path}"
            )
            run_cmd([cmd])
            os.replace(tmp_path, in_path)

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


def qiita_visualization(fastq_path,metadata_path,data_type, verbose_print):

    verbose_print("\n")
    verbose_print(datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S'))
    dir_path = os.path.commonpath([os.path.dirname(fastq_path), os.path.dirname(metadata_path)])

    verbose_print("\n")
    verbose_print('Checking metadata...', end= " ")
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
        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start 'Demultiplex & trim reads (cutadapt)' (1/2)")
        demux_and_trim_qiita(dir_path, fastq_path, metadata_path)
        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish 'Demultiplex & trim reads (cutadapt)' (1/2)")
        verbose_print("\n")

        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Start 'QC visualization (FastQC + MultiQC)' (2/2)")
        vis_file_path = run_qc_visualization(dir_path)
        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish 'QC visualization (FastQC + MultiQC)' (2/2)")

        reads_data = detect_reads_data(dir_path)
        verbose_print(f"{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')} -- Finish creating visualization\n")

        pickle.dump(reads_data, open(os.path.join(reads_data.dir_path, "reads_data.pkl"), "wb"))

        print(f"QC report is located in {vis_file_path}\n")
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