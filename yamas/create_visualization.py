from __future__ import annotations
import os
import csv
import os.path
import pickle
import datetime
from tqdm import tqdm
from metaphlan.utils.merge_metaphlan_tables import merge
from .utilities import run_cmd, ReadsData, check_conda_qiime2
import json
import shutil
from collections import Counter

# Correct import
from .generate_pathways import run_humann_pipeline
from pathlib import Path

from .generate_dehost import run_dehost_pipeline  
from typing import Union

CONDA_PREFIX = os.environ.get("CONDA_PREFIX", None)

def check_input(acc_list: str):
    print("Input path:", acc_list, "... ", end=" ")
    if not (os.path.exists(acc_list)):
        print("Invalid. File does not exist.")
    elif not (os.path.isfile(acc_list)):
        print("Invalid. Not a file.")
    else:
        print("Valid.")
        
def create_dir(dir_name, specific_location):
    if specific_location:
        dir_path = os.path.join(os.path.abspath(specific_location), dir_name)
    else:
        dir_path = os.path.abspath(dir_name)

    os.makedirs(dir_path, exist_ok=True)
    print(f"{dir_path} created.")

    for sub in ['sra', 'fastq', 'qza', 'vis', 'humann_results', 'export']:
        os.makedirs(os.path.join(dir_path, sub), exist_ok=True)

    return dir_path

def download_data_from_sra(dir_path: str, acc_list: str = ""):
    sra_dir = os.path.join(dir_path, 'sra')
    run_cmd(['prefetch',
             "--option-file", acc_list,
             "--output-directory", sra_dir,
             "--max-size", "100G"])
    
    repo_root = Path(os.environ.get("NCBI_VDB_REPOSITORY_ROOT", Path.home() / "ncbi"))
    src_dir   = repo_root / "public" / "sra"
    dst_dir   = Path(dir_path) / "sra"
    
    if src_dir.exists():
        for sra_file in src_dir.glob("*.sra"):
            shutil.move(str(sra_file), str(dst_dir / sra_file.name))
            print(f"moved {sra_file} -> {dst_dir/sra_file.name}")

def sra_to_fastq(dir_path: str, as_single):
    print(f"converting files from .sra to .fastq.")
    sra_root = os.path.join(dir_path, "sra")
    fastq_path = os.path.join(dir_path, "fastq")
    
    if not os.path.exists(sra_root):
        if os.path.exists(fastq_path) and os.listdir(fastq_path):
             pass 
        else:
             print("Warning: SRA directory missing and no FastQ files found.")
        return ReadsData(dir_path, fwd=True, rev=False)

    for item in tqdm(os.listdir(sra_root), desc="converted files"):
        full_path = os.path.join(sra_root, item)
        if os.path.isdir(full_path):
             nested_files = os.listdir(full_path)
             if not nested_files: continue
             target_sra = os.path.join(full_path, nested_files[0])
        else:
             target_sra = full_path
             
        run_cmd(["fasterq-dump", "--split-files", target_sra, "-O", fastq_path])

    fastqs = sorted([f for f in os.listdir(fastq_path) if f.endswith(".fastq") or f.endswith(".fq")])[:3]
    if not fastqs:
        return ReadsData(dir_path, fwd=True, rev=False)

    prefixes = [f.split("_")[0] for f in fastqs]
    prefix_counts = Counter(prefixes)
    all_have_two = any(count == 2 for count in prefix_counts.values())

    if all_have_two:
        if as_single:
            run_cmd([f"rm {os.path.join(fastq_path, '*_2.fastq')}"])
            print("Single reads requested - reverse reads deleted.")
            return ReadsData(dir_path, fwd=True, rev=False)
        else:
            return ReadsData(dir_path, fwd=True, rev=True)

    return ReadsData(dir_path, fwd=True, rev=False)

def create_manifest(reads_data: ReadsData):
    base_dir = os.path.abspath(reads_data.dir_path)
    fastq_path = os.path.join(base_dir, "fastq")
    manifest_path = os.path.join(base_dir, 'manifest.tsv')
    
    if os.path.exists(os.path.join(base_dir, "sra")) and os.listdir(os.path.join(base_dir, "sra")):
        names = sorted([os.path.splitext(os.path.basename(f))[0] for f in os.listdir(os.path.join(base_dir, "sra"))])
    else:
        fq_files = sorted([f for f in os.listdir(fastq_path) if f.endswith(".fastq")])
        names = sorted(list(set([f.split("_")[0] for f in fq_files])))

    with open(manifest_path, 'w', newline='') as manifest:
        tsv_writer = csv.writer(manifest, delimiter='\t')
        if not reads_data.rev:
            tsv_writer.writerow(["SampleID", "absolute-filepath"])
            for n in names:
                f_path = os.path.join(fastq_path, f"{n}.fastq")
                if not os.path.exists(f_path): f_path = os.path.join(fastq_path, f"{n}_1.fastq")
                if os.path.exists(f_path): tsv_writer.writerow([n, os.path.abspath(f_path)])
        else:
            tsv_writer.writerow(["SampleID", "forward-absolute-filepath", "reverse-absolute-filepath"])
            for n in names:
                f1 = os.path.join(fastq_path, f"{n}_1.fastq")
                f2 = os.path.join(fastq_path, f"{n}_2.fastq")
                if os.path.exists(f1) and os.path.exists(f2):
                    tsv_writer.writerow([n, os.path.abspath(f1), os.path.abspath(f2)])

def qiime_import(reads_data: ReadsData):
    paired = reads_data.rev and reads_data.fwd
    out_path = os.path.join(reads_data.dir_path, "qza", f"demux-{'paired' if paired else 'single'}-end.qza")
    command = [
        "qiime", "tools", "import",
        "--type", f"SampleData[{'PairedEndSequencesWithQuality' if paired else 'SequencesWithQuality'}]",
        "--input-path", os.path.join(reads_data.dir_path, 'manifest.tsv'),
        "--input-format", "PairedEndFastqManifestPhred33V2" if paired else "SingleEndFastqManifestPhred33V2",
        "--output-path", out_path,
    ]
    run_cmd(command)
    return out_path

def qiime_demux(reads_data: ReadsData, qza_file_path: str, dataset_id):
    vis_path = os.path.join(reads_data.dir_path, "vis", dataset_id + ".qzv")
    run_cmd(["qiime", "demux", "summarize", "--i-data", qza_file_path, "--o-visualization", vis_path])
    return vis_path

def metaphlan_extraction(reads_data, dataset_id, threads=8):
    paired = reads_data.rev and reads_data.fwd
    fastq_path = os.path.join(reads_data.dir_path, "fastq")
    export_path = os.path.join(reads_data.dir_path, "export")
    final_output_path = os.path.join(export_path, f'{dataset_id}_final.txt')
    
    os.makedirs(export_path, exist_ok=True)
    open(final_output_path, 'a').close()

    fastq_files = sorted([a for a in os.listdir(fastq_path) if a.endswith(".fastq")])
    qza_dir = os.path.join(reads_data.dir_path, 'qza')

    # FIXED: Enforce the 'Jun23' database which is compatible with HUMAnN 3.x/4.x
    target_index = "mpa_vJun23_CHOCOPhlAnSGB_202307"

    if paired:
        r1_files = [f for f in fastq_files if "_1.fastq" in f]
        for f1 in tqdm(r1_files):
            fname = f1.replace("_1.fastq", "")
            f2 = f"{fname}_2.fastq"
            p1 = os.path.join(fastq_path, f1)
            p2 = os.path.join(fastq_path, f2)
            if not os.path.exists(p2): continue

            map_out = os.path.join(fastq_path, f"{fname}.bowtie2.bz2")
            profile_out = os.path.join(qza_dir, f'{fname}_profile.txt')

            cmd = (
                f"metaphlan {p1},{p2} --input_type fastq "
                f"--nproc {threads} --bowtie2out {map_out} -o {profile_out} "
                f"--index {target_index}"
            )
            run_cmd([cmd])
    else:
        for f in tqdm(fastq_files):
            p = os.path.join(fastq_path, f)
            profile_out = os.path.join(qza_dir, f'{f}_profile.txt')
            map_out = os.path.join(fastq_path, f"{f}.bowtie2.bz2")
            
            cmd = (
                f"metaphlan {p} --input_type fastq "
                f"--nproc {threads} --bowtie2out {map_out} -o {profile_out} "
                f"--index {target_index}"
            )
            run_cmd([cmd])

    profile_files = [os.path.join(qza_dir, f) for f in os.listdir(qza_dir) if f.endswith("_profile.txt")]
    if profile_files:
        with open(final_output_path, 'w') as out:
            merge(profile_files, out, False)

def metaphlan_txt_csv(reads_data, dataset_id):
    export_path = os.path.join(reads_data.dir_path, "export")
    input_file = os.path.join(export_path, f"{dataset_id}_final.txt")
    output_file = os.path.join(export_path, f"{dataset_id}_final_table.csv")
    
    if not os.path.exists(input_file):
        print(f"Error: {input_file} not found.")
        return

    with open(input_file, 'r') as txt:
        lines = txt.readlines()
    if not lines: return

    headers = lines[0].strip().split('\t')
    data = [line.strip().split('\t') for line in lines[1:]]
    headers = [h.replace('|', ',') for h in headers]
    data = [[e.replace('|', ',') for e in row] for row in data]
    transposed = list(map(list, zip(*data)))

    with open(output_file, 'w', newline='') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(headers)
        writer.writerows(transposed)

def run_cleaning_pipeline(dir_path, threads=8, backup=True):
    run_dehost_pipeline(dir_path, threads=threads)
    
    base = Path(dir_path)
    fastq_dir = base / "fastq"
    fastq_clean = base / "fastq_clean"
    fastq_raw = base / "fastq_raw"
    
    if not fastq_clean.exists() or not list(fastq_clean.iterdir()):
        print("Warning: Cleaning pipeline produced no output. Using raw reads.")
        return

    if backup:
        fastq_raw.mkdir(exist_ok=True)
        for f in fastq_dir.iterdir():
            if f.is_file():
                shutil.move(str(f), str(fastq_raw / f.name))
    else:
        for f in fastq_dir.iterdir():
            if f.is_file():
                os.remove(f)
    
    swapped_count = 0
    for f in fastq_clean.iterdir():
        if f.suffix not in ['.fastq', '.fq']: continue
        
        if "paired_1" in f.name:
            clean_name = f.name.split("_kneaddata")[0] + "_1.fastq"
            shutil.copy2(str(f), str(fastq_dir / clean_name))
            swapped_count += 1
        elif "paired_2" in f.name:
            clean_name = f.name.split("_kneaddata")[0] + "_2.fastq"
            shutil.copy2(str(f), str(fastq_dir / clean_name))
            swapped_count += 1
        elif "_1.fastq" in f.name and "paired" not in f.name: 
             shutil.copy2(str(f), str(fastq_dir / f.name))
             swapped_count +=1
            
    print(f"Swapped {swapped_count} cleaned paired files into active fastq folder.")

# --- Main Logic ---

def visualization(acc_list, dataset_id, data_type, verbose_print, specific_location, as_single, 
                  threads: int = 8, pathways: str = "no", clean: bool = False):
    
    verbose_print("\n" + datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S'))
    dir_name = f"{dataset_id}-{datetime.datetime.now().strftime('%d-%m-%Y_%H-%M-%S')}"

    check_conda_qiime2() 

    verbose_print("Checking inputs...")
    check_input(acc_list)

    verbose_print(f"Creating directory: {dir_name}")
    dir_path = create_dir(dir_name, specific_location)
    json_file_path = f"{dir_path}/metadata.json"

    verbose_print("Start prefetch (1/6)")
    download_data_from_sra(dir_path, acc_list)
    
    data_json = {"dir_path": dir_path, "dataset_id": dataset_id}
    with open(json_file_path, "w") as jf: json.dump(data_json, jf)

    verbose_print("Start conversion (2/6)")
    reads_data = sra_to_fastq(dir_path, as_single)

    if clean:
        verbose_print("Running KneadData cleaning...")
        run_cleaning_pipeline(dir_path, threads=threads, backup=True)

    data_json.update({
        "type": data_type,
        "read_data_fwd": reads_data.fwd,
        "read_data_rev": reads_data.rev
    })
    with open(json_file_path, "w") as jf: json.dump(data_json, jf)

    if data_type == '16S' or data_type == '18S':
        verbose_print("Start 16S flow...")
        create_manifest(reads_data)
        qza_path = qiime_import(reads_data)
        vis_path = qiime_demux(reads_data, qza_path, dataset_id)
        pickle.dump(reads_data, open(os.path.join(reads_data.dir_path, "reads_data.pkl"), "wb"))
        return reads_data.dir_path
    else:
        verbose_print("Start Shotgun flow...")
        metaphlan_extraction(reads_data, dataset_id, threads)
        metaphlan_txt_csv(reads_data, dataset_id)
        if pathways == "yes":
            run_humann_pipeline(dir_path, dataset_id, threads)
        print("Shotgun analysis finished successfully.")

def visualization_continue_fastq(dataset_id, continue_path, data_type, verbose_print, specific_location, 
                                  threads, pathways, clean):
    check_conda_qiime2()
    continue_path = Path(continue_path)
    
    if not (continue_path / "fastq").exists():
        print(f"Warning: Fastq directory not found in {continue_path}")

    reads_data = sra_to_fastq(str(continue_path), as_single=False)

    if clean:
        run_cleaning_pipeline(str(continue_path), threads=threads, backup=True)

    if data_type in ['16S', '18S']:
         create_manifest(reads_data)
         qza_path = qiime_import(reads_data)
         qiime_demux(reads_data, qza_path, dataset_id)
    else:
         metaphlan_extraction(reads_data, dataset_id, threads)
         metaphlan_txt_csv(reads_data, dataset_id)
         if pathways == "yes":
             run_humann_pipeline(str(continue_path), dataset_id, threads)

def visualization_continue(dataset_id, continue_path, data_type, verbose_print, specific_location, threads, pathways, clean):
    check_conda_qiime2()
    
    if clean:
        run_cleaning_pipeline(str(continue_path), threads=threads, backup=True)
    
    try:
        with open(os.path.join(continue_path, 'metadata.json'), 'r') as jf:
            meta = json.load(jf)
            reads_data = ReadsData(str(continue_path), fwd=meta.get("read_data_fwd", True), rev=meta.get("read_data_rev", False))
    except:
        reads_data = ReadsData(str(continue_path), fwd=True, rev=False)

    if data_type in ['16S', '18S']:
         create_manifest(reads_data)
         qza_path = qiime_import(reads_data)
         qiime_demux(reads_data, qza_path, dataset_id)
    else:
         metaphlan_extraction(reads_data, dataset_id, threads)
         metaphlan_txt_csv(reads_data, dataset_id)
         if pathways == "yes":
             run_humann_pipeline(str(continue_path), dataset_id, threads)