from __future__ import annotations

import csv
import datetime
import os
import pickle
import pandas as pd

# Removed top-level imports of biom, qiime2, skbio, Bio to prevent crashes in Shotgun env
from .utilities import ReadsData, run_cmd, download_classifier_url, check_conda_qiime2

nodes_names = []

def _lazy_import_16s_deps():
    """Helper to import 16S dependencies only when strictly needed."""
    try:
        from biom import Table, load_table
        import qiime2 as q2
        from skbio import TreeNode
        from Bio import Phylo
        return Table, load_table, q2, TreeNode, Phylo
    except ImportError as e:
        raise ImportError(
            "Could not import QIIME2 or BioPython libraries. "
            "These are required for 16S export but NOT for Shotgun analysis. "
            "Please install them or use a QIIME2 environment."
        ) from e

def trim_trunc_check(reads_data: ReadsData, trim: int | tuple[int, int], trunc: int | tuple[int, int]):
    if reads_data.fwd and reads_data.rev:
        if not isinstance(trim, tuple) or not isinstance(trunc, tuple):
            raise TypeError(f"Trim/Trunc must be tuples for paired reads.")
        if len(trim) != 2 or len(trunc) != 2:
            raise ValueError("Trim/Trunc must be tuples of length 2.")
    if not isinstance(trim, int) or not isinstance(trunc, int):
        raise TypeError("Trim/Trunc must be integers for single reads.")

def classifier_exists(classifier_path: str):
    if not (os.path.exists(classifier_path) and os.path.isfile(classifier_path)):
        raise FileNotFoundError(f"Classifier not found! Download it from: {download_classifier_url()}\n")

def qiime_dada2(reads_data: ReadsData, input_path: str, left, right, threads: int = 12):
    # This uses shell commands via run_cmd, so it is safe without imports
    paired = reads_data.fwd and reads_data.rev
    trim_range = ["--p-trim-left-f", str(left.split(',')[0]), "--p-trim-left-r", str(left.split(',')[1])] if paired \
        else ["--p-trim-left", str(left)]
    trunc_range = ["--p-trunc-len-f", str(right.split(',')[0]), "--p-trunc-len-r", str(right.split(',')[1])] if paired \
        else ["--p-trunc-len", str(right)]
    
    command = [
        "qiime", "dada2", "denoise-paired" if paired else "denoise-single",
        "--i-demultiplexed-seqs", input_path,
    ] + trim_range + trunc_range + [
        "--o-table", os.path.join(reads_data.dir_path, "qza", "dada2_table.qza"),
        "--p-n-threads", str(threads),
        "--p-chimera-method", "consensus",
        "--o-representative-sequences", os.path.join(reads_data.dir_path, "qza", "dada2_rep-seqs.qza"),
        "--o-denoising-stats", os.path.join(reads_data.dir_path, "qza", "dada2_denoising-stats.qza"),
        "--verbose"
    ]
    run_cmd(command)

def cluster_features(reads_data: ReadsData):
    qza_path = lambda filename: os.path.join(reads_data.dir_path, "qza", filename)
    command = [
        "qiime", "vsearch", "cluster-features-de-novo",
        "--i-table", qza_path("dada2_table.qza"),
        "--i-sequences", qza_path("dada2_rep-seqs.qza"),
        "--p-perc-identity", "0.99",
        "--o-clustered-table", qza_path("table-dn-99.qza"),
        "--o-clustered-sequences", qza_path("rep-seqs-dn-99.qza")
    ]
    run_cmd(command)

def assign_taxonomy(reads_data: ReadsData, data_type, classifier_path: str):
    qza_path = lambda filename: os.path.join(reads_data.dir_path, "qza", filename)
    out_name = "gg-13-8-99-nb-classified.qza" if data_type == '16S' else "silva-132-99-nb-classifier.qza"
    command = [
        "qiime", "feature-classifier", "classify-sklearn",
        "--i-reads", qza_path("rep-seqs-dn-99.qza"),
        "--i-classifier", classifier_path,
        "--o-classification", qza_path(out_name)
    ]
    run_cmd(command)

def clean_taxonomy1(reads_data: ReadsData, data_type):
    qza_path = lambda filename: os.path.join(reads_data.dir_path, "qza", filename)
    tax_name = "gg-13-8-99-nb-classified.qza" if data_type == '16S' else "silva-132-99-nb-classifier.qza"
    command = [
        "qiime", "taxa", "filter-table",
        "--i-table", qza_path("table-dn-99.qza"),
        "--i-taxonomy", qza_path(tax_name),
        "--p-exclude", "mitochondria,chloroplast",
        "--o-filtered-table", qza_path("clean_table.qza")
    ]
    run_cmd(command)

def clean_taxonomy2(reads_data: ReadsData):
    qza_path = lambda filename: os.path.join(reads_data.dir_path, "qza", filename)
    command = [
        "qiime", "feature-table", "filter-features",
        "--i-table", qza_path("clean_table.qza"),
        "--p-min-samples", "3",
        "--p-min-frequency", "10",
        "--o-filtered-table", qza_path("feature-frequency-filtered-table.qza")
    ]
    run_cmd(command)

def export_otu(reads_data: ReadsData):
    output_file = os.path.join(reads_data.dir_path, "exports", "otu.tsv")
    command = [
        "qiime", "tools", "export",
        "--input-path", os.path.join(reads_data.dir_path, "qza", "feature-frequency-filtered-table.qza"),
        "--output-path", os.path.join(reads_data.dir_path, "exports")
    ]
    run_cmd(command)
    command = ["biom", "convert", "-i", os.path.join(reads_data.dir_path, "exports", "feature-table.biom"), "-o", output_file, "--to-tsv"]
    run_cmd(command)

def export_taxonomy(reads_data: ReadsData, data_type, classifier_file_path):
    output_file = os.path.join(reads_data.dir_path, "exports", "tax.tsv")
    tax_name = "gg-13-8-99-nb-classified.qza" if data_type == '16S' else "silva-132-99-nb-classifier.qza"
    command = [
        "qiime", "tools", "export",
        "--input-path", os.path.join(reads_data.dir_path, "qza", tax_name),
        "--output-path", output_file
    ]
    run_cmd(command)

def export_phylogeny(reads_data: ReadsData):
    input_file_path = os.path.join(reads_data.dir_path, 'qza', 'rep-seqs-dn-99.qza')
    output_file_path = os.path.join(reads_data.dir_path, 'qza', 'aligned-rep-seqs.qza')
    run_cmd(["qiime", "alignment", "mafft", "--i-sequences", input_file_path, "--o-alignment", output_file_path])

    input_file_path = output_file_path
    output_file_path = os.path.join(reads_data.dir_path, 'exports', 'fasttree-tree.qza')
    run_cmd(["qiime", "phylogeny", "fasttree", "--i-alignment", input_file_path, "--o-tree", output_file_path, "--verbose"])

    input_file_path = output_file_path
    output_file_path = os.path.join(reads_data.dir_path, 'exports', "fasttree-tree-rooted.qza")
    run_cmd(["qiime", "phylogeny", "midpoint-root", "--i-tree", input_file_path, "--o-rooted-tree", output_file_path])

def export_tree(reads_data: ReadsData):
    # LAZY IMPORT - Only loads if this specific 16S function is called
    Table, load_table, q2, TreeNode, Phylo = _lazy_import_16s_deps()
    
    input_file_path = os.path.join(reads_data.dir_path, 'exports', 'fasttree-tree-rooted.qza')
    output_file_path = os.path.join(reads_data.dir_path, 'exports', 'tree.nwk')
    tree = q2.Artifact.load(input_file_path)
    tntree = tree.view(TreeNode)
    tntree.prune()
    tntree.write(output_file_path)

def convert_to_csv(reads_data: ReadsData):
    otu_tsv = os.path.join(reads_data.dir_path, 'exports', 'otu.tsv')
    otu_csv = os.path.splitext(otu_tsv)[0] + '.csv'
    with open(otu_tsv, 'r', newline='') as tsv, open(otu_csv, 'w', newline='') as csv_out:
        writer = csv.writer(csv_out)
        reader = csv.reader(tsv, delimiter='\t')
        next(reader)
        for row in reader: writer.writerow(row)

    tax_tsv = os.path.join(reads_data.dir_path, 'exports', 'tax.tsv', 'taxonomy.tsv')
    tax_csv = os.path.splitext(tax_tsv)[0] + '.csv'
    with open(tax_tsv, 'r', newline='') as tsv, open(tax_csv, 'w', newline='') as csv_out:
        writer = csv.writer(csv_out)
        reader = csv.reader(tsv, delimiter='\t')
        for row in reader: writer.writerow(row)

def export_otu_padding_for_tree(reads_data: ReadsData):
    # LAZY IMPORT
    Table, load_table, q2, TreeNode, Phylo = _lazy_import_16s_deps()

    tree_file = os.path.join(reads_data.dir_path, 'exports', 'tree.nwk')
    otu_path = os.path.join(reads_data.dir_path, 'exports', 'otu.csv')
    otu_padding_path = os.path.join(reads_data.dir_path, 'exports', 'otu_padding.csv')
    otu = pd.read_csv(otu_path)
    asv_list = otu['#OTU ID'].tolist()
    tree = Phylo.read(tree_file, "newick")

    append_nodes_names(tree.root)
    
    in_tree_not_asv = [i for i in nodes_names if i not in asv_list]
    print(f"Adding {len(in_tree_not_asv)} ASVs from tree to CSV...")
    
    with open(otu_path, 'r', newline='') as orig, open(otu_padding_path, 'w', newline='') as pad:
        reader = csv.reader(orig)
        writer = csv.writer(pad)
        for row in reader: writer.writerow(row)
    
    with open(otu_padding_path, 'a', newline='') as pad:
        writer = csv.writer(pad)
        for i in in_tree_not_asv:
            new_row = [i] + [0] * len(otu.columns[1:])
            writer.writerow(new_row)

def append_nodes_names(clade):
    if clade.is_terminal():
        nodes_names.append(clade.name)
    else:
        for sub_clade in clade.clades:
            append_nodes_names(sub_clade)

def export(output_dir: str, data_type, trim, trunc, classifier_file_path: str, threads: int = 12):
    print(f"\n{datetime.datetime.now().strftime('%d/%m/%Y %H:%M:%S')}")
    print(f"### Exporting {data_type} ###")
    
    # Non-blocking check (from our updated utilities.py)
    check_conda_qiime2()
    
    reads_data: ReadsData = pickle.load(open(os.path.join(output_dir, "reads_data.pkl"), "rb"))
    classifier_exists(classifier_file_path)

    paired = reads_data.rev and reads_data.fwd
    output_path = os.path.join(reads_data.dir_path, "qza", f"demux-{'paired' if paired else 'single'}-end.qza")

    print("Running DADA2...")
    qiime_dada2(reads_data, output_path, left=trim, right=trunc, threads=threads)
    
    print("Clustering features...")
    cluster_features(reads_data)

    print("Assigning taxonomy...")
    assign_taxonomy(reads_data, data_type, classifier_file_path)

    run_cmd(["mkdir", "-p", os.path.join(reads_data.dir_path, "exports")])

    print("Cleaning taxonomy...")
    clean_taxonomy1(reads_data, data_type)
    clean_taxonomy2(reads_data)

    print("Exporting OTU & Taxonomy...")
    export_otu(reads_data)
    export_taxonomy(reads_data, data_type, classifier_file_path)

    print("Exporting Phylogeny & Tree...")
    export_phylogeny(reads_data)
    export_tree(reads_data)
    convert_to_csv(reads_data)
    export_otu_padding_for_tree(reads_data)
    print("Export finished.")