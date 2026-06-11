# YaMAS (YOLO Microbiome Analysis System)

### Table of Contents
- [Pipeline Overview](#pipeline-overview)
- [Installation & Dependencies](#installation--dependencies)
- [Flags & Options](#flags--options)
- [Output Structure](#output-structure)
- [Downloading a project](#downloading-a-project)
    - [Download from NCBI SRA](#download-from-ncbi-sra)
    - [Continue data downloading](#continue-data-downloading)
    - [Download from ENA](#download-from-ena)
    - [Download using fastq files](#download-using-fastq-files)
- [Exporting a project (only for 16S/18S)](#exporting-a-project-only-for-16s18s)
    - [The --export-as-single Flag](#the---export-as-single-flag)
- [Troubleshooting](#troubleshooting)
- [Arguments and configurations](#arguments-and-configurations)


YaMAS is a standalone microbiome analysis pipeline for downloading and processing DNA datasets from NCBI SRA and ENA. It is developed by the [YOLO lab team](https://yolo.math.biu.ac.il), and is designed to be simple, efficient, and easy to use for non-programmers. The updated 16S/18S module eliminates QIIME2 dependencies, using **VSEARCH** for rapid amplicon processing, **cutadapt** for demultiplexing and trimming, **SINTAX** for taxonomy classification, and **FastTree** for phylogeny — all within a single Python 3.13 conda environment.

## Pipeline Overview

The YaMAS pipeline consists of several stages, depending on the sequencing type:

**Supported input sources:**
- **Project ID** from SRA/ENA (automatic download)
- **Local FASTQ file** (process without download)
- **Existing FASTQ folder** (`--continue_from_fastq`)
- **Existing project folder** (`--continue_from`)

**For Shotgun datasets:**
1. **Download** dataset from SRA/ENA
2. **KneadData** *(optional, `--clean`)* – Host removal and quality control
3. **MetaPhlAn4** – Taxonomic profiling
4. **HUMAnN3** *(if `--pathways yes` is set)* – Functional profiling and pathway analysis
5. **Export & Visualization** – Merged abundance tables and QC reports

**For 16S/18S datasets:**
1. **Download** dataset
2. **Demultiplexing** – cutadapt (EMP / Qiita / manifest formats)
3. **Quality Control** – FastQC + MultiQC HTML reports
4. **Trimming** – cutadapt 5′ trim and read truncation
5. **Denoising** – VSEARCH UNOISE3 (ASV detection, chimera removal, 99% OTU clustering)
6. **Taxonomy** – VSEARCH SINTAX with 0.8 bootstrap confidence cutoff
7. **Phylogeny** – MAFFT alignment + FastTree (GTR model, midpoint-rooted)
8. **Export** – OTU table, taxonomy table, phylogenetic tree

> **Note:** HUMAnN3 integration is available **only** for Shotgun datasets and runs immediately after MetaPhlAn4.

---

## Installation & Dependencies

Before proceeding with the installation and execution of YaMAS, please ensure that you have a clean environment set up on your system. Follow the steps below:

### Step 1: Install Conda with miniforge3

If you don't already have Conda installed, download and install [miniforge3](https://github.com/conda-forge/miniforge):

```bash
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh -b -p $HOME/miniforge3
$HOME/miniforge3/bin/conda init
source ~/.bashrc
```

### Step 2: Clone the Repository

```bash
git clone -b main https://github.com/louzounlab/YaMAS.git
cd YaMAS
```

### Step 3: Create the Conda Environment and Install

```bash
conda create -n yamas_env python=3.13 -y
conda activate yamas_env

conda install -y -c bioconda -c conda-forge \
    sra-tools=3.2.1 \
    entrez-direct=24.0 \
    bowtie2=2.5.5 \
    diamond=2.1.23 \
    glpk=5.0 \
    fastqc=0.12.1 \
    vsearch=2.30.5 \
    mafft=7.525 \
    fasttree=2.2.0

pip install -e .
```

`pip install -e .` installs all Python dependencies from `setup.py`, including:

| Package | Version | Purpose |
|---------|---------|---------|
| `pandas` | ≥1.0.0 | Data manipulation |
| `tqdm` | ≥4.0.0 | Progress bars |
| `metaphlan` | 4.0.6 | Shotgun taxonomic profiling |
| `humann` | 3.9 | Functional pathway analysis |
| `kneaddata` | 0.12.4 | QC and host decontamination |
| `biopython` | 1.86 | Phylogeny tree I/O |
| `cutadapt` | 5.2 | Demultiplexing and adapter trimming |
| `multiqc` | 1.33 | QC report aggregation |

### Step 4: Set Your Database Root Directory

```bash
export DB_ROOT="/path/to/databases"  # Change this to your desired location
```

### Step 5a: Download Shotgun Databases (Optional)

**If you plan to use the Shotgun pathway, download these databases:**

```bash
# KneadData: host removal
kneaddata_database --download human_genome bowtie2 $DB_ROOT/kneaddata_db

# HUMAnN3: functional profiling
humann_databases --download chocophlan full $DB_ROOT/humann_db --update-config yes
humann_databases --download uniref uniref90_diamond $DB_ROOT/humann_db --update-config yes
humann_databases --download utility_mapping full $DB_ROOT/humann_db --update-config yes

# MetaPhlAn4: taxonomic profiling
metaphlan --install --index mpa_vJun23_CHOCOPhlAnSGB_202307 --bowtie2db $DB_ROOT/metaphlan_db
```

> **Note:** These databases are large (~20+ GB total). If you have already downloaded them, skip to Step 7.

### Step 5b: Download the 16S Reference Classifier

Download a VSEARCH `--sintax` compatible FASTA reference:

1. Visit: https://www.drive5.com/usearch/manual/sintax_downloads.html
2. Download `gg_16s_13.5.fa.gz` (GreenGenes 13.5, 99% OTUs)
3. Create the classifier directory and extract:

```bash
mkdir -p $DB_ROOT/16S_classifiers
cp /path/to/gg_16s_13.5.fa.gz $DB_ROOT/16S_classifiers/
gunzip $DB_ROOT/16S_classifiers/gg_16s_13.5.fa.gz
```

### Step 6: Persist Environment Variables

```bash
conda env config vars set YAMAS_HOST_DB=$DB_ROOT/kneaddata_db
conda env config vars set HUMANN_DB_DIR=$DB_ROOT/humann_db
conda env config vars set METAPHLAN_DB_DIR=$DB_ROOT/metaphlan_db
```

### Step 7: Configure HUMAnN Paths (If Reusing Existing Databases)

If you have already downloaded HUMAnN databases, ensure they are configured:

```bash
humann_config --update database_folders nucleotide $DB_ROOT/humann_db/chocophlan
humann_config --update database_folders protein $DB_ROOT/humann_db/uniref
humann_config --update database_folders utility_mapping $DB_ROOT/humann_db/utility_mapping
```

### Step 8: Refresh Environment

```bash
conda deactivate
conda activate yamas_env
```

**You're all set!** Verify the installation:

```bash
yamas --version
vsearch --version
mafft --version
fastqc -v
```

---

## Flags & Options

| Flag | Description |
|------|-------------|
| `--download <ACCESSIONS...>` | Download dataset(s) from SRA/ENA |
| `--type <16S/18S/Shotgun>` | Type of sequencing data |
| `--export <PATH> <TYPE> <TRIM> <TRUNC> <CLASSIFIER> <THREADS>` | Export 16S/18S analysis (OTU, taxonomy, tree) |
| `--export-as-single` | Force single-end processing when PE merge fails (use with `--export`) |
| `--as_single` | Treat paired-end reads as single-end at the download stage |
| `--pathways <yes/no>` | Enable HUMAnN3 for pathway profiling (Shotgun only) |
| `--clean` | Run KneadData for host removal and QC (Shotgun only) |
| `--threads <N>` | Number of threads to use |
| `--acc_list <FILE>` | Path to a file with one accession per line |
| `--verbose` | Print detailed progress information |
| `--continue_from_fastq <ID> <PATH> <TYPE>` | Continue processing from an existing FASTQ folder |
| `--continue_from <ID> <PATH> <TYPE>` | Continue processing from an existing dataset folder |
| `--config <FILE>` | Load configuration from a JSON file |
| `--ready <OS_TYPE>` | Initialize the YaMAS environment (Ubuntu/CentOS) |

---

## Output Structure

After running YaMAS, the project folder will contain:

```
<PROJECT_ID>-<DATE>_<TIME>/
│
├── sra/                      # Raw SRA files (from prefetch)
├── fastq_raw/                # Original FASTQ files (before any QC or host removal)
├── fastq_clean/              # Clean FASTQ files (after KneadData, if --clean is used)
├── fastq/                    # Active FASTQ folder used for downstream steps
│
├── knead_out/                # KneadData logs and QC reports (if --clean is set)
│   ├── *_kneaddata.log
│   ├── *.fastqc.html
│   └── multiqc_report.html
│
├── qza/                      # Intermediate files for 16S/18S processing
│   ├── trimmed/              # Quality-trimmed FASTQ files
│   ├── derep.fasta           # Deduplicated sequences
│   ├── denoised.fasta        # Denoised ASVs (UNOISE3 output)
│   └── rep-seqs-dn-99.fasta  # 99% OTU representative sequences
│
├── vis/                      # QC visualization reports
│   └── multiqc_report.html   # FastQC summary (read counts, qualities, lengths)
│
├── export/                   # Final analysis outputs (16S/18S, created with --export)
│   ├── otu.csv               # OTU abundance table (OTUs × samples)
│   ├── taxonomy.csv          # Taxonomy assignments with confidence scores
│   ├── tree.nwk              # Phylogenetic tree (Newick, midpoint-rooted)
│   └── otu_padding.csv       # OTU table padded with tree-only tips
│
├── metaphlan_results/        # Taxonomic profiles (Shotgun only)
│   ├── *_profile.tsv
│   └── merged_abundance.tsv
│
└── humann_results/           # Functional profiles (Shotgun only, if --pathways yes)
    ├── *_pathabundance.tsv
    ├── *_pathcoverage.tsv
    ├── *_pathabundance_stratified.tsv
    └── merged_pathways.tsv
```

**Key output files:**
- **16S/18S**: `export/otu.csv`, `export/taxonomy.csv`, `export/tree.nwk`, `vis/multiqc_report.html`
- **Shotgun**: `metaphlan_results/merged_abundance.tsv`, `humann_results/merged_pathways.tsv` *(if `--pathways yes`)*

---

## Downloading a project

To download a project from **NCBI SRA** or **ENA**, use one of the following templates:

### Download from NCBI SRA

```bash
yamas --download <dataset_id> --type <data_type> --pathways <pathways>
```

Arguments:
- `dataset_id`: the accession ID from NCBI SRA. For example: `PRJEB01234` or `SRR11415443`
- `data_type`: choose one of: `16S` / `18S` / `Shotgun`
- `pathways`: generate HUMAnN3 pathway tables (Shotgun only). Choose: `yes` / `no`

To download multiple accessions using a list file:
```bash
yamas --download --type 16S --acc_list /path/to/accessions.txt --verbose
```

### Continue data downloading

1. Continue from **after** SRA download, **before** FASTQ conversion:
```bash
yamas --continue_from_fastq <dataset_id> <project_path> <data_type> --pathways <pathways>
```

Arguments:
- `dataset_id`: the accession ID. For example: `PRJEB01234`
- `project_path`: path to the existing project directory created by YaMAS
- `data_type`: choose one of: `16S` / `18S` / `Shotgun`
- `pathways`: choose: `yes` / `no`

2. Continue from **after** FASTQ conversion (skip download and SRA conversion):
```bash
yamas --continue_from <dataset_id> <project_path> <data_type> --pathways <pathways>
```

Arguments:
- `dataset_id`: the accession ID. For example: `PRJEB01234`
- `project_path`: path to the existing project directory created by YaMAS
- `data_type`: choose one of: `16S` / `18S` / `Shotgun`
- `pathways`: choose: `yes` / `no`

### Download from ENA

```bash
yamas --qiita <preprocessed_fastq_path> <metadata_path> <data_type>
```

Arguments:
- `preprocessed_fastq_path`: path to the preprocessed FASTQ file. Rename the file to `forward.fastq.gz`
- `metadata_path`: path to the metadata file. Rename the file to `metadata.tsv`
- `data_type`: choose one of: `16S` / `18S`

The output will be created in the same folder as the FASTQ and metadata files; it is recommended to organise inputs in a dedicated directory.

### Download using fastq files

```bash
yamas --fastq <preprocessed_fastq_path> <barcode_path> <metadata_path> <data_type>
```

Arguments:
- `preprocessed_fastq_path`: path to the preprocessed FASTQ file, renamed to `preprocessed_fastq.gz`
- `barcode_path`: path to the barcode file, renamed to `barcodes.fastq.gz`
- `metadata_path`: path to the metadata TSV file (must contain a `barcode` column), renamed to `metadata.tsv`
- `data_type`: choose one of: `16S` / `18S` / `Shotgun`

---

## Exporting a project (only for 16S/18S)

After downloading a 16S/18S dataset, run `--export` to produce the OTU table, taxonomy assignments, and phylogenetic tree.

#### Paired-end (reads successfully merge)

```bash
yamas --export <project_path> <data_type> <trim_fwd>,<trim_rev> <trunc_fwd>,<trunc_rev> <classifier_file> <threads>
```

> Example: `yamas --export /PATH/PRJNA633959-30-09-2025_17-55-54 16S 13,13 175,175 $DB_ROOT/16S_classifiers/gg_16s_13.5.fa 24`

#### Single-end, or paired-end fallback with `--export-as-single`

```bash
yamas --export <project_path> <data_type> <trim> <trunc> <classifier_file> <threads> --export-as-single
```

> Example: `yamas --export /PATH/PRJNA633959-30-09-2025_17-55-54 16S 13 175 $DB_ROOT/16S_classifiers/gg_16s_13.5.fa 24 --export-as-single`

| Argument | Format | Description |
|----------|--------|-------------|
| `project_path` | string | Path to the project directory created by YaMAS |
| `data_type` | `16S` \| `18S` | Analysis type |
| `trim` | `"13,13"` (PE) or `"13"` (SE/fallback) | Bases to remove from the 5′ end of each read |
| `trunc` | `"175,175"` (PE) or `"175"` (SE/fallback) | Final read length after truncation |
| `classifier_file` | path | Path to the VSEARCH SINTAX-compatible FASTA reference database |
| `threads` | int | Number of CPU threads (default: 12) |
| `--export-as-single` | flag | Process only forward reads — see below |

### The `--export-as-single` Flag

By default, YaMAS attempts to **merge Paired-End (PE) reads** using VSEARCH. This is optimal when forward and reverse reads overlap sufficiently. Two scenarios require the single-end fallback:

**Scenario 1: Long Amplicon + Short Read Length**

> Example: The V3-V4 16S region is ~420 bp, but 150 bp PE reads produce only 300 bp combined — no physical overlap is possible.
>
> Result: VSEARCH `mergepairs` drops **99.9% of reads** with the error `"too few kmers found on same diagonal"`, making the dataset unusable.
>
> Solution: Use `--export-as-single` to process forward reads only.

**Scenario 2: Poor Reverse Read Quality**

> Reverse reads have high error rates due to sequencing issues. PE merging may "succeed" but produces degraded consensus sequences that inflate the OTU table with artifacts.
>
> Solution: Use `--export-as-single` to discard the low-quality reverse reads entirely.

**What `--export-as-single` does:**

1. **Ignores reverse reads completely** – R2 FASTQ files (`_2.fastq`) are not processed
2. **Processes only forward reads** – `_1.fastq` files are retained in canonical orientation
3. **Prevents OTU table inflation** – processing `_1` and `_2` independently without this flag would generate reverse-complement duplicates, inflating the feature table by ~2×

**Decision tree:**

```
Do you have both _1 and _2 FASTQ files?
├─ YES, and PE merge succeeds?            → Use comma-separated params (default PE mode)
├─ YES, but 99%+ reads dropped at merge? → Add --export-as-single flag
├─ YES, but R2 quality is very poor?      → Add --export-as-single flag
└─ NO (native single-end data)?           → Use single integers (no flag needed)
```

---

## Troubleshooting

**Error: `"missing 5 arguments"`**
> Cause: Paired-end data detected but single-integer trim/trunc values passed without `--export-as-single`.
> ```bash
> # Paired-end (comma-separated):
> yamas --export DIR 16S 15,15 150,150 classifier.fa 24
>
> # Single-end fallback:
> yamas --export DIR 16S 15 150 classifier.fa 24 --export-as-single
> ```

**Error: `"too few kmers found on same diagonal"` (99.9% read dropout)**
> Cause: Amplicon is too long for the PE reads to overlap.
> Fix: Use `--export-as-single` to process forward reads only.

**Error: `"Reference database not found at..."`**
> Cause: The classifier FASTA path is incorrect.
> Fix: Verify the file exists: `ls -lh $DB_ROOT/16S_classifiers/gg_16s_13.5.fa`

**Very few OTUs / empty `taxonomy.csv`**
> Status: Fixed in the current version. The prevalence filter is now capped to the actual number of samples, preventing all OTUs from being removed in small or single-sample datasets.

**MultiQC report shows very low read counts after trimming**
> Cause: Trim/trunc parameters are too aggressive.
> Fix: Relax the `trunc` cutoff (e.g., increase from `150` to `175`):
> ```bash
> yamas --export DIR 16S 13 175 classifier.fa 24 --export-as-single
> ```

---

## Arguments and configurations

1. **config**: You can provide a configuration file (`--config <file>`) to save data to a different folder and adjust other pipeline settings.
2. **verbose**: Use `--verbose` to print detailed progress during download and export operations (highly recommended).
3. **Multiple accessions**: Listing more than one project ID, or providing `--acc_list`, will download them one by one into separate timestamped folders.

---

## Cite us

If you are using our package, YaMAS for **any** purpose, please cite us; Shtossel Oshrit, Sondra Turjeman, Alona Riumin, Michael R. Goldberg, Arnon Elizur, Yarin Bekor, Hadar Mor, Omry Koren, and Yoram Louzoun. "Recipient-independent, high-accuracy FMT-response prediction and optimization in mice and humans." Microbiome 11, no. 1 (2023): 181. https://link.springer.com/article/10.1186/s40168-023-01623-w
