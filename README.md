# YaMAS: YOLO Microbiome Analysis System

A fast, standalone 16S/18S rRNA and Shotgun metagenomics pipeline for microbiome analysis. The updated 16S module eliminates QIIME2 dependencies, using **VSEARCH** for rapid amplicon processing, **cutadapt** for demultiplexing and trimming, **SINTAX** for taxonomy classification, and **FastTree** for phylogeny.

---

## Table of Contents

- [Project Overview](#project-overview)
- [Key Features](#key-features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Flags & Options](#flags--options)
- [Commands Reference](#commands-reference)
  - [Data Download](#1-data-download)
  - [Export & Analysis](#2-export--analysis)
  - [**The --export-as-single Flag (CRITICAL)**](#3-the---export-as-single-flag-critical)
- [Output Structure](#output-structure)
- [Troubleshooting](#troubleshooting)
- [Pipeline Architecture](#pipeline-architecture)

---

## Project Overview

YaMAS is a modular microbiome analysis framework supporting:

### 16S/18S rRNA Processing
- **Demultiplexing**: EMP, Qiita, or manifest-based formats
- **Quality Control**: FastQC + MultiQC for comprehensive QC reports
- **Denoising**: VSEARCH UNOISE3 algorithm (equivalent to DADA2 ASV detection)
- **Taxonomy**: SINTAX k-mer voting for fast, accurate classification
- **Phylogeny**: MAFFT alignment + FastTree tree construction

### Shotgun Metagenomics (Unchanged)
- MetaPHLAn4 for community profiling
- HUMAnN3 for functional pathway analysis
- KneadData for host-read removal

### Design Philosophy
- **No QIIME2**: Standalone tools provide speed, flexibility, and reproducibility
- **Single Python 3.13 environment**: All dependencies (vsearch, mafft, fasttree, cutadapt, multiqc) installed via conda
- **Data preservation**: Forward reads used when PE merging fails (new `--export-as-single` flag)

---

## Key Features

| Feature | Benefit |
|---------|---------|
| **VSEARCH denoising** | UNOISE3 algorithm: single-nucleotide ASVs comparable to DADA2 |
| **Direct FASTA processing** | No `.qza` container overhead; plain text for transparency |
| **SINTAX taxonomy** | Bootstrap voting; 0.8 confidence cutoff (conservative) |
| **FastQC + MultiQC** | Beautiful HTML QC reports (replaces QIIME2 `.qzv` files) |
| **Bioinformatically-sound fallback** | `--export-as-single` flag salvages PE data when merging fails |
| **Parallel execution** | Full multi-threading support via `--threads` parameter |

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
git clone -b main https://github.com/KfirPinto/Updated_yamas_3.13.git
cd Updated_yamas_3.13
```

### Step 3: Create Conda Environment & Install with Conda

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

This installs Python dependencies via `setup.py`:
- `pandas>=1.0.0` – data manipulation
- `tqdm>=4.0.0` – progress bars
- `metaphlan==4.0.6` – Shotgun taxonomic profiling
- `humann==3.9` – functional pathway analysis
- `kneaddata==0.12.4` – QC and host decontamination
- `biopython==1.86` – phylogeny tree I/O
- `cutadapt==5.2` – demultiplexing and adapter trimming
- `multiqc==1.33` – QC report aggregation

### Step 4: Download & Configure Reference Databases

Set your database root directory:

```bash
export DB_ROOT="/path/to/databases"  # Change this to your desired location
```

### Step 5a: Download Shotgun Databases (Optional)

**If you plan to use the Shotgun pathway, download these databases:**

```bash
# KneadData databases for host removal
kneaddata_database --download human_genome bowtie2 $DB_ROOT/kneaddata_db

# HUMAnN functional profiling databases
humann_databases --download chocophlan full $DB_ROOT/humann_db --update-config yes
humann_databases --download uniref uniref90_diamond $DB_ROOT/humann_db --update-config yes
humann_databases --download utility_mapping full $DB_ROOT/humann_db --update-config yes

# MetaPHLAn taxonomic profiling database
metaphlan --install --index mpa_vJun23_CHOCOPhlAnSGB_202307 --bowtie2db $DB_ROOT/metaphlan_db
```

> **Note**: These databases are large (~20+ GB total). If you've already downloaded them, skip this step.

### Step 5b: Download 16S Reference Classifier

Download a VSEARCH --sintax compatible FASTA reference:

1. Visit: https://www.drive5.com/usearch/manual/sintax_downloads.html
2. Download `gg_16s_13.5.fa.gz` (GreenGenes 13.5, 99% OTUs)
3. Create the classifier directory and extract:

```bash
mkdir -p $DB_ROOT/16S_classifiers
cp /path/to/gg_16s_13.5.fa.gz $DB_ROOT/16S_classifiers/
gunzip $DB_ROOT/16S_classifiers/gg_16s_13.5.fa.gz
```

### Step 6: Update Environment Variables

Configure Conda environment variables to persist across sessions:

```bash
conda env config vars set YAMAS_HOST_DB=$DB_ROOT/kneaddata_db
conda env config vars set HUMANN_DB_DIR=$DB_ROOT/humann_db
conda env config vars set METAPHLAN_DB_DIR=$DB_ROOT/metaphlan_db
```

### Step 7: Configure HUMAnN Database Paths (If Reusing Existing Databases)

If you've already downloaded HUMAnN databases, ensure they're configured:

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

**You're all set!** Verify installation:

```bash
yamas --version
vsearch --version
mafft --version
fastqc -v
```

---

## Quick Start

### Scenario 1: Download & Process 16S Data (Paired-End, Successful Merge)

```bash
# Download SRR accessions
yamas --download SRR11415443 SRR11415445 --type 16S \
  --acc_list accessions.txt --verbose

# Analyze (paired-end with comma-separated trim/trunc)
yamas --export /path/SRR11415443-... 16S 13,13 150,150 \
  /path/gg_16s_13.5.fa 24
```

### Scenario 2: Long Amplicon, PE Merge Fails → Use --export-as-single

```bash
# Same download step...

# Export with SINGLE-END fallback
yamas --export /path/PRJEB30615-... 16S 15 150 \
  /path/gg_16s_13.5.fa 24 --export-as-single
```

### Scenario 3: Shotgun Metagenomics

```bash
yamas --download SRR11415443 --type Shotgun --verbose

# No --export needed; taxonomy and HUMAnN tables auto-generated
```

---

## Flags & Options

| Flag | Description | Context |
|------|-------------|---------|
| `--download ACCESSIONS...` | Download dataset(s) from SRA/ENA/Qiita | All types |
| `--type {16S,18S,Shotgun}` | Specify sequencing type | All operations |
| `--export PATH DIR TYPE` | Export 16S/18S analysis (OTU, taxonomy, tree) | 16S/18S only |
| `--export-as-single` | Force single-end processing when PE merge fails | 16S/18S with --export |
| `--clean` | Remove host reads with KneadData before analysis | Download & Shotgun |
| `--pathways {yes,no}` | Generate HUMAnN pathway analysis (default: no) | Shotgun only |
| `--threads INT` | Number of parallel threads (default: 8) | All operations |
| `--acc_list FILE` | Path to file with one accession per line | With --download |
| `--verbose` | Print detailed progress information | Download & export |
| `--as_single` | Process as single-end (keeps only forward reads) | Legacy; see `--export-as-single` |
| `--continue_from PATH` | Resume processing from existing project folder | Download only |
| `--continue_from_fastq PATH` | Resume from existing FASTQ folder | Download only |
| `--config FILE` | Load configuration from JSON file | Optional |
| `--ready OS_TYPE` | Initialize YaMAS environment (Ubuntu/CentOS) | Setup only |

---

## Commands Reference

### 1. Data Download

**Download 16S/18S reads from SRA:**
```bash
yamas --download ACCESSION1 ACCESSION2 --type 16S \
  --acc_list accessions.txt \
  --verbose \
  --threads 12
```

| Argument | Required? | Description |
|----------|-----------|-------------|
| `--download ACCESSIONS...` | Yes | Space-separated SRA accessions |
| `--type {16S,18S,Shotgun}` | Yes | Data type |
| `--acc_list FILE` | No | File with accessions (one per line) |
| `--verbose` | No | Print progress messages |
| `--threads INT` | No | Parallel threads (default: 8) |
| `--clean` | No | Remove human-associated reads with KneadData |
| `--pathways {yes,no}` | No | Generate HUMAnN pathway analysis (Shotgun only) |

**Output:** A timestamped directory containing:
- `fastq/` – demultiplexed FASTQ files
- `qza/` – intermediate files
- `vis/` – QC reports
- `export/` – final OTU tables, taxonomy, tree (16S/18S after --export phase)

---

### 2. Export & Analysis

**Paired-End (reads successfully merge):**
```bash
yamas --export /path/dataset-... 16S \
  13,13 \
  150,150 \
  /path/gg_16s_13.5.fa \
  24
```

**Single-End (or PE fallback with --export-as-single flag):**
```bash
yamas --export /path/dataset-... 16S \
  15 \
  150 \
  /path/gg_16s_13.5.fa \
  24 \
  --export-as-single
```

| Argument | Format | Description |
|----------|--------|-------------|
| `origin_dir_path` | string | Directory with reads_data.pkl |
| `data_type` | 16S\|18S\|Shotgun | Analysis type |
| `trim` | **"13,13"** (PE) or **"13"** (SE) | Bases trimmed from 5' end |
| `trunc` | **"150,150"** (PE) or **"150"** (SE) | Read length to truncate to |
| `classifier_file` | path | FASTA reference DB (VSEARCH --sintax format) |
| `threads` | int | CPU threads (default: 12) |
| `--export-as-single` | flag | Force single-end processing (see below) |

#### ⚠️ **The --export-as-single Flag: When & Why to Use It**

##### Background

By default, YaMAS attempts to **merge Paired-End (PE) reads**. This is optimal when:
- Forward and Reverse reads overlap significantly
- Both reads are high quality
- The targeted amplicon is short enough to allow physical overlap

##### When PE Merging Fails

Two common scenarios cause PE merging to fail:

**Scenario 1: Long Amplicon + Short Read Length**
- Example: V3-V4 16S region (~420 bp) with 150 bp PE reads
- Result: VSEARCH mergepairs drops **99.9% of reads** with error "too few kmers found on same diagonal"
- Biological consequence: Dataset becomes unusable

**Scenario 2: Poor Reverse Read Quality**
- Example: Reverse reads have high error rates from sequencing issues
- Result: PE merging succeeds but produces low-quality consensus from bad R2 sequences
- Biological consequence: OTU table inflated with sequence artifacts

##### The Solution: --export-as-single Flag

When PE merging fails or produces poor results, use:

```bash
yamas --export /path/dataset-... 16S 15 150 /path/classifier.fa 24 --export-as-single
```

**What this flag does:**

1. **Ignores Reverse reads completely** – R2 FASTQ files are not processed
2. **Processes ONLY Forward reads** – F reads (_1.fastq) are retained
3. **Prevents data duplication** – Reverse-complement sequences won't inflate the OTU table as "ghost" ASVs

##### Why This Matters Bioinformatically

Without selective filtering, if you processed both `_1.fastq` and `_2.fastq` independently as single-end reads:

```
_1.fastq (forward reads)     → Treated as novel sequences → OTU set A
_2.fastq (reverse reads)     → Reverse-complements of originals → OTU set B (duplicates)
                                                               Total: ~2x feature inflation
```

By using `--export-as-single`, the pipeline:
- Preserves only the canonical forward orientation
- Eliminates duplicate OTUs from reverse complements
- Maintains biological integrity of the OTU table

##### Usage Examples

**Standard paired-end (default behavior):**
```bash
yamas --export /path/SRR11415443-... 16S 13,13 150,150 /path/classifier.fa 24
# Uses comma-separated parameters for forward and reverse trim/trunc
```

**Long amplicon fallback:**
```bash
yamas --export /path/PRJEB30615-... 16S 15 150 /path/classifier.fa 24 --export-as-single
# Uses single integers; ignores R2 files; processes F1 only
```

**Native single-end data (rare):**
```bash
yamas --export /path/SRR041654-... 16S 13 150 /path/classifier.fa 24
# Uses single integers; no flag needed if data has no _2 files
```

##### Decision Tree

```
Do you have both _1 and _2 FASTQ files?
├─ YES, and merge succeeds? → Use comma-separated params (default PE mode)
├─ YES, but merge fails (99% dropout)? → Add --export-as-single flag
├─ YES, but R2 quality terrible? → Add --export-as-single flag
└─ NO (native SE data)? → Use single integers (no flag needed)
```

---

### 3. The --export-as-single Flag (CRITICAL)

**When to use this flag:**

```bash
yamas --export /path/dataset-... 16S 15 150 \
  /path/gg_16s_13.5.fa 24 --export-as-single
```

#### Why You Might Need This

##### Scenario A: Long Amplicon + Paired-End Mismatch
- **Situation**: Your targeted 16S region (e.g., V3-V4) is naturally longer than the combined read length allows for overlap
- **Example**: 150 bp forward + 150 bp reverse = 300 bp total, but V3-V4 is ~420 bp
- **Result**: VSEARCH mergepairs drops **99.9% of reads** with error "too few kmers found on same diagonal"
- **Solution**: Use `--export-as-single` to process only forward reads

##### Scenario B: Poor Reverse Read Quality
- **Situation**: Reverse reads have extremely high error rates due to sequencing issues
- **Result**: PE merging succeeds but produces garbage sequences; OTU table inflated with artifacts
- **Solution**: Use `--export-as-single` to discard low-quality reverse reads

#### How It Works Bioinformatically

**Without flag (standard paired-end):**
```
Forward reads (_1.fastq)  ─┐
                           ├─→ VSEARCH mergepairs ─→ Merged reads
Reverse reads (_2.fastq)  ─┘
                               (99.9% dropout if no overlap)
```

**With --export-as-single flag:**
```
Forward reads (_1.fastq) ──→ Process as single-end only ✓
Reverse reads (_2.fastq) ──→ IGNORED (data discarded safely)
```

#### Critical Biological Point

When you use `--export-as-single`, the pipeline:
1. **Ignores reverse reads completely** – no attempt to merge or process them
2. **Processes only forward (_1) reads** – maintains single-strand orientation
3. **Prevents ghost ASVs** – reverse-complement sequences would inflate the OTU table if processed independently

**Without this selective filtering**, processing both `_1` and `_2` as separate single-end reads would create:
- Duplicate OTUs from reverse-complement sequences
- ~2x inflation of your feature table
- Biological interpretation errors

#### Command Syntax

**Paired-end (default, requires comma-separated params):**
```bash
yamas --export DIR 16S 13,13 150,150 classifier.fa 24
```

**Single-end fallback (requires flag + single integer params):**
```bash
yamas --export DIR 16S 13 150 classifier.fa 24 --export-as-single
```

#### Example Workflow

```bash
# Step 1: Download PRJEB30615 (V3-V4, long amplicon)
yamas --download PRJEB30615 --type 16S --acc_list acc.txt

# Step 2: First attempt with paired-end
yamas --export /path/PRJEB30615-... 16S 15,15 150,150 \
  /data/gg_16s_13.5.fa 24

# >>> ERROR: Merging fails, 99.9% dropout detected <<<

# Step 3: Fallback to single-end
yamas --export /path/PRJEB30615-... 16S 15 150 \
  /data/gg_16s_13.5.fa 24 --export-as-single

# >>> SUCCESS: Process only forward reads, salvage dataset <<<
```

---

## Output Structure

### Project Directory Hierarchy

After downloading a dataset, YaMAS creates a timestamped project directory:

```
<ACCESSION>-<DATE>_<TIME>/
│
├── sra/                      # Raw SRA files (from prefetch)
├── fastq_raw/                # Original FASTQ files (demultiplexed, not processed)
├── fastq_clean/              # Cleaned FASTQ (after KneadData, if --clean used)
├── fastq/                    # Active FASTQ folder (symlink: points to raw or clean)
│
├── knead_out/                # KneadData logs & QC (if --clean used)
│   ├── *_kneaddata.log
│   ├── *.fastqc.html
│   └── multiqc_report.html
│
├── qza/                      # Intermediate files for 16S/18S (can be deleted after export)
│   ├── trimmed/              # Quality-trimmed FASTQ files
│   ├── derep.fasta           # Deduplicated sequences
│   ├── denoised.fasta        # Denoised ASVs (UNOISE3 output)
│   └── rep-seqs-dn-99.fasta  # 99% OTU representative sequences
│
├── vis/                      # QC visualization reports
│   └── multiqc_report.html   # FastQC summary (read counts, qualities, lengths)
│
├── export/                   # Final analysis outputs (16S/18S only, created with --export)
│   ├── otu.tsv / otu.csv     # OTU abundance table (rows=OTUs, cols=samples)
│   ├── taxonomy.tsv / taxonomy.csv  # Taxonomy assignments with confidence
│   ├── tree.nwk              # Phylogenetic tree (Newick, midpoint-rooted)
│   └── otu_padding.csv       # OTU table padded with tree-only tips
│
├── metaphlan_results/        # Taxonomic profiles (Shotgun only)
│   ├── *_profile.tsv         # MetaPHLAn taxonomic abundance
│   └── merged_abundance.tsv  # Merged across all samples
│
└── humann_results/           # Functional profiles (Shotgun only, if --pathways yes)
    ├── *_pathabundance.tsv   # Pathway abundance per sample
    ├── *_pathcoverage.tsv    # Pathway presence/absence
    ├── *_pathabundance_stratified.tsv  # Abundance stratified by species
    └── merged_pathways.tsv   # Merged pathway results
```

### File Lifecycle

| File/Folder | Created | Purpose | Persistent After Processing |
|------------|---------|---------|-----|
| `sra/` | Download | Raw SRA downloads | ✓ (can be deleted after fastq_raw/ created) |
| `fastq_raw/` | Download | Original demultiplexed reads | ❌ Deleted after `--clean` (if used) |
| `fastq_clean/` | Download + `--clean` | Host-removed reads | ❌ Deleted after analysis (only fastq/ symlink kept) |
| `fastq/` | Download | Active processing source | ✓ Symlink only (data elsewhere) |
| `knead_out/` | Download + `--clean` | QC logs | ✓ (optional reference) |
| `qza/` | Download | Processing intermediates | ❌ Can be safely deleted after export |
| `vis/` | Download | QC HTML reports | ✓ (for review) |
| `export/` | `--export` (16S/18S) | Final OTU/taxonomy/tree | ✓ Primary output |
| `metaphlan_results/` | Download (Shotgun) | Taxonomy tables | ✓ Primary output |
| `humann_results/` | Download + `--pathways yes` | Pathway analysis | ✓ Primary output |

### Key Output Files by Analysis Type

**16S/18S Analysis (after `--export`):**
- `export/otu.csv` – OTU abundance (feature × sample matrix)
- `export/taxonomy.csv` – Taxon assignments with confidence scores
- `export/tree.nwk` – Rooted phylogenetic tree
- `vis/multiqc_report.html` – QC summary

**Shotgun Analysis (automatic):**
- `metaphlan_results/merged_abundance.tsv` – Taxonomic composition
- `humann_results/merged_pathways.tsv` – Functional pathways (if `--pathways yes`)
- `vis/multiqc_report.html` – QC summary
- `knead_out/multiqc_report.html` – Host-removal stats (if `--clean`)

---

## Troubleshooting

### Error: "missing 5 arguments"

**Cause**: Paired-end data, but single integers passed without `--export-as-single` flag.

**Fix**:
```bash
# WRONG (causes error):
yamas --export DIR 16S 15 150 classifier.fa 24

# RIGHT (paired-end):
yamas --export DIR 16S 15,15 150,150 classifier.fa 24

# RIGHT (single-end fallback):
yamas --export DIR 16S 15 150 classifier.fa 24 --export-as-single
```

---

### Error: "Reference database not found at..."

**Cause**: Classifier file doesn't exist or path is wrong.

**Fix**: Verify FASTA reference database exists:
```bash
ls -lh /path/gg_16s_13.5.fa
```

Ensure it's formatted for VSEARCH --sintax (FASTA with proper headers).

---

### Error: "too few kmers found on same diagonal" (99.9% dropout during merge)

**Cause**: Amplicon too long; PE reads don't overlap.

**Fix**: Use the `--export-as-single` flag:
```bash
yamas --export DIR 16S 15 150 classifier.fa 24 --export-as-single
```

---

### Very few OTUs / Empty taxonomy.csv

**Cause**: Over-aggressive filtering (min_samples=3 removed all OTUs in single-sample datasets).

**Status**: FIXED in current version. Filter is now capped to the number of actual samples.

---

### MultiQC report shows very low read counts after trimming

**Cause**: Trim/trunc parameters too strict.

**Solution**: Adjust parameters (especially `trunc` / read length cutoff):
```bash
# More permissive:
yamas --export DIR 16S 13 175 classifier.fa 24 --export-as-single
```

---

## Pipeline Architecture

### Stage 0: Download & SRA Conversion
- `prefetch` downloads `.sra` files
- `fasterq-dump --split-files` converts to paired FASTQ

### Stage 1: Demultiplexing & QC
- **EMP format**: Pure Python streaming (barcode + sequence files in lockstep)
- **Qiita format**: Inline barcode extraction + cutadapt trim
- **Manifest format**: Direct FASTQ import
- **QC**: FastQC + MultiQC for HTML reports

### Stage 2: Trimming
- `cutadapt` removes `trim` bases from 5' end, truncates to `trunc` bp

### Stage 3: Denoising & Clustering
- `vsearch --fastq_mergepairs` (PE only, skipped with `--export-as-single`)
- `vsearch --fastq_filter` (quality filtering)
- `vsearch --derep_fulllength` (dereplication)
- `vsearch --cluster_unoise` (UNOISE3 denoising; ASV detection)
- `vsearch --uchime3_denovo` (de novo chimera removal)
- `vsearch --cluster_size` (99% clustering into OTUs)

### Stage 4: Taxonomy
- `vsearch --sintax` with reference FASTA (0.8 confidence cutoff)

### Stage 5: Filtering
- Remove mitochondria/chloroplast by taxonomic name
- Filter by prevalence (min 3 samples or fewer if SE)
- Filter by abundance (min 10 total counts)

### Stage 6: Export
- TSV→CSV conversion for OTU and taxonomy tables

### Stage 7: Phylogeny
- `mafft --auto` (multiple sequence alignment)
- `FastTree -gtr -nt` (GTR model phylogenetic tree)
- BioPython `root_at_midpoint()` (midpoint rooting)
- Padding to include tree tips not in OTU table

---

## Algorithmic Equivalence to QIIME2 Pipeline

| Step | Old (QIIME2) | New (Standalone) | Equivalence |
|------|-------------|-----------------|------------|
| Denoising | DADA2 (error model) | UNOISE3 (fixed threshold) | ✓ High concordance |
| Taxonomy | sklearn Naive Bayes | VSEARCH SINTAX | ✓ Comparable accuracy |
| Chimera | DADA2 consensus | UCHIME3 de novo | ✓ Equivalent |
| Alignment | MAFFT via QIIME2 | Direct MAFFT | ✓ Identical |
| Tree | FastTree via QIIME2 | Direct FastTree | ✓ Identical |
| Rooting | skbio midpoint-root | BioPython midpoint-root | ✓ Identical |

**Key Difference**: DADA2 learns per-sample error models; UNOISE3 uses fixed alpha threshold (2.0). For most well-sequenced datasets, results are indistinguishable.

---

## Citation

If you use YaMAS in your research, please cite:

```
YaMAS: YOLO Microbiome Analysis System (Updated)
https://github.com/[your-repo]
VSEARCH: Rognes T, et al. VSEARCH: a versatile open source tool for metagenomics. PeerJ 2016
SINTAX: Edgar R. SINTAX, a simple non-Bayesian taxonomy classifier for 16S and ITS sequences
FastTree: Price MN, et al. FastTree 2. PLoS ONE 2010
```

---

## License

MIT License – see LICENSE file for details.

---

## Support & Issues

For bugs, feature requests, or questions:
- Check the [Troubleshooting](#troubleshooting) section
- Review the [Algorithmic Equivalence](#algorithmic-equivalence-to-qiime2-pipeline) table
- Create an issue on GitHub with:
  - Error message (full traceback)
  - Command used
  - Output directory structure
  - Log file (if available)
