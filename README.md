# 🧬 piRAT – piRNA Annotation Tool

**piRAT** (v0.1.0) is a Python-based tool for annotation, visualization, and analysis of **PIWI-interacting RNAs (piRNAs)** from `.bam` files. Designed for ease-of-use and scientific rigor, piRAT enables researchers to detect **primary piRNA clusters** and **ping-pong amplification signatures** directly from raw NGS data.

> ✅ No pre-processing required &nbsp; | &nbsp; 📊 HTML reports &nbsp; | &nbsp; 🚀 Auto-tuned clustering &nbsp; | &nbsp; 🧵 Multi-threaded
---
## 🔧 Features

### 📂 Seamless Input Handling
- Direct analysis of `.bam` files — no pre-processing required
- Automatic `.bai` index generation if missing
- Multi-file support for batch analysis

### 🧬 Smart piRNA Annotation
- Auto-detection of piRNA size range (based on read length distribution)
- Adaptive DBSCAN clustering optimized for piRNA data
- Annotation of:
  - **Primary piRNA clusters**
  - **Ping-pong amplification signatures**

### 📊 Rich, Interactive Outputs
- Fully featured HTML reports (one per sample) including:
  - Cluster statistics
  - Read length histograms
  - 5′–5′ overlap distributions
  - SeqLogo visualizations
  - Heatmaps of ping-pong interactions
  - Venn diagrams of pathway overlaps
- Export of structured output files:
  - `.gff` cluster and ping-pong annotations
  - SeqLogo images
  - Heatmaps and plots

### ⚡ Performance & Usability
- Multi-threaded execution for large-scale datasets
- Auto-tuned clustering parameters (`eps`, `k`)
- Interactive or fully automated mode (`-a` flag)
- Integrated troubleshooting and error messages

---
## 🚀 Installation

Install from PyPI:

```bash
pip install pirat
```
Or from the source
```bash
git clone https://github.com/ylla-lab/piRAT.git
cd piRAT
pip install .
```

---
## ⚡ Quick Start

Analysis of `.bam` files in a folder:

```bash
pirat -p <path_directory_with_bam_files>
```
Full automatic run with 16 threads

```bash
pirat -p <path_directory_with_bam_files> -a -t 16
```
---
## 🧪 Example Use Cases
### Cluster detection only
```bash
pirat -p <path_directory_with_bam_files> -m primary
```
### Cluster detection + cluster plots generation
```bash
pirat -p <path_directory_with_bam_files> -m primary -d
```
### Ping-pong signature detection only
```bash
pirat -p <path_directory_with_bam_files> -m secondary
```
### Full pipeline: clusters + ping-pong signatures
```bash
pirat -p <path_directory_with_bam_files> -m both
```
### Full run with custom parameters
```bash
pirat -p <input_path> -o <output_path> -t 16 -a -d -r 28,30 -k 10 -e 800
```
---
## 🛠️ Command-line Options

Below are the available flags for running piRAT from the command line:

### 📂 Input/Output
- `-p`, `--path`  
  **Required.** Path to the input directory containing `.bam` files.

- `-o`, `--output-path`  
  Path to the directory for output files. If not specified, defaults to `./<name_of_the_input_directory>`.

### ⚙️ Analysis Control
- `-m`, `--module`  
  Which module(s) to run:
  - `primary` – Detect primary piRNA clusters
  - `secondary` – Detect ping-pong signatures
  - `both` *(default)* – Run full analysis

- `-r`, `--range_of_size`  
  Set custom piRNA size range. Format: `min,max` (e.g., `26,32`)

- `-k`, `--minreads`  
  Minimum number of reads required to form a cluster seed. *(Integer)*

- `-e`, `--eps`  
  Maximum distance between reads in a cluster. *(Integer)*

- `-v`, `--variation_threshold`  
  Threshold for read position variation filtering. *(Default: 3)*

### 🧵 Performance
- `-t`, `--threads`  
  Number of threads to use. *(Integer)*

- `--plot_iter`  
  Number of plots generated per iteration (affects RAM usage). *(Default: 16)*

### 🧪 Execution Mode
- `-a`  
  Run in automatic (non-interactive) mode. Useful for pipelines or scripts.

- `-d`  
  Generate cluster plots in addition to annotations.

### ℹ️ Miscellaneous
- `--version`  
  Display version info and exit.

---
## 📄 Output Files
Each run (of the full pipeline) generates:
- Cluster and ping-pong annotation files (.gff format)
- HTML reports for each part of the analysis
- Optional cluster plots
- Comprehensive analysis of each file
- SeqLogo visualizations
- Heatmaps of ping-pong pairs interactions
- Venn diagrams of primary/secondary reads
---
## ❓ Troubleshooting

- **Corrupted BAM file?**  
  Make sure your `.bam` files are BGZF-compressed and indexed. If missing, piRAT will try to index them using `samtools`.

- **Out of memory?**  
  Reduce thread count (`-t`) or number of clusters analyzed during statistical analysis (`--plot_iter`).

- **Reruns?**  
  Avoid modifying the output directory during execution. Re-running on the same folder will overwrite previous results.
---
## 📚 Citation

If you use piRAT in your research, please cite:

> **Dominik Robak, Guillem Ylla (2025)**  
> *piRAT: piRNA Annotation Tool for annotation, analyzing, and visualizing piRNAs*
---
## 📬 Contact

- **Guillem Ylla**  
  Laboratory of Bioinformatics and Genome Biology  
  Jagiellonian University, Kraków, Poland  
  📧 [guillem.ylla@uj.edu.pl](mailto:guillem.ylla@uj.edu.pl)  
  🌐 [Ylla Lab](https://ylla-lab.github.io/)

- **Dominik Robak**  
  📧 [dominikrobak03@gmail.com](mailto:dominikrobak03@gmail.com)
