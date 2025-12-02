<p align="center">

  <img src="https://camo.githubusercontent.com/87d2f7858eb2f36b36097cf257faaa40d3146b708d0afb3ebca17f0f2f99ab11/68747470733a2f2f692e696d6775722e636f6d2f31744a76696f442e706e67" alt="TransRate2" width="1000">

</p>



<p align="right">EC Bretz</p>



**TransRate2** is a tool for quality analysis of de-novo transcriptome assemblies. It evaluates assemblies using a wide range of metrics, leveraging read mapping, reference-based statistics, and contig-level features to provide a detailed assessment of transcriptome assembly quality. TransRate2 supports both single and paired-end reads, multiple aligners, and reference-based evaluation.



**TransRate2** is a reimplementation of the original bioinformatics tool, TransRate (Smith-Unna et al. 2015).

https://github.com/blahah/transrate





> [!CAUTION]
> $\Huge\textcolor[RGB]{248, 82, 73}{\textsf{TransRate2 is currently under development}}$



<h2><img src="https://i.imgur.com/3UA4xwp.png" width="20" align="top">&ensp;Dependencies</h2>



- **Python 3.8+**

- **pandas**, **numpy**, **pysam**, **psutil** (Python libraries)

- **Bowtie2** or **Hisat2** - Read aligners

- **Salmon** - Quantification

- **Samtools** - BAM/SAM processing

- **Diamond** (optional, for reference-based analysis)



#### Ubuntu/Debian:

```bash

sudo apt-get install bowtie2 hisat2 salmon samtools diamond

pip install pandas numpy pysam psutil

```



<h2><img src="https://i.imgur.com/3UA4xwp.png" width="20" align="top">&ensp;Installation</h2>



1. Clone the repository:

```bash

git clone https://github.com/ericbretz/transrate2.git

cd TransRate2

```



2. (Optional) Set up a virtual environment and install Python dependencies:

```bash

python3 -m venv venv

source venv/bin/activate

pip install pandas numpy pysam psutil

```



<h2><img src="https://i.imgur.com/3UA4xwp.png" width="20" align="top">&ensp;Usage</h2>



### Basic Usage



```bash

python transrate2.py --assembly assembly.fasta --left reads_1.fq --right reads_2.fq

```



### Command Line Options




| Option             | Short | Description                                      | Default         |
|--------------------|-------|--------------------------------------------------|-----------------|
| `--assembly`       | `-a`  | Assembly file(s), comma separated                | *Required*      |
| `--left`           | `-l`  | Left reads file (FASTQ)                          |                 |
| `--right`          | `-r`  | Right reads file (FASTQ)                         |                 |
| `--reference`      | `-f`  | Reference file (FASTA, optional)                 |                 |
| `--bam`            | `-x`  | BAM file for pre-aligned reads                   |                 |
| `--output-dir`     | `-o`  | Output directory                                 | Current dir     |
| `--threads`        | `-t`  | Number of threads                                | `4`             |
| `--hisat2`         | `-s`  | Use Hisat2 aligner (default: Bowtie2)            | `False`         |
| `--bowtie2`        | `-b`  | Use Bowtie2 aligner                              | `True`          |
| `--clutter`        | `-c`  | Remove intermediate files                        | `False`         |
| `--quiet`          | `-q`  | Suppress terminal output                         | `False`         |
| `--debug`          |       | Enable debug mode                                | `False`         |
| `--nocolor`        |       | Disable colored terminal output                  | `False`         |
| `--config`         |       | Path to configuration file                       |                 |
| `--config-create`  |       | Create configuration template                    | `config.yaml`   |
| `--config-save`    |       | Save current arguments to config file            | `config.yaml`   |
| `--version`        | `-v`  | Print version                                    |                 |
| `--help`           | `-h`  | Display help message                             |                 |



> **Note:** TransRate2 supports multiple input modes:
> - **Assembly-only mode**: Provide only `--assembly`
> - **Assembly + reads mode**: Provide `--assembly` with `--left` and/or `--right`
> - **BAM + assembly mode**: Provide `--bam` with `--assembly`
> - **BAM-only mode**: Provide only `--bam` with `--left` and/or `--right`
> - **Multi-assembly mode**: Provide multiple assemblies comma-separated to `--assembly`
> - **Reference-based analysis**: Add `--reference` to any mode for enhanced evaluation



<h2><img src="https://i.imgur.com/3UA4xwp.png" width="20" align="top">&ensp;Examples</h2>



#### Basic run with paired-end reads:
```bash
python transrate2.py -a assembly.fasta -l reads_1.fq -r reads_2.fq
```

#### Single-end reads:
```bash
python transrate2.py -a assembly.fasta -l reads.fq
```

#### Assembly-only analysis:
```bash
python transrate2.py -a assembly.fasta
```

#### Multi-assembly processing:
```bash
python transrate2.py -a assembly1.fasta,assembly2.fasta,assembly3.fasta -l reads_1.fq -r reads_2.fq
```

#### Using pre-aligned BAM file:
```bash
python transrate2.py -a assembly.fasta -x aligned_reads.bam
```

#### BAM-only mode (no assembly file):
```bash
python transrate2.py -x aligned_reads.bam -l reads_1.fq -r reads_2.fq
```

#### Specify output directory and use Hisat2:
```bash
python transrate2.py -a assembly.fasta -l reads_1.fq -r reads_2.fq -o results -s
```

#### Reference-based evaluation:
```bash
python transrate2.py -a assembly.fasta -l reads_1.fq -r reads_2.fq -f reference.fasta
```

#### Using configuration files:
```bash
# Create a config template
python transrate2.py --config-create my_config.yaml

# Run with config file
python transrate2.py --config my_config.yaml

# Save current arguments to config
python transrate2.py -a assembly.fasta -l reads_1.fq -r reads_2.fq --config-save my_settings.yaml
```

#### Advanced options:
```bash
# Remove intermediate files, use 8 threads, quiet mode
python transrate2.py -a assembly.fasta -l reads_1.fq -r reads_2.fq -c -t 8 -q

# Debug mode with no color output
python transrate2.py -a assembly.fasta -l reads_1.fq -r reads_2.fq --debug --nocolor
```

#### Show help and version:
```bash
python transrate2.py -h
python transrate2.py -v
```



<h2><img src="https://i.imgur.com/3UA4xwp.png" width="20" align="top">&ensp;Overview</h2>


1. **Input Validation** - Comprehensive validation of input files and formats with support for multiple input modes.

2. **Multi-Assembly Processing** - Process multiple assemblies in a single run for comparative analysis.

3. **Flexible Input Support** - Support for FASTA assemblies, FASTQ reads, and pre-aligned BAM files.

4. **Read Mapping** - Maps reads to assemblies using Bowtie2 or Hisat2 aligners.

5. **Quantification** - Uses Salmon for transcript quantification and expression analysis.

6. **Contig Metrics** - Calculates comprehensive contig-level and assembly-level metrics (length, GC content, N50, coverage, mapping statistics, fragmentation, segmentation, etc.).

7. **Reference-based Evaluation** (optional) - Uses Diamond to compare assemblies to reference sequences and computes reciprocal best hits and coverage statistics.

8. **Scoring & Classification** - Computes overall and optimal assembly scores, and classifies contigs as "good" or "bad" based on multiple quality criteria.

9. **Advanced Logging** - Comprehensive logging system with stage tracking and performance monitoring.

10. **Configuration Management** - YAML-based configuration system for reproducible analyses.

11. **Reporting** - Generates detailed CSV reports, summary statistics, and filtered FASTA files.



<h2><img src="https://i.imgur.com/3UA4xwp.png" width="20" align="top">&ensp;Output</h2>



TransRate2 creates an output directory structure containing:


```
TransRate2/
├── results/                                        # Main analysis results
│   ├── <assembly_name>.contigs.csv                 # Per-contig metrics
│   ├── good.<assembly_name>.fa                     # Good contigs (FASTA)
│   ├── bad.<assembly_name>.fa                      # Bad contigs (FASTA)
│   ├── assembly_score_optimisation.csv             # Score optimization data
│   └── assembly.csv                                # Combined summary for all assemblies
├── logs/                                           # Comprehensive log files with stage tracking
│   ├── transrate2_<assembly_name>_<timestamp>.log  # Main analysis log
│   ├── <aligner>_<assembly_name>_stdout.log        # Aligner output logs
│   ├── <aligner>_<assembly_name>_stderr.log        # Aligner error logs
│   ├── salmon_<assembly_name>_stdout.log           # Salmon output logs
│   ├── salmon_<assembly_name>_stderr.log           # Salmon error logs
│   └── samtools_<assembly_name>_stdout.log         # Samtools logs
└── temp/                                           # Temporary files (removed with --clutter)
    ├── alignments/
    │   └── <aligner>/
    │       ├── index/                              # Aligner index files
    │       └── <assembly_name>_aligned.bam         # Alignment results
    ├── bam/
    │   ├── <assembly_name>_sorted.bam              # Sorted BAM file
    │   └── <assembly_name>_sorted.bam.bai          # BAM index
    ├── salmon/                                     # Salmon quantification files
    │   ├── aux_info/
    │   ├── logs/
    │   ├── quant.sf
    │   └── ...
    ├── analysis/
    │   └── <assembly_name>.assembly.json           # Assembly analysis data
    └── reference/                                  # Reference-based results
        └── reciprocal_hits.csv
```

<h2><img src="https://i.imgur.com/3UA4xwp.png" width="20" align="top">&ensp;Output Files</h2>



#### Main Results
- **`<assembly_name>.contigs.csv`** - Per-contig metrics and quality scores
- **`good.<assembly_name>.fa`** - High-quality contigs (FASTA)
- **`bad.<assembly_name>.fa`** - Low-quality contigs (FASTA)
- **`assembly_score_optimisation.csv`** - Score optimization data
- **`assembly.csv`** - Combined summary for all assemblies

#### Logs
- **`transrate2_<assembly_name>_<timestamp>.log`** - Main analysis log with stage tracking
- **`<aligner>_<assembly_name>_stdout.log`** - Aligner output logs
- **`salmon_<assembly_name>_stdout.log`** - Salmon quantification logs

#### Temporary Files
- **`alignments/<aligner>/`** - Aligner indices and BAM files
- **`salmon/quant.sf`** - Salmon quantification results
- **`bam/<assembly_name>_sorted.bam`** - Sorted alignment file
- **`analysis/<assembly_name>.assembly.json`** - Assembly analysis data
- **`reference/reciprocal_hits.csv`** - Reference-based reciprocal best hits

<h2><img src="https://i.imgur.com/3UA4xwp.png" width="20" align="top">&ensp;Configuration Management</h2>

TransRate2 includes a YAML configuration system for reproducible analyses:

#### Creating Configuration Files
```bash
# Create a default configuration template
python transrate2.py --config-create

# Create with custom name
python transrate2.py --config-create my_analysis.yaml
```

#### Using Configuration Files
```bash
# Run analysis with config file
python transrate2.py --config my_analysis.yaml

# Override config values with command line arguments
python transrate2.py --config my_analysis.yaml -t 8 --debug
```

#### Saving Current Settings
```bash
# Save current command line arguments to config
python transrate2.py -a assembly.fasta -l reads_1.fq -r reads_2.fq -t 8 --config-save my_settings.yaml
```
