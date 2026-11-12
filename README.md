# LncKit

LncKit is an automated, configuration-driven pipeline for identifying and classifying long non-coding RNAs (lncRNAs) from RNA-seq data. It accepts sequencing data from NCBI SRA or local FASTQ files and uses a reference genome and annotation to produce candidate lncRNAs, intermediate results, and a summary report.

The pipeline runs from the command line and can be resumed if an analysis is interrupted. Intermediate files are retained so that users can inspect the results of individual stages.

## Workflow overview

The default workflow includes sequencing-data retrieval, quality control, trimming, optional replicate merging, reference indexing, read alignment, transcript assembly, comparison with the reference annotation, candidate filtering, expression filtering, ORF filtering, coding-potential assessment, protein-homology screening, known ncRNA filtering, classification, and report generation.

## Installation

Unpack the release archive and move into the directory:

```bash
tar xJf LncKit.tar.xz
cd LncKit
```

Run the setup script from the unpacked LncKit directory:

```bash
bash setup.sh
```

The setup script creates a self-contained micromamba/Conda installation and the environments used by LncKit: `lnckit`, `lnckit_feelnc`, `lnckit_lncfinder`, and `lnckit_py2`.

Download the reference databases required for your analysis:

```bash
./download_databases.sh
```

### Database download options

| Option | Database content | Approximate size |
|:---|:---|---:|
| `--uniprot` | UniProt Swiss-Prot protein database | 0.3 GB |
| `--rfam` | Rfam covariance models for Infernal | 1.5 GB |
| `--cog` | COG protein FASTA database | 1.6 GB |
| `--pfam` or `--pfam-fasta` | Pfam-A protein FASTA database | 20 GB |
| `--homology` | UniProt, Pfam, and COG | 22 GB |
| `--all` | UniProt, Rfam, Pfam, and COG | 24 GB |
| `--cpat=human` | One CPAT model; supported examples include human, mouse, fly, and zebrafish | Small |
| `--cpat` | All available animal CPAT models | Small |

Database sizes are approximate and may change. The tools and databases are kept within the LncKit installation or configured database directory rather than being installed system-wide.

## Prepare the input files

Create a project folder for the input files and results. The project sheet is usually named `samples.txt`. Edit the supplied template and provide the paths to the reference genome, annotation, and samples.

### Required inputs

| Input | Description |
|:---|:---|
| `samples.txt` | Project sheet containing reference paths and sample or group information |
| Reference genome | Genome sequence in FASTA format |
| Reference annotation | Gene annotation in GTF format; GFF may also be supplied where required |
| FASTQ files or SRA accessions | Sequencing reads for each sample |

### Example `samples.txt` for NCBI SRA input

```text
genome_fasta,/path/to/YourGenome.fasta
annotation_gtf,/path/to/YourAnnotation.gtf
# or supply a GFF instead of the GTF above (LncKit converts it automatically):
# annotation_gff,/path/to/YourAnnotation.gff
accession,group
SRRXXXXXXXX,disease
SRRXXXXXXXX,disease
SRRXXXXXXXX,normal
SRRXXXXXXXX,normal
```

Give the paths to your reference genome (`genome_fasta`) and annotation (`annotation_gtf`), then one `SRR/ERR/DRR` run accession per sample row under the `accession,group` header. Replace the example paths, accessions, and group labels with your own values.

No reference genome is bundled - you supply your own, for whatever organism you study. The file names are arbitrary: LncKit uses whatever path you give and recognises each file's role from the row key (`genome_fasta`, `annotation_gtf`), not from the filename. You may give a **URL** instead of a local path (`http`/`https`/`ftp`) and it is downloaded automatically (`.gz` is unpacked). You supply just ONE annotation: either a GTF (`annotation_gtf`) or a GFF (`annotation_gff`). If you give a GFF, LncKit converts it to GTF automatically and uses it everywhere the annotation is needed, so there is no need to supply both.

### Input options

For SRA input, provide the SRR, ERR, or DRR accession in the sample rows.

For local FASTQ input, add a `reads_dir,<folder>` row giving the path to your FASTQ files, and list sample names that match the FASTQ filename prefixes. LncKit attempts to detect paired-end or single-end reads from the filenames. You can also provide a local FASTQ folder directly through the `--input` option.

### Example `samples.txt` for local FASTQ files

```text
genome_fasta,/path/to/YourGenome.fasta
annotation_gtf,/path/to/YourAnnotation.gtf
# or supply a GFF instead of the GTF above (LncKit converts it automatically):
# annotation_gff,/path/to/YourAnnotation.gff
reads_dir,/path/to/fastq-folder
accession,group
sample1,disease
sample2,disease
sample3,normal
sample4,normal
```

Give the same reference paths as above (any filename, or a URL), a `reads_dir` row with the path to the folder holding your FASTQ files, and one row per sample under `accession,group`. Each sample name must match its FASTQ prefix (`sample1_1.fastq.gz` / `sample1_2.fastq.gz` for paired-end, or `sample1.fastq.gz` for single-end). Both compressed and uncompressed reads are accepted - `.fastq.gz`, `.fastq`, `.fq.gz` and `.fq` all work.

Example:

```bash
./lnckit.sh /path/to/run-folder --input /path/to/fastq-folder
```

## Run LncKit

LncKit is controlled through command-line options. If no options are supplied for a step, the default configuration is used.

### Run the default workflow

```bash
./lnckit.sh /path/to/run-folder
```

### Change selected tools or parameters

```bash
./lnckit.sh /path/to/run-folder --aligner star
./lnckit.sh /path/to/run-folder --organism 1 --coding 1,2,5 --strandness RF
```

The first command selects STAR as the aligner. The second selects the plant organism profile, combines the specified coding-potential tools, and sets the library strandness explicitly to RF (overriding auto-detection).

### Preview the workflow

Use `--dry-run` to display the resolved execution plan without running the pipeline:

```bash
./lnckit.sh /path/to/run-folder --dry-run
```

### Run on a compute node

```bash
nohup ./lnckit.sh /data/project --organism 1 --threads 16 > lnckit.log 2>&1 &
```

Example SLURM submission:

```bash
sbatch --cpus-per-task=16 --wrap \
  "/path/to/LncKit/lnckit.sh /path/to/run-folder --organism 1 --threads 16"
```

Example PBS submission:

```bash
echo "/path/to/LncKit/lnckit.sh /path/to/run-folder --organism 1 --threads 16" \
  | qsub -l nodes=1:ppn=16 -l walltime=48:00:00 -N lnckit
```

## Containers (Docker and Singularity)

LncKit can be run from a container image with every tool pre-installed, so no separate `setup.sh` or Conda installation is needed on the host. The reference databases and your project data are mounted at run time; they are not stored inside the image.

Usage is identical to the normal command line: whatever you pass after the image name is handed to `lnckit.sh`, so every option in this document works unchanged.

First, on the host, download the reference databases once (this needs internet access, so run it on a login node). The directory you choose is mounted into the container at `/db`:

```bash
./download_databases.sh                 # into the databases/ folder, by default
```

### Docker

```bash
# Build the image once, from the unpacked LncKit directory
docker build -t lnckit:latest .

# See all options
docker run --rm lnckit:latest --help

# Run: mount the project folder at /work and the databases at /db
docker run --rm -v /path/to/project:/work -v /path/to/databases:/db \
  lnckit:latest /work --organism plant --threads 16
```

### Singularity / Apptainer (HPC)

```bash
# Build the image (no root needed with --fakeroot)
apptainer build --fakeroot lnckit.sif Singularity.def

# See all options
apptainer run lnckit.sif --help

# Run: bind the project folder to /work and the databases to /db
apptainer run -B /path/to/project:/work -B /path/to/databases:/db \
  lnckit.sif /work --organism plant --threads 16
```

Containers run without an interactive terminal, so LncKit uses its non-interactive mode automatically, which makes the images suitable for batch schedulers. Example PBS job using the image:

```bash
#!/bin/bash
#PBS -N lnckit
#PBS -l nodes=1:ppn=16
#PBS -l walltime=48:00:00
#PBS -j oe

cd "$PBS_O_WORKDIR"
apptainer run -B "$PBS_O_WORKDIR":/work -B /path/to/databases:/db \
  /path/to/lnckit.sif /work --organism plant --threads 16
```

Build the image and download the databases on a node with internet access (usually the login node). If the compute node has no internet, fetch the SRA reads beforehand or use local FASTQ files with the `--input` option, because the sequencing-data retrieval step needs network access.

## Command-line options

Unless stated otherwise, each option is optional and the default setting is used when it is not supplied. Tool selections can generally be provided by name or by number.

### Input and read layout

| Option | Description |
|:---|:---|
| `--input <path>` | Path to `samples.txt` or a local FASTQ directory |
| `--paired` | Force paired-end read handling |
| `--single` | Force single-end read handling |
| `--merge` | Pool replicates by group before alignment; enabled by default |
| `--no-merge` | Disable replicate merging |

### Tool selection

| Option | Available values |
|:---|:---|
| `--organism` | `1=plant`, `2=human`, `3=mouse`, `4=fly`, `5=zebrafish`, `6=other` |
| `--strandness` | `1=auto`, `2=RF`, `3=FR`, `4=none` (auto is the default) |
| `--qc` | `1=fastqc`, `2=fastp`, `3=multiqc` |
| `--trim` | `1=fastp`, `2=trimmomatic`, `3=cutadapt` |
| `--aligner` | `1=hisat2`, `2=star`, `3=bowtie2`, `4=bwa` |
| `--assembler` | `1=stringtie`, `2=cufflinks` |
| `--coding` | `1=cpat`, `2=cpc2`, `3=cnci`, `4=plek`, `5=lncfinder`, `6=rnaplonc` |
| `--consensus` | `all`, `majority`, or a minimum number of supporting tools |
| `--protein` | `1=diamond`, `2=blastx` |
| `--housekeeping` | `1=infernal`, `2=blastn`, `3=none` |
| `--classifier` | `1=feelnc`, `2=gffcompare` |

Multiple coding-potential tools can be selected as a comma-separated list. For example:

```bash
./lnckit.sh /path/to/run-folder --coding 1,2,5
```

The `--consensus` option controls how the coding-potential results are combined:

| Setting | Meaning |
|:---|:---|
| `all` | Retain a candidate only when every selected tool predicts it to be non-coding |
| `majority` | Retain a candidate when more than half of the selected tools predict it to be non-coding |
| `N` | Retain a candidate when at least N tools predict it to be non-coding |

For example, `--consensus 3` requires at least three selected tools to classify a candidate as non-coding.

### Run control

| Option | Description |
|:---|:---|
| `-c, --config <file.yaml>` | Use a custom configuration file instead of `pipeline.yaml` |
| `-o, --output <dir>` | Set the output directory |
| `-p, --threads <N>` | Set the number of CPU threads |
| `--from <step>` | Start execution from a specified step |
| `--until <step>` | Stop execution after a specified step |
| `--dry-run` | Print the resolved execution plan without running the pipeline |
| `--help` | Display available command-line options |

Short forms include `-i` for `--input`, `-org` for `--organism`, `-str` for `--strandness`, `-a` for `--aligner`, `-o` for `--output`, `-c` for `--config`, and `-p` for `--threads`.

## Default settings

| Workflow setting | Default |
|:---|:---|
| Organism | `other` (works for any organism; selects species-specific coding-potential models when set to plant/human/mouse/fly/zebrafish) |
| Quality control | `fastp` |
| Trimming | `fastp` |
| Aligner | `hisat2` |
| Assembler | `stringtie` |
| Strandness | `auto` (detects RF/FR/unstranded from the aligned reads; falls back to unstranded if the signal is weak or ambiguous) |
| Protein-homology tool | `diamond` |
| Known ncRNA filtering | `infernal` |
| Classifier | `feelnc` |
| Candidate filtering | FEELnc; minimum transcript length of 200 nt, exon-structure rules, and mRNA-overlap filtering |
| Coding-potential assessment | CPAT, CPC2, CNCI, PLEK, LncFinder, and RNAplonc (RNAplonc is applied for plant organisms only) |
| Coding-potential consensus | Majority |
| Input source | SRA accessions listed in `samples.txt` |

## Library strandedness

By default (`strandness: auto`), LncKit infers the library type from the aligned reads before assembly: it compares the strand that the splice-aware aligner assigns to each spliced read (the `XS` tag) against the read's own mapped strand, the way RSeQC's `infer_experiment` does. It commits to `RF` (dUTP first-strand, the most common) or `FR` (second-strand) only when the evidence is overwhelming and based on enough reads; otherwise it falls back to unstranded, so it never mis-assigns strand. The decision is written to the run log.

This matters for lncRNA discovery specifically: running genuinely stranded data as unstranded silently loses antisense and mono-exonic lncRNAs (StringTie cannot strand them, and they are then dropped). If you already know your library type, set it explicitly with `--strandness RF` (or `FR`, or `none`), or in `pipeline.yaml`; an explicit value always skips detection.

## Choosing tools

HISAT2 is the default aligner and is splice-aware, making it suitable for RNA-seq. STAR is also splice-aware and may require substantially more memory. Bowtie2 and BWA are not splice-aware and are intended only when reproducing a protocol that specifically requires them; they are generally not recommended for transcript assembly in lncRNA discovery.

`fastp` performs quality filtering and trimming and produces a quality report. FastQC generates quality-control reports, while MultiQC combines reports from multiple samples. Trimmomatic and Cutadapt are available as alternative trimming tools.

The default coding-potential analysis uses five tools (CPAT, CPC2, CNCI, PLEK, and LncFinder), with RNAplonc added for plant organisms, combined under a majority-consensus rule. The `--consensus` option can be changed to make the selection stricter or more permissive.

## Protein-homology databases

UniProt, COG, and Pfam are used by the protein-homology screening stage; they are not coding-potential predictors. Download the databases you intend to use with `download_databases.sh`. Databases found in the configured reference directory may be detected automatically. You can also specify database paths in the `reference.protein_db` section of `pipeline.yaml`.

## Output files

By default, LncKit writes results to `lnckit_output/`. Intermediate files are retained to support inspection and troubleshooting.

| Directory | Contents |
|:---|:---|
| `raw_reads/` | Downloaded or supplied FASTQ files |
| `trimmed_reads/` | Adapter- and quality-trimmed reads |
| `merged_reads/` | Replicates pooled by group, when merging is enabled |
| `qc/` | Quality-control reports |
| `alignment/` | Alignment files, including sorted BAM files |
| `assembly/` | Per-sample and merged transcript assemblies |
| `coding_prediction/` | Coding-potential prediction results |
| Other filtering directories | Candidate-filtering, protein-homology, known ncRNA, and classification results |
| `final/` | Final lncRNA outputs |
| `reports/` | Summary reports |
| `logs/` | Execution logs |

The final output includes the candidate lncRNA set and a summary report. Exact files and directories may depend on the selected tools and configuration.

## Configuration

Use `pipeline.yaml` to change default tools, thresholds, expression settings, ORF settings, coding-potential consensus rules, and other supported parameters.

Configuration precedence is:

1. Packaged defaults
2. User YAML configuration
3. Command-line options

Command-line options override YAML settings, and YAML settings override packaged defaults. For reference paths (genome and annotation), the project sheet (`samples.txt`) takes precedence over `pipeline.yaml`. There is no filename-based auto-detection: you always give an explicit path or URL.

## Resuming an analysis

LncKit retains intermediate files and supports resuming an interrupted run. Re-run the command with the same project folder and configuration to continue from the available checkpoint. Review the execution plan and logs if the previous run stopped unexpectedly.

## Notes

Use a splice-aware aligner, such as HISAT2 or STAR, for standard RNA-seq transcript assembly. Bowtie2 and BWA are included as selectable options for compatibility with specific workflows, but they are not splice-aware.

Review the execution plan before running a large analysis, particularly when downloading large reference databases or using STAR with limited memory.

For the complete option list, run:

```bash
./lnckit.sh --help
```
