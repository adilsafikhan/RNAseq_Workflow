# Transcriptome Analysis Workflow (Maize RNA‑seq: 2024–2026)

This document describes the complete workflow used to process maize RNA‑seq data across **three seasons** (2024, 2025, 2026): environment setup, directory structure, QC, trimming, alignment with STAR, transcript-level quantification with Salmon, gene-level TPM/count aggregation with `tximport`, strandness inference, and read counting with featureCounts.

---

## 1) Conda Environments & Tooling

```bash
# Create primary environment
conda create -n bioinfo -y
conda activate bioinfo

# Core tools
conda install -y trimmomatic fastqc multiqc samtools STAR
conda install -y -c bioconda bedops
conda install -y agat
conda install -y -c bioconda subread   # provides featureCounts

# Separate environment for Salmon (recommended)
conda create -n salmon-env -c conda-forge -c bioconda salmon -y
```

---

## 2) Project Structure

```bash
mkdir -p data feature_count final_results logs mapping qc reference salmon_quant scripts
```

Within `data/`:

```
data/
├── fastq_files_2024/             # raw FASTQ (2024)
├── fastq_files_2025/             # raw FASTQ (2025)
├── fastq_files_2026/             # raw FASTQ (2026; may include *_merged.fastq.gz)
├── trimmed_fastq_file_2024/      # trimmed FASTQ outputs (2024)
├── trimmed_fastq_file_2025/      # trimmed FASTQ outputs (2025)
└── trimmed_fastq_file_2026/      # trimmed FASTQ outputs (2026)
```

**Note:** For 2026, some samples were sequenced on two lanes and merged (e.g., `*_R1_merged.fastq.gz`, `*_R2_merged.fastq.gz`). For 2024, one sample’s R2 is missing; alignment/trimming skipped for that sample pending re-supply. Later I found the sample. 

---

## 3) FASTQ QC (Raw & Trimmed)

If the facility did not provide FastQC:

```bash
# Raw QC (optional if facility already provided reports)
fastqc -o qc data/fastq_files_202*/*.fastq.gz
multiqc qc -o qc/multiqc_raw
```

After trimming (section 5), run:

```bash
# Simple, parallel FastQC on all paired trimmed reads
mkdir -p qc logs
parallel -j 50 fastqc -o qc ::: data/trimmed_fastq_file_202*/*_paired.fq.gz
multiqc qc -o qc/multiqc_trimmed
```

---

## 4) Reference: Genome & Annotation (Ensembl Plants, Zea mays V5)

```bash
# GTF (annotation), primary assembly FASTA (genome), and cDNA (transcripts for Salmon)
wget -P reference https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-62/plants/gtf/zea_mays/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.62.gtf.gz
wget -P reference https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-62/plants/fasta/zea_mays/dna/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa.gz
wget -P reference https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-62/plants/fasta/zea_mays/cdna/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.cdna.all.fa.gz
```

---

## 5) Trimming with Trimmomatic

**Typical parameters for Illumina PE150 RNA‑seq:**

- `LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36`
- Adapters: `TruSeq3-PE.fa`

**Example (single sample):**
```bash
# If using conda 'trimmomatic' wrapper and large files, give Java heap:
export JAVA_TOOL_OPTIONS="-Xmx32g -Xms4g"

trimmomatic PE -threads 4 -phred33 \
  R1.fastq.gz R2.fastq.gz \
  sample_R1_paired.fq.gz sample_R1_unpaired.fq.gz \
  sample_R2_paired.fq.gz sample_R2_unpaired.fq.gz \
  ILLUMINACLIP:$CONDA_PREFIX/share/trimmomatic*/adapters/TruSeq3-PE.fa:2:30:10:2:True \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

> For parallel batch trimming, keep **per-job threads low (2–4)** and **set heap per job (8–16 GB)**; scale number of jobs to your RAM/cores.

---

## 6) STAR Genome Index (V5)

```bash
STAR --runMode genomeGenerate \
  --runThreadN 104 \
  --genomeDir reference/Zea_mays_RefGen_V5_star_index \
  --genomeFastaFiles reference/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa.gz \
  --sjdbGTFfile reference/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.62.gtf.gz \
  --sjdbOverhang 149
```

---

## 7) STAR Mapping (Parallel on All Years)

**scripts/run_star.sh**
```bash
#!/usr/bin/env bash
set -euo pipefail

GENOME_DIR="reference/Zea_mays_RefGen_V5_star_index"
OUT_DIR="mapping"
LOG_DIR="logs"
mkdir -p "$OUT_DIR" "$LOG_DIR"

READ_CMD="zcat"
PARALLEL_JOBS=12
THREADS_PER_JOB=10
BAMSORT_RAM=$((16 * 1024 * 1024 * 1024))   # 16 GB

export GENOME_DIR OUT_DIR LOG_DIR READ_CMD THREADS_PER_JOB BAMSORT_RAM

SAMPLE_LIST=$(mktemp)
ls data/trimmed_fastq_file_202*/*_R1_paired.fq.gz > "$SAMPLE_LIST"

cat "$SAMPLE_LIST" | parallel -j "$PARALLEL_JOBS" --halt now,fail=1 --eta '
  r1="{}"
  prefix=$(basename "$r1" "_R1_paired.fq.gz")
  r2="${r1/_R1_paired.fq.gz/_R2_paired.fq.gz}"

  STAR \
    --genomeDir "$GENOME_DIR" \
    --readFilesCommand "$READ_CMD" \
    --readFilesIn "$r1" "$r2" \
    --runThreadN "$THREADS_PER_JOB" \
    --outFileNamePrefix "$OUT_DIR/${prefix}_" \
    --quantMode GeneCounts \
    --outSAMtype BAM SortedByCoordinate \
    --limitBAMsortRAM "$BAMSORT_RAM" \
    --outReadsUnmapped Fastx \
    > "$LOG_DIR/${prefix}.out" 2> "$LOG_DIR/${prefix}.err"
'
```

**Run:**
```bash
chmod +x scripts/run_star.sh
bash scripts/run_star.sh
```

Outputs per sample:
```
mapping/<sample>_Aligned.sortedByCoord.out.bam
mapping/<sample>_ReadsPerGene.out.tab
mapping/<sample>_Log.final.out
```

---

## 8) Salmon Index (V5 cDNA)

```bash
salmon index \
  -t reference/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.cdna.all.fa.gz \
  -i reference/salmon_index_Zea_mays_RefGen_V5 \
  -k 31
```

---

## 9) Salmon Quantification (Parallel on All Years)

**scripts/run_salmon_quant.sh**
```bash
#!/usr/bin/env bash
set -euo pipefail

INDEX_DIR="reference/salmon_index_Zea_mays_RefGen_V5"
OUT_DIR="salmon_quant"
LOG_DIR="logs/salmon"
mkdir -p "$OUT_DIR" "$LOG_DIR"

PARALLEL_JOBS=16
THREADS_PER_JOB=8

export INDEX_DIR OUT_DIR LOG_DIR THREADS_PER_JOB

SAMPLE_LIST=$(mktemp)
ls data/trimmed_fastq_file_202*/*_R1_paired.fq.gz > "$SAMPLE_LIST"

cat "$SAMPLE_LIST" | parallel -j "$PARALLEL_JOBS" --halt now,fail=1 --eta '
  r1="{}"
  sample=$(basename "$r1" "_R1_paired.fq.gz")
  r2="${r1/_R1_paired.fq.gz/_R2_paired.fq.gz}"

  salmon quant \
    -i "$INDEX_DIR" \
    -l A \
    -1 "$r1" \
    -2 "$r2" \
    -p "$THREADS_PER_JOB" \
    --validateMappings \
    --seqBias --gcBias \
    --writeUnmappedNames \
    -o "$OUT_DIR/$sample" \
    > "$LOG_DIR/${sample}.out" 2> "$LOG_DIR/${sample}.err"
'
```

**Run:**
```bash
chmod +x scripts/run_salmon_quant.sh
bash scripts/run_salmon_quant.sh
```

Outputs per sample:
```
salmon_quant/<sample>/quant.sf
salmon_quant/<sample>/meta_info.json
salmon_quant/<sample>/lib_format_counts.json
```

---

## 10) Build `tx2gene` and Aggregate Gene-level TPM/Counts (R / tximport)

Create `tx2gene` from the GTF:

```bash
awk '$3=="transcript"{
  split($0,a,"transcript_id \"");
  split(a[2],b,"\"");
  split($0,c,"gene_id \"");
  split(c[2],d,"\"");
  print b[1]"\t"d[1]
}' Zea_mays.Zm-B73-REFERENCE-NAM-5.0.62.gtf > tx2gene.tsv

```

**scripts/tximport_salmon.R**
```r
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tximport)
  library(readr)
})

salmon_parent <- "salmon_quant"
tx2gene_file  <- "reference/Zea_mays.tx2gene.tsv"
out_dir       <- "final_results"
out_prefix    <- "maize_salmon_gene"

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

quant_files <- list.files(salmon_parent, pattern="^quant\\.sf$", recursive=TRUE, full.names=TRUE)
if (length(quant_files) == 0L) stop("No quant.sf files found under: ", salmon_parent)

sample_names <- basename(dirname(quant_files))
names(quant_files) <- sample_names

tx2gene <- read_tsv(tx2gene_file, col_names=c("TX","GENE"), show_col_types = FALSE)

txi <- tximport(quant_files, type="salmon", tx2gene=tx2gene, countsFromAbundance="no")

write.table(round(txi$abundance, 6),
            file=file.path(out_dir, paste0(out_prefix, "_TPM.tsv")),
            sep="\t", quote=FALSE, col.names=NA)

write.table(round(txi$counts),
            file=file.path(out_dir, paste0(out_prefix, "_counts.tsv")),
            sep="\t", quote=FALSE, col.names=NA)

write.table(round(txi$length, 2),
            file=file.path(out_dir, paste0(out_prefix, "_length.tsv")),
            sep="\t", quote=FALSE, col.names=NA)
```

---

## 11) Strandness Inference (for featureCounts)

Convert GTF → BED:

```bash
gtf2bed < reference/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.62.gtf > reference/Zea_mays.bed
```

Infer strandness (RSeQC):

```bash
infer_experiment.py -i mapping/<sample>_Aligned.sortedByCoord.out.bam -r reference/Zea_mays.bed
```

Interpretation:
- Majority of reads in pattern **`1+-,1-+,2++,2--`** → **reverse-stranded** → use `-s 2`
- Majority in **`1++,1--,2+-,2-+`** → **forward-stranded** → use `-s 1`
- If ambiguous or unstranded → `-s 0`

**Empirical result in this project:** **reverse-stranded** (use `-s 2`).

---

## 12) Gene-level Counts with featureCounts

```bash
featureCounts \
  -T 16 \
  -p \
  -s 2 \
  -t exon \
  -g gene_id \
  -a reference/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.62.gtf.gz \
  -o feature_count/gene_counts.txt \
  mapping/*.bam
```

Output:
```
feature_count/gene_counts.txt
```

---

## 13) Summary of Key Outputs

- **STAR BAMs** → `mapping/*_Aligned.sortedByCoord.out.bam`  
- **STAR per-sample gene counts** → `mapping/*_ReadsPerGene.out.tab`  
- **Salmon quant (per sample)** → `salmon_quant/<sample>/quant.sf`  
- **Gene-level TPM (Salmon)** → `final_results/maize_salmon_gene_TPM.tsv`  
- **Gene-level counts (Salmon)** → `final_results/maize_salmon_gene_counts.tsv`  
- **featureCounts gene counts** → `feature_count/gene_counts.txt`  
- **QC reports:** raw (optional) → `qc/multiqc_raw`, trimmed → `qc/multiqc_trimmed`

