# summer_internship_2024
```
project_directory/
├── raw_data/
│   ├── A1/
│   │   ├── A1_1.fastq.gz
│   │   ├── A1_2.fastq.gz
│   ├── A2/
│   │   ├── A2_1.fastq.gz
│   │   ├── A2_2.fastq.gz
│   ├── A3/
│   │   ├── A3_1.fastq.gz
│   │   ├── A3_2.fastq.gz
│   ├── B1/
│   │   ├── B1_1.fastq.gz
│   │   ├── B1_2.fastq.gz
│   ├── B2/
│   │   ├── B2_1.fastq.gz
│   │   ├── B2_2.fastq.gz
│   ├── B3/
│   │   ├── B3_1.fastq.gz
│   │   ├── B3_2.fastq.gz
│   ├── C1/
│   │   ├── C1_1.fastq.gz
│   │   ├── C1_2.fastq.gz
│   ├── C2/
│   │   ├── C2_1.fastq.gz
│   │   ├── C2_2.fastq.gz
│   ├── C3/
│   │   ├── C3_1.fastq.gz
│   │   ├── C3_2.fastq.gz
│   ├── D1/
│   │   ├── D1_1.fastq.gz
│   │   ├── D1_2.fastq.gz
│   ├── D2/
│   │   ├── D2_1.fastq.gz
│   │   ├── D2_2.fastq.gz
│   ├── D3/
│   │   ├── D3_1.fastq.gz
│   │   ├── D3_2.fastq.gz
│   ├── E1/
│   │   ├── E1_1.fastq.gz
│   │   ├── E1_2.fastq.gz
│   ├── E2/
│   │   ├── E2_1.fastq.gz
│   │   ├── E2_2.fastq.gz
│   ├── E3/
│   │   ├── E3_1.fastq.gz
│   │   ├── E3_2.fastq.gz
│   ├── F1/
│   │   ├── F1_1.fastq.gz
│   │   ├── F1_2.fastq.gz
│   ├── F2/
│   │   ├── F2_1.fastq.gz
│   │   ├── F2_2.fastq.gz
│   ├── F3/
│   │   ├── F3_1.fastq.gz
│   │   ├── F3_2.fastq.gz
├── reference/
│   ├── GRCh38.primary_assembly.genome.fa
│   ├── gencode.v38.annotation.gtf
├── results/
│   ├── alignments/
│   ├── sorted_bams/
│   ├── counts/
├── scripts/
│   ├── run_cutadapt.sh
│   ├── Snakefile

```

```
conda install STAR
```
```
pip install multiqc
```
```
pip install cutadapt
```
``` 
conda install -c bioconda subread
```
```
conda install -c bioconda snakemake
```


## 1. Quality Control and Trimming
Run FastQC and Cutadapt for quality control and trimming of the raw sequencing reads.

## 2. Generate STAR Genome Index
Generate a new STAR genome index using STAR version 2.7.11b:
```
STAR --runThreadN 4 --runMode genomeGenerate --genomeDir samples/hg38_125_gencodeV38_STAR2.7.11b/ --genomeFastaFiles samples/GRCh38.primary_assembly.genome.fa --sjdbGTFfile samples/gencode.v38.annotation.gtf --sjdbOverhang 100
```
