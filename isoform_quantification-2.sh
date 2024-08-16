#!/bin/bash -l
#$ -P langchip
#$ -N STAR_Lang
#$ -j y
#$ -M luluj@bu.edu
#$ -m ae
#$ -l cpu_arch=!bulldozer
#$ -o qlog/
#$ -cwd
#$ -l h_rt=12:00:00
#$ -pe omp 8
#$ -l mem_per_core=16G


module load star/2.7.10b
module load python3/3.8.10
module load cutadapt/3.4
module load R/4.2.1
module load rsem/1.3.3

STARREFGENOME="/projectnb/langchip/Data/Genomes/hg38_125_gencodeV38_STAR2.7.10b/"
GENCODEGTF="/projectnb/langchip/Data/Genomes/gencode.v38.annotation.gtf"
TEMPLOC="/projectnb/langchip/PersonalFolder/Lulu/temp3/"
ADDRESS=luluj@bu.edu

BASEDIR="/projectnb/langchip/Data/2405_Noah_NMD_NovoGene/data/01.RawData/"
sample=$1
CPU=8

ADAPTERS_FILE="/projectnb/langchip/PersonalFolder/Simon/adapters.fasta"
REF="/projectnb/langchip/PersonalFolder/Lulu/ref/"


# Removing adapter sequences
echo "Running cutadapt to remove adapter sequences for ${sample}..."
cutadapt -j ${CPU} -a file:${ADAPTERS_FILE} -A file:${ADAPTERS_FILE} --minimum-length 1 \
  -o ${TEMPLOC}/${sample}_trimmed.R1.fastq -p ${TEMPLOC}/${sample}_trimmed.R2.fastq \
  ${BASEDIR}/${sample}/${sample}_1.fq.gz ${BASEDIR}/${sample}/${sample}_2.fq.gz

# Alignment with STAR
echo "Running STAR alignment for ${sample}..."
STAR --genomeDir ${STARREFGENOME} --quantMode TranscriptomeSAM --readFilesIn \
  ${TEMPLOC}/${sample}_trimmed.R1.fastq ${TEMPLOC}/${sample}_trimmed.R2.fastq \
  --runThreadN ${CPU} --outFileNamePrefix ${TEMPLOC}/${sample}_ --outSAMunmapped Within

# 
# rsem-calculate-expression -p ${CPU} --paired-end --bam --estimate-rspd --append-names --output-genome-bam \
#   ${TEMPLOC}/${sample}_Aligned.toTranscriptome.out.bam ${REF}/human_gencode ${sample}

rsem-calculate-expression -p ${CPU} --paired-end --bam --estimate-rspd --append-names --output-genome-bam \
  --keep-intermediate-files ${TEMPLOC}/${sample}_Aligned.toTranscriptome.out.bam ${REF}/human_gencode ${sample}





