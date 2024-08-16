#!/bin/bash -l
#$ -P adscdata
#$ -N isoformQuant_Lang
#$ -j y
#$ -M slianglu@bu.edu
#$ -l cpu_arch=!bulldozer
#$ -o qlog/
#$ -cwd
#$ -l h_rt=12:00:00
#$ -pe omp 4
#$ -l mem_per_core=4G


#qsub isoform_quantification.sh siPAX3-1 temp
#qsub isoform_quantification.sh siPAX3-2 
#qsub isoform_quantification.sh siPAX3-3 
#qsub isoform_quantification.sh siScramble-1
#qsub isoform_quantification.sh siScramble-2 temp2
#qsub isoform_quantification.sh siScramble-3 

module load star/2.7.2b
module load python3/3.8.10
module load cutadapt/3.4
module load rsem/1.3.3

STARREFGENOME="/projectnb/adscdata/CZ/Genomes/Human/hg38_125_gencodeV38_STAR2.7.2b/"
GENCODEGTF="/projectnb/adscdata/CZ/Genomes/Human/GTF/gencode.v38.annotation.gtf"
TEMPLOC="/projectnb/czlab/A00/simon/Lang"
ADDRESS=slianglu@bu.edu

BASEDIR="/projectnb/adscdata/A04_Lang_RNASeq/azenta"
# sample="siScramble-2"
# sample="siPAX3-2"
#BASEDIR=$1
sample=$1
TEMP=$2
CPU=4

## Removing adaptor
echo "cutadapt"
cutadapt -j ${CPU} -a file:adapters.fasta -A file:adapters.fasta --minimum-length 1 -o ${TEMPLOC}/${TEMP}/${sample}_trimmed.R1.fastq -p ${TEMPLOC}/${TEMP}/${sample}_trimmed.R2.fastq ${BASEDIR}/00_fastq/${sample}_R1_001.fastq.gz ${BASEDIR}/00_fastq/${sample}_R2_001.fastq.gz


## Alignment with transcriptome quantification
echo "Alignment"
echo "Running Star for mRNA..."
# STAR --genomeDir ${STARREFGENOME} --readFilesIn ${TEMPLOC}/${TEMP}/${sample}_R1_001.fq ${TEMPLOC}/${TEMP}/${sample}_R2_001.fq --runThreadN ${CPU} --outFileNamePrefix ${TEMPLOC}/${TEMP}/ --outSAMunmapped Within
STAR --genomeDir ${STARREFGENOME} --quantMode TranscriptomeSAM --readFilesIn ${TEMPLOC}/${TEMP}/${sample}_trimmed.R1.fastq ${TEMPLOC}/${TEMP}/${sample}_trimmed.R2.fastq --runThreadN ${CPU} --outFileNamePrefix ${TEMPLOC}/${TEMP}/ --outSAMunmapped Within

# Prepare reference for rsem-calculate-expression
# GENOME="/projectnb/adscdata/CZ/Genomes/Human/GTF/hg38.p13.chronly.fa"
# rsem-prepare-reference --gtf ${GENCODEGTF} --star ${GENOME} -p 4 ref/human_gencode


# RSEM isoform quantification
echo "Running RSEM for isoform quantification..."
rsem-calculate-expression -p 4 \
                          --paired-end \
                          --bam \
                          --estimate-rspd \
                          --append-names \
                          --output-genome-bam \
                          # ${TEMPLOC}/star/${sample}/isoform/Aligned.toTranscriptome.out.bam \
                          ${TEMPLOC}/${TEMP}/Aligned.toTranscriptome.out.bam \
                          ref/human_gencode \
                          ${sample}
