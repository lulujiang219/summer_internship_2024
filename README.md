# summer_internship_2024



**Sample Information**: Sample data and information are stored in [sample_info](sample_info)


software: 
running on SCC, but if run on local:
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
![image](https://github.com/lulujiang219/summer_internship_2024/assets/52995244/033cdf41-1663-4fa8-9f97-036716d099ec)


# Workflow
## 1. Quality Control and Trimming
Run FastQC and Cutadapt for quality control and trimming of the raw sequencing reads.

## 2. Generate STAR Genome Index
* Generate a new STAR genome index using STAR version 2.7.11b:
```
STAR --runThreadN 4 --runMode genomeGenerate --genomeDir samples/hg38_125_gencodeV38_STAR2.7.11b/ --genomeFastaFiles samples/GRCh38.primary_assembly.genome.fa --sjdbGTFfile samples/gencode.v38.annotation.gtf --sjdbOverhang 100
```

* Create a bash script (generate_star_genome_index.sh) that submits this job through qsub with the -P projectlang option. 

* Define variables(file location) using "$" within the script.



## 3. Run STAR Alignment for Trimmed FASTQ Files
Create a script (run_star.sh)

* defined a temporary folder temp3/ with each $sample having a individual folder within to prevent overide
* delete the intermediate/temp folder after.

* Run STAR alignment --> Run samtools to index the BAM file --> Run featureCounts to generate counts --> Remove intermediate files 

* Test the script on a single sample before submitting the full job using qsub. Load the necessary modules and run the following commands:
```
# Load required modules
module load python/3.8.10
module load cutadapt/3.4
module load samtools/1.10
module load subread/2.0.1
module load star/2.7.10b

# Define the paths for the test sample
STARREFGENOME="/projectnb/langchip/Data/Genomes/hg38_125_gencodeV38_STAR2.7.10b/"
GTF_FILE="/projectnb/langchip/Data/Genomes/gencode.v38.annotation.gtf"
SAMPLE="C1"
TEMPOUT="/projectnb/langchip/PersonalFolder/Lulu/temp3/${SAMPLE}/"
FINALOUT="/projectnb/langchip/PersonalFolder/Lulu/star/"
R1_FILE="/projectnb/langchip/PersonalFolder/Lulu/trimmed/${SAMPLE}_trimmed.R1.fastq.gz"
R2_FILE="/projectnb/langchip/PersonalFolder/Lulu/trimmed/${SAMPLE}_trimmed.R2.fastq.gz"

# Create necessary directories
mkdir -p $TEMPOUT
mkdir -p $FINALOUT

echo "Running STAR alignment for sample: $SAMPLE"
STAR --genomeDir $STARREFGENOME --readFilesIn $R1_FILE $R2_FILE --runThreadN 4 --readFilesCommand zcat --outFileNamePrefix ${TEMPOUT}${SAMPLE}_ --outSAMtype BAM SortedByCoordinate

echo "Indexing BAM file for sample: $SAMPLE"
samtools index ${TEMPOUT}${SAMPLE}_Aligned.sortedByCoord.out.bam

echo "Running featureCounts for sample: $SAMPLE"
featureCounts -T 4 -a $GTF_FILE -o ${TEMPOUT}${SAMPLE}_counts.txt ${TEMPOUT}${SAMPLE}_Aligned.sortedByCoord.out.bam

```
