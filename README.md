# summer_internship_2024
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


## generate new STAR genome index for STAR2.7.11b
```
STAR --runThreadN 4 --runMode genomeGenerate --genomeDir samples/hg38_125_gencodeV38_STAR2.7.11b/ --genomeFastaFiles samples/GRCh38.primary_assembly.genome.fa --sjdbGTFfile samples/gencode.v38.annotation.gtf --sjdbOverhang 100
```
