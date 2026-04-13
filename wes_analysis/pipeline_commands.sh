#!/bin/bash
# WES Pipeline — LRRK2 S1647T Variant Discovery
# Dataset: SRR19500742

# Step 1: Download raw data
prefetch SRR19500742
fastq-dump --split-files SRR19500742

# Step 2: Quality control
fastqc SRR19500742_1.fastq SRR19500742_2.fastq

# Step 3: Trimming
trimmomatic PE SRR19500742_1.fastq SRR19500742_2.fastq \
    clean_R1.fastq unpaired_R1.fastq \
    clean_R2.fastq unpaired_R2.fastq \
    ILLUMINACLIP:adapters.fa:2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# Step 4: Alignment
bwa mem -t 8 GRCh38.fa clean_R1.fastq clean_R2.fastq > aligned.sam
samtools sort aligned.sam -o sorted.bam
samtools index sorted.bam

# Step 5: Mark duplicates
picard MarkDuplicates I=sorted.bam O=dedup.bam M=metrics.txt
samtools index dedup.bam

# Step 6: Variant calling
gatk HaplotypeCaller -R GRCh38.fa -I dedup.bam -O raw_variants.g.vcf -ERC GVCF
gatk GenotypeGVCFs -R GRCh38.fa -V raw_variants.g.vcf -O genotyped.vcf

# Step 7: Hard filtering
gatk VariantFiltration -R GRCh38.fa -V genotyped.vcf \
    --filter-expression "QD < 2.0" --filter-name "QD2" \
    --filter-expression "FS > 60.0" --filter-name "FS60" \
    --filter-expression "MQ < 40.0" --filter-name "MQ40" \
    -O filtered.vcf

# Step 8: Annotation
perl annovar/table_annovar.pl filtered.vcf annovar/humandb/ \
    -buildver hg38 -out annotated \
    -remove -protocol refGene,clinvar_20221231,gnomad30_genome \
    -operation g,f,f -nastring . -vcfinput
