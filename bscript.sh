#!/bin/bash
# Pre Alignment QC
# 1.0 fastqc on raw data
fastqc /home/ubuntu/assessment/dnaseq/data/untrimmed_fastq/NGS0001.R1.fastq.gz \
       /home/ubuntu/assessment/dnaseq/data/untrimmed_fastq/NGS0001.R2.fastq.gz
# 2.0 trimmomatic
trimmomatic PE -threads 4 -phred33 /home/ubuntu/assessment/dnaseq/data/untrimmed_fastq/NGS0001.R1.fastq.gz /home/ubuntu/assessment/dnaseq/data/untrimmed_fastq/NGS0001.R2.fastq.gz -baseout /home/ubuntu/assessment/dnaseq/data/trimmed_fastq/NGS0001_trimmed_R ILLUMINACLIP:/home/ubuntu/anaconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa:2:30:10 TRAILING:25 MINLEN:50
# 3.0 fastqc on paired trimmed fastq data
fastqc /home/ubuntu/assessment/dnaseq/data/trimmed_fastq/NGS0001_trimmed_R_1P /home/ubuntu/assessment/dnaseq/data/trimmed_fastq/NGS0001_trimmed_R_2P
# Alignment
# 0.0 Index the reference genome
bwa index /home/ubuntu/assessment/dnaseq/data/reference/hg19.fa.gz
# 1.0 alignment with BWA on trimmed data
bwa mem -t 4 -v 1 -R '@RG\tID:11V6WR1:111:D1375ACXX:1:NGS0001\tSM:NGS0001\tPL:ILLUMINA\tLB:nextera-NGS0001-blood\tDT:2025-03-17\tPU:11V6WR1 -I 250,50  ~/assessment/dnaseq/data/reference/hg19.fa.gz ~/assessment/dnaseq/data/trimmed_fastq/NGS0001_trimmed_R_1P ~/assessment/dnaseq/data/trimmed_fastq/NGS0001_trimmed_R_2P > ~/assessment/dnaseq/data/aligned_data/NGS0001.sam
# 2.0 convert sam to bam, sort bam
samtools view -h -b NGS0001.sam > NGS0001.bam
samtools sort NGS0001.bam > NGS0001_sorted.bam
# 2.1 index sorted bam
samtools index NGS0001_sorted.bam
# Basic Alignment post processing
# 1.0 Mark Duplicate
picard MarkDuplicates I=NGS0001_sorted.bam O=NGS0001_sorted_marked.bam M=marked_dup_metrics.txt
# 2.0 Filter BAM file
samtools view -F 1796  -q 20 -o NGS0001_sorted_filtered.bam NGS0001_sorted_marked.bam
# 2.1 Index filtered bam
samtools index NGS0001_sorted_filtered.bam
# 2.2 Alignment statistics
# 2.2.1 Flagstats
samtools flagstat NGS0001_sorted_filtered.bam > NGS0001_sorted_filtered_flagstat.txt
# 2.2.2 Idxstats
samtools idxstats NGS0001_sorted_filtered.bam > NGS0001_sorted_filtered_idxstats.txt
# 2.2.3 Depth of coverage
samtools depth NGS0001_sorted_filtered.bam > NGS0001_sorted_filtered_depth.txt
# 2.2.3 Insert size
samtools stats NGS0001_sorted_filtered.bam > NGS0001_sorted_filtered_stats.txt
# Basic Variant Calling
# 1.0 Call Variants with Freebayes
freebayes --bam ~/assessment/dnaseq/data/aligned_data/NGS0001_sorted_filtered.bam --fasta-reference ~/assessment/dnaseq/data/reference/hg19.fa --vcf ~/assessment/dnaseq/results/NGS0001.vcf
# 2.0 Hard Filter Variants
zcat ~/assessment/dnaseq/results/NGS0001.vcf.gz | vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" > ~/assessment/dnaseq/results/NGS0001_filtered.vcf
# ANNOVAR 
# 1.0 convert VCF to ANNOVAR input format 
./convert2annovar.pl -format vcf4 
~/assessment/dnaseq/results/NGS0001_filtered_annotation.vcf.gz > ~/assessment/dnaseq/results/NGS0001_filtered_annotation.avinput
# 2.0 Basic Variant Annotation
perl /home/ubuntu/annovar/table_annovar.pl ~/assessment/dnaseq/results/NGS0001_filtered_annotation.avinput /home/ubuntu/annovar/humandb -buildver hg19 -out ~/assessment/dnaseq/results/NGS0001_filtered_annotation -remove -protocol refGene,ensGene,clinvar_20180603,exac03,dbnsfp31a_interpro -operation g,g,f,f,f -otherinfo -nastring . –csvout
# 3.0 Basic Variant Prioritization
awk -F NR==1 || ($6 ~ /exonic/ && $11 == ".") ~/assessment/dnaseq/results/NGS0001_filtered_annotation_dbsnp.hg19_multianno.csv > ~/assessment/dnaseq/results/exonic_variants_not_in_dbSNP.csv
# Exit the script
exit 0