# VariantCalling
Danya Hassan GENE 5120 12/9/2024

# HUMAN VARIANT CALL
## In this project, we want to take an input of two FASTQ files of an individual's genomic sequence reads, and create a VCF file with it to analyze their human genomic variants
### Necessary softwares, and versions that I used:
> Software versions:
> bwa version (0.7.17-6).
> samtools version (1.13-4).
> bcftools is already the newest version (1.13-1).
> Pandas version 2.2.2
> cyvcf2 version 0.31.1
> matplotlib.pyplot

#### Overall steps: Download the files, download the hg38 reference genome and unzip it to get your database. Then utilize bwatools, samtools, and cyvcf2 to index the reference genome, create bam files for compressed alignment, and then make the VCF files.

```
!bwa index -p Human Human38/Human.fasta
```
```
#bwa alignment
!bwa aln Human HumanWGS.SUBSET.R1.fastq > HumanWGS.R1.sai
!bwa aln Human HumanWGS.SUBSET.R2.fastq > HumanWGS.R2.sai
!bwa sampe Human HumanWGS.R1.sai HumanWGS.R2.sai HumanWGS.SUBSET.R1.fastq HumanWGS.SUBSET.R2.fastq > Humanbwa.bam
```
```
!samtools view -bS -q 12 Humanbwa.bam | samtools sort -o Humanbwa.bam.sorted
#THIS ALIGNS THE READS.
## The pipe (|) says that instead of printing the output out, input it
##into the next operator. Sort it into the output file SARS.bam.sorted
```
```
!bcftools mpileup -f Human38/Human.fasta --max-depth 2000 Humanbwa.bam.sorted  | bcftools call --multiallelic-caller --variants-only --ploidy 2 -mv -Oz -o Human.Variants.vcf.gz
```
```
!gzip -df Human.Variants.vcf.gz
```
