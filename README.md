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
```
!pip install cyvcf2
vcf_path = 'Human.Variants.vcf'

## This is a new python library that is a high-performance Python library for reading and manipulating Variant Call Format (VCF) files.
#It is a Cython wrapper around the htslib C library, providing a fast and efficient way to access VCF data in Python.
```
```
import cyvcf2

# Create a VCF reader object
vcf_reader = cyvcf2.VCF(vcf_path)
#using the cyvcf2.VCF function on vcf_path

with open('humanvariants.tsv', 'w') as file:
# Loop through each variant in the VCF file
    for variant in vcf_reader:
    # Access various properties of the variant
        genotype = variant.gt_types[0]
        data_string = f"{variant.CHROM}\t{variant.POS}\t{variant.ID}\t{variant.REF}\t{','.join(variant.ALT)}\t{variant.QUAL}\t{variant.FILTER}\t{genotype}\n"
        #For every row/variant, make the genotype and data_string
        # Write the data string to the file
        file.write(data_string)

#Here, we are making a text file called humanvariants.tsv
```
```
import pandas
columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'GENOTYPE']
df1 = pandas.read_csv("humanvariants.tsv", delimiter = "\t", header=None, names=columns)
```

### Now we have a final VCF file, in the form of a pandas dataframe, and we can see all of the genomic variants in this sample.
