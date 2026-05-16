# LIFE750 Cycle 2 — Command Record

## Software and Versions
- BWA v0.7.19
- SAMtools v1.21
- FreeBayes v1.3.10
- BCFtools v1.21
- BEDtools v2.31.1
- DESeq2 (R package)
- IGV (Integrative Genomics Viewer)

---

## Setup

Download dataset:
```bash
wget -R "index.html*" -r -nd -nH -np https://cgr.liv.ac.uk/454/acdarby/LIFE750/Life750_datasets/dataset_hliabba2/
```

Create and activate conda environment:
```bash
conda create -n life750 -c bioconda -c conda-forge bwa samtools freebayes bcftools bedtools -y
conda activate life750
```

---

## Variant Calling

### 1. Index the reference sequence
```bash
bwa index GeneX_reference.fa
```

### 2. Map paired-end reads to the reference
```bash
bwa mem -R '@RG\tID:normal\tSM:normal\tPL:ILLUMINA\tLB:WGS\tPU:HiSeq2500' \
GeneX_reference.fa normal_GeneX_R1.fastq normal_GeneX_R2.fastq > normal_GeneX.sam

bwa mem -R '@RG\tID:mutant\tSM:mutant\tPL:ILLUMINA\tLB:WGS\tPU:HiSeq2500' \
GeneX_reference.fa mutant_GeneX_R1.fastq mutant_GeneX_R2.fastq > mutant_GeneX.sam
```

### 3. Convert SAM to BAM, sort, and index
```bash
samtools view -b normal_GeneX.sam | samtools sort -o normal_GeneX_sorted.bam
samtools index normal_GeneX_sorted.bam

samtools view -b mutant_GeneX.sam | samtools sort -o mutant_GeneX_sorted.bam
samtools index mutant_GeneX_sorted.bam
```

### 4. Call variants with FreeBayes
Gene X was cloned from a plasmid (haploid), so --ploidy 1 is used.
```bash
freebayes --ploidy 1 --min-mapping-quality 30 --min-base-quality 25 \
--fasta-reference GeneX_reference.fa normal_GeneX_sorted.bam > normal_GeneX.vcf

freebayes --ploidy 1 --min-mapping-quality 30 --min-base-quality 25 \
--fasta-reference GeneX_reference.fa mutant_GeneX_sorted.bam > mutant_GeneX.vcf
```

### 5. Filter variants
```bash
bcftools view -e 'QUAL < 10 || INFO/DP < 10 || INFO/AF = 0' \
mutant_GeneX.vcf -O v -o mutant_GeneX_filtered.vcf
```

### 6. Visualise alignments and variants in IGV
1. Genomes > Load Genome from File > GeneX_reference.fa
2. File > Load from File > normal_GeneX_sorted.bam, mutant_GeneX_sorted.bam, mutant_GeneX_filtered.vcf
3. Navigate to positions 456, 557, and 638 to visualise each variant

---

## Differential Gene Expression Analysis

Open RStudio and set working directory to the dataset folder.

### 1. Load libraries and data
```r
library(DESeq2)
library(tidyverse)

counts_data <- read.table("gene_counts.tsv", header=TRUE, row.names=1, sep="\t")
colData <- read.table("sample_metadata.tsv", header=TRUE, row.names=1, sep="\t")
```

### 2. Prepare count data
```r
# Reorder counts columns to match metadata
counts_data <- counts_data[, rownames(colData)]

# Check columns match rows
all(colnames(counts_data) %in% rownames(colData))
all(colnames(counts_data) == rownames(colData))
```

### 3. Construct DESeq2 object
```r
dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = colData,
                              design = ~ condition)

# Pre-filter low count genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Set reference level
dds$condition <- relevel(dds$condition, ref = "normal")
```

### 4. Run DESeq2
```r
dds <- DESeq(dds)
res <- results(dds)
res <- res[order(res$padj),]
summary(res)
```

### 5. Apply stricter threshold and explore results
```r
res0.01 <- results(dds, alpha = 0.01)
res0.01 <- res0.01[order(res0.01$padj),]
summary(res0.01)

# Most upregulated and downregulated genes
sig_genes <- as.data.frame(res0.01[!is.na(res0.01$padj) & res0.01$padj < 0.01,])
most_up <- sig_genes[which.max(sig_genes$log2FoldChange),]
most_down <- sig_genes[which.min(sig_genes$log2FoldChange),]

# Save results
write.csv(sig_genes, file="DESeq2_significant_0.01.csv")
```

### 6. MA plot
```r
plotMA(res, alpha=0.01,
       main="Differential Expression: Mutant vs Normal Gene X",
       ylab="Log2 Fold Change",
       xlab="Mean of Normalised Counts")
abline(h=0, col="black", lty=1)
```

---

## Intersecting DEGs with Gene X Binding Sites

### 1. Extract gene coordinates from GFF into BED format
```bash
awk '$3=="gene" {print $1"\t"$4-1"\t"$5"\t"$9}' genes.gff3 | \
sed 's/ID=//;s/;.*//' > all_genes.bed
```

### 2. Extract significant DEG IDs
```bash
tail -n +2 DESeq2_significant_0.01.csv | cut -d',' -f1 | sed 's/"//g' > sig_gene_ids.txt
```

### 3. Filter BED file to significant DEGs only
```bash
grep -Fw -f sig_gene_ids.txt all_genes.bed > sig_genes.bed
```

### 4. Intersect DEGs with Gene X binding sites
```bash
bedtools intersect -a sig_genes.bed -b gene_x_binding_sites.bed -wa -wb > intersect_results.bed
```

### 5. Merge intersection results with DESeq2 output
```r
intersect <- read.table("intersect_results.bed", sep="\t",
    col.names=c("chr","start","end","gene_id","peak_chr",
                "peak_start","peak_end","peak_name"))

res_df <- as.data.frame(res0.01)
res_df$gene_id <- rownames(res_df)

final_table <- merge(intersect, res_df, by="gene_id")
final_table <- final_table[order(final_table$padj),
    c("gene_id","log2FoldChange","padj","peak_name")]

write.csv(final_table, "GeneX_direct_targets.csv", row.names=FALSE)
```
