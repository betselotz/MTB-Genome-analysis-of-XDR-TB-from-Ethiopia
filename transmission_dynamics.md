🧬 MTB SNP-based Transmission Analysis

This repository contains a workflow to analyze Mycobacterium tuberculosis (MTB) whole genome sequences using SNP-based clustering to infer transmission dynamics.

The workflow starts from consensus genome sequences and produces SNP alignments, distance matrices, phylogenetic trees, and clusters.
```bash
.
├── consensus_sequences/        # Contains all consensus FASTA files or concatenated all_consensus.fasta
├── core_results/               # Output directory for alignments, SNPs, distances, and trees
├── run_snp_clustering.sh       # Optional bash script to automate the workflow
└── README.md
```

1. Prerequisites

You need the following tools installed (via Conda or system):

MAFFT → multiple sequence alignment

snp-sites → extract SNP sites

snp-dists → compute pairwise SNP distances

IQ-TREE2 → phylogenetic tree construction

R (optional) → cluster analysis and visualization


Example Conda environment:
```bash
conda create -n snp_env -c bioconda -c conda-forge snp-sites snp-dists -y
conda activate snp_env
```
Install R via Conda (recommended for bioinformatics)

```bash
conda create -n r_env -c conda-forge r-base -y
conda activate r_env
```
We can also install additional R packages at the same time:
```bash
conda install -n r_env -c conda-forge r-ape r-dendextend
```




prepare Input

Place your concatenated FASTA file here:
```bash
consensus_sequences/all_consensus.fasta
```

Alternatively, we can place individual sample FASTAs in consensus_sequences/ and concatenate them:
```bash
cat consensus_sequences/*.fasta > consensus_sequences/all_consensus.fasta
```

Run Core SNP Workflow
Manual commands:
```bash
# Create output folder
mkdir -p core_results

# Align all sequences
mafft --auto consensus_sequences/all_consensus.fasta > core_results/core.aln

# Extract SNPs only
snp-sites -o core_results/core.snps.aln core_results/core.aln

# Compute pairwise SNP distances
snp-dists core_results/core.snps.aln > core_results/snp_distances.tsv

# Build phylogenetic tree
iqtree2 -s core_results/core.snps.aln -m GTR+G -nt AUTO -bb 1000 -pre core_results/core

```

Optional: R code to take your snp_distances.tsv (P distance matrix) and identify clusters with ≤5 SNP differences.

Start an interactive R session
```bash
conda activate r_env
```

```bash
R
```
Then in the R prompt:

```bash
# Load libraries
library(ape)
library(dendextend)

# Read the SNP distance matrix correctly
dist_matrix <- read.table("core_results/snp_distances.tsv",
                          header=TRUE, 
                          sep="\t",       # tab-separated
                          row.names=1,    # first column = row names
                          check.names=FALSE,
                          stringsAsFactors=FALSE)

# Convert to numeric matrix (in case values are read as character)
dist_matrix <- as.matrix(dist_matrix)
mode(dist_matrix) <- "numeric"

# Check first few rows and dimensions
head(dist_matrix)
dim(dist_matrix)

# Convert to distance object
dist_obj <- as.dist(dist_matrix)

# Hierarchical clustering (average linkage)
hc <- hclust(dist_obj, method="average")

# Cut clusters at ≤5 SNPs
clusters <- cutree(hc, h=5)

# Save clusters
write.table(clusters, file="core_results/cluster_assignments.tsv",
            sep="\t", quote=FALSE, col.names=NA)

# Plot dendrogram with colored branches
dend <- as.dendrogram(hc)
dend <- color_branches(dend, h=5)
png("core_results/cluster_dendrogram.png", width=1200, height=800)
plot(dend, main="MTB SNP Clusters (≤5 SNPs)")
dev.off()

# Print clusters to console
print(clusters)

```


Optional Bash Script

You can automate all steps using run_snp_clustering.sh:

```bash
chmod +x run_snp_clustering.sh
./run_snp_clustering.sh
```
