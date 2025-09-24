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


Create and activate SNP environment:
```bash
conda create -n snp_env -c bioconda -c conda-forge snp-sites snp-dists -y
conda activate snp_env
```
Create and activate R environment:

```bash
conda create -n r_env -c conda-forge r-base -y
conda activate r_env
conda install -n r_env -c conda-forge r-ape r-dendextend -y
```


Step 2: Prepare Input

Concatenate all FASTA files:
```bash
cat consensus_sequences/*.fasta > consensus_sequences/all_consensus.fasta
```
Step 3: Extract SNPs
```bash
conda activate snp_env
snp-sites -o consensus_sequences/all_consensus_snps.fasta consensus_sequences/all_consensus.fasta
```
Step 4: Compute Pairwise SNP Distances
```bash
snp-dists consensus_sequences/all_consensus_snps.fasta > consensus_sequences/snp_distances.tsv
```
Step 5: Cluster Genomes in R
```bash
conda activate r_env
R
```

Inside R:
```bash
library(ape)
library(dendextend)

snp_dist <- read.table("consensus_sequences/snp_distances.tsv", header=TRUE, row.names=1, check.names=FALSE, sep="\t")
dist_matrix <- as.matrix(snp_dist)
dist_obj <- as.dist(dist_matrix)
hc <- hclust(dist_obj, method="average")
clusters <- cutree(hc, h=5)
write.table(clusters, file="consensus_sequences/cluster_assignments.tsv", sep="\t", quote=FALSE, col.names=NA)
dend <- as.dendrogram(hc)
dend <- color_branches(dend, h=5)
plot(dend, main="MTB SNP Clustering (≤5 SNPs)")
pdf("consensus_sequences/snp_clustering_dendrogram.pdf", width=10, height=8)
plot(dend, main="MTB SNP Clustering (≤5 SNPs)")
dev.off()
```


Step 6: Cluster and Visualize in R
```bash
conda activate r_env
R
```
Then run this R code:
```bash
library(ape)
library(dendextend)

snp_dist <- read.table("consensus_sequences/snp_distances.tsv", header=TRUE, row.names=1, check.names=FALSE, sep="\t")
dist_matrix <- as.matrix(snp_dist)
dist_obj <- as.dist(dist_matrix)
hc <- hclust(dist_obj, method="average")
clusters <- cutree(hc, h=5)
write.table(clusters, file="consensus_sequences/cluster_assignments.tsv", sep="\t", quote=FALSE, col.names=NA)
dend <- as.dendrogram(hc)
dend <- color_branches(dend, h=5)
plot(dend, main="MTB SNP Clustering (≤5 SNPs)")
pdf("consensus_sequences/snp_clustering_dendrogram.pdf", width=10, height=8)
plot(dend, main="MTB SNP Clustering (≤5 SNPs)")
dev.off()
```


