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
conda create -n mtb_snp_env mafft snp-sites snp-dists iqtree r-base -y
conda activate mtb_snp_env
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

Optional: Cluster isolates in R 
```bash
library(ape)

# Load SNP distance matrix
dist <- as.dist(read.table("core_results/snp_distances.tsv", header=TRUE, row.names=1))

# Hierarchical clustering
hc <- hclust(dist, method="average")
plot(hc, main="MTB SNP-based clustering")

# Cut clusters by SNP threshold (e.g., ≤12 SNPs)
clusters <- cutree(hc, h=12)
print(clusters)
```
≤5 SNPs → likely recent transmission

≤12 SNPs → possible epidemiological link




Optional Bash Script

You can automate all steps using run_snp_clustering.sh:

```bash
chmod +x run_snp_clustering.sh
./run_snp_clustering.sh
```
