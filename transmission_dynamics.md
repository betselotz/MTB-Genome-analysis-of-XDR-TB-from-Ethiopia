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
Step 3: Align FASTA sequences using MAFFT
```bash
conda activate mafft_env
mafft --auto consensus_sequences/all_consensus.fasta > consensus_sequences/all_consensus_aligned.fasta
```
Verify sequence lengths
```bash
conda activate snp_env
seqkit fx2tab -l consensus_sequences/all_consensus_aligned.fasta
```

Step 4: Extract SNPs
```bash
conda activate snp_env
snp-sites -o consensus_sequences/all_consensus_snps.fasta consensus_sequences/all_consensus_aligned.fasta
```
Step 5: Compute Pairwise SNP Distances
```bash
snp-dists consensus_sequences/all_consensus_snps.fasta > consensus_sequences/snp_distances.tsv
```
Step 6: Cluster and Visualize in R
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

Step 7: Annotate Clusters with Metadata

1. Prepare a metadata file (metadata.tsv) with at least these columns:
SampleID    CollectionDate    Location
SRR12345    2025-01-10        AddisAbaba
SRR12346    2025-01-12        BahirDar
...
2. Merge cluster assignments with metadata in R:
```bash
library(tidyverse)

clusters <- read.table("consensus_sequences/cluster_assignments.tsv", header=FALSE, row.names=1, sep="\t")
metadata <- read.table("consensus_sequences/metadata.tsv", header=TRUE, sep="\t")
metadata <- metadata %>% left_join(data.frame(SampleID=rownames(clusters), Cluster=clusters[,1]), by="SampleID")
write.table(metadata, "consensus_sequences/metadata_clusters.tsv", sep="\t", quote=FALSE, row.names=FALSE)
```

Step 8: Build Phylogenetic Tree
```bash
library(ape)

snp_alignment <- read.dna("consensus_sequences/all_consensus_snps.fasta", format="fasta")
tree <- nj(dist.dna(snp_alignment, model="N"))
plot(tree, main="NJ Tree for MTB Transmission")
write.tree(tree, file="consensus_sequences/mtb_tree.nwk")
```
Optional: Color tree tips by clusters:
```bash
tip_colors <- rainbow(length(unique(metadata$Cluster)))[metadata$Cluster]
plot(tree, tip.color=tip_colors, main="NJ Tree with Clusters")
```
Step 9: Build Transmission Network
```bash
library(igraph)

snp_dist <- read.table("consensus_sequences/snp_distances.tsv", header=TRUE, row.names=1, check.names=FALSE, sep="\t")
edges <- which(as.matrix(snp_dist) <= 5, arr.ind=TRUE)
edges_df <- data.frame(from=rownames(snp_dist)[edges[,1]], to=colnames(snp_dist)[edges[,2]])
g <- graph_from_data_frame(edges_df, directed=FALSE)
plot(g, vertex.label=V(g)$name, vertex.size=15, main="MTB Transmission Network (≤5 SNPs)")
```



