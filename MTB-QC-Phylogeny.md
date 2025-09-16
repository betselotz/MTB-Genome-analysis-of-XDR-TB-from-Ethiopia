
## 🧪 Methodology: MTB Genome Analysis of XDR-TB from Ethiopia

This workflow performs **high-confidence genomic analysis** of **extensively drug-resistant Mycobacterium tuberculosis (XDR-TB)** strains from Ethiopia, combining raw data QC, variant calling, and phylogenetic inference.  

### 1️⃣ Explore Raw FASTQ Files
- Inspect **paired-end reads**: count, read length distribution, and base composition (A/T/G/C).  
- Balanced bases ✅ → good-quality sequencing data.

### 2️⃣ Quality Control & Trimming with **FASTP**
- Remove adapters, trim low-quality bases, and filter short reads.  
- Evaluate **pre- vs post-trimming metrics** to measure effectiveness.  
- Generate **per-sample HTML & JSON QC reports**.

### 3️⃣ Aggregate QC Reports with **MultiQC**
- Combine all FASTQ and FASTP QC reports into a **single interactive HTML report**.  
- Quickly identify samples with potential issues.

### 4️⃣ Drug Resistance Screening with **TB-Profiler**
- Check **raw FASTQ quality** and remove extremely low-quality samples.  
- Generate initial **drug-resistance predictions** per isolate.

### 5️⃣ Variant Calling with **Snippy**
- Align reads to **H37Rv reference genome**.  
- Produce **BAM files** and preliminary **VCFs**.

### 6️⃣ BAM Quality Check with **Qualimap**
- Assess **coverage, mapping quality, and read distribution** of BAM files.  
- Aggregate all **Qualimap QC reports** into a **single summary** to easily review per-sample metrics.  
- Ensure reliability for downstream variant filtering.

### 7️⃣ High-Confidence Variant Filtering with **tb_variant_filter**
- Mask **Refined Low Confidence (RLC) regions**.  
- Retain **only high-confidence variants** for accurate analysis.

### 8️⃣ Consensus Genome Generation with **BCFtools** & Outgroup Inclusion
- Generate **per-sample consensus FASTA sequences** using **BCFtools**, incorporating only high-confidence variants.  
- Include **SRR10828835 as an outgroup** to root the phylogenetic tree and provide directionality for evolutionary analysis.

### 9️⃣ Multiple Sequence Alignment using **MAFFT** 
- Align all consensus sequences using **MAFFT** for consistent comparison.  

### 🔟 Phylogenetic Tree Construction with **IQ-TREE**
- Build **maximum-likelihood trees** with **ultrafast bootstrap (1000 replicates)**.  
- Visualize with **TB-Profiler ITOL outputs** to include resistance and metadata.

---

## 🛠 Tool Installation (Separate Conda Environments)
We use separate Conda environments for each tool to ensure reproducibility, avoid software conflicts, and simplify maintenance. Different tools may require different versions of Python or libraries, and installing them together can cause one tool to break when another requires a different version. By isolating each tool in its own environment, we can safely update or reinstall software without affecting the rest of the pipeline, keep the base system clean, and allow other researchers to reproduce the exact setup by exporting and sharing environment files. This approach ensures a robust, conflict-free, and easily maintainable MTB WGS workflow.

### Prerequisites
```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

FASTP
```bash
conda create -n fastp_env fastp=0.23.3 -y
conda activate fastp_env
fastp --version
```

MultiQC
```bash
conda create -n multiqc_env multiqc=1.16 -y
conda activate multiqc_env
multiqc --version
```
TB-Profiler
```bash
conda create -n tbprofiler_env python=3.9 -y
conda activate tbprofiler_env
pip install tb-profiler
tb-profiler --version
```
Snippy
```bash
conda create -n snippy_env snippy=4.6.0 -y
conda activate snippy_env
snippy --version
```
Qualimap
```bash
conda create -n qualimap_env qualimap=2.2.2d -y
conda activate qualimap_env
qualimap --help
```
tb_variant_filter (COMBAT-TB)
```bash
conda create -n tb_variant_filter_env python=3.9 -y
conda activate tb_variant_filter_env
pip install tb-variant-filter
tb_variant_filter --help
```
BCFtools / SAMtools
```bash
conda create -n bcftools_env bcftools=1.17 samtools=1.17 -y
conda activate bcftools_env
bcftools --version
samtools --version
```
MAFFT
```bash
conda create -n mafft_env mafft=7.505 -y
conda activate mafft_env
mafft --version
```
IQ-TREE
```bash
conda create -n iqtree_env iqtree=2.2.2.3 -y
conda activate iqtree_env
iqtree2 --version
```
iTOL
Web-based tool: https://itol.embl.de



# 1️⃣ Explore Raw FASTQ Files

### FASTQ summary

Paired-end FASTQ files are first inspected using simple command-line tools.

### 1. Peek at the first few reads
Before analyzing sequencing data, it’s useful to inspect the first few reads in each FASTQ file. This allows a quick check of file format, sequence identifiers, and the general structure of paired-end reads:
```bash
zcat raw_data/SRR28821350_1.fastq.gz | head -n 16
zcat raw_data/SRR28821350_2.fastq.gz | head -n 16
```
<details>
<summary>🔍 Inspect First Reads of FASTQ Files</summary>
- `zcat` is a Linux/Unix command used to **view the contents of compressed files** without manually uncompressing them.  
- Works primarily with `.gz` files (gzip-compressed).  
- Unlike `gunzip`, it **prints the uncompressed data to standard output** instead of creating a new file.  
- Example workflow:  
  - `zcat raw_data/SRR28821350_1.fastq.gz | head -n 16` →  Show the first 16 lines of R1 (forward) FASTQ file of a compressed file without extracting it. 
- Each read in FASTQ format consists of 4 lines:  
  1. Header line (`@`) with read ID  
  2. Sequence line  
  3. Separator line (`+`)  
  4. Quality scores line  
  - Useful for:  
  - Checking file format  
  - Inspecting sequences and quality  
  - Ensuring R1 and R2 files are correctly paired  
</details>


### 2. Count total reads
To determine the total number of reads in each FASTQ file, divide the total number of lines by 4, since each read consists of four lines (identifier, sequence, optional identifier, and quality scores). This ensures paired-end files are consistent and provides an initial check on sequencing depth:
For single sample
```bash
echo $(( $(zcat raw_data/SRR28821350_1.fastq.gz | wc -l) / 4 ))
echo $(( $(zcat raw_data/SRR28821350_2.fastq.gz | wc -l) / 4 ))
```
<details>
<summary>📊 Counting Reads in FASTQ Files</summary>

- `zcat raw_data/SRR28821350_1.fastq.gz` → Decompresses the R1 FASTQ file and sends its content to standard output.  
- `| wc -l` → Counts the total number of lines in the decompressed file.  
- `$(( ... / 4 ))` → **Arithmetic expansion** in Bash: evaluates the expression inside `$(( ))`.  
  - Divides the total number of lines by 4 because **each sequencing read occupies 4 lines** in a FASTQ file:  
    1. Header line (`@`)  
    2. Sequence line  
    3. Separator line (`+`)  
    4. Quality scores line  
- `echo` → Prints the resulting number of reads.  
- The second command works the same for the R2 (reverse) FASTQ file.  
- Overall, these commands give a **quick read count** for paired-end sequencing files.  

</details>


For batch processing 
##### Step 1: Open a new script
```bash
nano count_reads.sh
```
##### Step 2: Paste the following code
This script counts reads in paired-end FASTQ files and saves results to a CSV.

```bash
#!/bin/bash
set -euo pipefail

INDIR="raw_data"
OUTDIR="csv_output"
mkdir -p "$OUTDIR"

OUTFILE="$OUTDIR/fastq_read_counts.csv"
echo "Sample,R1_reads,R2_reads" > "$OUTFILE"
echo "📊 Counting reads in FASTQ files from '$INDIR'..."

for R1 in "$INDIR"/*_1.fastq.gz "$INDIR"/*_R1.fastq.gz; do
    [[ -f "$R1" ]] || continue
    SAMPLE=$(basename "$R1" | sed -E 's/_R?1.*\.fastq\.gz//')
    R2=""
    for suffix in "_2.fastq.gz" "_R2.fastq.gz" "_R2_*.fastq.gz"; do
        [[ -f "$INDIR/${SAMPLE}${suffix}" ]] && R2="$INDIR/${SAMPLE}${suffix}" && break
    done
    R1_COUNT=$(( $(zcat "$R1" | wc -l) / 4 ))
    R2_COUNT=$([[ -n "$R2" ]] && echo $(( $(zcat "$R2" | wc -l) / 4 )) || echo "NA")
    echo "$SAMPLE,$R1_COUNT,$R2_COUNT" >> "$OUTFILE"
    echo "✅ $SAMPLE → R1: $R1_COUNT | R2: $R2_COUNT"
done

echo "🎉 All done! Read counts saved to '$OUTFILE'"

```

<details>
<summary>📊 FASTQ Read Count Script Overview</summary>

- `INDIR="raw_data"` → Directory containing input FASTQ files.  
- `OUTDIR="csv_output"` → Directory where the CSV file will be saved; created automatically if it doesn’t exist.  
- `OUTFILE="$OUTDIR/fastq_read_counts.csv"` → CSV file to store read counts.  
- `echo "Sample,R1_reads,R2_reads" > "$OUTFILE"` → Creates the CSV header.  
- `echo "📊 Counting reads in FASTQ files from '$INDIR'..."` → Prints starting message.  
- `for R1 in "$INDIR"/*_1.fastq.gz "$INDIR"/*_R1.fastq.gz; do ... done` → Loops through all R1 FASTQ files.  
- `[[ -f "$R1" ]] || continue` → Skips if the R1 file does not exist.  
- `SAMPLE=$(basename "$R1" | sed -E 's/_R?1.*\.fastq\.gz//')` → Extracts sample name from the file name.  
- `for suffix in "_2.fastq.gz" "_R2.fastq.gz" "_R2_*.fastq.gz"; do ... done` → Finds the corresponding R2 file if it exists.  
- `R1_COUNT=$(( $(zcat "$R1" | wc -l) / 4 ))` → Counts reads in R1 by dividing total lines by 4.  
- `R2_COUNT=$([[ -n "$R2" ]] && echo $(( $(zcat "$R2" | wc -l) / 4 )) || echo "NA")` → Counts reads in R2 if present; otherwise outputs "NA".  
- `echo "$SAMPLE,$R1_COUNT,$R2_COUNT" >> "$OUTFILE"` → Appends counts to the CSV file.  
- `echo "✅ $SAMPLE → R1: $R1_COUNT | R2: $R2_COUNT"` → Prints progress for each sample.  
- `echo "🎉 All done! Read counts saved to '$OUTFILE'"` → Prints completion message.  

</details>unsaved

### 3. Base composition
Checking the nucleotide composition of each FASTQ file helps assess sequencing quality. Balanced proportions of A, T, G, and C indicate high-quality data with minimal bias. This script counts the occurrence of each base in both paired-end files:

```bash
#!/bin/bash
for fq in raw_data/SRR28821350_1.fastq.gz raw_data/SRR28821350_2.fastq.gz; do
    [ -f "$fq" ] || continue
    echo "Counting bases in $fq..."
    zcat "$fq" | awk 'NR%4==2 { for(i=1;i<=length($0);i++) b[substr($0,i,1)]++ } 
    END { for(base in b) print base, b[base] }'
    echo "----------------------"
done
```
<details>
<summary>🧬 Nucleotide Counting Script Explanation</summary>

- `for fq in raw_data/SRR28821350_1.fastq.gz raw_data/SRR28821350_2.fastq.gz; do ... done` → Loops over the two specified FASTQ files (R1 and R2).  
- `[ -f "$fq" ] || continue` → Skips the iteration if the file does not exist.  
- `echo "Counting bases in $fq..."` → Prints which file is being processed.  
- `zcat "$fq"` → Decompresses the FASTQ file and streams its content to standard output.  
- `awk 'NR%4==2 { for(i=1;i<=length($0);i++) b[substr($0,i,1)]++ } END { for(base in b) print base, b[base] }'` → Counts nucleotides:  
  - `NR%4==2` → Only processes the **sequence line** of each FASTQ read.  
  - `for(i=1;i<=length($0);i++)` → Iterates over each nucleotide in the sequence line.  
  - `b[substr($0,i,1)]++` → Increments a counter for each base (A, T, G, C, N, or other).  
  - `END { for(base in b) print base, b[base] }` → Prints the total counts for each base after processing the file.  
- `echo "----------------------"` → Adds a visual separator between files for readability.  

</details>


### 4. Quality score summary
FASTQ files encode base quality scores on the 4th line of every read. Checking these scores provides an initial assessment of sequencing quality before trimming or downstream analysis:
First 10 quality lines
```bash
zcat raw_data/SRR28821350_1.fastq.gz | sed -n '4~4p' | head -n 10
zcat raw_data/SRR28821350_2.fastq.gz | sed -n '4~4p' | head -n 10
```
<details>
<summary>🔍 View FASTQ Quality Scores</summary>

- `zcat raw_data/SRR28821350_1.fastq.gz | sed -n '4~4p' | head -n 10`  
  - Decompresses R1 FASTQ.  
  - `sed -n '4~4p'` → Prints every 4th line starting from line 4 (the **quality score line** for each read).  
  - `head -n 10` → Shows only the first 10 quality lines for quick inspection.  

- `zcat raw_data/SRR28821350_2.fastq.gz | sed -n '4~4p' | head -n 10`  
  - Same as above, but for R2 FASTQ.  
</details>

FASTQ quality scores are encoded as ASCII characters. Counting the occurrence of each character provides a quantitative overview of base quality across the reads:
```bash
zcat raw_data/SRR28821350_1.fastq.gz | sed -n '4~4p' | awk '{for(i=1;i<=length($0);i++){q[substr($0,i,1)]++}} END{for (k in q) print k,q[k]}'
zcat raw_data/SRR28821350_2.fastq.gz | sed -n '4~4p' | awk '{for(i=1;i<=length($0);i++){q[substr($0,i,1)]++}} END{for (k in q) print k,q[k]}'
```
<details>
<summary>🔢 Count Quality Score Frequencies</summary>

- `zcat raw_data/SRR28821350_1.fastq.gz | sed -n '4~4p' | awk '{for(i=1;i<=length($0);i++){q[substr($0,i,1)]++}} END{for (k in q) print k,q[k]}'`  
  - Decompresses R1 FASTQ.  
  - `sed -n '4~4p'` → Selects every 4th line (the **quality line**).  
  - `awk '{for(i=1;i<=length($0);i++){q[substr($0,i,1)]++}} END{for (k in q) print k,q[k]}'` → Counts occurrences of each quality score character.  

- `zcat raw_data/SRR28821350_2.fastq.gz | sed -n '4~4p' | awk '{for(i=1;i<=length($0);i++){q[substr($0,i,1)]++}} END{for (k in q) print k,q[k]}'`  
  - Same as above, but for R2 FASTQ.  
</details>

### 5.   Checking FASTQ Pairing 

We ensured all our FASTQ files are correctly paired before running any bioinformatics analysis.

##### Step 1: Create the script
```bash
nano check_fastq_pairs.sh
```
##### Step 2: Paste the following into `check_fastq_pairs.sh`
```bash
#!/bin/bash
set -euo pipefail

INDIR="raw_data"
[[ "$(basename "$PWD")" != "raw_data" ]] && cd "$INDIR" || { echo "❌ raw_data directory not found"; exit 1; }

echo "🔍 Checking FASTQ pairings in $PWD ..."

MISSING=false
PAIRED_COUNT=0
TOTAL_COUNT=0

for R1 in *_1.fastq.gz *_R1.fastq.gz *_R1_*.fastq.gz *_001.fastq.gz; do
    [[ -f "$R1" ]] || continue
    TOTAL_COUNT=$((TOTAL_COUNT+1))
    SAMPLE=${R1%_1.fastq.gz}; SAMPLE=${SAMPLE%_R1.fastq.gz}; SAMPLE=${SAMPLE%_R1_*.fastq.gz}; SAMPLE=${SAMPLE%_001.fastq.gz}; SAMPLE=${SAMPLE%_R1_001.fastq.gz}

    if [[ -f "${SAMPLE}_2.fastq.gz" || -f "${SAMPLE}_R2.fastq.gz" || -f "${SAMPLE}_R2_*.fastq.gz" || -f "${SAMPLE}_002.fastq.gz" ]]; then
        echo "✅ $SAMPLE — paired"
        PAIRED_COUNT=$((PAIRED_COUNT+1))
    else
        echo "❌ $SAMPLE — missing R2 file"
        MISSING=true
    fi
done

echo -e "\nTotal samples checked: $TOTAL_COUNT"
echo "Correctly paired samples: $PAIRED_COUNT"
$MISSING && echo "⚠ Some samples are missing pairs. Fix before running fastp." || echo "✅ All FASTQ files are correctly paired."
```
<details>
<summary>🔗 FASTQ Pairing Check Script Explanation</summary>

- `#!/bin/bash` → Runs the script using Bash.  
- `set -euo pipefail` → Exits on errors, unset variables, or failed commands.  
- `INDIR="raw_data"` → Directory containing the FASTQ files.  
- `[[ "$(basename "$PWD")" != "raw_data" ]] && cd "$INDIR" ...` → Changes to `raw_data` directory if not already there.  
- `MISSING=false` → Flag to track if any R2 files are missing.  
- `PAIRED_COUNT=0` / `TOTAL_COUNT=0` → Counters for paired samples and total samples checked.  
- `for R1 in *_1.fastq.gz *_R1.fastq.gz *_R1_*.fastq.gz *_001.fastq.gz; do ... done` → Loops over all R1 FASTQ files.  
- `[[ -f "$R1" ]] || continue` → Skips if the file does not exist.  
- `SAMPLE=...` → Strips common R1 suffixes to extract the sample name.  
- `if [[ -f "${SAMPLE}_2.fastq.gz" || ... ]]; then ... fi` → Checks if a corresponding R2 file exists.  
- `echo "✅ $SAMPLE — paired"` → Prints a message if the pair is found.  
- `echo "❌ $SAMPLE — missing R2 file"` → Prints a message if the pair is missing and sets `MISSING=true`.  
- `TOTAL_COUNT` and `PAIRED_COUNT` → Track the total and successfully paired samples.  
- Final messages:  
  - `⚠ Some samples are missing pairs` → Warns user if any R2 files are missing.  
  - `✅ All FASTQ files are correctly paired` → Confirms all samples are paired.  

</details>

##### Step 4: Make the script executable
```bash
chmod +x check_fastq_pairs.sh
```
##### Step 5: Run the script
```bash
./check_fastq_pairs.sh
```
> **Tip:** Ensure all R1/R2 naming conventions in your directory match the patterns used in the script.  
> You can adjust the patterns (`*_1.fastq.gz`, `*_R1.fastq.gz`, etc.) if needed.

Calculating Minimum, Maximum, and Average Read Lengths for Paired-End Reads

### 6.   Checking raw FASTQ Read Length Summary

<details>
<summary>📏 Read Length Summary – Importance and Benefits</summary>

Before performing any downstream bioinformatics analysis, it is important to understand the quality and characteristics of your sequencing data. One key metric is the **read length** of FASTQ files.  

## 🔹 Why read length matters  
- **Minimum read length:** Identifies very short reads that may result from sequencing errors or trimming. Extremely short reads can cause mapping errors or low-quality variant calls.  
- **Maximum read length:** Confirms whether reads were sequenced to the expected length and detects unusually long reads that may indicate adapter contamination or sequencing artifacts.  
- **Average read length:** Provides an overall measure of sequencing quality and consistency across the dataset.  

## 🔹 Importance in paired-end sequencing  
Calculating these metrics for **both R1 and R2 reads** is crucial:  
- Ensures both reads in a pair are of comparable lengths → essential for accurate alignment and variant calling.  
- Detects discrepancies between forward and reverse reads that may indicate technical issues during sequencing or library preparation.  
- Allows early filtering of problematic samples before computationally intensive steps such as mapping, variant calling, or assembly.  

## 🔹 Benefits of summarizing into a CSV  
By compiling read lengths into a **CSV file**, you can:  
- Quickly inspect and compare samples.  
- Identify outliers or problematic datasets.  
- Make informed decisions on trimming, filtering, or quality control.  
- Improve reliability and reproducibility of downstream analyses.  

</details>

---

##### Step 1:  Open nano to create a new script
```bash
nano fastq_read_length_summary.sh
```
##### Step 2: Paste the following code into nano

```bash
#!/bin/bash
set -euo pipefail

FASTQ_DIR="raw_data"
OUTDIR="csv_output"
OUTPUT_CSV="${OUTDIR}/read_length_summary.csv"

mkdir -p "$OUTDIR"

echo "Sample,R1_min,R1_max,R1_avg,R2_min,R2_max,R2_avg" > "$OUTPUT_CSV"

for R1 in "$FASTQ_DIR"/*_1.fastq.gz "$FASTQ_DIR"/*_R1.fastq.gz; do
    [[ -f "$R1" ]] || continue
    SAMPLE=$(basename "$R1" | sed -E 's/_R?1.*\.fastq\.gz//')
    R2="$FASTQ_DIR/${SAMPLE}_2.fastq.gz"
    [[ -f "$R2" ]] || R2="$FASTQ_DIR/${SAMPLE}_R2.fastq.gz"

    if [[ -f "$R2" ]]; then
        echo "Processing sample $SAMPLE"

        calc_stats() {
            zcat "$1" | awk 'NR%4==2 {len=length($0); sum+=len; if(min==""){min=len}; if(len<min){min=len}; if(len>max){max=len}; count++} END{avg=sum/count; printf "%d,%d,%.2f", min, max, avg}'
        }

        STATS_R1=$(calc_stats "$R1")
        STATS_R2=$(calc_stats "$R2")

        echo "$SAMPLE,$STATS_R1,$STATS_R2" >> "$OUTPUT_CSV"
    else
        echo "⚠ Missing R2 for $SAMPLE, skipping."
    fi
done

echo "✅ Read length summary saved to $OUTPUT_CSV"

```
<details>
<summary>📊 Read Length Summary Script Explanation</summary>

- `FASTQ_DIR="raw_data"` → Directory containing FASTQ files.  
- `OUTDIR="csv_output"` → Directory to save CSV output; created automatically if missing.  
- `OUTPUT_CSV="${OUTDIR}/read_length_summary.csv"` → Output CSV file path.  
- `mkdir -p "$OUTDIR"` → Ensure output directory exists.  
- `echo "Sample,R1_min,R1_max,R1_avg,R2_min,R2_max,R2_avg" > "$OUTPUT_CSV"` → CSV header.  
- `for R1 in "$FASTQ_DIR"/*_1.fastq.gz "$FASTQ_DIR"/*_R1.fastq.gz; do ...` → Loop over all R1 FASTQ files.  
- `[[ -f "$R1" ]] || continue` → Skip if R1 file does not exist.  
- `SAMPLE=$(basename "$R1" | sed -E 's/_R?1.*\.fastq\.gz//')` → Extract sample name.  
- `R2="$FASTQ_DIR/${SAMPLE}_2.fastq.gz"` → Get paired R2 filename; also checks `_R2` naming.  
- `if [[ -f "$R2" ]]; then ... else ... fi` → Skip sample if R2 is missing.  
- `calc_stats() { ... }` → Function to calculate min, max, avg read lengths for a FASTQ file.  
- `STATS_R1=$(calc_stats "$R1")` → Stats for R1.  
- `STATS_R2=$(calc_stats "$R2")` → Stats for R2.  
- `echo "$SAMPLE,$STATS_R1,$STATS_R2" >> "$OUTPUT_CSV"` → Append sample stats to CSV.  
- `echo "⚠ Missing R2 for $SAMPLE, skipping."` → Warning if R2 missing.  
- `echo "✅ Read length summary saved to $OUTPUT_CSV"` → Final confirmation message.

</details>

##### Step 3: Save and exit nano
Press Ctrl + O → Enter (to write the file)
Press Ctrl + X → Exit nano
##### Step 4: Make the script executable
```bash
chmod +x fastq_read_length_summary.sh
```
##### Step 5: Run the script
```bash
./fastq_read_length_summary.sh
```

# 2️⃣ Quality Control & Trimming with **FASTP**

<details>
<summary>⚡ FASTQ Preprocessing with fastp</summary>

[`fastp`](https://github.com/OpenGene/fastp) is a widely used **FASTQ preprocessor** for quality control (QC), read trimming, and adapter removal. It is efficient, multithreaded, and provides both **JSON and HTML reports** for each sample.  

In this pipeline, `fastp` ensures that only **high-quality reads** are retained before mapping and variant calling. High-quality read preprocessing is crucial in *Mycobacterium tuberculosis* (TB) WGS analysis because poorly trimmed or unfiltered reads can lead to:
- False-positive SNP calls  
- Mapping errors, especially in repetitive regions (e.g., PE/PPE genes)  
- Biased coverage, affecting downstream variant interpretation  

### Advantages of fastp over Trimmomatic and other tools
- **All-in-one solution**: Handles trimming, adapter detection, filtering, and QC in a single step.  
- **Automatic adapter detection**: Reduces human error, especially for large TB projects with mixed sequencing batches.  
- **Speed and multithreading**: Written in C++, much faster than Java-based Trimmomatic.  
- **Comprehensive QC output**: HTML (interactive) and JSON (machine-readable) reports for quality distribution, duplication rates, adapter content, and polyG/polyX tails.  
- **Better polyG/polyX handling**: Important for Illumina NovaSeq/NextSeq data.  
- **UMI support**: Useful for advanced TB sequencing protocols using Unique Molecular Identifiers.  
- **Minimal parameter tuning**: Default settings are optimized, reducing manual adjustments.  

### Why this matters for TB analysis
- **Accurate SNP calling**: TB drug-resistance prediction relies on high-confidence SNPs.  
- **Low genetic diversity detection**: TB isolates often differ by only a few SNPs. Filtering errors while retaining true variants is critical.  
- **Scalability for large cohorts**: Efficient, reproducible preprocessing is essential for thousands of public TB isolates from ENA/NCBI.  

</details>

---

### Steps to Run FASTP
##### Step 1: **Open nano to create the script**
```bash
nano run_fastp.sh
```
##### Step 2: Paste the following code into nano
```bash
#!/bin/bash
set -euo pipefail

INDIR="raw_data"
OUTDIR="fastp_results_min_50"
mkdir -p "$OUTDIR"

SAMPLES=()

for R1 in "$INDIR"/*_1.fastq.gz "$INDIR"/*_R1.fastq.gz "$INDIR"/*_001.fastq.gz "$INDIR"/*_R1_001.fastq.gz; do
    [[ -f "$R1" ]] || continue

    SAMPLE=$(basename "$R1")
    SAMPLE=${SAMPLE%%_1.fastq.gz}
    SAMPLE=${SAMPLE%%_R1.fastq.gz}
    SAMPLE=${SAMPLE%%_001.fastq.gz}
    SAMPLE=${SAMPLE%%_R1_001.fastq.gz}

    if   [[ -f "$INDIR/${SAMPLE}_2.fastq.gz" ]]; then R2="$INDIR/${SAMPLE}_2.fastq.gz"
    elif [[ -f "$INDIR/${SAMPLE}_R2.fastq.gz" ]]; then R2="$INDIR/${SAMPLE}_R2.fastq.gz"
    elif [[ -f "$INDIR/${SAMPLE}_002.fastq.gz" ]]; then R2="$INDIR/${SAMPLE}_002.fastq.gz"
    elif [[ -f "$INDIR/${SAMPLE}_R2_001.fastq.gz" ]]; then R2="$INDIR/${SAMPLE}_R2_001.fastq.gz"
    else
        echo "⚠ No R2 file found for $SAMPLE — skipping."
        continue
    fi

    if [[ -f "$OUTDIR/${SAMPLE}_1.trim.fastq.gz" && -f "$OUTDIR/${SAMPLE}_2.trim.fastq.gz" ]]; then
        echo "⏩ Skipping $SAMPLE (already processed)."
        continue
    fi

    SAMPLES+=("$SAMPLE,$R1,$R2")
done

if [[ ${#SAMPLES[@]} -eq 0 ]]; then
    echo "❌ No paired FASTQ files found in $INDIR"
    exit 1
fi

THREADS=$(nproc)
FASTP_THREADS=$(( THREADS / 2 ))

run_fastp() {
    SAMPLE=$1
    R1=$2
    R2=$3
    echo "✅ Processing sample: $SAMPLE"
    fastp \
        -i "$R1" \
        -I "$R2" \
        -o "$OUTDIR/${SAMPLE}_1.trim.fastq.gz" \
        -O "$OUTDIR/${SAMPLE}_2.trim.fastq.gz" \
        -h "$OUTDIR/${SAMPLE}_fastp.html" \
        -j "$OUTDIR/${SAMPLE}_fastp.json" \
        --length_required 50 \
        --qualified_quality_phred 20 \
        --detect_adapter_for_pe \
        --thread $FASTP_THREADS \
        &> "$OUTDIR/${SAMPLE}_fastp.log"
}

export -f run_fastp
export OUTDIR FASTP_THREADS

printf "%s\n" "${SAMPLES[@]}" | parallel -j 3 --colsep ',' run_fastp {1} {2} {3}

echo "🎉 Completed fastp for $(ls "$OUTDIR"/*_fastp.json | wc -l) samples."

```
<details>
<summary>🧹 fastp Trimming Script Explanation</summary>

- `#!/bin/bash` → Run script with Bash.  
- `set -euo pipefail` → Exit on errors, undefined variables, or failed pipelines.  
- `INDIR="raw_data"` → Raw FASTQ directory.  
- `OUTDIR="fastp_results_min_50"` → Directory for trimmed FASTQs.  
- `mkdir -p "$OUTDIR"` → Create output directory.  
- `SAMPLES=()` → Initialize array to store sample info.  
- `for R1 in ...` → Loop over R1 files with common naming patterns.  
- `SAMPLE=...` → Extract sample name from R1 filename.  
- `if ... elif ... else` → Detect corresponding R2 under multiple naming conventions.  
- `if [[ -f "$OUTDIR/${SAMPLE}_1.trim.fastq.gz" && ... ]]` → Skip already processed samples.  
- `SAMPLES+=("$SAMPLE,$R1,$R2")` → Store sample info for parallel execution.  
- `THREADS=$(nproc)` → Detect total CPU cores.  
- `FASTP_THREADS=$(( THREADS / 2 ))` → Allocate threads per fastp process.  
- `run_fastp() { ... }` → Function to run fastp per sample:  
  - `-i / -I` → Input R1/R2  
  - `-o / -O` → Output trimmed FASTQs  
  - `-h / -j` → HTML and JSON reports  
  - `--length_required 50` → Minimum read length  
  - `--qualified_quality_phred 20` → Quality threshold  
  - `--detect_adapter_for_pe` → Auto adapter trimming  
  - `--thread` → Threads for fastp  
- `export -f run_fastp` → Make function available to GNU Parallel.  
- `printf "%s\n" "${SAMPLES[@]}" | parallel -j 3 --colsep ',' run_fastp {1} {2} {3}` → Run 3 fastp jobs in parallel.  
- `echo "🎉 Completed fastp ..."` → Display completion message.

</details>


##### Step 3: Save & exit nano
Press CTRL+O, Enter (save)
Press CTRL+X (exit)
##### Step 4: Make the script executable
```bash
chmod +x run_fastp.sh
```
##### Step 5: Activate your conda env and run
```bash
conda activate fastp_env
./run_fastp.sh
```

Count R1 trimmed files
```bash
ls -lth fastp_results_min_50/*_1.trim.fastq.gz | wc -l
```
Count R2 trimmed files
```bash
ls -lth fastp_results_min_50/*_2.trim.fastq.gz | wc -l
```

View first 10 quality lines in trimmed FASTQ

```bash
 Show first 10 quality lines from R1
```bash
echo "🔹 First 10 quality lines from ET3_S55_1 (R1):"
zcat fastp_results_min_50/ET3_S55_1.trim.fastq.gz \
| sed -n '4~4p' \
| head -n 10 \
| awk '{print "✅ " $0}'
```
 Show first 10 quality lines from R2
```bash
echo "🔹 First 10 quality lines from ET3_S55_2 (R2):"
zcat fastp_results_min_50/ET3_S55_2.trim.fastq.gz \
| sed -n '4~4p' \
| head -n 10 \
| awk '{print "✅ " $0}'

```
<details>
  <summary>🔍 How it works</summary>

- `zcat` → Decompresses the trimmed FASTQ.  
- `sed -n '4~4p'` → Prints every 4th line starting from line 4 (the quality score line of each read).  
- `head -n 10` → Shows the first 10 quality lines for quick inspection.  

</details>

Count ASCII characters in quality lines:
```bash
    Count base composition in R1
```bash
zcat fastp_results_min_50/ET3_S55_1.trim.fastq.gz \
| sed -n '4~4p' \
| awk '{
    for(i=1;i<=length($0);i++){ q[substr($0,i,1)]++ }
} END {
    for (k in q) print k, q[k]
}' \
| awk '{print "✅ Base " $1 ": " $2 " occurrences"}'
```
   Count base composition in R2
```bash
zcat fastp_results_min_50/ET3_S55_2.trim.fastq.gz \
| sed -n '4~4p' \
| awk '{
    for(i=1;i<=length($0);i++){ q[substr($0,i,1)]++ }
} END {
    for (k in q) print k, q[k]
}' \
| awk '{print "✅ Base " $1 ": " $2 " occurrences"}'
```
<details>
  <summary>🔢 How it works</summary>

- `zcat` → Decompresses trimmed FASTQ.  
- `sed -n '4~4p'` → Selects every 4th line (quality score line).  
- `awk '{for(i=1;i<=length($0);i++){q[substr($0,i,1)]++}} END{for (k in q) print k,q[k]}'` → Counts occurrences of each ASCII character in the quality scores.  

This helps quickly identify if the quality encoding is correct (usually Phred+33 for Illumina) and whether trimming improved overall quality.

</details>

For batch processing trimmed FASTQ
##### Step 1: Open a new script
```bash
nano count_trimmed_reads.sh
```
##### Step 2: Paste the following code
This script counts reads in trimmed paired-end FASTQ files and saves results to a CSV.
```bash
#!/bin/bash
set -euo pipefail

INDIR="fastp_results_min_50"
OUTDIR="csv_output"
mkdir -p "$OUTDIR"

OUTFILE="$OUTDIR/trimmed_read_counts.csv"
echo "Sample,R1_reads,R2_reads" > "$OUTFILE"
echo "📊 Counting reads in trimmed FASTQ files from '$INDIR'..."

for R1 in "$INDIR"/*_1.trim.fastq.gz "$INDIR"/*_R1.trim.fastq.gz; do
    [[ -f "$R1" ]] || continue
    SAMPLE=$(basename "$R1" | sed -E 's/_R?1.*\.trim\.fastq\.gz//')
    R2=""
    for suffix in "_2.trim.fastq.gz" "_R2.trim.fastq.gz" "_R2_*.trim.fastq.gz"; do
        [[ -f "$INDIR/${SAMPLE}${suffix}" ]] && R2="$INDIR/${SAMPLE}${suffix}" && break
    done
    R1_COUNT=$(( $(zcat "$R1" | wc -l) / 4 ))
    R2_COUNT=$([[ -n "$R2" ]] && echo $(( $(zcat "$R2" | wc -l) / 4 )) || echo "NA")
    echo "$SAMPLE,$R1_COUNT,$R2_COUNT" >> "$OUTFILE"
    echo "✅ $SAMPLE → R1: $R1_COUNT | R2: $R2_COUNT"
done

echo "🎉 All done! Read counts saved to '$OUTFILE'"


```
<details>
  <summary>📊 Trimmed FASTQ Read Count Script Explanation</summary>

- `#!/bin/bash` → Runs the script using Bash.  
- `set -euo pipefail` → Exits on errors, unset variables, or failed commands.  
- `INDIR="fastp_results_min_50"` → Directory containing trimmed FASTQ files.  
- `OUTDIR="csv_output"` → Directory to save the output CSV file.  
- `OUTFILE="$OUTDIR/trimmed_read_counts.csv"` → Path of the output CSV file.  
- `mkdir -p "$OUTDIR"` → Creates the output directory if it doesn’t exist.  
- `echo "Sample,R1_reads,R2_reads" > "$OUTFILE"` → Writes the CSV header.  
- `echo "📊 Counting reads in trimmed FASTQ files from '$INDIR'..."` → Prints starting message.  
- `for R1 in "$INDIR"/*_1.trim.fastq.gz "$INDIR"/*_R1.trim.fastq.gz; do ... done` → Iterates over all R1 trimmed FASTQ files.  
- `[[ -f "$R1" ]] || continue` → Skips if the R1 file does not exist.  
- `SAMPLE=$(basename "$R1" | sed -E 's/_R?1.*\.trim\.fastq\.gz//')` → Extracts the sample name from the filename.  
- `for suffix in "_2.trim.fastq.gz" "_R2.trim.fastq.gz" "_R2_*.trim.fastq.gz"; do ... done` → Searches for the corresponding R2 file.  
- `R1_COUNT=$(( $(zcat "$R1" | wc -l) / 4 ))` → Counts number of reads in R1 (lines divided by 4).  
- `R2_COUNT=$([[ -n "$R2" ]] && echo $(( $(zcat "$R2" | wc -l) / 4 )) || echo "NA")` → Counts number of reads in R2 if present; outputs "NA" for single-end samples.  
- `echo "$SAMPLE,$R1_COUNT,$R2_COUNT" >> "$OUTFILE"` → Appends sample read counts to the CSV file.  
- `echo "✅ $SAMPLE → R1: $R1_COUNT | R2: $R2_COUNT"` → Prints per-sample progress.  
- `echo "🎉 All done! Read counts saved to '$OUTFILE'"` → Prints final completion message.  

</details>


##### Step 3: Save and exit nano
Press Ctrl + O → Enter (to write the file)
Press Ctrl + X → Exit nano

##### Step 4: Make the script executable
```bash
chmod +x count_trimmed_reads.sh
```
##### Step 5: Run the script
```bash
./count_trimmed_reads.sh
```
Using zcat and wc (line counts)
Each read in FASTQ has 4 lines, so dividing by 4 gives read count. 
```bash
RAW_R1_COUNT=$(( $(zcat raw_data/ET3_S55_1.fastq.gz | wc -l) / 4 ))
RAW_R2_COUNT=$(( $(zcat raw_data/ET3_S55_2.fastq.gz | wc -l) / 4 ))

TRIM_R1_COUNT=$(( $(zcat fastp_results_min_50/ET3_S55_1.trim.fastq.gz | wc -l) / 4 ))
TRIM_R2_COUNT=$(( $(zcat fastp_results_min_50/ET3_S55_2.trim.fastq.gz | wc -l) / 4 ))

echo "R1 trimmed reads: $(( RAW_R1_COUNT - TRIM_R1_COUNT ))"
echo "R2 trimmed reads: $(( RAW_R2_COUNT - TRIM_R2_COUNT ))"
```
read length summary on trimmed FASTQ files
##### Step 1: Open nano to create a new script
```bash
nano trimmed_fastq_read_length_summary.sh
```
##### Step 2: Paste the following code into nano
```bash
#!/bin/bash
set -euo pipefail

FASTQ_DIR="fastp_results_min_50"
OUTDIR="csv_output"
OUTPUT_CSV="${OUTDIR}/trimmed_read_length_summary.csv"

mkdir -p "$OUTDIR"

echo "Sample,R1_min,R1_max,R1_avg,R2_min,R2_max,R2_avg" > "$OUTPUT_CSV"

for R1 in "$FASTQ_DIR"/*_1.trim.fastq.gz; do
    SAMPLE=$(basename "$R1" _1.trim.fastq.gz)
    R2="${FASTQ_DIR}/${SAMPLE}_2.trim.fastq.gz"

    if [[ -f "$R2" ]]; then
        echo "Processing sample $SAMPLE"

        calc_stats() {
            zcat "$1" | awk 'NR%4==2 {
                len=length($0)
                sum+=len
                if(min==""){min=len}
                if(len<min){min=len}
                if(len>max){max=len}
                count++
            } END {
                avg=sum/count
                printf "%d,%d,%.2f", min, max, avg
            }'
        }

        STATS_R1=$(calc_stats "$R1")
        STATS_R2=$(calc_stats "$R2")

        echo "$SAMPLE,$STATS_R1,$STATS_R2" >> "$OUTPUT_CSV"
    else
        echo "⚠ Missing R2 for $SAMPLE, skipping."
    fi
done

echo "✅ Trimmed read length summary saved to $OUTPUT_CSV"

```
<details>
  <summary>📊 Trimmed Read Length Summary Script Explanation</summary>

- `#!/bin/bash` → Run script with Bash.  
- `FASTQ_DIR="fastp_results_min_50"` → Directory containing trimmed FASTQ files.  
- `OUTDIR="csv_output"` → Directory to save CSV output.  
- `OUTPUT_CSV="${OUTDIR}/trimmed_read_length_summary.csv"` → Output CSV file path.  
- `mkdir -p "$OUTDIR"` → Create output directory if missing.  
- `echo "Sample,R1_min,R1_max,R1_avg,R2_min,R2_max,R2_avg" > "$OUTPUT_CSV"` → Write CSV header.  
- `for R1 in "$FASTQ_DIR"/*_1.trim.fastq.gz; do ...` → Loop over all R1 trimmed FASTQ files.  
- `SAMPLE=$(basename "$R1" _1.trim.fastq.gz)` → Extract sample name.  
- `R2="${FASTQ_DIR}/${SAMPLE}_2.trim.fastq.gz"` → Find corresponding R2 file.  
- `if [[ -f "$R2" ]]; then ... else ... fi` → Skip sample if R2 is missing.  
- `calc_stats() { ... }` → Function to calculate min, max, and average read lengths.  
- `STATS_R1=$(calc_stats "$R1")` → Compute stats for R1.  
- `STATS_R2=$(calc_stats "$R2")` → Compute stats for R2.  
- `echo "$SAMPLE,$STATS_R1,$STATS_R2" >> "$OUTPUT_CSV"` → Append results to CSV.  
- `echo "⚠ Missing R2 for $SAMPLE, skipping."` → Print warning if R2 not found.  
- `echo "✅ Trimmed read length summary saved to $OUTPUT_CSV"` → Final confirmation message.  

</details>

##### Step 3: Save and exit nano
Press Ctrl + O → Enter
Press Ctrl + X → Exit
##### Step 4: Make the script executable
```
chmod +x trimmed_fastq_read_length_summary.sh
```
##### Step 5: Run the script
```bash
./trimmed_fastq_read_length_summary.sh
```
# 3️⃣ Aggregate QC Reports with **MultiQC**
[`MultiQC`](https://github.com/MultiQC/MultiQC) 

<details>
<summary>📊 Aggregating QC with MultiQC</summary>

After preprocessing with `fastp`, we often generate dozens or hundreds of per-sample QC reports (`.html` and `.json`). Instead of checking each report manually, **MultiQC** aggregates all results into a single interactive HTML report.  

### Why we use MultiQC in TB analysis
- **Aggregated QC overview**: Summarizes all `fastp` results in one place.  
- **Consistency check**: Quickly detects outlier samples (e.g., unusually short reads, poor quality, or failed trimming).  
- **Scalable**: Handles hundreds or thousands of TB isolates efficiently.  
- **Standardized reporting**: Facilitates sharing results across teams or for publications.  

</details>

---

### Script: Run MultiQC
##### Step 1: **Open nano to create the script `run_multiqc.sh`
```bash
nano run_multiqc.sh
```
##### Step 2: Paste the following code into nano
```bash
#!/bin/bash
set -euo pipefail

INPUT_DIR="fastp_results_min_50"
OUTPUT_DIR="multiqc/fastp_multiqc"

mkdir -p "$OUTPUT_DIR"

if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Input directory '$INPUT_DIR' does not exist!"
    exit 1
fi

multiqc "$INPUT_DIR" -o "$OUTPUT_DIR"

echo "MultiQC report generated in '$OUTPUT_DIR'."

```
<details>
<summary>📊 MultiQC Script Explanation</summary>

- `#!/bin/bash` → Run script with Bash.  
- `set -euo pipefail` → Exit on errors, undefined variables, or failed pipelines.  
- `INPUT_DIR="fastp_results_min_50"` → Directory containing fastp JSON/HTML outputs.  
- `OUTPUT_DIR="multiqc/fastp_multiqc"` → Directory where the aggregated MultiQC report will be saved.  
- `mkdir -p "$OUTPUT_DIR"` → Create `multiqc` and `fastp_multiqc` directories if they don’t exist.  
- `if [ ! -d "$INPUT_DIR" ]; then ... fi` → Check that the input directory exists; exit with an error if not.  
- `multiqc "$INPUT_DIR" -o "$OUTPUT_DIR"` → Run MultiQC on all files in `INPUT_DIR` and save the combined report in `OUTPUT_DIR`.  

</details>

##### Step 3: Save & exit nano
Press CTRL+O, Enter (save)
Press CTRL+X (exit)
##### Step 4: Make the script executable
```bash
chmod +x run_multiqc.sh
```
##### Step 5: Activate your conda env and run
```bash
conda activate multiqc_env
./run_multiqc.sh
```
# 4️⃣ Drug Resistance Screening with **TB-Profiler**
Before running [`TB-Profiler`](https://github.com/jodyphelan/TBProfiler), we perform quality checks on raw FASTQ files and exclude samples with extremely low-quality reads. However, since TB-Profiler internally uses Trimmomatic to trim adapters and low-quality bases, it is not necessary to pre-trim the reads. Therefore, only quality-checked FASTQ files are provided as input to TB-Profiler, allowing it to handle trimming and variant calling internally.

<details>
<summary>🧬 TB-Profiler: Variant Calling, Lineage, and Drug Resistance</summary>

**TB-Profiler** is a specialized tool for *Mycobacterium tuberculosis* whole-genome sequencing (WGS) data. It performs **variant calling, lineage determination, and drug resistance prediction** in a single pipeline.

### Why we use TB-Profiler
- 🧪 **Drug Resistance Prediction** → Detects known resistance mutations for first- and second-line TB drugs.  
- 🌍 **Lineage Typing** → Classifies isolates into recognized TB lineages (e.g., Lineage 1–7).  
- 🔄 **Flexible Input** → Accepts FASTQ, BAM, or VCF files.  
- 📊 **Clear Outputs** → Produces human-readable (`.txt`) and machine-readable (`.json`) reports.  
- ⚡ **Speed & Integration** → Efficient and easily incorporated into TB genomics pipelines.  

### Why it matters
- Provides **actionable insights** for public health, including drug resistance and outbreak tracking.  
- Avoids manual cross-referencing with multiple TB resistance databases.  
- Ensures **standardized results** comparable across studies.  

</details>
 
---

## Steps
##### Step 1: Create or edit the script
```bash
nano run_tbprofiler.sh
```
##### Step 2: Paste the following code
```bash
#!/bin/bash
set -euo pipefail

FASTQ_DIR="raw_data"

echo "📊 Starting TBProfiler runs for all samples in $FASTQ_DIR ..."

for R1 in "$FASTQ_DIR"/*_1.fastq.gz; do
    SAMPLE=$(basename "$R1" _1.fastq.gz)
    R2="$FASTQ_DIR/${SAMPLE}_2.fastq.gz"

    if [[ ! -f "$R2" ]]; then
        echo "❌ Warning: missing paired file for $SAMPLE, skipping."
        continue
    fi

    echo "▶️ Processing sample: $SAMPLE"

    tb-profiler profile \
        -1 "$R1" \
        -2 "$R2" \
        --threads 8 \
        --prefix "$SAMPLE" \
        --txt \
        --spoligotype

    echo "✅ Finished $SAMPLE"
done

echo "📌 All samples processed!"
```
<details>
<summary>🧪 TB-Profiler Script Explanation</summary>

- `#!/bin/bash` → Run script with Bash.  
- `set -euo pipefail` → Exit on errors or undefined variables.  
- `FASTQ_DIR="raw_data"` → Folder containing paired-end FASTQ files.  
- `for R1 in "$FASTQ_DIR"/*_1.fastq.gz; do ... done` → Loop through all R1 files.  
- `SAMPLE=$(basename "$R1" _1.fastq.gz)` → Extract sample name from filename.  
- `R2="$FASTQ_DIR/${SAMPLE}_2.fastq.gz"` → Construct path for paired R2 file.  
- `if [[ ! -f "$R2" ]]; then ... fi` → Skip sample if paired R2 file is missing.  
- `tb-profiler profile -1 "$R1" -2 "$R2" --threads 8` → Run TBProfiler on paired reads using 8 threads; outputs are saved automatically in `results/tbprofiler_results/$SAMPLE/`.  
- `echo "✅ Finished $SAMPLE"` → Completion message per sample.  
- `echo "📌 All samples processed!"` → Final message after all samples are run.

</details>

##### Step 3: Save and exit nano
Press Ctrl + O → Enter (to write the file)
Press Ctrl + X → Exit nano

##### Step 4: Make the script executable
```bash
chmod +x run_tbprofiler.sh
```
##### Step 5: Activate environment and run
```bash
conda activate tbprofiler_env
./run_tbprofiler.sh
```
##### Step 6: move the output from tbprofiler into new directory 
make directory tbprofiler_results
```bash
mkdir -p tbprofiler_results
```
then move the tbprofiler output `bam`, `vcf` and `results` directories from current directory into tbprofiler_results directory 
```bash
mv bam vcf results tbprofiler_results/
echo "📌 All TBProfiler outputs moved to tbprofiler_results/"
```
##### Step 6: change directory to tbprofiler_results
```bash
cd ./tbprofiler_results
```
##### Step 7: Collate results
The results from numerous runs can be collated into one table using the following command:
```bash
tb-profiler collate
```
##### Step 8: Generate iTOL config for phylogeny visualization
If we want to visualize your phylogeny alongside drug-resistance types, lineages, and sublineages, run:
```bash
tb-profiler collate --itol
```
##### Step 8: Combine all tbprofiler.txt into one CSV
```bash
find . -name "tbprofiler.txt" \
| exec awk 'FNR==1 && NR!=1{next} {print}' {} + \
| sed 's/\t/,/g' \
> tbprofiler_collated.csv

echo "✅ Collated all tbprofiler.txt files into tbprofiler_collated.csv"

```
<details>
<summary>📌 Explanation of the TBProfiler collate command</summary>

- `find . -name "tbprofiler.txt"`  
  Searches the current directory (`.`) and all subdirectories for files named `tbprofiler.txt`.

- `-exec awk 'FNR==1 && NR!=1{next} {print}' {} +`  
  For each `tbprofiler.txt` found:  
  - `FNR==1 && NR!=1{next}` → skips the header line of every file except the first one.  
  - `{print}` → prints all other lines (the data).  
  - `NR` = total number of lines processed so far; `FNR` = line number in current file.

- `| sed 's/\t/,/g'`  
  Replaces tabs (`\t`) with commas, converting the tab-delimited text into CSV format.

- `> tbprofiler_collated.csv`  
  Redirects the final output into a file called `tbprofiler_collated.csv`.

</details>

# 5️⃣ Variant Calling with **Snippy**

[`Snippy`](https://github.com/tseemann/snippy) 
<details>
<summary>🧬 Detailed Overview: Variant Calling with Snippy</summary>

`Snippy` is a rapid and reproducible pipeline designed for bacterial genome variant calling. It is particularly well-suited for **Mycobacterium tuberculosis (TB) WGS analysis** due to its efficiency and accuracy.

### Key Features
- **Reference-based mapping**: Maps raw sequencing reads directly to a reference genome (commonly *M. tuberculosis* H37Rv), ensuring accurate alignment even in repetitive regions.  
- **High-confidence SNP calling**: Detects single nucleotide polymorphisms (SNPs) with stringent filtering criteria to minimize false positives.  
- **Consensus sequence generation**: Produces high-quality consensus FASTA files per sample, which can be used for phylogenetic analysis or comparative genomics.  
- **Reproducibility**: Standardized workflow ensures consistent results across multiple samples and datasets.  
- **Lightweight and scalable**: Efficient for processing large cohorts of TB isolates without requiring extensive computational resources.  
- **Standardized output files**: Outputs include VCF files for SNPs, consensus FASTA sequences, and optional alignment summaries, facilitating downstream analyses such as:
  - Phylogenetic tree reconstruction  
  - Drug resistance prediction using TBProfiler  
  - Comparative genomics across multiple TB strains  
- **Easy integration**: Works seamlessly with other bioinformatics tools and pipelines, allowing automated WGS analysis workflows.  

### Importance in TB Analysis
- Ensures accurate variant detection for **drug resistance prediction**.  
- Allows monitoring of **microevolution** within TB outbreaks.  
- Provides **reliable consensus sequences** for large-scale phylogenetic studies.  

</details>

---

##### Step 1: Create the script
```bash
nano run_snippy.sh
```
##### Step 2: Paste the following into `run_snippy.sh`
```bash
#!/bin/bash
set -euo pipefail

REF="H37Rv.fasta"
FASTP_DIR="fastp_results_min_50"
OUTDIR="snippy_results"
THREADS=8
BWA_THREADS=30
JOBS=4

mkdir -p "$OUTDIR"

run_snippy_sample() {
    SAMPLE="$1"
    R1="${FASTP_DIR}/${SAMPLE}_1.trim.fastq.gz"
    R2="${FASTP_DIR}/${SAMPLE}_2.trim.fastq.gz"

    if [[ ! -f "$R1" || ! -f "$R2" ]]; then
        echo "⚠ Missing R1/R2 for $SAMPLE"
        return
    fi

    echo "Running Snippy on sample: $SAMPLE"
    TMP_DIR="${OUTDIR}/${SAMPLE}_tmp"
    mkdir -p "$TMP_DIR"

    snippy --cpus "$THREADS" --outdir "$TMP_DIR" --ref "$REF" \
           --R1 "$R1" --R2 "$R2" --force --bwaopt "-T $BWA_THREADS"

    if [[ -f "$TMP_DIR/snps.vcf" ]]; then
        mv "$TMP_DIR/snps.vcf" "${OUTDIR}/${SAMPLE}.vcf"
    fi

    for f in "$TMP_DIR"/*; do
        base=$(basename "$f")
        case "$base" in
            *.consensus.fa) mv "$f" "${OUTDIR}/${SAMPLE}.consensus.fa" ;;
            *.bam) mv "$f" "${OUTDIR}/${SAMPLE}.bam" ;;
            *.bam.bai) mv "$f" "${OUTDIR}/${SAMPLE}.bam.bai" ;;
            *.tab) mv "$f" "${OUTDIR}/${SAMPLE}.snps.tab" ;;
        esac
    done

    rm -rf "$TMP_DIR"

    [[ -f "${OUTDIR}/${SAMPLE}.vcf" ]] && echo "✅ Full VCF generated for $SAMPLE" || echo "⚠ No VCF produced for $SAMPLE"
}

export -f run_snippy_sample
export REF FASTP_DIR OUTDIR THREADS BWA_THREADS

ls "${FASTP_DIR}"/*_1.trim.fastq.gz \
    | sed 's|.*/||; s/_1\.trim\.fastq\.gz//' \
    | parallel -j "$JOBS" run_snippy_sample {}

ls "${FASTP_DIR}"/*_1.trim.fastq.gz \
    | sed 's|.*/||; s/_1\.trim\.fastq\.gz//' | sort > fastq_samples.txt
ls "${OUTDIR}"/*.vcf 2>/dev/null \
    | sed 's|.*/||; s/\.vcf//' | sort > snippy_samples.txt

echo "FASTQ pairs count: $(wc -l < fastq_samples.txt)"
echo "Snippy outputs count: $(wc -l < snippy_samples.txt)"

if diff fastq_samples.txt snippy_samples.txt >/dev/null; then
    echo "✅ All FASTQ pairs have corresponding Snippy results."
else
    echo "⚠ Missing samples detected:"
    diff fastq_samples.txt snippy_samples.txt || true
fi

rm -f fastq_samples.txt snippy_samples.txt

echo "🎯 All steps completed!"
echo "Snippy results are in: ${OUTDIR}/"

```
<details>
<summary>🌳 Snippy Pipeline Script Explanation</summary>

- `#!/bin/bash` → Run script with Bash.  
- `set -euo pipefail` → Exit on errors, undefined variables, or pipeline failures.  
- `REF="H37Rv.fasta"` → Reference genome.  
- `FASTP_DIR="fastp_results_min_50"` → Directory containing trimmed FASTQ files.  
- `OUTDIR="snippy_results"` → Directory to store Snippy outputs.  
- `THREADS=8` & `BWA_THREADS=30` → Threads for Snippy and BWA alignment.  
- `JOBS=4` → Number of samples to run in parallel.  
- `run_snippy_sample() { ... }` → Function for processing a single sample:  
  - Checks if paired FASTQ files exist.  
  - Runs Snippy with specified threads and BWA options.  
  - Moves key outputs (`.vcf`, `.consensus.fa`, `.bam`, `.bam.bai`, `.snps.tab`) to final directory.  
  - Deletes temporary Snippy directory.  
  - Prints confirmation if full VCF is generated.  
- `export -f run_snippy_sample` → Makes function available for GNU Parallel.  
- `ls ... | parallel -j "$JOBS" run_snippy_sample {}` → Runs multiple samples in parallel.  
- Verification section:  
  - Compares FASTQ sample list vs VCF output list.  
  - Prints warnings if any sample is missing.  
- `rm -f fastq_samples.txt snippy_samples.txt` → Cleans temporary lists.  
- `echo "🎯 All steps completed!"` → Final message indicating pipeline completion.  

</details>

##### Step 3: Save and exit nano
Press Ctrl + O, then Enter (save)
Press Ctrl + X (exit)

##### Step 4: Make the script executable
```bash
chmod +x run_snippy.sh
```
##### Step 5: Activate environment and run
```bash
conda activate snippy_env
./run_snippy.sh
```
Check Snippy VCFs Before Variant Filtering
- `tb_variant_filter` relies on a correctly formatted VCF with the `#CHROM` line to parse variants.  
- VCFs missing this line will **fail** during filtering, causing errors like:  
  `Missing line starting with "#CHROM"`.  
- Running this check before variant filtering saves time and ensures downstream analysis runs smoothly.

This script ensures that all VCF files generated by **Snippy** contain the required `#CHROM` header line.
```bash
#!/bin/bash
OUTDIR="snippy_results"

echo "Checking Snippy VCFs in $OUTDIR ..."

for vcf in "$OUTDIR"/*.vcf; do
    SAMPLE=$(basename "$vcf")
    if grep -q "^#CHROM" "$vcf"; then
        echo "✅ $SAMPLE contains #CHROM line"
    else
        echo "⚠ $SAMPLE is missing #CHROM line"
    fi
done
```

# 6️⃣ BAM Quality Check with **Qualimap**
[`Qualimap`](http://qualimap.conesalab.org/) 
<details>
<summary>📈 BAM Quality Assessment with Qualimap</summary>

After generating BAM files with Snippy, it is crucial to evaluate their quality before downstream analyses. **Qualimap** is widely used for this purpose.

### Why we use Qualimap in TB genomics
- **Mapping quality evaluation**: Assesses alignment accuracy and read placement.  
- **Coverage distribution and GC content**: Detects uneven coverage or biases across the genome.  
- **Visual reports**: Generates **HTML** and **PDF** reports with clear QC metrics.  
- **Identify problematic samples**: Flags samples with low coverage, poor alignment, or uneven depth, which may affect variant calling.  
- **Integration with MultiQC**: QC results can be aggregated for large TB cohorts for batch reporting.  

</details>

---

##### Step 1: Create the script
```bash
nano run_qualimap.sh
```
#####  Step 2: Paste the following into `run_qualimap.sh`
```bash
#!/bin/bash
set -euo pipefail

SNIPPY_DIR="snippy_results"
QUALIMAP_OUT="qualimap_reports"
mkdir -p "$QUALIMAP_OUT"

for bam in "$SNIPPY_DIR"/*.bam; do
    sample=$(basename "$bam" .bam)
    echo "Running Qualimap BAM QC for sample: $sample"

    outdir="${QUALIMAP_OUT}/${sample}"
    mkdir -p "$outdir"

    qualimap bamqc \
        -bam "$bam" \
        -outdir "$outdir" \
        -outformat pdf:html \
        --java-mem-size=4G
done
```
<details>
<summary>📊 Qualimap BAM QC Script Explanation</summary>

- `#!/bin/bash` → Run script with Bash.  
- `set -euo pipefail` → Exit on errors or undefined variables.  
- `SNIPPY_DIR="all_bams"` → Directory with Snippy BAM files.  
- `QUALIMAP_OUT="qualimap_reports"` → Output directory for QC reports.  
- `mkdir -p "$QUALIMAP_OUT"` → Ensure output directory exists.  
- `for bam in "$SNIPPY_DIR"/*.bam; do ... done` → Loop over all BAM files.  
- `sample=$(basename "$bam" .bam)` → Extract sample name.  
- `outdir="${QUALIMAP_OUT}/${sample}"` → Unique folder per sample.  
- `qualimap bamqc -bam "$bam" -outdir "$outdir" -outformat pdf:html --java-mem-size=4G` → Run QC, generate PDF & HTML, allocate 4 GB memory.  
- `echo "Running Qualimap BAM QC for sample: $sample"` → Print progress.

</details>

##### Step 3: Save and exit nano

Press Ctrl + O, then Enter (save)
Press Ctrl + X (exit)
##### Step 4: Make the script executable
```bash
chmod +x run_qualimap.sh
```
##### Step 5: Activate environment and install GNU Parallel into your `qualimap_env`:
```bash
conda activate qualimap_env
```
##### Step 6: run:
```bash
./run_qualimap.sh
```

### Run MultiQC on Qualimap outputs
[`MultiQC`](https://github.com/MultiQC/MultiQC) 
<details>
<summary>📊 Aggregating BAM QC with MultiQC</summary>

After running **Qualimap BAM QC**, each sample produces individual HTML and PDF reports. For large TB cohorts, opening these reports one by one is inefficient. **MultiQC** solves this by aggregating all results.

### Why we use MultiQC after Qualimap
- **Scans all Qualimap output folders**: Automatically detects per-sample QC reports.  
- **Aggregates QC metrics**: Combines mapping quality, depth, and coverage statistics across all samples.  
- **Quick overview of problematic samples**: Highlights low coverage, uneven depth, or poor alignment.  
- **Standardized reporting**: Ensures results are consistent, shareable, and suitable for team collaboration or publications.  

</details>

---

##### Step 1: **Open nano to create the script `run_multiqc_qualimap.sh`
```bash
nano run_multiqc_qualimap.sh
```
##### Step 2: Paste the following code into nano
```bash
#!/bin/bash
set -euo pipefail

INPUT_DIR="qualimap_reports"
OUTPUT_DIR="multiqc/qualimap_multiqc"

mkdir -p "$OUTPUT_DIR"

if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Input directory '$INPUT_DIR' does not exist!"
    exit 1
fi

multiqc "$INPUT_DIR" -o "$OUTPUT_DIR"

echo "MultiQC report generated in '$OUTPUT_DIR'."

```
<details>
<summary>📊 MultiQC for Qualimap Reports – Script Explanation</summary>

- `#!/bin/bash` → Run the script with Bash.  
- `set -euo pipefail` → Exit on errors, undefined variables, or failed pipelines.  
- `INPUT_DIR="qualimap_reports"` → Directory containing Qualimap BAM QC reports.  
- `OUTPUT_DIR="multiqc/qualimap_multiqc"` → Directory where the aggregated MultiQC report will be saved.  
- `mkdir -p "$OUTPUT_DIR"` → Create `multiqc` and `qualimap_multiqc` directories if they don’t exist.  
- `if [ ! -d "$INPUT_DIR" ]; then ... fi` → Check that the input directory exists; exit with an error if not.  
- `multiqc "$INPUT_DIR" -o "$OUTPUT_DIR"` → Run MultiQC on all files in `INPUT_DIR` and save the combined report in `OUTPUT_DIR`.  
- `echo "MultiQC report generated in '$OUTPUT_DIR'."` → Confirmation message after successful completion.

</details>

##### Step 3: Save & exit nano
Press CTRL+O, Enter (save)
Press CTRL+X (exit)
##### Step 4: Make the script executable
```bash
chmod +x run_multiqc_qualimap.sh
```
##### Step 5: Activate your conda env and run
```bash
conda activate multiqc_env
./run_multiqc_qualimap.sh
```

# 7️⃣ High-Confidence Variant Filtering with **tb_variant_filter**
[`tb_variant_filter`](https://github.com/COMBAT-TB/tb_variant_filter) 
<details>
<summary>🧬 TB-Specific Variant Filtering with <code>tb_variant_filter</code></summary>

The **tb_variant_filter** tool is designed specifically for **Mycobacterium tuberculosis (M. tb)** sequencing data. Unlike generic variant filters, it takes into account TB-specific genomic features and problematic regions in the **H37Rv reference genome**, ensuring that only **high-confidence variants** are retained.

### 🛠 Key Filtering Options

**tb_variant_filter** offers several ways to refine your VCF files:

1. **Region-based filtering**  
   Mask out variants in defined genomic regions. Region lists include:  
   - **RLC (Refined Low Confidence) regions** – Marin et al 2022 (default)  
   - **RLC + Low Mappability regions** – Marin et al 2022  
   - **PE/PPE genes** – Fishbein et al 2015  
   - **Antibiotic resistance genes** – TBProfiler and MTBseq lists  
   - **Repetitive loci** – UVP list  

   > ⚠️ **Default:** Use RLC regions. These are parts of the H37Rv genome where Illumina reads map poorly. For reads shorter than 100 bp or single-ended reads, consider using the **RLC + Low Mappability filter**. PE/PPE and UVP filters are mainly for backward compatibility, but they may exclude too much of the genome.

2. **Window around indels**  
   Masks variants within a set distance (default 5 bases) of insertions or deletions.

3. **Alternate allele percentage**  
   Removes variants with fewer than the minimum percentage (default 90%) of alternative alleles.

4. **Depth of aligned reads**  
   Filters variants based on sequencing depth to remove low-confidence calls.

5. **SNV-only filtering**  
   Optionally discard all variants that are not single nucleotide variants.

### ⚡ How Filters Work Together

When multiple filters are applied, they **stack**: a variant will be masked if **any filter** flags it. This ensures a conservative, high-confidence set of variants for downstream analyses.

</details>

---
## Steps

##### Step 1: Activate environment
```bash
conda activate tb_variant_filter_env
```
##### Step 2:  We need this directory to store all BED files containing genomic regions to be masked during variant filtering
```bash
mkdir -p ./masking_regions
```
#####  Step 3: Generate each BED file in that folder
<details>
<summary>🧬 BED Files for TB Variant Filtering</summary>

We need **Refined Low Confidence (RLC) regions** by default, but we may also need other regions. You can download and store the BED files in your directory to mask genomic regions during variant filtering, such as:

- **Refined Low Confidence (RLC) regions**  
- **Low mappability regions**  
- **PE/PPE genes**  
- **Antibiotic resistance gene lists** (TBProfiler, MTBseq)  
- **Repetitive loci** (UVP)

</details>

```bash
tb_region_list_to_bed --chromosome_name H37Rv farhat_rlc ./masking_regions/RLC_Marin2022.bed
tb_region_list_to_bed --chromosome_name H37Rv farhat_rlc_lowmap ./masking_regions/RLC_and_LowMappability_Marin2022.bed
tb_region_list_to_bed --chromosome_name H37Rv pe_ppe ./masking_regions/PE_PPE_Fishbein2015.bed
tb_region_list_to_bed --chromosome_name H37Rv tbprofiler ./masking_regions/TBProfiler_resistance_genes.bed
tb_region_list_to_bed --chromosome_name H37Rv mtbseq ./masking_regions/MTBseq_resistance_genes.bed
tb_region_list_to_bed --chromosome_name H37Rv uvp ./masking_regions/UVP_repetitive_loci.bed

```
<details>
<summary>🧬 tb_region_list_to_bed Commands Explanation</summary>

- `tb_region_list_to_bed` → Command from `tb_variant_filter` to export predefined genomic regions as BED files.  
- `--chromosome_name H37Rv` → Specifies the reference genome (M. tuberculosis H37Rv).  
- `{farhat_rlc, farhat_rlc_lowmap, pe_ppe, tbprofiler, mtbseq, uvp}` → Region list names available in `tb_variant_filter`:  
  - `farhat_rlc` → Refined Low Confidence regions (Marin et al., 2022).  
  - `farhat_rlc_lowmap` → RLC + Low Mappability regions (Marin et al., 2022).  
  - `pe_ppe` → PE/PPE gene regions (Fishbein et al., 2015).  
  - `tbprofiler` → TBProfiler antibiotic resistance genes.  
  - `mtbseq` → MTBseq antibiotic resistance genes.  
  - `uvp` → UVP repetitive loci in the genome.  
- Last argument → Output BED file path where the exported regions will be saved.  

</details>

##### Step 4: Create or edit the script  
```bash
nano run_tb_variant_filter.sh
```
#####  Step 5: Paste the following into `run_tb_variant_filter.sh`
```bash
#!/bin/bash
set -euo pipefail

CURDIR=$(pwd)
SNIPPY_DIR="$CURDIR/snippy_results"
OUTDIR="$CURDIR/tb_variant_filter_results"
mkdir -p "$OUTDIR"

REGION_FILTER="farhat_rlc"

for vcf in "$SNIPPY_DIR"/*.vcf; do
    sample=$(basename "$vcf")
    echo "Filtering $sample ..."
    
    tb_variant_filter --region_filter "$REGION_FILTER" \
        "$vcf" \
        "$OUTDIR/${sample%.vcf}.filtered.vcf"

        if [ ! -s "$OUTDIR/${sample%.vcf}.filtered.vcf" ]; then
        echo "⚠ $sample filtered VCF is empty."
    fi
done

echo "✅ All VCFs filtered using $REGION_FILTER and saved in $OUTDIR"

```
<details>
<summary>🧬 TB Variant Filter Script Explanation</summary>

- `#!/bin/bash` → Run script with Bash.  
- `set -euo pipefail` → Exit on errors or undefined variables.  
- `CURDIR=$(pwd)` → Save current working directory.  
- `SNIPPY_DIR="$CURDIR/snippy_results"` → Folder containing Snippy VCFs.  
- `OUTDIR="$CURDIR/tb_variant_filter_results"` → Output folder for filtered VCFs.  
- `mkdir -p "$OUTDIR"` → Ensure output directory exists.  
- `REGION_FILTER="farhat_rlc"` → Predefined region filter for TB variant filtering.  
- `for vcf in "$SNIPPY_DIR"/*.vcf; do ... done` → Loop through all Snippy VCFs.  
- `sample=$(basename "$vcf")` → Extract filename for naming filtered outputs.  
- `tb_variant_filter --region_filter "$REGION_FILTER" "$vcf" "$OUTDIR/${sample%.vcf}.filtered.vcf"` → Filter each VCF using the specified region filter and save result.  
- `echo "✅ All VCFs filtered using $REGION_FILTER and saved in $OUTDIR"` → Prints confirmation when all VCFs are filtered.

</details>

##### Step 6: Save and exit nano
Press Ctrl + O → Enter (to write the file)
Press Ctrl + X → Exit nano

##### Step 7: Make the script executable
```bash
chmod +x run_tb_variant_filter.sh
```
##### Step 8: run
```bash
./run_tb_variant_filter.sh
```
<details>
<summary>📊 VCF QC Script Step-by-Step Guide</summary>

### Purpose

This script is designed to compare unfiltered Snippy VCFs with filtered VCFs produced by `tb_variant_filter`. The main goals are:

- Count the total number of variants in each VCF.
- Identify variants that meet a custom “PASS” criterion, based on user-defined quality thresholds.
- Calculate the PASS retention ratio, i.e., the fraction of high-quality variants retained after filtering.

This QC step ensures that your filtered VCFs retain high-confidence variants and helps detect potential over-filtering or loss of important variants before downstream analysis.

---

### Key Parameters

- `MIN_DP=20` → Minimum read depth required for a variant to be considered “high quality.”  
  - Read depth (DP in the VCF INFO field) indicates the number of reads supporting a variant.  
  - Variants with fewer than 20 supporting reads are considered low-confidence and are excluded from the PASS count.

- `MIN_QUAL=30` → Minimum variant quality score (QUAL in the VCF) to consider a variant as PASS.  
  - The QUAL score represents the confidence that the variant is real.  
  - Variants with QUAL < 30 are considered low-confidence and are excluded from the PASS count.

> These thresholds can be adjusted depending on your experimental design and sequencing quality.

---

### What the Script Does
For each sample:
**Unfiltered VCF analysis**
- Counts all variants (`Unfiltered_total`).
- Counts variants meeting `DP ≥ 20` and `QUAL ≥ 30` (`Unfiltered_PASS`).
**Filtered VCF analysis**
- Counts all variants in the filtered VCF (`Filtered_total`).
- Counts variants meeting `DP ≥ 20` and `QUAL ≥ 30` (`Filtered_PASS`).
**Calculate PASS retention ratio**
- `PASS_retention_ratio = Filtered_PASS / Unfiltered_PASS`  
- Measures how many high-quality variants remain after filtering.
---

### Output Table Columns

| Column               | Description                                                                 |
|----------------------|-----------------------------------------------------------------------------|
| Sample               | Name of the sample (derived from the VCF filename)                           |
| Unfiltered_total     | Total variants in the original Snippy VCF                                    |
| Unfiltered_PASS      | Variants meeting `DP ≥ 20` and `QUAL ≥ 30` in unfiltered VCF                |
| Filtered_total       | Total variants in tb_variant_filter output                                   |
| Filtered_PASS        | Variants meeting `DP ≥ 20` and `QUAL ≥ 30` in filtered VCF                  |
| PASS_retention_ratio | Fraction of high-quality variants retained after filtering                  |

</details>

##### Step 1: Open a new file in nano
```bash
nano compare_vcf_qc.sh
```
##### Step 2: Paste the script
```bash
#!/bin/bash
set -euo pipefail

SNIPPY_DIR="snippy_results"
FILTERED_DIR="tb_variant_filter_results"
OUTDIR="csv_output"
OUTFILE="$OUTDIR/variant_filter_summary.csv"

MIN_DP=20
MIN_QUAL=30

mkdir -p "$OUTDIR"
echo "Sample,Unfiltered_total,Unfiltered_PASS,Filtered_total,Filtered_PASS,PASS_retention_ratio" | tee "$OUTFILE"

count_pass_variants() {
    local vcf=$1
    awk -v min_dp=$MIN_DP -v min_qual=$MIN_QUAL '
        BEGIN{count=0}
        /^#/ {next}
        {
            qual=$6
            dp=0
            split($8, info_arr, ";")
            for(i in info_arr){
                if(info_arr[i] ~ /^DP=/){
                    split(info_arr[i], kv, "=")
                    dp=kv[2]
                }
            }
            if(qual >= min_qual && dp >= min_dp) count++
        }
        END{print count}
    ' "$vcf"
}

for vcf in "$SNIPPY_DIR"/*.vcf; do
    sample=$(basename "$vcf" .vcf)

    # safety: skip if VCF is empty
    if [[ ! -s "$vcf" ]]; then
        echo "$sample,NA,NA,NA,NA,NA" | tee -a "$OUTFILE"
        continue
    fi

    unfiltered_total=$(grep -v "^#" "$vcf" | wc -l)
    unfiltered_pass=$(count_pass_variants "$vcf")

    if [[ ! -f "$FILTERED_DIR/$sample.filtered.vcf" ]]; then
        echo "$sample,$unfiltered_total,$unfiltered_pass,NA,NA,NA" | tee -a "$OUTFILE"
        continue
    fi

    filtered_total=$(grep -v "^#" "$FILTERED_DIR/$sample.filtered.vcf" | wc -l)
    filtered_pass=$(count_pass_variants "$FILTERED_DIR/$sample.filtered.vcf")

    if [[ $unfiltered_pass -eq 0 ]]; then
        ratio="NA"
    else
        ratio=$(awk -v pass=$filtered_pass -v total=$unfiltered_pass 'BEGIN{printf "%.2f", pass/total}')
    fi

    echo "$sample,$unfiltered_total,$unfiltered_pass,$filtered_total,$filtered_pass,$ratio" | tee -a "$OUTFILE"
done
```

##### Step  3: Save and exit nano
   Press Ctrl+O → Enter to save.
   Press Ctrl+X → Exit nano.
##### Step  4:Make the script executable
```bash
chmod +x compare_vcf_qc.sh
```
##### Step 5: Run the script
```bash
./compare_vcf_qc.sh
```

# 8️⃣ Consensus Genome Generation with **BCFtools / SAMtools** & Outgroup Inclusion
[`BCFtools`](https://github.com/samtools/bcftools) / [`SAMtools`](https://github.com/samtools/samtools)
<details>
<summary>🧬 Generate Sample-Specific Consensus Sequences</summary>

After filtering VCFs with **tb_variant_filter**, we generate **consensus FASTA sequences** for each sample. These sequences represent the **full genome of each isolate**, including only **high-confidence variants** relative to the reference genome (*H37Rv*).

### Why consensus sequences matter
- Provide a **single representative genome** per sample for downstream analyses.  
- Used in **phylogenetic reconstruction**, outbreak investigation, and comparative genomics.  
- Incorporate **only reliable SNPs and indels**, minimizing noise from sequencing errors.  
- Standardize genome representations across multiple isolates, ensuring comparability.  

</details>
---

##### Step 1: Compress and index each filtered VCF
```bash
for vcf in tb_variant_filter_results/*.vcf; do
   
    if [ ! -s "$vcf" ]; then
        echo "Skipping empty VCF: $vcf"
        continue
    fi

    gz_file="${vcf}.gz"
    echo "Compressing $vcf → $gz_file"
    bgzip -c "$vcf" > "$gz_file"

    echo "Indexing $gz_file"
    bcftools index "$gz_file"
done
```
##### Step 2: Create a script to generate consensus sequences
```bash
nano generate_consensus_all.sh
```
##### Step 3: Paste this code
```bash
#!/bin/bash
set -euo pipefail

CURDIR=$(pwd)
VCFDIR="$CURDIR/tb_variant_filter_results"
OUTDIR="$CURDIR/consensus_sequences"
mkdir -p "$OUTDIR"

for vcf in "$VCFDIR"/*.vcf; do
    sample=$(basename "$vcf" .vcf)

    if [ $(grep -v '^#' "$vcf" | wc -l) -eq 0 ]; then
        echo "$sample VCF is empty. Skipping."
        continue
    fi

    gz_file="${vcf}.gz"
    bgzip -c "$vcf" > "$gz_file"
    bcftools index "$gz_file"
    bcftools consensus -f "$CURDIR/H37Rv.fasta" "$gz_file" | sed "1s/.*/>$sample/" > "$OUTDIR/${sample}.consensus.fasta"
done

```
<details>
<summary>🧬 VCF-to-Consensus Script Explanation</summary>

- `#!/bin/bash` → Run script with Bash.  
- `set -euo pipefail` → Exit on errors or undefined variables.  
- `CURDIR=$(pwd)` → Save current working directory.  
- `VCFDIR="$CURDIR/tb_variant_filter_results"` → Folder with filtered VCFs.  
- `OUTDIR="$CURDIR/consensus_sequences"` → Folder for consensus FASTA sequences.  
- `mkdir -p "$OUTDIR"` → Ensure output directory exists.  
- `for vcf in "$VCFDIR"/*.vcf; do ... done` → Loop through all filtered VCF files.  
- `sample=$(basename "$vcf" .vcf)` → Extract sample name.  
- `bgzip -c "$vcf" > "$vcf.gz"` → Compress VCF with bgzip.  
- `bcftools index "$vcf.gz"` → Index compressed VCF.  
- `bcftools consensus -f "$CURDIR/H37Rv.fasta" "$vcf.gz" | sed "1s/.*/>$sample/" > "$OUTDIR/${sample}.consensus.fasta"` → Generate consensus FASTA and replace header with sample name.  
- `echo "✅ $sample consensus generated"` → Confirmation per sample.  
- `echo "🎉 All consensus sequences saved in $OUTDIR"` → Final message.  

**⚠ Note:** Activate the `tb_consensus_env` before running this script.

</details>

##### Step 4: Save and exit nano
Press Ctrl + O → Enter (to write the file)
Press Ctrl + X → Exit nano

##### Step 5: Make the script executable
```bash
chmod +x generate_consensus_all.sh
```
##### Step 6: Run the script
```bash
conda activate tb_consensus_env
./generate_consensus_all.sh
```
Check Consensus FASTA Lengths

After generating consensus sequences, it's important to **verify the genome length** for each sample.  
This ensures no sequences are truncated or incomplete due to missing coverage or filtering.

---

### 📏 Calculating Consensus Genome Lengths

We can check the length of each consensus FASTA sequence to ensure completeness and consistency.  
This helps verify that consensus sequences cover the full *M. tuberculosis* genome (~4.4 Mbp) and can reveal missing regions.

### Rename the FASTA files

```bash
#!/bin/bash
FASTA_DIR="consensus_sequences"

for f in "$FASTA_DIR"/*.filtered.consensus.fasta; do
    mv "$f" "${f/.filtered.consensus/}"
done

echo "✅ All consensus FASTA files have been renamed to .fasta."

```

<details>
<summary>📝 Rename Consensus FASTA Files</summary>

- `FASTA_DIR="consensus_sequences"` → Directory containing consensus FASTA files.  
- `for f in "$FASTA_DIR"/*.filtered.consensus.fasta; do ... done` → Loop over all FASTA files ending with `.filtered.consensus.fasta`.  
- `mv "$f" "${f/.filtered.consensus/}"` → Rename each file by removing `.filtered.consensus` from filename.  
- `echo "✅ All consensus FASTA files have been renamed to .fasta."` → Confirmation message after renaming.

</details>

###  Update headers inside the FASTA files
```bash
#!/bin/bash
FASTA_DIR="consensus_sequences"

for f in "$FASTA_DIR"/*.fasta; do
    sample=$(basename "$f" .fasta)
    awk -v s="$sample" '/^>/{print ">" s; next} {print}' "$f" > "${f}.tmp" && mv "${f}.tmp" "$f"
    echo "✅ Updated header in: $(basename "$f")"
done
echo "🎉 All FASTA headers have been successfully updated."
```
<details>
<summary>📝 Update FASTA Headers with Sample Names</summary>

- `FASTA_DIR="consensus_sequences"` → Directory containing FASTA files.  
- `for f in "$FASTA_DIR"/*.fasta; do ... done` → Loop through all FASTA files.  
- `sample=$(basename "$f" .fasta)` → Extract sample name from filename.  
- `awk -v s="$sample" '/^>/{print ">" s; next} {print}' "$f" > "${f}.tmp" && mv "${f}.tmp" "$f"` → Replace FASTA header with `>sample`, keep sequence lines unchanged.  
- `echo "✅ Updated header in: $(basename "$f")"` → Logs each updated file.  
- `echo "🎉 All FASTA headers have been successfully updated."` → Completion message after all files processed.

</details>


### Using `grep` and `wc`
We remove the FASTA headers (`>` lines) and count the remaining nucleotides to get the total genome length:

```bash
for f in consensus_sequences/*.fasta; do
    sample=$(basename "$f")
    length=$(grep -v ">" "$f" | tr -d '\n' | wc -c)
    echo "$sample : $length bp"
done
```

to save the result in csv file 
```bash
#!/bin/bash
FASTA_DIR="consensus_sequences"
OUTDIR="csv_output"
mkdir -p "$OUTDIR"  
OUTPUT_CSV="${OUTDIR}/consensus_lengths.csv"

echo "Sample,Length_bp" > "$OUTPUT_CSV"

for f in "$FASTA_DIR"/*.fasta; do
    sample=$(basename "$f" .fasta)
    length=$(grep -v ">" "$f" | tr -d '\n' | wc -c)
    echo "$sample,$length" >> "$OUTPUT_CSV"
done

echo "✅ Consensus genome lengths saved to $OUTPUT_CSV"

```
<details>
<summary>📖 Explanation of calculating consensus genome lengths and saving the result in CSV</summary>

- `FASTA_DIR="consensus_sequences"` → sets the directory containing consensus FASTA files.  
- `OUTDIR="csv_output"` → defines the directory where the CSV will be saved.  
- `mkdir -p "$OUTDIR"` → ensures the output directory exists before writing the file.  
- `OUTPUT_CSV="${OUTDIR}/consensus_lengths.csv"` → defines the CSV file path inside `OUTDIR`.  
- `echo "Sample,Length_bp" > "$OUTPUT_CSV"` → creates the CSV file and writes the header line.  
- `for f in "$FASTA_DIR"/*.fasta; do ... done` → loops over all FASTA files in the directory.  
- `sample=$(basename "$f" .fasta)` → extracts the sample name from the FASTA filename.  
- `length=$(grep -v ">" "$f" | tr -d '\n' | wc -c)` → removes header lines (`>`), concatenates sequences into one line, and counts nucleotides.  
- `echo "$sample,$length" >> "$OUTPUT_CSV"` → appends the sample name and its sequence length to the CSV file.  
- `echo "✅ Consensus genome lengths saved to $OUTPUT_CSV"` → prints a completion message when finished.  

</details>

# 9️⃣ Multiple Sequence Alignment

[`MAFFT`](https://github.com/GSLBiotech/mafft) requires **a single FASTA file** as input.  
It **cannot take multiple FASTA files** on the command line directly, otherwise it interprets filenames as options.

---

Before performing phylogenetic analysis, it is crucial to include an **outgroup sequence**. The outgroup serves as a reference point to root the phylogenetic tree, providing directionality for evolutionary relationships.  

For our analysis, we searched the NCBI database for sequences belonging to **Lineage 8** of *Mycobacterium tuberculosis*. We selected the sequence `SRR10828835`, which originates from Rwanda, as our outgroup. This choice was intentional because Lineage 8 is **not present in Ethiopia**, ensuring it is sufficiently divergent from our study isolates and providing a robust root for the tree.  

The selected outgroup sequence was downloaded and saved in the `consensus_sequences/` directory alongside our aligned TB consensus sequences.

---
Perform multiple sequence alignment on the merged file:
Use a faster algorithm than --auto

--auto lets MAFFT decide, but for hundreds/thousands of sequences, it often picks a slower iterative refinement.
You can explicitly pick faster modes:

Since Mycobacterium tuberculosis genomes are highly conserved and we care more about speed than very small accuracy gains, we can safely drop the expensive iterative refinement steps in MAFFT.
fast and TB-suitable command

##### Step 1: Merge all consensus FASTAs
We combine all individual consensus sequences into one multi-FASTA file:
- **Aligning sequences**  
  - All sequences, including the outgroup, must be aligned together.  
  - This ensures that **homologous positions are matched across all sequences**, which is essential for correct tree inference.  
  - Example workflow:  
    1. Combine all consensus sequences and the outgroup FASTA into a single folder.  
    2. Run MAFFT on the combined sequences to produce a single aligned file (e.g., `aligned_consensus.fasta`).  
  - The outgroup sequence must already be included in the alignment.  
  - When running IQ-TREE, we need to specify the outgroup using the `-o` option with the **FASTA header name** of the outgroup:  
  - IQ-TREE will then **root the phylogenetic tree using this outgroup**.
    
```bash
cat consensus_sequences/*.fasta > consensus_sequences/all_consensus.fasta
```

##### Step 2: Paste the script into nano
```bash
#!/bin/bash
mkdir -p mafft_results
mafft --auto --parttree consensus_sequences/all_consensus.fasta > mafft_results/aligned_consensus.fasta
```
##### Step 3: Verify the alignment

A. Quickly inspect the top of the aligned FASTA:
```bash
head mafft_results/aligned_consensus.fast
```
B. Quick visual inspection with `less`
```bash
fold -w 100 mafft_results/aligned_consensus.fasta | less
```
C. Checking MAFFT alignment output to 
- Ensure the alignment file **exists** and contains all intended sequences.  
- Verify the **number of sequences** matches expectations.  
- Check **sequence lengths** (total, min, max, average) to confirm proper alignment.  
- Detect if sequences have **varying lengths**, which may indicate misalignment or excessive gaps.  
- A **good alignment** is essential for accurate phylogenetic tree inference with IQ-TREE or other downstream analyses.

```bash
#!/bin/bash
set -euo pipefail

FILE="mafft_results/aligned_consensus.fasta"

if [[ ! -f "$FILE" ]]; then
    echo "❌ $FILE not found!"
    exit 1
fi

SEQ_COUNT=$(grep -c ">" "$FILE")
LENGTHS=($(grep -v ">" "$FILE" | awk 'BEGIN{RS=">"} NR>1{print length($0)}'))
TOTAL_LENGTH=$(IFS=+; echo "$((${LENGTHS[*]}))")
MIN_LENGTH=$(printf "%s\n" "${LENGTHS[@]}" | sort -n | head -n1)
MAX_LENGTH=$(printf "%s\n" "${LENGTHS[@]}" | sort -n | tail -n1)
AVG_LENGTH=$(( TOTAL_LENGTH / SEQ_COUNT ))

echo "✅ $FILE exists"
echo "Sequences: $SEQ_COUNT"
echo "Total length (bp): $TOTAL_LENGTH"
echo "Min length: $MIN_LENGTH"
echo "Max length: $MAX_LENGTH"
echo "Avg length: $AVG_LENGTH"

if [[ $MIN_LENGTH -eq $MAX_LENGTH ]]; then
    echo "✅ All sequences have equal length, alignment looks good."
else
    echo "⚠️ Sequence lengths vary; check for gaps or misalignment."
fi

```
<details>
<summary>📖 Explanation of MAFFT output check script</summary>

- `FILE="mafft_results/aligned_consensus.fasta"` → sets the path to the MAFFT alignment output file.  

- `if [[ ! -f "$FILE" ]]; then ... fi` → checks if the file exists; exits with a message if it does not.  

- `SEQ_COUNT=$(grep -c ">" "$FILE")` → counts the number of sequences in the alignment.  
  - In FASTA format, each sequence starts with a `>` header.  

- `LENGTHS=($(grep -v ">" "$FILE" | awk 'BEGIN{RS=">"} NR>1{print length($0)}'))` → creates an array of sequence lengths.  
  - Removes header lines and calculates the length of each sequence.  

- `TOTAL_LENGTH=$(IFS=+; echo "$((${LENGTHS[*]}))")` → sums all sequence lengths to get the total number of base pairs.  

- `MIN_LENGTH=$(printf "%s\n" "${LENGTHS[@]}" | sort -n | head -n1)` → finds the shortest sequence length.  

- `MAX_LENGTH=$(printf "%s\n" "${LENGTHS[@]}" | sort -n | tail -n1)` → finds the longest sequence length.  

- `AVG_LENGTH=$(( TOTAL_LENGTH / SEQ_COUNT ))` → calculates the average sequence length.  

- `echo ...` → prints summary statistics for the alignment: number of sequences, total, min, max, and average lengths.  

- `if [[ $MIN_LENGTH -eq $MAX_LENGTH ]]; then ... else ... fi` → checks if all sequences are the same length.  
  - If yes → alignment looks good.  
  - If no → prints a warning that sequence lengths vary, indicating possible gaps or misalignment.  

</details>


D. Using `seqkit` stats (recommended)
seqkit is a fast toolkit for FASTA/Q file summaries. It gives a detailed report of sequences in a file:
```bash
conda activate seqkit_env
seqkit stats mafft_results/aligned_consensus.fasta
```
E. Check for gaps / alignment columns
<details>
<summary>💡 Understanding Gaps in Multiple Sequence Alignment (MSA)</summary>

When you do a **multiple sequence alignment (MSA):**

- All sequences should line up **column by column**.
- **MAFFT** will add **gaps (`-`)** where needed so that homologous positions align.
- A good alignment usually has **all sequences the same length**, including gaps.
- If sequence lengths differ, it often indicates:
  - Truncated sequences
  - Missing data
  - Alignment problems

So, checking the **length of each sequence** in the alignment is a quick way to see if **gaps and sequences are consistent**.

</details>

```bash
grep -v ">" mafft_results/aligned_consensus.fasta \
| awk '{print length($0)}' \
| sort -n \
| uniq -c \
| awk -v L=60 '{print ($2<L?"⚠️ ":"✅ ") $1 " sequences of length " $2 " bp"}'
```
<details>
<summary>🔍 Line-by-Line Explanation of Alignment Check Command</summary>

- `grep -v ">" mafft_results/aligned_consensus.fasta` → removes FASTA headers, leaving only sequence lines.  
- `awk '{print length($0)}'` → prints the length of each sequence line.  
- `sort -n` → sorts sequence lengths numerically.  
- `uniq -c` → counts how many sequences have each unique length.  
- `awk -v L=60 '{print ($2<L?"⚠️ ":"✅ ") $1 " sequences of length " $2 " bp"}'` → formats output: ⚠️ for lengths <60, ✅ for ≥60, showing count and length.

</details>


F. Compute pairwise identity
```bash
awk '/^>/{if(seqlen){print seqlen}; seqlen=0; next} {seqlen+=length($0)} END{print seqlen}' mafft_results/aligned_consensus.fasta \
| sort -n \
| uniq -c \
| awk -v L=60 '{print ($2<L?"⚠️ ":"✅ ") $1 " sequences of length " $2 " bp"}'
```
G. Use AMAS (Alignment Manipulation and Summary)
AMAS is a Python tool to summarize alignments:
```bash
conda activate amas_env
AMAS.py summary -f fasta -d dna -i mafft_results/aligned_consensus.fasta
```
H. Use aliview or MEGA for GUI inspection
Load the FASTA alignment in AliView, MEGA, or Geneious.
Advantages:
  Can visually check gaps, conserved regions, and misaligned sequences.
  Highlight sequences that differ significantly.


# 🔟 Phylogenetic Tree Construction with **IQ-TREE**
[`IQ-TREE`](https://github.com/iqtree/iqtree2) 
<details>
<summary>📖 Overview of IQ-TREE and important concepts</summary>

- **What is IQ-TREE?**  
  - IQ-TREE is a **phylogenetic tree inference software** used to construct evolutionary trees from aligned sequences.  
  - It implements **maximum likelihood (ML) methods**, which estimate the tree that best explains the observed sequences under a substitution model.  

- **Input**  
  - Requires an **aligned sequence file** (FASTA, PHYLIP, or NEXUS).  
  - Sequences must be aligned so that homologous positions are in the same columns.  

- **Substitution models (-m)**  
  - Models describe how nucleotides (or amino acids) change over time.  
  - Common nucleotide models:  
    - **GTR**: General Time Reversible model, allows different rates for each nucleotide substitution.  
    - **+G**: Gamma distribution to account for rate variation across sites.  
    - **+I**: Proportion of invariant sites.  
  - IQ-TREE can **automatically select the best model** using `-m MFP` (ModelFinder Plus).  

- **Bootstrap support (-bb)**  
  - Ultrafast bootstrap replicates measure **confidence in the tree branches**.  
  - Example: `-bb 1000` runs 1000 replicates.  

- **Outgroup rooting (-o)**  
  - Specify the outgroup sequence (by its FASTA header) to **root the tree** correctly.  
  - The outgroup must be included in the alignment.  

- **Multithreading (-nt)**  
  - IQ-TREE can use multiple CPU threads for faster computation, e.g., `-nt 4`.  

- **Output files**  
  - `.treefile` → final tree in Newick format.  
  - `.log` → details of the run.  
  - `.iqtree` → summary of models, likelihoods, and bootstrap values.  

- **Key points**  
  - Input sequences must be properly aligned.  
  - The choice of substitution model affects tree accuracy.  
  - Use the outgroup for rooting but include it in the alignment.  
  - Bootstrapping provides branch support values for interpreting the tree.  

</details>

 Steps 
##### Step 1: activate iqtree environment
```bash
conda activate iqtree_env
```
##### Step 2: Run the following script
```bash
mkdir -p iqtree_results

iqtree2 -s mafft_results/aligned_consensus.fasta \
        -m GTR+G \
        -bb 1000 \
        -nt 4 \
        -o SRR10828835 \
        -pre iqtree_results/aligned_consensus
```
<details>
<summary>📖 Explanation of IQ-TREE command for aligned consensus sequences</summary>

- `mkdir -p iqtree_results` → creates the directory to store IQ-TREE output files if it doesn’t exist.  
- `iqtree2` → runs IQ-TREE version 2, a program for phylogenetic tree inference.  
- `-s mafft_results/aligned_consensus.fasta` → specifies the input alignment file generated by MAFFT.  
- `-m GTR+G` → sets the substitution model to GTR (General Time Reversible) with Gamma rate heterogeneity.  
- `-bb 1000` → performs 1000 ultrafast bootstrap replicates to assess branch support.  
- `-nt 4` → uses 4 CPU threads for faster computation.  
- `-pre iqtree_results/aligned_consensus` → sets the output file prefix and saves all IQ-TREE results in `iqtree_results/` with this prefix.  

</details>

### 🌳 Visualization with TB-Profiler + iTOL

TB-Profiler outputs can be integrated into [iTOL](https://itol.embl.de/) for rich phylogenetic visualization, combining resistance, lineage, and metadata.

* **Drug resistance types**:  
  - **HR-TB**: Isoniazid-resistant TB  
  - **RR-TB**: Rifampicin-resistant TB  
  - **MDR-TB**: Multidrug-resistant TB (at least INH + RIF resistant)  
  - **XDR-TB**: Extensively drug-resistant TB (MDR + fluoroquinolone + second-line injectable resistant)  
  - **Other**: Cases that do not fit the above categories  

  These categories can be represented in iTOL as colored rings around the phylogeny for quick classification.  

* **Individual drug resistance profiles**:  
  Each anti-TB drug is represented as a binary dot in iTOL:  
  - **Black dot** = Resistant  
  - **White dot** = Susceptible  

  This scheme provides a clear and compact way to compare resistance across isolates without cluttering the tree.  

* **Lineages & Sublineages**:  
  - Major M. tuberculosis lineages (1–10) can be highlighted with branch colors or background shading.  
  - Sublineages (e.g. 4.2.2.2, 2.2.1) can be shown as labels or an extra ring.  



# 📖 References

1. WHO. *Catalogue of mutations in MTBC and their association with drug resistance*, 2nd ed, 2023.  

2. Chen, S., Zhou, Y., Chen, Y., & Gu, J. (2018). [*fastp: an ultra-fast all-in-one FASTQ preprocessor*](https://doi.org/10.1093/bioinformatics/bty560). **Bioinformatics**, 34(17), i884–i890.  

3. Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). [*MultiQC: summarize analysis results for multiple tools and samples in a single report*](https://doi.org/10.1093/bioinformatics/btw354). **Bioinformatics**, 32(19), 3047–3048.  

4. Coll, F., McNerney, R., Preston, M. D., Guerra-Assunção, J. A., Warry, A., Hill-Cawthorne, G., … Clark, T. G. (2015). [*Rapid determination of anti-tuberculosis drug resistance from whole-genome sequences*](https://doi.org/10.1186/s13073-015-0164-0). **Genome Medicine**, 7, 51.  

5. Seemann, T. (2015). *Snippy: rapid haploid variant calling and core genome alignment*. Available at: [https://github.com/tseemann/snippy](https://github.com/tseemann/snippy)  

6. Okonechnikov, K., Conesa, A., & García-Alcalde, F. (2016). [*Qualimap 2: advanced multi-sample quality control for high-throughput sequencing data*](https://doi.org/10.1093/bioinformatics/btv566). **Bioinformatics**, 32(2), 292–294.  

7. Walker, T. M., Merker, M., Knoblauch, A. M., Helbling, P., Schoch, O. D., van der Werf, T. S., … Niemann, S. (2022). *A cluster of multidrug-resistant Mycobacterium tuberculosis among patients arriving in Europe from the Horn of Africa: a molecular epidemiological study*. **Lancet Microbe**, 3(9), e672–e681. (Supplementary methods describe tb_variant_filter). Available at: [https://github.com/iqbal-lab-org/tb_variant_filter](https://github.com/iqbal-lab-org/tb_variant_filter)  

8. Danecek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O., … Li, H. (2021). [*Twelve years of SAMtools and BCFtools*](https://doi.org/10.1093/gigascience/giab008). **GigaScience**, 10(2), giab008.  

9. Katoh, K., & Standley, D. M. (2013). [*MAFFT multiple sequence alignment software version 7: improvements in performance and usability*](https://doi.org/10.1093/molbev/mst010). **Molecular Biology and Evolution**, 30(4), 772–780.  

10. Minh, B. Q., Schmidt, H. A., Chernomor, O., Schrempf, D., Woodhams, M. D., von Haeseler, A., & Lanfear, R. (2020). [*IQ-TREE 2: new models and efficient methods for phylogenetic inference in the genomic era*](https://doi.org/10.1093/molbev/msaa015). **Molecular Biology and Evolution**, 37(5), 1530–1534.  

11. Letunic, I., & Bork, P. (2021). [*Interactive Tree Of Life (iTOL) v5: an online tool for phylogenetic tree display and annotation*](https://doi.org/10.1093/nar/gkab301). **Nucleic Acids Research**, 49(W1), W293–W296.  

