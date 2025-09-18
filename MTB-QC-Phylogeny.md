
## 🧪 MTB Genome Analysis of XDR-TB from Ethiopia

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
Every time we install tools in a new environment, we deactivate the conda 
```bash
conda deactivate 
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
conda deactivate
conda install -n base -c conda-forge mamba -y
mamba create -n tbprofiler_env -c bioconda -c conda-forge python=3.9 tb-profiler -y
conda activate tbprofiler_env
tb-profiler --version
```

Snippy
```bash
conda deactivate
conda create -n snippy_env -c bioconda -c conda-forge snippy -y
conda activate snippy_env
snippy --version
```
Qualimap
```bash
conda create -n qualimap_env -c bioconda -c conda-forge qualimap=2.2.2d -y
conda activate qualimap_env
qualimap --help

```
tb_variant_filter (COMBAT-TB)
```bash
conda install mamba -c conda-forge
mamba create -n tb_variant_filter_env python=3.9
mamba activate tb_variant_filter_env
mamba install -c bioconda tb_variant_filter

```
BCFtools / SAMtools
```bash
conda create -n bcftools_env -c bioconda -c conda-forge bcftools=1.17 samtools=1.17 -y
conda activate bcftools_env
bcftools --version
samtools --version

```
MAFFT
```bash
conda create -n mafft_env -c bioconda -c conda-forge mafft=7.505 -y
conda activate mafft_env
mafft --version

```
IQ-TREE
```bash
conda create -n iqtree_env -c bioconda -c conda-forge iqtree=2.2.2.3 -y
conda activate iqtree_env
iqtree2 --version
```
iTOL
Web-based tool: https://itol.embl.de



# 1️⃣ Explore Raw FASTQ Files

### FASTQ summary


Before starting any job we have changed to our working directory
```bash
cd /home/betselot/Project_data/TB/NTRL
```
We renamed all paired-end sequencing files so that forward reads end with `_1.fastq.gz` and reverse reads end with `_2.fastq.gz`, following the standard naming convention expected by most bioinformatics tools.

```bash
for f in raw_data/*_R1.fastq.gz; do
    base=${f%_R1.fastq.gz}
    mv "$f" "${base}_1.fastq.gz"
done

for f in raw_data/*_R2.fastq.gz; do
    base=${f%_R2.fastq.gz}
    mv "$f" "${base}_2.fastq.gz"
done
```
Paired-end FASTQ files are first inspected using simple command-line tools.

### 1. Peek at the first few reads
Before analyzing sequencing data, it’s useful to inspect the first few reads in each FASTQ file. This allows a quick check of file format, sequence identifiers, and the general structure of paired-end reads:

```bash
zcat raw_data/ETRS-003_1.fastq.gz | head -n 16
zcat raw_data/ETRS-003_2.fastq.gz | head -n 16
```
<details>
<summary>🔍 Inspect First Reads of FASTQ Files</summary>
- `zcat` is a Linux/Unix command used to **view the contents of compressed files** without manually uncompressing them.  
- Works primarily with `.gz` files (gzip-compressed).  
- Unlike `gunzip`, it **prints the uncompressed data to standard output** instead of creating a new file.  
- Example workflow:  
  - `zcat raw_data/ETRS-003_1.fastq.gz | head -n 16` →  Show the first 16 lines of R1 (forward) FASTQ file of a compressed file without extracting it. 
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
echo $(( $(zcat raw_data/ETRS-003_1.fastq.gz | wc -l) / 4 ))
echo $(( $(zcat raw_data/ETRS-003_2.fastq.gz | wc -l) / 4 ))
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

for R1 in "$INDIR"/*_1.fastq.gz; do
    [[ -f "$R1" ]] || continue
    SAMPLE=$(basename "$R1" | sed -E 's/_1\.fastq\.gz//')
    R2="$INDIR/${SAMPLE}_2.fastq.gz"

    R1_COUNT=$(( $(zcat "$R1" | wc -l) / 4 ))
    if [[ -f "$R2" ]]; then
        R2_COUNT=$(( $(zcat "$R2" | wc -l) / 4 ))
    else
        R2_COUNT="NA"
    fi

    echo "$SAMPLE,$R1_COUNT,$R2_COUNT" >> "$OUTFILE"
    echo "✅ $SAMPLE → R1: $R1_COUNT | R2: $R2_COUNT"
done

echo "🎉 All done! Read counts saved to '$OUTFILE'"

```

<details>
<summary>📄 Script Overview</summary>

Counts reads in paired-end FASTQ files and outputs a CSV.  

- 🛡 `set -euo pipefail` → safe script execution  
- 📂 `INDIR="raw_data"` → input FASTQ folder  
- 📂 `OUTDIR="csv_output"` → output CSV folder (auto-created)  
- 📝 `OUTFILE="$OUTDIR/fastq_read_counts.csv"` → CSV to store counts  
- 🔄 Loop through R1 files: `for R1 in "$INDIR"/*_1.fastq.gz`  
- ⛔ Skip missing files: `[[ -f "$R1" ]] || continue`  
- 🏷 Extract sample name: `basename "$R1" | sed -E 's/_1\.fastq\.gz//'`  
- 🔗 Find R2: `R2="$INDIR/${SAMPLE}_2.fastq.gz"`  
- 📊 Count reads: `R1_COUNT=$(zcat "$R1" | wc -l /4)`, `R2_COUNT=…` or `"NA"`  
- ➕ Append to CSV: `echo "$SAMPLE,$R1_COUNT,$R2_COUNT" >> "$OUTFILE"`  
- ✅ Progress message per sample  
- 🎉 Completion message  

</details>


### 3. Base composition
Checking the nucleotide composition of each FASTQ file helps assess sequencing quality. Balanced proportions of A, T, G, and C indicate high-quality data with minimal bias. This script counts the occurrence of each base in both paired-end files:

```bash
#!/bin/bash
for fq in raw_data/ETRS-003_1.fastq.gz raw_data/ETRS-003_2.fastq.gz; do
    [ -f "$fq" ] || continue
    echo "Counting bases in $fq..."
    zcat "$fq" | awk 'NR%4==2 { for(i=1;i<=length($0);i++) b[substr($0,i,1)]++ } 
    END { for(base in b) print base, b[base] }'
    echo "----------------------"
done
```
<details>
<summary>🧬 Base Counting Script Overview</summary>

Counts nucleotide bases in specific FASTQ files.  

- 🔄 Loops through specified FASTQ files: `raw_data/ETRS-003_1.fastq.gz` and `_2.fastq.gz`  
- ⛔ Skips missing files: `[ -f "$fq" ] || continue`  
- 🖨 Prints message: `"Counting bases in $fq..."`  
- 🧮 Counts bases using `awk`:  
  - `NR%4==2` → selects sequence lines  
  - Loops over each character to tally `A, C, G, T, N`  
  - Prints counts per base at the end  
- ➖ Prints separator `"----------------------"` between files  

</details>



### 4. Quality score summary
FASTQ files encode base quality scores on the 4th line of every read. Checking these scores provides an initial assessment of sequencing quality before trimming or downstream analysis:
First 10 quality lines
```bash
zcat raw_data/ETRS-003_1.fastq.gz | sed -n '4~4p' | head -n 10
zcat raw_data/ETRS-003_2.fastq.gz | sed -n '4~4p' | head -n 10
```
<details>
<summary>🔍 FASTQ Quality Lines Preview</summary>

Displays the first 10 quality lines from paired-end FASTQ files.  

- 🔄 Reads R1 and R2 FASTQ files using `zcat`  
- 📝 `sed -n '4~4p'` → selects every 4th line (quality scores)  
- 👀 `head -n 10` → shows only the first 10 quality lines  
- 🖨 Prints output directly to the terminal  

</details>


FASTQ quality scores are encoded as ASCII characters. Counting the occurrence of each character provides a quantitative overview of base quality across the reads:
```bash
zcat raw_data/ETRS-003_1.fastq.gz | sed -n '4~4p' | awk '{for(i=1;i<=length($0);i++){q[substr($0,i,1)]++}} END{for (k in q) print k,q[k]}'
zcat raw_data/ETRS-003_2.fastq.gz | sed -n '4~4p' | awk '{for(i=1;i<=length($0);i++){q[substr($0,i,1)]++}} END{for (k in q) print k,q[k]}'
```
<details>
<summary>📊 Quality Score Base Counting</summary>

Counts the occurrence of each quality score character in FASTQ files.  

- 🔄 Reads R1 and R2 FASTQ files using `zcat`  
- 📝 `sed -n '4~4p'` → selects every 4th line (quality score lines)  
- 🧮 `awk` → loops over each character in quality lines to tally counts  
- 🖨 Prints counts per quality score character for each file  

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

for R1 in *_1.fastq.gz; do
    [[ -f "$R1" ]] || continue
    TOTAL_COUNT=$((TOTAL_COUNT+1))
    SAMPLE=${R1%_1.fastq.gz}

    if [[ -f "${SAMPLE}_2.fastq.gz" ]]; then
        echo "✅ $SAMPLE — paired"
        PAIRED_COUNT=$((PAIRED_COUNT+1))
    else
        echo "❌ $SAMPLE — missing _2.fastq.gz file"
        MISSING=true
    fi
done

echo -e "\nTotal samples checked: $TOTAL_COUNT"
echo "Correctly paired samples: $PAIRED_COUNT"
$MISSING && echo "⚠ Some samples are missing pairs. Fix before running fastp." || echo "✅ All FASTQ files are correctly paired."

```
<details>
<summary>🔗 FASTQ Pairing Check Script</summary>

Verifies that all R1 FASTQ files have corresponding R2 files.  

- 🛡 `set -euo pipefail` → safe script execution  
- 📂 `INDIR="raw_data"` → input FASTQ folder  
- 🔄 Checks if current directory is `raw_data`; cd if not  
- 🔍 Loops through all `_1.fastq.gz` files  
- 🏷 Extracts sample name: `${R1%_1.fastq.gz}`  
- ✅ Prints `"paired"` if R2 exists; ❌ prints `"missing"` otherwise  
- 🔢 Tracks total and paired sample counts  
- ⚠ Summarizes results: missing pairs warning or all paired confirmation  

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

for R1 in "$FASTQ_DIR"/*_1.fastq.gz; do
    [[ -f "$R1" ]] || continue
    SAMPLE=$(basename "$R1" | sed -E 's/_1\.fastq\.gz//')
    R2="$FASTQ_DIR/${SAMPLE}_2.fastq.gz"

    if [[ -f "$R2" ]]; then
        echo "Processing sample $SAMPLE"

        calc_stats() {
            zcat "$1" | awk 'NR%4==2 {len=length($0); sum+=len; if(min==""){min=len}; if(len<min){min=len}; if(len>max){max=len}; count++} END{avg=sum/count; printf "%d,%d,%.2f", min, max, avg}'
        }

        STATS_R1=$(calc_stats "$R1")
        STATS_R2=$(calc_stats "$R2")

        echo "$SAMPLE,$STATS_R1,$STATS_R2" >> "$OUTPUT_CSV"
    else
        echo "⚠ Missing _2.fastq.gz for $SAMPLE, skipping."
    fi
done

echo "✅ Read length summary saved to $OUTPUT_CSV"

```
<details>
<summary>📏 FASTQ Read Length Summary Script</summary>

Calculates min, max, and average read lengths for paired-end FASTQ files.  

- 🛡 `set -euo pipefail` → safe script execution  
- 📂 `FASTQ_DIR="raw_data"` → input FASTQ folder  
- 📂 `OUTDIR="csv_output"` → output folder; auto-created  
- 📝 `OUTPUT_CSV="read_length_summary.csv"` → stores read length stats  
- 🔄 Loops through all `_1.fastq.gz` files  
- 🏷 Extracts sample name: `basename ... | sed -E 's/_1\.fastq\.gz//'`  
- 🔗 Finds corresponding R2 file; skips if missing  
- 🧮 `calc_stats()` → calculates min, max, average read lengths using `awk`  
- ➕ Appends results to CSV: `Sample,R1_min,R1_max,R1_avg,R2_min,R2_max,R2_avg`  
- ✅ Prints completion message when finished  

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

for R1 in "$INDIR"/*_1.fastq.gz; do
    [[ -f "$R1" ]] || continue

    SAMPLE=$(basename "$R1" | sed -E 's/_1\.fastq\.gz//')
    R2="$INDIR/${SAMPLE}_2.fastq.gz"

    if [[ ! -f "$R2" ]]; then
        echo "⚠ No _2.fastq.gz file found for $SAMPLE — skipping."
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
<summary>✂️ FASTQ Trimming with fastp Script</summary>

Trims paired-end FASTQ files using `fastp` and generates JSON/HTML reports.  

- 🛡 `set -euo pipefail` → safe script execution  
- 📂 `INDIR="raw_data"` → input FASTQ folder  
- 📂 `OUTDIR="fastp_results_min_50"` → trimmed output folder; auto-created  
- 🔄 Loops through `_1.fastq.gz` files  
- 🏷 Extracts sample name and finds corresponding `_2.fastq.gz`; skips missing or already processed samples  
- 🧵 Determines threads: `FASTP_THREADS=$((nproc / 2))`  
- ✅ `run_fastp()` → runs `fastp` with:  
  - min read length 50 (`--length_required 50`)  
  - quality Phred ≥20 (`--qualified_quality_phred 20`)  
  - paired-end adapter detection (`--detect_adapter_for_pe`)  
  - outputs: trimmed FASTQ, JSON, HTML, and log files  
- ⚡ Uses GNU `parallel` to process multiple samples concurrently  
- 🎉 Prints completion message with total processed samples  

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
echo "🔹 First 10 quality lines from ETRS-003_1 (R1):"
zcat fastp_results_min_50/ETRS-003_1.trim.fastq.gz \
| sed -n '4~4p' \
| head -n 10 \
| awk '{print "✅ " $0}'
```
 Show first 10 quality lines from R2
```bash
echo "🔹 First 10 quality lines from ETRS-003_2 (R2):"
zcat fastp_results_min_50/ETRS-003_2.trim.fastq.gz \
| sed -n '4~4p' \
| head -n 10 \
| awk '{print "✅ " $0}'

```
<details>
<summary>🔍 Preview Trimmed FASTQ Quality Lines</summary>

Displays the first 10 quality score lines from trimmed FASTQ files.  

- 🔹 Reads trimmed R1 and R2 FASTQ files (`fastp_results_min_50/*.trim.fastq.gz`)  
- 📝 `sed -n '4~4p'` → selects every 4th line (quality scores)  
- 👀 `head -n 10` → shows only the first 10 lines  
- ✅ `awk '{print "✅ " $0}'` → adds checkmark for visualization  
- 🖨 Prints output to terminal with clear labels for R1 and R2  

</details>


Count ASCII characters in quality lines:
```bash
    Count base composition in R1
```bash
zcat fastp_results_min_50/ETRS-003_1.trim.fastq.gz \
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
zcat fastp_results_min_50/ETRS-003_2.trim.fastq.gz \
| sed -n '4~4p' \
| awk '{
    for(i=1;i<=length($0);i++){ q[substr($0,i,1)]++ }
} END {
    for (k in q) print k, q[k]
}' \
| awk '{print "✅ Base " $1 ": " $2 " occurrences"}'
```
<details>
<summary>📊 Quality Score Base Counts (Trimmed FASTQ)</summary>

Counts occurrences of each quality score character in trimmed FASTQ files.  

- 🔹 Reads trimmed R1 and R2 FASTQ files (`fastp_results_min_50/*.trim.fastq.gz`)  
- 📝 `sed -n '4~4p'` → selects every 4th line (quality score lines)  
- 🧮 `awk` → loops over each character to tally occurrences  
- ✅ `awk '{print "✅ Base " $1 ": " $2 " occurrences"}'` → formats counts for clear visualization  
- 🖨 Prints counts per quality score for R1 and R2  

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

for R1 in "$INDIR"/*_1.trim.fastq.gz; do
    [[ -f "$R1" ]] || continue
    SAMPLE=$(basename "$R1" | sed -E 's/_1\.trim\.fastq\.gz//')
    R2="$INDIR/${SAMPLE}_2.trim.fastq.gz"

    if [[ -f "$R2" ]]; then
        R1_COUNT=$(( $(zcat "$R1" | wc -l) / 4 ))
        R2_COUNT=$(( $(zcat "$R2" | wc -l) / 4 ))
    else
        R1_COUNT=$(( $(zcat "$R1" | wc -l) / 4 ))
        R2_COUNT="NA"
        echo "⚠ Missing _2 file for $SAMPLE"
    fi

    echo "$SAMPLE,$R1_COUNT,$R2_COUNT" >> "$OUTFILE"
    echo "✅ $SAMPLE → R1: $R1_COUNT | R2: $R2_COUNT"
done

echo "🎉 All done! Read counts saved to '$OUTFILE'"

```
<details>
<summary>📊 Trimmed FASTQ Read Counts</summary>

Counts reads in paired trimmed FASTQ files and outputs a CSV.  

- 🛡 `set -euo pipefail` → safe script execution  
- 📂 `INDIR="fastp_results_min_50"` → folder with trimmed FASTQ files  
- 📂 `OUTDIR="csv_output"` → CSV output folder; auto-created  
- 📝 `OUTFILE="trimmed_read_counts.csv"` → stores read counts per sample  
- 🔄 Loops through `_1.trim.fastq.gz` files  
- 🏷 Extracts sample name: `basename ... | sed -E 's/_1\.trim\.fastq\.gz//'`  
- 🔗 Finds corresponding R2 file; if missing, R2 count = "NA"  
- 🧮 Counts reads: `$(zcat file | wc -l) / 4`  
- ➕ Appends counts to CSV: `Sample,R1_reads,R2_reads`  
- ✅ Prints progress for each sample  
- 🎉 Prints completion message when done  

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
<details>
<summary>✂️ Read Trimming Summary</summary>

Calculates the number of reads removed during trimming for a sample.  

- 🔹 Counts reads in raw FASTQ files: `_1.fastq.gz` and `_2.fastq.gz`  
- 🔹 Counts reads in trimmed FASTQ files: `_1.trim.fastq.gz` and `_2.trim.fastq.gz`  
- 🧮 Calculates trimmed reads: `RAW_COUNT - TRIM_COUNT` for R1 and R2  
- 🖨 Prints the number of reads removed during trimming  

</details>


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
    [[ -f "$R1" ]] || continue
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
        echo "⚠ Missing _2 file for $SAMPLE, skipping."
    fi
done

echo "✅ Trimmed read length summary saved to $OUTPUT

```
<details>
<summary>📏 Trimmed FASTQ Read Length Summary</summary>

Calculates min, max, and average read lengths for paired trimmed FASTQ files.  

- 🛡 `set -euo pipefail` → safe script execution  
- 📂 `FASTQ_DIR="fastp_results_min_50"` → folder with trimmed FASTQ files  
- 📂 `OUTDIR="csv_output"` → output folder; auto-created  
- 📝 `OUTPUT_CSV="trimmed_read_length_summary.csv"` → stores read length stats  
- 🔄 Loops through `_1.trim.fastq.gz` files  
- 🏷 Extracts sample name: `basename "$R1" _1.trim.fastq.gz`  
- 🔗 Finds corresponding R2 file; skips sample if missing  
- 🧮 `calc_stats()` → calculates min, max, average read lengths using `awk`  
- ➕ Appends results to CSV: `Sample,R1_min,R1_max,R1_avg,R2_min,R2_max,R2_avg`  
- ✅ Prints completion message when done  

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
<summary>📑 Generate MultiQC Report</summary>

Generates a MultiQC report from trimmed FASTQ files.  

- 🛡 `set -euo pipefail` → safe script execution  
- 📂 `INPUT_DIR="fastp_results_min_50"` → folder with trimmed FASTQ files  
- 📂 `OUTPUT_DIR="multiqc/fastp_multiqc"` → folder to store MultiQC report; auto-created  
- ⚠ Checks if input directory exists; exits if not  
- 🔄 Runs `multiqc` on the input folder and outputs to the specified directory  
- 🖨 Prints confirmation message with report location  

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
<summary>🦠 TBProfiler Sample Profiling</summary>

Runs TBProfiler on paired-end FASTQ files to profile Mycobacterium tuberculosis samples.  

- 🛡 `set -euo pipefail` → safe script execution  
- 📂 `FASTQ_DIR="raw_data"` → folder with input FASTQ files  
- 🔄 Loops through `_1.fastq.gz` files and finds corresponding `_2.fastq.gz`  
- ❌ Skips samples if paired R2 file is missing  
- ▶️ Runs `tb-profiler profile` with:  
  - 8 threads (`--threads 8`)  
  - output prefix as sample name (`--prefix`)  
  - TXT report (`--txt`)  
  - spoligotype analysis (`--spoligotype`)  
- ✅ Prints progress for each sample  
- 📌 Prints completion message when all samples are processed  

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
<summary>📑 Collate TBProfiler Results</summary>

Merges all `tbprofiler.txt` outputs into a single CSV file.  

- 🔍 Finds all `tbprofiler.txt` files recursively in the current directory  
- 🧮 Uses `awk` to skip repeated headers and combine contents  
- 🔄 Replaces tabs with commas (`sed 's/\t/,/g'`)  
- ➕ Writes output to `tbprofiler_collated.csv`  
- ✅ Prints confirmation message after collation  

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
<summary>🧬 Snippy Variant Calling Pipeline</summary>

Performs variant calling on paired trimmed FASTQ files against a reference genome using Snippy.  

- 🛡 `set -euo pipefail` → safe script execution  
- 📂 Input/output directories:  
  - `FASTP_DIR="fastp_results_min_50"` → trimmed FASTQ files  
  - `OUTDIR="snippy_results"` → stores Snippy outputs  
- 🧩 Reference genome: `REF="H37Rv.fasta"`  
- 🧵 Threads: `THREADS=8`, `BWA_THREADS=30`; parallel jobs: `JOBS=4`  
- 🔄 `run_snippy_sample()` → for each sample:  
  - Checks paired R1/R2 files  
  - Runs Snippy with specified threads and BWA options  
  - Moves outputs (`.vcf`, `.consensus.fa`, `.bam`, `.bam.bai`, `.snps.tab`) to output folder  
  - Cleans temporary directories  
  - Confirms if full VCF generated  
- ⚡ Uses GNU `parallel` to process multiple samples concurrently  
- 🔍 Compares FASTQ sample list with Snippy outputs and reports missing samples  
- 🎯 Prints completion message and output directory  

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
<details>
<summary>📝 Snippy VCF Integrity Check</summary>

Verifies that Snippy-generated VCF files contain the mandatory `#CHROM` header line.  

- 📂 `OUTDIR="snippy_results"` → folder with Snippy VCF files  
- 🔄 Loops through all `.vcf` files in the folder  
- ✅ Prints confirmation if `#CHROM` header line exists  
- ⚠ Warns if VCF is missing the header line  

</details>

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
<summary>📊 Qualimap BAM QC</summary>

Performs quality control on BAM files using Qualimap and generates PDF/HTML reports.  

- 🛡 `set -euo pipefail` → safe script execution  
- 📂 `SNIPPY_DIR="snippy_results"` → folder with BAM files  
- 📂 `QUALIMAP_OUT="qualimap_reports"` → output folder for QC reports; auto-created  
- 🔄 Loops through all `.bam` files  
- 🏷 Extracts sample name from BAM file  
- 🧩 Creates per-sample output directory  
- ⚡ Runs `qualimap bamqc` with:  
  - PDF and HTML output (`-outformat pdf:html`)  
  - Java memory allocation (`--java-mem-size=4G`)  
- 🖨 Prints progress for each sample  

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
<summary>📑 Generate MultiQC Report for Qualimap</summary>

Generates a MultiQC report from per-sample Qualimap BAM QC results.  

- 🛡 `set -euo pipefail` → safe script execution  
- 📂 `INPUT_DIR="qualimap_reports"` → folder with Qualimap per-sample reports  
- 📂 `OUTPUT_DIR="multiqc/qualimap_multiqc"` → folder to store MultiQC report; auto-created  
- ⚠ Checks if input directory exists; exits if not  
- 🔄 Runs `multiqc` on input folder and outputs to the specified directory  
- 🖨 Prints confirmation message with report location  

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
<summary>🧬 TB Variant Filtering</summary>

Filters Snippy-generated VCF files using `tb_variant_filter` with a specified region filter.  

- 🛡 `set -euo pipefail` → safe script execution  
- 📂 `SNIPPY_DIR="snippy_results"` → folder with Snippy VCF files  
- 📂 `OUTDIR="tb_variant_filter_results"` → output folder for filtered VCFs; auto-created  
- 🏷 `REGION_FILTER="farhat_rlc"` → region filter applied to variants  
- 🔄 Loops through all `.vcf` files  
- ⚡ Runs `tb_variant_filter --region_filter` for each VCF and saves as `.filtered.vcf`  
- ⚠ Warns if filtered VCF is empty  
- ✅ Prints completion message with output directory  

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
<details>
<summary>📊 Variant Filter Summary</summary>

Generates a CSV summary comparing unfiltered and filtered TB VCF variants.  

- 🛡 `set -euo pipefail` → safe script execution  
- 📂 Input directories:  
  - `SNIPPY_DIR="snippy_results"` → unfiltered VCFs  
  - `FILTERED_DIR="tb_variant_filter_results"` → filtered VCFs  
- 📂 `OUTDIR="csv_output"` → stores summary CSV (`variant_filter_summary.csv`)  
- ⚡ Filters variants based on:  
  - Minimum depth: `MIN_DP=20`  
  - Minimum quality: `MIN_QUAL=30`  
- 🧮 `count_pass_variants()` → counts variants meeting depth & quality thresholds  
- 🔄 Loops through samples to calculate:  
  - Total variants (unfiltered/filtered)  
  - PASS variants (unfiltered/filtered)  
  - PASS retention ratio  
- 🖨 Appends results to CSV and prints per sample  

</details>

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
<details>
<summary>🗄️ Compress & Index VCFs</summary>

Compresses and indexes filtered VCF files for downstream analysis.  

- 🔄 Loops through all `.vcf` files in `tb_variant_filter_results/`  
- ⚠ Skips empty VCF files  
- 🧩 Compresses VCFs using `bgzip` → creates `.vcf.gz`  
- 📂 Indexes compressed VCFs with `bcftools index`  
- 🖨 Prints progress messages for each file  

</details>

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
<summary>🧬 Generate Consensus Sequences</summary>

Generates per-sample consensus FASTA sequences from filtered VCFs.  

- 🛡 `set -euo pipefail` → safe script execution  
- 📂 `VCFDIR="tb_variant_filter_results"` → folder with filtered VCFs  
- 📂 `OUTDIR="consensus_sequences"` → stores consensus FASTA files; auto-created  
- 🔄 Loops through all `.vcf` files  
- ⚠ Skips empty VCFs  
- 🗄 Compresses VCFs using `bgzip` and indexes with `bcftools index`  
- 🧩 Generates consensus sequence using `bcftools consensus` and reference `H37Rv.fasta`  
- 🏷 Replaces FASTA header with sample name and saves output  
- 🖨 Prints progress messages for each sample  

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
<summary>✏️ Rename Consensus FASTA Files</summary>

Simplifies filenames of consensus sequences by removing `.filtered.consensus` from the name.  

- 📂 `FASTA_DIR="consensus_sequences"` → folder with consensus FASTA files  
- 🔄 Loops through all `.filtered.consensus.fasta` files  
- ✂ Renames files to remove `.filtered.consensus` → results in `.fasta`  
- 🖨 Prints confirmation message after renaming  

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
<summary>📝 Update FASTA Headers</summary>

Updates FASTA headers so each sequence header matches its sample name.  

- 📂 `FASTA_DIR="consensus_sequences"` → folder with consensus FASTA files  
- 🔄 Loops through all `.fasta` files  
- ✂ Replaces header (`>`) with the sample name extracted from the filename  
- 🗄 Saves changes by overwriting original files  
- 🖨 Prints confirmation for each file and a final success message  

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
<details>
<summary>📏 Report Consensus Sequence Lengths</summary>

Calculates and prints the length (in base pairs) of each consensus FASTA file.  

- 📂 `consensus_sequences/*.fasta` → folder with consensus FASTA files  
- 🔄 Loops through all FASTA files  
- ✂ Strips headers (`>`) and concatenates sequence lines  
- 🧮 Counts total bases using `wc -c`  
- 🖨 Prints sample name and sequence length in bp  

</details>

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
<summary>🧮 Export Consensus Genome Lengths</summary>

Calculates the length of each consensus FASTA sequence and saves results to a CSV.  

- 📂 `FASTA_DIR="consensus_sequences"` → folder with consensus FASTA files  
- 📂 `OUTDIR="csv_output"` → folder for CSV output; auto-created  
- 🗄 `OUTPUT_CSV="csv_output/consensus_lengths.csv"` → output CSV file  
- 🔄 Loops through all FASTA files  
- ✂ Strips headers (`>`) and counts bases  
- 🖨 Appends sample name and length (bp) to CSV  
- ✅ Prints confirmation message with CSV location  

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
<details>
<summary>🧬 Align Consensus Sequences with MAFFT</summary>

Performs multiple sequence alignment of consensus genomes using MAFFT.  

- 📂 `consensus_sequences/all_consensus.fasta` → input FASTA with all consensus sequences  
- 📂 `mafft_results` → output folder for aligned sequences; auto-created  
- ⚡ Runs `mafft --auto --parttree` for optimized alignment  
- 🗄 Saves aligned sequences to `mafft_results/aligned_consensus.fasta`  

</details>

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
<summary>🔍 Check MAFFT Alignment Quality</summary>

Performs basic QC on the aligned consensus sequences to assess length consistency.  

- 📂 `FILE="mafft_results/aligned_consensus.fasta"` → input aligned FASTA  
- ⚠ Exits if file does not exist  
- 🧮 Counts sequences, total bases, minimum, maximum, and average sequence lengths  
- 🖨 Prints summary stats per alignment  
- ✅ Checks if all sequences have equal length; flags potential misalignment if lengths vary  

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
<details>
<summary>📏 Aligned Sequence Length Summary</summary>

Summarizes the length distribution of sequences in an aligned FASTA.  

- 📂 `mafft_results/aligned_consensus.fasta` → input aligned FASTA  
- 🔄 Computes the length of each sequence  
- 📊 Counts unique lengths and number of sequences per length  
- ⚠ Marks sequences shorter than 60 bp with a warning  
- ✅ Prints summary with sequence counts and lengths  

</details>

G. Use AMAS (Alignment Manipulation and Summary)
AMAS is a Python tool to summarize alignments:
```bash
conda activate amas_env
AMAS.py summary -f fasta -d dna -i mafft_results/aligned_consensus.fasta
```
<details>
<summary>📊 AMAS Alignment Summary</summary>

Generates a summary of the aligned consensus sequences using AMAS.  

- 🐍 `conda activate amas_env` → activates the AMAS environment  
- ⚡ `AMAS.py summary` → summarizes alignment statistics  
- 🔹 Options:  
  - `-f fasta` → input format is FASTA  
  - `-d dna` → sequence type is DNA  
  - `-i mafft_results/aligned_consensus.fasta` → input aligned FASTA file  
- 🖨 Outputs alignment summary including sequence counts, lengths, and base composition  

</details>

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
<summary>🌳 Build Phylogenetic Tree with IQ-TREE</summary>

Constructs a maximum likelihood phylogenetic tree from aligned consensus sequences.  

- 📂 `mafft_results/aligned_consensus.fasta` → input aligned FASTA  
- 📂 `iqtree_results` → output folder; auto-created  
- ⚡ `iqtree2` options:  
  - `-m GTR+G` → substitution model (GTR with gamma rate heterogeneity)  
  - `-bb 1000` → 1000 ultrafast bootstrap replicates  
  - `-nt 4` → use 4 threads  
  - `-o SRR10828835` → designate outgroup  
  - `-pre iqtree_results/aligned_consensus` → prefix for output files  
- 🖨 Generates tree files and bootstrap support values  

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

