#!/bin/bash

# eboVar.sh - EBOV Variant Analysis Pipeline with full QC and custom trimming
# Usage: ./eboVar.sh -i rawreads/ -o results/ -r ebov_ref.fa -t 4

set -euo pipefail

threads=4

usage() {
  echo "Usage: $0 -i <input_dir> -o <output_dir> -r <reference.fasta> [-t <threads>]"
  echo ""
  echo "Required arguments:"
  echo "  -i, --input      Directory containing raw FASTQ files"
  echo "  -o, --output     Directory to store results"
  echo "  -r, --ref        EBOV reference genome in FASTA format"
  echo ""
  echo "Optional:"
  echo "  -t, --threads    Number of threads to use (default: 4)"
  echo "  -h, --help       Show this message"
  exit 1
}

# Parse arguments
while [[ "$#" -gt 0 ]]; do
  case $1 in
    -i|--input) input_dir="$2"; shift ;;
    -o|--output) output_dir="$2"; shift ;;
    -r|--ref) ref_fa="$2"; shift ;;
    -t|--threads) threads="$2"; shift ;;
    -h|--help) usage ;;
    *) echo "❌ Unknown parameter: $1"; usage ;;
  esac
  shift
done

# Validate input
[[ -z "${input_dir:-}" || -z "${output_dir:-}" || -z "${ref_fa:-}" ]] && usage
[[ ! -d "$input_dir" ]] && echo "❌ Input directory does not exist!" && exit 1
[[ ! -f "$ref_fa" ]] && echo "❌ Reference genome not found!" && exit 1

if ! ls "$input_dir"/*_1.fastq.gz 1> /dev/null 2>&1; then
  echo "❌ No FASTQ files matching *_1.fastq.gz found in $input_dir"
  exit 1
fi

# Define output folders
mkdir -p "$output_dir"/{qc_raw,qc_trimmed,trimmed,bam,vcf,logs,multiqc}
RAW_QC_DIR="$output_dir/qc_raw"
TRIM_QC_DIR="$output_dir/qc_trimmed"
FASTP_REPORTS="$output_dir/qc_trimmed"
TRIMMED_DIR="$output_dir/trimmed"

# Index reference genome if needed
if [ ! -f "$ref_fa.bwt" ]; then
  echo "🔄 Indexing reference genome..."
  bwa index "$ref_fa"
fi

# Process each sample
for R1 in "$input_dir"/*_1.fastq.gz; do
  sample=$(basename "$R1" | sed 's/_1\.fastq\.gz//')
  R2="$input_dir/${sample}_2.fastq.gz"
  log_file="$output_dir/logs/${sample}.log"

  echo "🚀 Processing sample: $sample" | tee "$log_file"

  # Step 1: FastQC on raw reads
  echo "📊 Step 1: FastQC on raw reads..." | tee -a "$log_file"
  fastqc -t "$threads" "$R1" "$R2" -o "$RAW_QC_DIR" >> "$log_file" 2>&1

  # Step 2: Trimming with fastp
  echo "✂️  Step 2: Trimming with fastp..." | tee -a "$log_file"
  fastp \
    -i "$R1" -I "$R2" \
    -o "$TRIMMED_DIR/${sample}_trimmed_R1.fastq.gz" \
    -O "$TRIMMED_DIR/${sample}_trimmed_R2.fastq.gz" \
    --thread "$threads" \
    --length_required 50 \
    --trim_tail1 1 \
    --trim_tail2 1 \
    --n_base_limit 0 \
    --detect_adapter_for_pe \
    -h "$FASTP_REPORTS/${sample}_fastp.html" \
    -j "$FASTP_REPORTS/${sample}_fastp.json" \
    --report_title "$sample fastp report" \
    >> "$log_file" 2>&1

  # Step 3: FastQC on trimmed reads
  echo "📊 Step 3: FastQC on trimmed reads..." | tee -a "$log_file"
  fastqc -t "$threads" \
    "$TRIMMED_DIR/${sample}_trimmed_R1.fastq.gz" \
    "$TRIMMED_DIR/${sample}_trimmed_R2.fastq.gz" \
    -o "$TRIM_QC_DIR" >> "$log_file" 2>&1

  # Step 4: Alignment with BWA
  echo "🧬 Step 4: Aligning with BWA..." | tee -a "$log_file"
  bwa mem -t "$threads" "$ref_fa" \
    "$TRIMMED_DIR/${sample}_trimmed_R1.fastq.gz" \
    "$TRIMMED_DIR/${sample}_trimmed_R2.fastq.gz" \
    | samtools view -bS - \
    | samtools sort -@ "$threads" -o "$output_dir/bam/${sample}.bam" >> "$log_file" 2>&1

  # Step 5: Index BAM
  echo "📌 Step 5: Indexing BAM..." | tee -a "$log_file"
  samtools index "$output_dir/bam/${sample}.bam" >> "$log_file" 2>&1

  # Step 6: Variant calling
  echo "🧪 Step 6: Variant calling..." | tee -a "$log_file"
  bcftools mpileup -Ou -f "$ref_fa" "$output_dir/bam/${sample}.bam" \
    | bcftools call -mv -Oz -o "$output_dir/vcf/${sample}.raw.vcf.gz" >> "$log_file" 2>&1

  # Step 7: Sort + index VCF
  echo "🗂️  Step 7: Sorting + indexing VCF..." | tee -a "$log_file"
  bcftools sort "$output_dir/vcf/${sample}.raw.vcf.gz" -Oz -o "$output_dir/vcf/${sample}.vcf.gz" >> "$log_file" 2>&1
  tabix -p vcf "$output_dir/vcf/${sample}.vcf.gz" >> "$log_file" 2>&1
  rm "$output_dir/vcf/${sample}.raw.vcf.gz" >> "$log_file" 2>&1

  echo "✅ Sample $sample complete!" | tee -a "$log_file"
done

# Step 8: Run MultiQC
echo "📊 Step 8: Running MultiQC..."
multiqc "$output_dir/qc_raw" "$output_dir/qc_trimmed" -o "$output_dir/multiqc"

echo "🎉 All samples processed successfully! Results saved to $output_dir"
