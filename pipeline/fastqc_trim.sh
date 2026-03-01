#!/bin/bash
set -euo pipefail

# -------------------------
# Paths and variables
# -------------------------
PIPELINE_DIR="$(pwd)"
PROJECT_DIR="$(dirname "$PIPELINE_DIR")"

ADAPTERS="$PIPELINE_DIR/Trimmomatic-0.39/adapters/TruSeq3-PE.fa"
TRIMMOMATIC_JAR="$PIPELINE_DIR/Trimmomatic-0.39/trimmomatic-0.39.jar"

THREADS=${THREADS:-4}
MINLEN=36

# -------------------------
# Sanity checks
# -------------------------
[[ -f "$ADAPTERS" ]] || { echo "❌ Adapter file missing"; exit 1; }
[[ -f "$TRIMMOMATIC_JAR" ]] || { echo "❌ Trimmomatic JAR missing"; exit 1; }

echo "🔍 Pipeline directory: $PIPELINE_DIR"
echo "🔍 Project directory:  $PROJECT_DIR"
echo "🔍 Threads:            $THREADS"
echo

# -------------------------
# Function: detect R1/R2
# -------------------------
find_pair() {
  local dir="$1"
  R1=$(ls "$dir"/*_1.fastq.gz 2>/dev/null | head -n1 || true)
  R2=$(ls "$dir"/*_2.fastq.gz 2>/dev/null | head -n1 || true)
  printf "%s\n%s\n" "$R1" "$R2"
}

# -------------------------
# Loop through samples
# -------------------------
for sample_dir in "$PROJECT_DIR"/*/; do
  sample=$(basename "$sample_dir")

  case "$sample" in
    pipeline|hisat2_index|genome_data|Trimmomatic-0.39)
      continue ;;
  esac

  echo "----------------------------------------------"
  echo "📁 Sample: $sample"

  readpair=$(find_pair "$sample_dir")
  R1=$(echo "$readpair" | sed -n '1p')
  R2=$(echo "$readpair" | sed -n '2p')

  if [[ -z "$R1" || -z "$R2" ]]; then
    echo "❌ No FASTQ pair found — skipping"
    continue
  fi

  echo "📌 Raw R1: $R1"
  echo "📌 Raw R2: $R2"

  # -------------------------
  # STEP 1: FastQC on RAW
  # -------------------------
  raw_qc_dir="$sample_dir/fastqc_raw"
  mkdir -p "$raw_qc_dir"

  raw_qc_html="$raw_qc_dir/$(basename "$R1" .fastq.gz)_fastqc.html"

  if [[ -f "$raw_qc_html" ]]; then
    echo "⏭ Raw FastQC already present — skipping FastQC"
  else
    echo "🔍 Running FastQC on RAW reads..."
    fastqc -t "$THREADS" "$R1" "$R2" --outdir "$raw_qc_dir"
    echo "✔ Raw FastQC completed"
  fi

  # -------------------------
  # STEP 2: Trimming
  # -------------------------
  trimmed_dir="$sample_dir/trimmed"
  mkdir -p "$trimmed_dir"

  out_paired_R1="$trimmed_dir/${sample}_1_paired.fastq.gz"
  out_paired_R2="$trimmed_dir/${sample}_2_paired.fastq.gz"

  if [[ -f "$out_paired_R1" && -f "$out_paired_R2" ]]; then
    echo "⏭ Trimmed files already present — skipping trimming"
  else
    echo "✂️ Running Trimmomatic..."

    java -jar "$TRIMMOMATIC_JAR" PE -phred33 \
      -threads "$THREADS" \
      "$R1" "$R2" \
      "$trimmed_dir/${sample}_1_paired.fastq.gz" "$trimmed_dir/${sample}_1_unpaired.fastq.gz" \
      "$trimmed_dir/${sample}_2_paired.fastq.gz" "$trimmed_dir/${sample}_2_unpaired.fastq.gz" \
      ILLUMINACLIP:"$ADAPTERS":2:30:10 \
      LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:$MINLEN

    echo "✔ Trimming completed"
  fi

  # -------------------------
  # STEP 3: FastQC on TRIMMED
  # -------------------------
  trimmed_qc_html="$trimmed_dir/$(basename "$out_paired_R1" .fastq.gz)_fastqc.html"

  if [[ -f "$trimmed_qc_html" ]]; then
    echo "⏭ Trimmed FastQC already present — skipping"
  else
    echo "🔍 Running FastQC on trimmed reads..."
    fastqc -t "$THREADS" "$out_paired_R1" "$out_paired_R2" --outdir "$trimmed_dir"
    echo "✔ Trimmed FastQC completed"
  fi

done

echo
echo "🎉 Pipeline finished (resume-safe mode)"
