#!/bin/bash
set -uo pipefail   # no -e → do not exit on error

PIPELINE_DIR="$(pwd)"
PROJECT_DIR="$(dirname "$PIPELINE_DIR")"

ADAPTERS="$PIPELINE_DIR/Trimmomatic-0.39/adapters/TruSeq3-PE.fa"
TRIMMOMATIC_JAR="$PIPELINE_DIR/Trimmomatic-0.39/trimmomatic-0.39.jar"

THREADS=${THREADS:-4}
MINLEN=36

[[ -f "$ADAPTERS" ]] || { echo "❌ Adapter file missing"; exit 1; }
[[ -f "$TRIMMOMATIC_JAR" ]] || { echo "❌ Trimmomatic JAR missing"; exit 1; }

echo "🔍 Project directory:  $PROJECT_DIR"
echo

# ----------- Helper functions -----------

find_pair() {
  local dir="$1"
  R1=$(ls "$dir"/*_1.fastq.gz 2>/dev/null | head -n1 || true)
  R2=$(ls "$dir"/*_2.fastq.gz 2>/dev/null | head -n1 || true)
  printf "%s\n%s\n" "$R1" "$R2"
}

check_gzip() {
  gzip -t "$1" 2>/dev/null
  return $?
}

# ----------- Main Loop -----------

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
    echo "❌ FASTQ pair missing — skipping"
    continue
  fi

  # -------- Corruption Check --------
  if ! check_gzip "$R1"; then
    echo "❌ Corrupt file: $R1 — skipping sample"
    continue
  fi

  if ! check_gzip "$R2"; then
    echo "❌ Corrupt file: $R2 — skipping sample"
    continue
  fi

  echo "✔ Raw FASTQ integrity OK"

  trimmed_dir="$sample_dir/trimmed"
  mkdir -p "$trimmed_dir"

  out_paired_R1="$trimmed_dir/${sample}_1_paired.fastq.gz"
  out_paired_R2="$trimmed_dir/${sample}_2_paired.fastq.gz"

  # =====================================================
  # 🔥 SMART SKIP CHECK (NEW)
  # If trimmed files exist AND are valid → skip sample
  # =====================================================

  if [[ -f "$out_paired_R1" && -f "$out_paired_R2" ]]; then
    if check_gzip "$out_paired_R1" && check_gzip "$out_paired_R2"; then
      echo "⏭ Trimmed files already exist and valid — skipping entire sample"
      continue
    else
      echo "⚠ Trimmed files corrupt — reprocessing sample"
      rm -f "$out_paired_R1" "$out_paired_R2"
    fi
  fi

  # -------- RAW FastQC --------
  raw_qc_dir="$sample_dir/fastqc_raw"
  mkdir -p "$raw_qc_dir"

  raw_qc_html="$raw_qc_dir/$(basename "$R1" .fastq.gz)_fastqc.html"

  if [[ ! -f "$raw_qc_html" ]]; then
    echo "🔍 Running FastQC on RAW..."
    fastqc -t "$THREADS" "$R1" "$R2" --outdir "$raw_qc_dir" || {
      echo "❌ FastQC failed — skipping sample"
      continue
    }
  else
    echo "⏭ Raw FastQC already present"
  fi

  # -------- Trimming --------
  echo "✂️ Running Trimmomatic..."

  java -jar "$TRIMMOMATIC_JAR" PE -phred33 \
    -threads "$THREADS" \
    "$R1" "$R2" \
    "$trimmed_dir/${sample}_1_paired.fastq.gz" "$trimmed_dir/${sample}_1_unpaired.fastq.gz" \
    "$trimmed_dir/${sample}_2_paired.fastq.gz" "$trimmed_dir/${sample}_2_unpaired.fastq.gz" \
    ILLUMINACLIP:"$ADAPTERS":2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:$MINLEN || {
      echo "❌ Trimming failed — skipping sample"
      continue
    }

  echo "✔ Trimming completed"

  # -------- Trimmed FastQC --------
  trimmed_qc_html="$trimmed_dir/${sample}_1_paired_fastqc.html"

  if [[ ! -f "$trimmed_qc_html" ]]; then
    echo "🔍 Running FastQC on trimmed reads..."
    fastqc -t "$THREADS" "$out_paired_R1" "$out_paired_R2" --outdir "$trimmed_dir" || {
      echo "❌ Trimmed FastQC failed — skipping sample"
      continue
    }
  else
    echo "⏭ Trimmed FastQC already present"
  fi

  echo "✅ Sample completed successfully"

done

echo
echo "🎉 Pipeline finished — resume-safe + corruption-safe mode"
