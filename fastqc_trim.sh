#!/bin/bash
set -euo pipefail

# Run this inside pipeline/
PIPELINE_DIR="$(pwd)"
PROJECT_DIR="$(dirname "$PIPELINE_DIR")"

ADAPTERS="$PIPELINE_DIR/Trimmomatic-0.39/adapters/TruSeq3-PE.fa"
TRIMMOMATIC_JAR="$PIPELINE_DIR/Trimmomatic-0.39/trimmomatic-0.39.jar"

THREADS=${THREADS:-4}
MINLEN=36

# Check paths
[[ -f "$ADAPTERS" ]] || { echo "❌ Adapter file missing"; exit 1; }
[[ -f "$TRIMMOMATIC_JAR" ]] || { echo "❌ Trimmomatic JAR missing"; exit 1; }

echo "🔍 Pipeline directory: $PIPELINE_DIR"
echo "🔍 Project directory:  $PROJECT_DIR"
echo "🔍 Threads:            $THREADS"
echo

###############################################
# Function: Detect read pair automatically
###############################################
find_pair() {
  local dir="$1"
  R1=$(find "$dir" -maxdepth 1 -type f \( -iname "*1*.fastq*" -o -iname "*r1*.fastq*" -o -iname "*1*.fq*" -o -iname "*r1*.fq*" \) | head -n1 || true)
  R2=$(find "$dir" -maxdepth 1 -type f \( -iname "*2*.fastq*" -o -iname "*r2*.fastq*" -o -iname "*2*.fq*" -o -iname "*r2*.fq*" \) | head -n1 || true)
  printf "%s\n%s\n" "$R1" "$R2"
}

###############################################
# LOOP THROUGH SAMPLES
###############################################
for sample_dir in "$PROJECT_DIR"/*/; do
  sample=$(basename "$sample_dir")

  case "$sample" in
    pipeline|hisat2_index|genome_data|Trimmomatic-0.39)
      continue ;;
  esac

  echo "----------------------------------------------"
  echo "📁 Sample: $sample"

  ###############################################
  # Detect raw read pairs
  ###############################################
  readpair=$(find_pair "$sample_dir")
  R1=$(echo "$readpair" | sed -n '1p')
  R2=$(echo "$readpair" | sed -n '2p')

  if [[ -z "$R1" || -z "$R2" ]]; then
    echo "❌ No FASTQ pair found in $sample_dir"
    continue
  fi

  echo "📌 Raw R1: $R1"
  echo "📌 Raw R2: $R2"

  ###############################################
  # STEP 1: FastQC on RAW reads
  ###############################################
  raw_qc_dir="$sample_dir/fastqc_raw"
  mkdir -p "$raw_qc_dir"

  if [[ -f "$raw_qc_dir/$(basename "$R1" .fastq.gz)_fastqc.html" ]]; then
      echo "⏭ Raw FastQC already exists — skipping"
  else
      echo "🔍 Running FastQC on RAW reads..."
      fastqc -t "$THREADS" "$R1" "$R2" --outdir "$raw_qc_dir"
      echo "✔ Raw FastQC completed"
  fi

  ###############################################
  # PAUSE: Ask user
  ###############################################
  echo
  read -p "▶ Press ENTER to run trimming for sample '$sample' (or type 'n' to skip): " choice
  if [[ "$choice" =~ ^[Nn]$ ]]; then
      echo "⏭ Skipping trimming for $sample"
      continue
  fi

  ###############################################
  # STEP 2: Trimming
  ###############################################
  trimmed_dir="$sample_dir/trimmed"
  mkdir -p "$trimmed_dir"

  out_paired_R1="$trimmed_dir/${sample}_1_paired.fastq.gz"
  out_unpaired_R1="$trimmed_dir/${sample}_1_unpaired.fastq.gz"
  out_paired_R2="$trimmed_dir/${sample}_2_paired.fastq.gz"
  out_unpaired_R2="$trimmed_dir/${sample}_2_unpaired.fastq.gz"

  echo "✂️ Running Trimmomatic..."
  java -jar "$TRIMMOMATIC_JAR" PE \
    -threads "$THREADS" \
    "$R1" "$R2" \
    "$out_paired_R1" "$out_unpaired_R1" \
    "$out_paired_R2" "$out_unpaired_R2" \
    ILLUMINACLIP:"$ADAPTERS":2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:$MINLEN

  echo "✔ Trimming completed"

  ###############################################
  # STEP 3: FastQC on TRIMMED reads
  ###############################################
  echo "🔍 Running FastQC on trimmed reads..."
  fastqc -t "$THREADS" \
    "$out_paired_R1" "$out_paired_R2" \
    --outdir "$trimmed_dir"

  echo "✔ Trimmed FastQC completed"

done

echo
echo "🎉 Pipeline step completed: FastQC → optional trimming → FastQC"

