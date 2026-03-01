#!/bin/bash
set -euo pipefail

###########################################
# AUTO DETECT DIRECTORIES
###########################################

PIPELINE_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(dirname "$PIPELINE_DIR")"

GENOME_DIR="$PIPELINE_DIR/genome_data"
INDEX_DIR="$PIPELINE_DIR/hisat2_index"
INDEX_PREFIX="$INDEX_DIR/genome_index"

THREADS=8

echo "🔍 Pipeline Dir : $PIPELINE_DIR"
echo "🔍 Project Dir  : $PROJECT_DIR"
echo

mkdir -p "$INDEX_DIR"

###########################################
# 1️⃣ DETECT GENOME FASTA
###########################################

FASTA=$(find "$GENOME_DIR" -maxdepth 1 -type f \
    \( -iname "*.fa" -o -iname "*.fasta" -o -iname "*.fna" \) | head -n1 || true)

if [[ -z "$FASTA" ]]; then
    echo "❌ No FASTA found in genome_data"
    exit 1
fi

echo "✔ Genome FASTA: $(basename "$FASTA")"

###########################################
# 2️⃣ DETECT GTF
###########################################

GTF=$(find "$GENOME_DIR" -maxdepth 1 -type f -iname "*.gtf" | head -n1 || true)

if [[ -z "$GTF" ]]; then
    echo "❌ No GTF found in genome_data"
    exit 1
fi

echo "✔ GTF: $(basename "$GTF")"
echo

###########################################
# 3️⃣ BUILD INDEX IF NEEDED
###########################################

if compgen -G "$INDEX_DIR/*.ht2" > /dev/null; then
    echo "⏭ HISAT2 index already exists — skipping build"
else
    echo "🔨 Building HISAT2 index..."
    hisat2-build -p "$THREADS" "$FASTA" "$INDEX_PREFIX"
    echo "✔ Index built"
fi

echo

###########################################
# 4️⃣ PROCESS EACH SAMPLE
###########################################

for sample_dir in "$PROJECT_DIR"/SRR*/; do

    [[ -d "$sample_dir" ]] || continue

    sample=$(basename "$sample_dir")

    echo "-----------------------------------------"
    echo "📁 SAMPLE: $sample"

    trimmed_dir="$sample_dir/trimmed"
    hisat_dir="$sample_dir/hisat2"
    count_dir="$sample_dir/htseq_counts"

    mkdir -p "$hisat_dir"
    mkdir -p "$count_dir"

    R1=$(ls "$trimmed_dir"/*_1_paired.fastq.gz 2>/dev/null | head -n1 || true)
    R2=$(ls "$trimmed_dir"/*_2_paired.fastq.gz 2>/dev/null | head -n1 || true)

    if [[ -z "$R1" || -z "$R2" ]]; then
        echo "❌ No trimmed paired reads — skipping"
        continue
    fi

    BAM="$hisat_dir/${sample}.sorted.bam"
    SAM="$hisat_dir/${sample}.sam"
    COUNT_FILE="$count_dir/${sample}_counts.txt"

    ###########################################
    # HISAT2 ALIGNMENT (SKIP IF BAM EXISTS)
    ###########################################

    if [[ -f "$BAM" ]]; then
        echo "⏭ BAM already exists — skipping HISAT2"
    else
        echo "🧬 Running HISAT2..."

        hisat2 -p "$THREADS" \
            -x "$INDEX_PREFIX" \
            -1 "$R1" \
            -2 "$R2" \
            -S "$SAM"

        echo "🔄 Sorting BAM..."
        samtools sort -@ "$THREADS" -o "$BAM" "$SAM"
        samtools index "$BAM"

        rm -f "$SAM"

        echo "✔ Alignment completed"
    fi

    ###########################################
    # HTSEQ COUNT (SKIP IF COUNT EXISTS)
    ###########################################

    if [[ -f "$COUNT_FILE" ]]; then
        echo "⏭ HTSeq-count already exists — skipping"
    else
        echo "🧮 Running HTSeq-count..."

        htseq-count \
            -f bam \
            -r pos \
            -s no \
            -t exon \
            -i gene_id \
            "$BAM" \
            "$GTF" > "$COUNT_FILE"

        echo "✔ Counts generated"
    fi

done

echo
echo "🎉 FULL PIPELINE COMPLETED SUCCESSFULLY!"
