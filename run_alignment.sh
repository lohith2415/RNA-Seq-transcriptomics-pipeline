#!/bin/bash
set -euo pipefail

PIPELINE_DIR="$(pwd)"
PROJECT_DIR="$(dirname "$PIPELINE_DIR")"
GENOME_DIR="$PIPELINE_DIR/genome_data"
INDEX_DIR="$PIPELINE_DIR/hisat2_index"

echo "🔍 Project Folder: $PROJECT_DIR"
echo "🔍 Genome Folder:  $GENOME_DIR"
echo

mkdir -p "$INDEX_DIR"

#############################################
# 1. AUTO-DETECT GENOME FASTA
#############################################
FASTA=$(find "$GENOME_DIR" -maxdepth 1 -type f \( -iname "*.fa" -o -iname "*.fasta" -o -iname "*.fna" \) | head -n1 || true)

if [[ -z "$FASTA" ]]; then
    GZ_FASTA=$(find "$GENOME_DIR" -maxdepth 1 -type f \( -iname "*.fa.gz" -o -iname "*.fasta.gz" -o -iname "*.fna.gz" \) | head -n1 || true)
    if [[ -n "$GZ_FASTA" ]]; then
        echo "📦 Unzipping FASTA: $(basename "$GZ_FASTA")"
        gunzip -k "$GZ_FASTA"
        FASTA="${GZ_FASTA%.gz}"
    fi
fi

if [[ -z "$FASTA" ]]; then
    echo "❌ No FASTA (.fa / .fa.gz) found inside genome_data!"
    exit 1
fi

echo "✔ Using Genome FASTA: $(basename "$FASTA")"

#############################################
# 2. AUTO-DETECT GTF
#############################################
GTF=$(find "$GENOME_DIR" -maxdepth 1 -type f -iname "*.gtf" | head -n1 || true)

if [[ -z "$GTF" ]]; then
    GTF_GZ=$(find "$GENOME_DIR" -maxdepth 1 -type f -iname "*.gtf.gz" | head -n1 || true)
    if [[ -n "$GTF_GZ" ]]; then
        echo "📦 Unzipping GTF: $(basename "$GTF_GZ")"
        gunzip -k "$GTF_GZ"
        GTF="${GTF_GZ%.gz}"
    fi
fi

if [[ -z "$GTF" ]]; then
    GTF_TAR=$(find "$GENOME_DIR" -maxdepth 1 -type f -iname "*.tar.gz" | head -n1 || true)
    if [[ -n "$GTF_TAR" ]]; then
        echo "📦 Extracting GTF TAR: $(basename "$GTF_TAR")"
        tar -xzf "$GTF_TAR" -C "$GENOME_DIR"
        GTF=$(find "$GENOME_DIR" -maxdepth 1 -type f -iname "*.gtf" | head -n1)
    fi
fi

if [[ -z "$GTF" ]]; then
    echo "❌ No GTF found inside genome_data!"
    exit 1
fi

echo "✔ Using GTF: $(basename "$GTF")"
echo

#############################################
# 3. BUILD HISAT2 INDEX IF NOT EXISTS
#############################################
if ls "$INDEX_DIR"/*ht2 &>/dev/null; then
    echo "⏭ HISAT2 index already exists — skipping"
else
    echo "🔨 Building HISAT2 genome index..."
    hisat2-build -p 8 "$FASTA" "$INDEX_DIR/genome_index"
    echo "✔ HISAT2 index complete."
fi

#############################################
# 4. PROCESS EACH SAMPLE
#############################################
for sample_dir in "$PROJECT_DIR"/*/; do
    sample=$(basename "$sample_dir")

    case "$sample" in
        pipeline|hisat2_index|genome_data) continue ;;
    esac

    echo "---------------------------------------"
    echo "📁 SAMPLE: $sample"

    trimmed="$sample_dir/trimmed"
    hisat_out="$sample_dir/hisat2"
    mkdir -p "$hisat_out"

    R1=$(ls "$trimmed"/*_1_paired.fastq.gz 2>/dev/null | head -n1 || true)
    R2=$(ls "$trimmed"/*_2_paired.fastq.gz 2>/dev/null | head -n1 || true)

    if [[ -z "$R1" || -z "$R2" ]]; then
        echo "❌ No trimmed paired reads found — skipping sample."
        continue
    fi

    BAM="$hisat_out/${sample}.sorted.bam"

    # Skip if BAM already exists
    if [[ -f "$BAM" ]]; then
        echo "⏭ BAM already exists — skipping alignment."
    else
        SAM="$hisat_out/${sample}.sam"

        echo "🧬 Running HISAT2 for $sample..."
        hisat2 -p 8 -x "$INDEX_DIR/genome_index" -1 "$R1" -2 "$R2" -S "$SAM"

        echo "🔄 Converting SAM → Sorted BAM"
        samtools sort "$SAM" -o "$BAM"
        samtools index "$BAM"

        echo "🧹 Removing SAM to save space"
        rm -f "$SAM"

        echo "✔ Alignment complete for $sample"
    fi

    #############################################
    # 5. RUN HTSEQ-COUNT
    #############################################
    HTSEQ_DIR="$sample_dir/htseq_counts"
    mkdir -p "$HTSEQ_DIR"
    COUNT_FILE="$HTSEQ_DIR/${sample}_counts.txt"

    if [[ -f "$COUNT_FILE" ]]; then
        echo "⏭ HTSeq-count already exists — skipping."
    else
        echo "🧮 Running HTSeq-count for $sample..."
        htseq-count -f bam -r pos -s no -t exon -i gene_id "$BAM" "$GTF" > "$COUNT_FILE"
        echo "✔ HTSeq-count saved at $COUNT_FILE"
    fi
done

echo "🎉 ALL ALIGNMENTS AND COUNTING FINISHED SUCCESSFULLY!"

