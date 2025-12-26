#!/bin/bash

set -e

echo "========================================"
echo " RNA-Seq Pipeline: Tool Setup"
echo "========================================"

# Always run from pipeline directory
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PIPELINE_DIR"

echo "📁 Pipeline directory: $PIPELINE_DIR"

# -------------------------
# Helper functions
# -------------------------
command_exists () {
    command -v "$1" &> /dev/null
}

python_pkg_exists () {
    python3 - <<EOF
import importlib.util
exit(0) if importlib.util.find_spec("$1") else exit(1)
EOF
}

# -------------------------
# System update
# -------------------------
echo "🔄 Updating system package list..."
sudo apt update -y

# -------------------------
# Core dependencies
# -------------------------
CORE_PKGS=(
    build-essential
    wget
    unzip
    git
    curl
    default-jre
    python3
    python3-pip
)

echo "📦 Checking core system dependencies..."
for pkg in "${CORE_PKGS[@]}"; do
    if dpkg -s "$pkg" &> /dev/null; then
        echo "✅ $pkg already installed"
    else
        echo "⬇️ Installing $pkg"
        sudo apt install -y "$pkg"
    fi
done

# -------------------------
# Bioinformatics tools
# -------------------------
echo "🧬 Checking bioinformatics tools..."

TOOLS=(
    hisat2
    samtools
    fastqc
)

for tool in "${TOOLS[@]}"; do
    if command_exists "$tool"; then
        echo "✅ $tool already installed"
    else
        echo "⬇️ Installing $tool"
        sudo apt install -y "$tool"
    fi
done

# -------------------------
# MultiQC (pip-based)
# -------------------------
if command_exists multiqc; then
    echo "✅ MultiQC already installed"
else
    echo "⬇️ Installing MultiQC"
    pip3 install --user multiqc
fi

# -------------------------
# Python libraries for pipeline
# -------------------------
PY_PKGS=(
    pandas
    numpy
)

echo "🐍 Checking Python libraries..."
for pkg in "${PY_PKGS[@]}"; do
    if python_pkg_exists "$pkg"; then
        echo "✅ Python package '$pkg' present"
    else
        echo "⬇️ Installing Python package '$pkg'"
        pip3 install --user "$pkg"
    fi
done

# -------------------------
# Final verification
# -------------------------
echo ""
echo "========================================"
echo " Tool Verification"
echo "========================================"

hisat2 --version | head -n 1
samtools --version | head -n 1
fastqc --version
multiqc --version
python3 --version
java -version

echo "========================================"
echo " ✅ All tools ready. Pipeline can run!"
echo "========================================"

