#!/usr/bin/env python3
import os
import pandas as pd

# Automatically detect the base project directory (one level up from this script)
BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

print(f"🔍 Base project directory detected: {BASE_DIR}")

# List only numeric sample folder names (757, 758, 759, 760...)
sample_dirs = [d for d in os.listdir(BASE_DIR) if d.isdigit()]

if not sample_dirs:
    print("❌ No sample folders detected.")
    exit(1)

count_tables = []

for sample in sample_dirs:
    sample_path = os.path.join(BASE_DIR, sample)
    
    # Automatically detect HTSeq folder (any folder starting with 'htseq')
    htseq_folder = None
    for item in os.listdir(sample_path):
        if os.path.isdir(os.path.join(sample_path, item)) and item.lower().startswith("htseq"):
            htseq_folder = os.path.join(sample_path, item)
            break
    
    if htseq_folder is None:
        print(f"⚠ No HTSeq folder found for sample {sample}")
        continue

    # Automatically detect count file (any txt inside the HTSeq folder)
    count_file = None
    for f in os.listdir(htseq_folder):
        if f.endswith(".txt") and sample in f:
            count_file = os.path.join(htseq_folder, f)
            break

    if count_file and os.path.exists(count_file):
        print(f"Reading {count_file}")
        df = pd.read_csv(count_file, sep="\t", header=None, names=["Gene", sample])
        count_tables.append(df)
    else:
        print(f"⚠ Count file missing for sample {sample}")

# Merge all count tables on "Gene"
if count_tables:
    merged_df = count_tables[0]
    for df in count_tables[1:]:
        merged_df = merged_df.merge(df, on="Gene", how="outer")

    # Output file
    output_file = os.path.join(BASE_DIR, "All_samples_counts.xlsx")
    merged_df.to_excel(output_file, index=False)
    print(f"✓ Merged count matrix saved as: {output_file}")
else:
    print("❌ No count tables found to merge.")

