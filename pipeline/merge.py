#!/usr/bin/env python3

import os
import pandas as pd

# Detect project directory (one level above pipeline)
BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

print(f"🔍 Base project directory detected: {BASE_DIR}")

# Detect SRR sample folders
sample_dirs = [
    d for d in os.listdir(BASE_DIR)
    if d.startswith("SRR") and os.path.isdir(os.path.join(BASE_DIR, d))
]

if not sample_dirs:
    print("❌ No SRR sample folders detected.")
    exit(1)

count_tables = []

for sample in sorted(sample_dirs):
    sample_path = os.path.join(BASE_DIR, sample)

    htseq_folder = os.path.join(sample_path, "htseq_counts")

    if not os.path.isdir(htseq_folder):
        print(f"⚠ No htseq_counts folder for {sample}")
        continue

    # Find count file
    count_file = None
    for f in os.listdir(htseq_folder):
        if f.endswith("_counts.txt"):
            count_file = os.path.join(htseq_folder, f)
            break

    if count_file is None:
        print(f"⚠ No count file found for {sample}")
        continue

    print(f"📖 Reading: {count_file}")

    df = pd.read_csv(count_file, sep="\t", header=None, names=["Gene", sample])

    # Remove HTSeq summary rows
    df = df[~df["Gene"].str.startswith("__")]

    count_tables.append(df)

# Merge tables
if not count_tables:
    print("❌ No count tables found.")
    exit(1)

merged_df = count_tables[0]

for df in count_tables[1:]:
    merged_df = merged_df.merge(df, on="Gene", how="outer")

merged_df = merged_df.fillna(0)

# Save Excel
output_file = os.path.join(BASE_DIR, "All_samples_counts.xlsx")
merged_df.to_excel(output_file, index=False)

print(f"\n✅ Merged count matrix saved as: {output_file}")
