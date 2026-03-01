#!/usr/bin/env python3

import os
import re
import shutil

# Project directory (one level above pipeline)
project_dir = os.path.abspath(os.path.join(os.getcwd(), ".."))

# Match SRR + read number
pattern = re.compile(r'(SRR\d+)_(\d)\.fastq\.gz')

print(f"\n🔍 Scanning project directory: {project_dir}")

# Walk through project directory and subfolders
for root, dirs, files in os.walk(project_dir):

    # Skip pipeline folder
    if "pipeline" in root:
        continue

    for file in files:

        if not file.endswith(".gz"):
            continue

        match = pattern.search(file)

        if match:
            srr = match.group(1)
            read = match.group(2)

            correct_folder = os.path.join(project_dir, srr)
            os.makedirs(correct_folder, exist_ok=True)

            old_path = os.path.join(root, file)
            new_name = f"{srr}_{read}.fastq.gz"
            new_path = os.path.join(correct_folder, new_name)

            # If file already correct and in correct place
            if old_path == new_path:
                continue

            # Remove duplicates like (1)
            if os.path.exists(new_path):
                print(f"⚠ Duplicate found — removing {old_path}")
                os.remove(old_path)
                continue

            shutil.move(old_path, new_path)
            print(f"✔ Moved & Renamed: {file} → {srr}/{new_name}")

print("\n✅ Organization + Renaming completed.")
