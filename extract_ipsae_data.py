import os
import csv
import re
import numpy as np

MAIN = os.getcwd()

OUTPUT_CSV = "ipsae_summary.csv"


def extract_mean_ipsae(txt_file):
    """
    Reads a *_10_10.txt ipsae output file and computes the mean IPSAE
    of rows containing 'max'. IPSAE is expected to be in column 6.
    Returns None if no valid numbers found.
    """

    ipsae_values = []

    with open(txt_file, "r") as f:
        for line in f:
            if "max" not in line:
                continue

            # Split on whitespace
            parts = line.strip().split()

            if len(parts) < 6:
                continue

            # Column 6 is index 5
            try:
                value = float(parts[5])
                ipsae_values.append(value)
            except ValueError:
                continue

    if len(ipsae_values) == 0:
        return None

    return float(np.mean(ipsae_values))


def find_ipsae_output_files(folder):
    """
    Return dict: rank -> .txt output file
    Matches files like:
        blah_rank_003_something_10_10.txt
    """

    rank_re = re.compile(r"rank_(\d{3}).*?_15_15\.txt$")
    results = {}

    for file in os.listdir(folder):
        match = rank_re.search(file)
        if match:
            rank = f"rank_{match.group(1)}"
            results[rank] = os.path.join(folder, file)

    return results


# ----------------------------------------------------------------------
# MAIN PROCESSING
# ----------------------------------------------------------------------

rows = []

for protein in sorted(os.listdir(MAIN)):
    protein_path = os.path.join(MAIN, protein)
    if not os.path.isdir(protein_path):
        continue

    for repeat in sorted(os.listdir(protein_path)):
#        if not repeat.startswith("repeat_"):
 #           continue

        pdb_folder = os.path.join(protein_path, repeat, "pdb")
        if not os.path.isdir(pdb_folder):
            continue

        # Find all ipa output files for this protein/repeat
        rank_files = find_ipsae_output_files(pdb_folder)
        if not rank_files:
            continue

        for rank, txt_path in sorted(rank_files.items()):
            mean_value = extract_mean_ipsae(txt_path)

            if mean_value is None:
                print(f"[WARNING] No IPSAE values found in {txt_path}")
                continue

            rows.append([protein, repeat, rank, mean_value])


# Write output CSV
with open(OUTPUT_CSV, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["protein", "repeat", "rank", "mean_ipsae"])
    writer.writerows(rows)

print(f"\n✨ IPSAE summary written to: {OUTPUT_CSV}\n")
print(f"Total entries: {len(rows)}")

"extract_ipsae_data.py" 105L, 2670C     