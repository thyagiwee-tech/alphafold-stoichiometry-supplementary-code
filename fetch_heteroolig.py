import csv
import os
import shutil
import requests

CHAIN_LABELS = ["A", "B", "C", "D", "E"]


def fetch_fasta(uniprot_id):
    """
    Fetch the FASTA sequence for a UniProt ID.
    """
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    response = requests.get(url, timeout=30)

    if response.status_code == 200 and response.text.strip():
        fasta_data = response.text.strip()
        fasta_sequence = "".join(
            line.strip() for line in fasta_data.splitlines() if not line.startswith(">")
        )
        return fasta_sequence

    print(f"Error fetching data for UniProt ID: {uniprot_id} (HTTP {response.status_code})")
    return None


def sanitize_folder_name(name):
    """
    Make folder names filesystem-safe.
    """
    bad_chars = '<>:"/\\|?*'
    cleaned = "".join("_" if c in bad_chars else c for c in name.strip())
    return cleaned.replace(" ", "_")


def get_populated_chain_ids(row):
    """
    Return populated chain IDs from columns A-E.
    """
    chain_ids = {}
    for label in CHAIN_LABELS:
        value = row.get(label, "").strip()
        if value:
            chain_ids[label] = value
    return chain_ids


def generate_combinations(labels, total_chains):
    """
    Recursively generate all non-negative combinations of counts
    for the given labels that sum to total_chains.

    Example:
      labels = ["A", "B"], total_chains = 3
      -> {"A":3,"B":0}, {"A":2,"B":1}, {"A":1,"B":2}, {"A":0,"B":3}
    """
    combos = []

    def recurse(idx, remaining, current):
        if idx == len(labels) - 1:
            current[labels[idx]] = remaining
            combos.append(current.copy())
            return

        label = labels[idx]
        for count in range(remaining, -1, -1):
            current[label] = count
            recurse(idx + 1, remaining - count, current)

    recurse(0, total_chains, {})
    return combos


def combo_folder_name(combo, labels):
    """
    Turn a combo dict into a folder name like 3A0B or 1A1B1C.
    """
    return "".join(f"{combo[label]}{label}" for label in labels)


def build_colon_separated_sequence(combo, labels, sequences):
    """
    Build one colon-separated sequence string in label order.
    Example for 1A2B:
      seqA:seqB:seqB
    """
    parts = []
    for label in labels:
        count = combo[label]
        seq = sequences[label]
        parts.extend([seq] * count)
    return ":".join(parts)


def copy_support_scripts(destination_folder):
    """
    Copy af1.sh and af2.sh into each destination folder.
    """
    for script in ["af1.sh", "af2.sh"]:
        if os.path.exists(script):
            shutil.copy(script, destination_folder)
            print(f"Copied {script} to {destination_folder}")
        else:
            print(f"Warning: {script} not found in current directory")

"fetch_heteroolig.py" 182L, 5257C                                                                                                                                                                                                                                                                                                                                   1,1           Top

