import requests
import random
import pandas as pd
import time

# ============================================================
# SETTINGS
# ============================================================

N_PER_CLASS = 1

OUTPUT_CSV = "heteromer_samples_summary.csv"

# ============================================================
# FUNCTIONS
# ============================================================

def search_rcsb(min_proteins=2, rows=3000):

    query = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_entry_info.polymer_entity_count_protein",
                "operator": "greater_or_equal",
                "value": min_proteins
            }
        },
        "return_type": "entry",
        "request_options": {
            "paginate": {
                "start": 0,
                "rows": rows
            }
        }
    }

    url = "https://search.rcsb.org/rcsbsearch/v2/query"

    r = requests.post(url, json=query)
    r.raise_for_status()

    results = r.json()["result_set"]

    return [x["identifier"] for x in results]


# ------------------------------------------------------------

def get_entry_details(pdb_id):

    graphql_url = "https://data.rcsb.org/graphql"

    graphql_query = """
    query($id: String!) {
      entry(entry_id: $id) {

        struct {
          title
        }

        polymer_entities {

          rcsb_polymer_entity_container_identifiers {
            asym_ids
            uniprot_ids
          }

          entity_poly {
            type
          }
        }
      }
    }
    """

    r = requests.post(
        graphql_url,
        json={
            "query": graphql_query,
            "variables": {"id": pdb_id}
        }
    )

    r.raise_for_status()

    return r.json()["data"]["entry"]


# ------------------------------------------------------------

def classify_complex(entry):

    if entry is None:
        return None

    entities = entry["polymer_entities"]

    total_subunits = 0

    composition = []

    seen = set()

    for ent in entities:

        poly_type = ent["entity_poly"]["type"]

        # proteins only
        if poly_type != "polypeptide(L)":
            continue

        ids = ent["rcsb_polymer_entity_container_identifiers"]

        chains = ids.get("asym_ids", [])
        uniprots = ids.get("uniprot_ids", [])

        # skip entities without UniProt mapping
        if not uniprots:
            continue

        uniprot = uniprots[0]

        stoich = len(chains)

        total_subunits += stoich

        # avoid duplicate UniProt entries
        if uniprot in seen:
            continue

        seen.add(uniprot)

        composition.append({
            "uniprot": uniprot,
            "stoich": stoich,
            "chains": ",".join(chains)
        })

    n_unique = len(composition)

    # must be heteromeric
    if n_unique < 2:
        return None

    # classify based on TOTAL number of chains
    if total_subunits == 2:
        oligomer = "heterodimer"

    elif total_subunits == 3:
        oligomer = "heterotrimer"

    elif total_subunits == 4:
        oligomer = "heterotetramer"

    else:
        return None

    return {
        "oligomeric_state": oligomer,
        "total_subunits": total_subunits,
        "unique_proteins": n_unique,
        "composition": composition
    }


# ============================================================
# MAIN
# ============================================================

print("Searching RCSB...")

all_pdbs = search_rcsb(
    min_proteins=2,
    rows=3000
)

print(f"Retrieved {len(all_pdbs)} candidate entries")

random.shuffle(all_pdbs)

results = []

counts = {
    "heterodimer": 0,
    "heterotrimer": 0,
    "heterotetramer": 0
}

needed = {
    "heterodimer": N_PER_CLASS,
    "heterotrimer": N_PER_CLASS,
    "heterotetramer": N_PER_CLASS
}

# ------------------------------------------------------------

for i, pdb_id in enumerate(all_pdbs):

    done = all(
        counts[k] >= needed[k]
        for k in counts
    )

    if done:
        break

    print(f"\n[{i+1}] Checking {pdb_id}")

    try:

        entry = get_entry_details(pdb_id)

        classified = classify_complex(entry)

        if classified is None:
            continue

        oligomer = classified["oligomeric_state"]

        # skip if already enough examples
        if counts[oligomer] >= needed[oligomer]:
            continue

        title = entry["struct"]["title"]

        print(f"  -> {oligomer}")

        row = {
            "pdb_id": pdb_id,
            "title": title,
            "oligomeric_state": oligomer,
            "total_subunits": classified["total_subunits"],
            "unique_proteins": classified["unique_proteins"]
        }

        # ----------------------------------------------------
        # Add UniProt + stoichiometry + chain columns
        # ----------------------------------------------------

        for j, comp in enumerate(classified["composition"], start=1):

            row[f"uniprot_{j}"] = comp["uniprot"]

            row[f"stoich_{j}"] = comp["stoich"]

            row[f"chains_{j}"] = comp["chains"]

        results.append(row)

        counts[oligomer] += 1

        print("  Counts:")
        print(counts)

        # be polite to API
        time.sleep(0.2)

    except Exception as e:

        print(f"  ERROR: {e}")

# ============================================================
# SAVE
# ============================================================

df = pd.DataFrame(results)

# ------------------------------------------------------------
# Sort columns nicely
# ------------------------------------------------------------

base_cols = [
    "pdb_id",
    "title",
    "oligomeric_state",
    "total_subunits",
    "unique_proteins"
]

uniprot_cols = sorted(
    [c for c in df.columns if c.startswith("uniprot_")]
)

stoich_cols = sorted(
    [c for c in df.columns if c.startswith("stoich_")]
)

chain_cols = sorted(
    [c for c in df.columns if c.startswith("chains_")]
)

ordered_cols = (
    base_cols
    + uniprot_cols
    + stoich_cols
    + chain_cols
)

df = df[ordered_cols]

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

df.to_csv(OUTPUT_CSV, index=False)

print("\n====================================")
print("FINAL COUNTS")
print("====================================")

for k, v in counts.items():
    print(f"{k}: {v}")

print(f"\nSaved to: {OUTPUT_CSV}")