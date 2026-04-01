#!/usr/bin/env python3
"""
01_bvbrc_query.py  (v5 - resumable + incremental save)
"""
import os, time, json
from urllib.parse import quote
import requests
import pandas as pd

os.makedirs("data", exist_ok=True)

BASE = "https://www.bv-brc.org/api"
CHECKPOINT = "data/amr_checkpoint.json"   # tracks last successful offset

def bvbrc_post_resumable(endpoint, rql, limit=2000, save_every=10000):
    """
    Resumable BV-BRC paginator.
    - Saves rows to disk every `save_every` records
    - Resumes from last checkpoint if interrupted
    - Retries with exponential backoff on connection errors
    """
    # Load checkpoint if exists
    if os.path.exists(CHECKPOINT):
        with open(CHECKPOINT) as f:
            ck = json.load(f)
        start = ck.get("offset", 0)
        rows  = ck.get("rows", [])
        print(f"  Resuming from offset {start} ({len(rows)} records already saved)")
    else:
        start, rows = 0, []

    url = f"{BASE}/{endpoint}/"
    consecutive_errors = 0

    while True:
        body = f"{rql}&limit({limit},{start})"
        try:
            r = requests.post(
                url, data=body,
                headers={
                    "Accept": "application/json",
                    "Content-Type": "application/rqlquery+x-www-form-urlencoded"
                },
                timeout=60
            )
            r.raise_for_status()
            chunk = r.json()
            consecutive_errors = 0

        except (requests.exceptions.ConnectionError,
                requests.exceptions.Timeout,
                requests.exceptions.ChunkedEncodingError) as e:
            consecutive_errors += 1
            wait = min(60, 5 * 2 ** consecutive_errors)
            print(f"  Connection error (attempt {consecutive_errors}): {e!s:.80}")
            print(f"  Retrying in {wait}s ...")
            time.sleep(wait)
            continue

        if not chunk:
            break

        rows.extend(chunk)
        start += len(chunk)
        print(f"  ... fetched {len(rows)} records total")

        # Save checkpoint every save_every records
        if len(rows) % save_every < limit:
            with open(CHECKPOINT, "w") as f:
                json.dump({"offset": start, "rows": rows}, f)
            print(f"  [checkpoint saved at offset {start}]")

        if len(chunk) < limit:
            break

        time.sleep(0.8)   # be polite to the server

    # Final save & clear checkpoint
    if os.path.exists(CHECKPOINT):
        os.remove(CHECKPOINT)
    return rows


# ── Discover antibiotic names ───────────────────────────────────────────────
print("Discovering antibiotic names for E. coli ...")
sample = bvbrc_post_resumable(
    "genome_amr",
    "eq(taxon_id,562)&select(antibiotic,resistant_phenotype)"
)
df_sample = pd.DataFrame(sample)
df_sample.to_csv("data/antibiotic_names_discovery.tsv", sep="\t", index=False)
print("\nTop antibiotic names:")
print(df_sample["antibiotic"].value_counts().head(20).to_string())

# ── Find TZP name ───────────────────────────────────────────────────────────
tzp_mask = df_sample["antibiotic"].str.lower().str.contains("piperacillin", na=False)
tzp_name = df_sample[tzp_mask]["antibiotic"].value_counts().idxmax()
print(f"\nUsing: '{tzp_name}'")

# ── Fetch resistant genomes ─────────────────────────────────────────────────
# Delete old checkpoint so this query starts fresh
if os.path.exists(CHECKPOINT):
    os.remove(CHECKPOINT)

print(f"\nFetching Resistant records for '{tzp_name}' ...")
amr_rows = bvbrc_post_resumable(
    "genome_amr",
    f"eq(taxon_id,562)&eq(antibiotic,{quote(tzp_name,safe='')})"
    f"&eq(resistant_phenotype,Resistant)"
    f"&select(genome_id,genome_name,antibiotic,resistant_phenotype,"
    f"laboratory_typing_method,measurement,measurement_sign,measurement_unit,evidence)"
)
df_amr = pd.DataFrame(amr_rows)
df_amr.to_csv("data/bvbrc_amr_metadata.tsv", sep="\t", index=False)
resistant_ids = df_amr["genome_id"].dropna().unique().tolist()
print(f"\nResistant records : {len(df_amr)}")
print(f"Unique genome IDs : {len(resistant_ids)}")

# ── Genome metadata ─────────────────────────────────────────────────────────
if os.path.exists(CHECKPOINT):
    os.remove(CHECKPOINT)

print("\nFetching genome metadata ...")
id_list = ",".join(str(i) for i in resistant_ids[:5000])
genome_rows = bvbrc_post_resumable(
    "genome",
    f"in(genome_id,({id_list}))"
    f"&select(genome_id,genome_name,strain,isolation_source,"
    f"assembly_accession,taxon_id,genome_length,contigs,"
    f"patric_cds,genbank_accessions)"
)
df_genome = pd.DataFrame(genome_rows)
df_genome.to_csv("data/bvbrc_genome_metadata.tsv", sep="\t", index=False)
print(f"Saved {len(df_genome)} genomes → data/bvbrc_genome_metadata.tsv")