#!/usr/bin/env python3
"""
04_filter_unexplained.py
Cross-reference phenotypic TZP resistance with ResFinder results.
Isolates with resistance but NO known beta-lactam ARGs = unexplained cohort.
Outputs: data/unexplained_cohort.tsv
"""
import pandas as pd
from pathlib import Path

df_amr    = pd.read_csv("data/bvbrc_amr_metadata.tsv",    sep="\t")
df_genome = pd.read_csv("data/bvbrc_genome_metadata.tsv", sep="\t")
df_rf     = pd.read_csv("data/resfinder_summary.tsv",     sep="\t")

# ── Step 1: phenotypically TZP-resistant genomes ─────────────────────────────
resistant = (
    df_amr[df_amr["resistant_phenotype"] == "Resistant"][["genome_id"]]
    .drop_duplicates()
)
print(f"Phenotypically TZP-resistant : {len(resistant)}")

# ── Step 2: merge with ResFinder summary ─────────────────────────────────────
merged = resistant.merge(df_rf, on="genome_id", how="inner")
print(f"After ResFinder merge        : {len(merged)}")

# ── Step 3: keep only those with NO known beta-lactam ARGs ───────────────────
unexplained = merged[merged["betalactam_genes"] == "NONE"].copy()
print(f"No known beta-lactam ARGs    : {len(unexplained)}  ← YOUR TARGET COHORT")

# ── Step 4: attach genome metadata ───────────────────────────────────────────
unexplained = unexplained.merge(df_genome, on="genome_id", how="left")

# ── Step 5: verify genome FASTA exists for each ──────────────────────────────
unexplained["fasta_exists"] = unexplained["genome_id"].apply(
    lambda g: Path(f"genomes/{g}.fna").exists()
)
missing = unexplained[~unexplained["fasta_exists"]]
if not missing.empty:
    print(f"\nWARN: {len(missing)} genomes missing FASTA — will be skipped in downstream steps")
    print(missing[["genome_id","genome_name"]].to_string())

unexplained.to_csv("data/unexplained_cohort.tsv", sep="\t", index=False)
print(f"\nSaved: data/unexplained_cohort.tsv")
print()
print(unexplained[["genome_id","genome_name","isolation_source",
                    "genome_length","fasta_exists"]].to_string())
