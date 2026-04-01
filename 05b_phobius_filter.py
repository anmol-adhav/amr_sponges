#!/usr/bin/env python3
"""
05b_phobius_filter.py
Predict signal peptides AND TM helices via Phobius web API.
Outputs: tmhmm_out/phobius_results.tsv
         candidate_proteins/filtered_sp_notm.faa
         data/candidate_sponges_final.faa
"""
import requests, time, subprocess
import pandas as pd
from pathlib import Path
from Bio import SeqIO

Path("tmhmm_out").mkdir(exist_ok=True)

recs  = list(SeqIO.parse("candidate_proteins/all_small.faa", "fasta"))
BATCH = 200
results = []

print(f"Submitting {len(recs)} sequences to Phobius ({BATCH}/batch) ...")

for i in range(0, len(recs), BATCH):
    batch     = recs[i:i+BATCH]
    fasta_str = "".join(f">{r.id}\n{r.seq}\n" for r in batch)

    for attempt in range(3):
        try:
            r = requests.post(
                "https://phobius.sbc.su.se/cgi-bin/predict.pl",
                data={"format": "short", "protseq": fasta_str},
                timeout=120
            )
            break
        except requests.exceptions.Timeout:
            print(f"  Timeout batch {i//BATCH}, retry {attempt+1}")
            time.sleep(10)

    # Parse short format output:
    # SEQNAME   TM  SP  PREDICTION
    for line in r.text.splitlines():
        line = line.strip()
        if not line or line.startswith("SEQNAME") or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) >= 4:
            results.append({
                "id":     parts[0],
                "n_tm":   int(parts[1]),
                "has_sp": parts[2] == "Y",
                "pred":   parts[3]
            })

    print(f"  [{min(i+BATCH, len(recs))}/{len(recs)}] processed")
    time.sleep(1.5)   # be polite to the server

# ── Save raw results ──────────────────────────────────────────────────────────
df = pd.DataFrame(results)
df.to_csv("tmhmm_out/phobius_results.tsv", sep="\t", index=False)
print(f"\nPhobius results: {len(df)} sequences")
print(f"  Has SP only (no TM) : {(df['has_sp'] & (df['n_tm']==0)).sum()}")
print(f"  TM only (no SP)     : {(~df['has_sp'] & (df['n_tm']>0)).sum()}")
print(f"  Cytoplasmic (no SP, no TM): {(~df['has_sp'] & (df['n_tm']==0)).sum()}")

# ── Filter: keep SP + no TM (periplasmic) ────────────────────────────────────
periplasmic_ids = set(df[df["has_sp"] & (df["n_tm"] == 0)]["id"])
print(f"\nPeriplasmic candidates (SP + no TM): {len(periplasmic_ids)}")

filtered = [r for r in recs if r.id in periplasmic_ids]
SeqIO.write(filtered, "candidate_proteins/filtered_sp_notm.faa", "fasta")

# ── Re-cluster at 50% identity ────────────────────────────────────────────────
print("\nRe-clustering filtered proteins at 50% identity ...")
subprocess.run([
    "mmseqs", "easy-cluster",
    "candidate_proteins/filtered_sp_notm.faa",
    "candidate_proteins/nr50_final",
    "tmp_mmseqs",
    "--min-seq-id", "0.50",
    "-c",           "0.80",
    "--cov-mode",   "1",
    "--threads",    "10",
    "-v",           "1"
])

rep = Path("candidate_proteins/nr50_final_rep_seq.fasta")
if rep.exists():
    final = list(SeqIO.parse(rep, "fasta"))
    SeqIO.write(final, "data/candidate_sponges_final.faa", "fasta")
    print(f"\nFinal non-redundant periplasmic candidates: {len(final)}")
    print("Saved: data/candidate_sponges_final.faa → ready for ColabFold")
else:
    SeqIO.write(filtered, "data/candidate_sponges_final.faa", "fasta")
    print(f"Saved unclusterd: {len(filtered)} proteins → data/candidate_sponges_final.faa")
