#!/usr/bin/env python3
"""
05_mine_candidate_sponges.py
For each unexplained isolate:
  1. Predict ORFs with Prodigal
  2. Filter: 50-250 aa, signal peptide (SignalP6), no TM helices (TMHMM/Phobius)
  3. Cluster at 50% identity with MMseqs2
Outputs: data/candidate_sponges.faa
Requires: prodigal, mmseqs2
          pip install biopython
"""
import subprocess, os
import pandas as pd
from pathlib import Path
from Bio import SeqIO

os.makedirs("prodigal_out",       exist_ok=True)
os.makedirs("candidate_proteins", exist_ok=True)
os.makedirs("tmp_mmseqs",         exist_ok=True)

df = pd.read_csv("data/unexplained_cohort.tsv", sep="\t")
df = df[df["fasta_exists"] == True]
print(f"Processing {len(df)} unexplained isolates\n")

# ── Step 1: Prodigal ORF prediction ─────────────────────────────────────────
print("Running Prodigal ...")
for _, row in df.iterrows():
    gid = row["genome_id"]
    fna = Path(f"genomes/{gid}.fna")
    faa = Path(f"prodigal_out/{gid}.faa")
    if faa.exists():
        print(f"  {gid}: already done")
        continue
    r = subprocess.run(
        ["prodigal", "-i", str(fna), "-a", str(faa), "-p", "single", "-q"],
        capture_output=True, text=True
    )
    n = len([1 for _ in SeqIO.parse(faa, "fasta")]) if faa.exists() else 0
    print(f"  {gid}: {n} ORFs predicted")

# ── Step 2: Filter small proteins (50-250 aa) ───────────────────────────────
print("\nFiltering small proteins (50-250 aa) ...")
all_small = []
for faa in Path("prodigal_out").glob("*.faa"):
    gid = faa.stem
    for rec in SeqIO.parse(faa, "fasta"):
        L = len(rec.seq.rstrip("*"))   # remove stop codon
        if 50 <= L <= 250:
            rec.id          = f"{gid}|{rec.id}"
            rec.description = ""
            rec.seq         = rec.seq.rstrip("*")
            all_small.append(rec)

SeqIO.write(all_small, "candidate_proteins/all_small.faa", "fasta")
print(f"Small proteins (50-250 aa): {len(all_small)}")

# ── Step 3: SignalP 6 instructions ──────────────────────────────────────────
print("""
── Next: run SignalP 6 (signal peptide prediction) ──────────────────────────
SignalP requires a separate academic licence from DTU.
If you have it installed:

  signalp6 --fastafile candidate_proteins/all_small.faa \\
           --output_dir signalp_out \\
           --format none \\
           --organism other \\
           --mode fast

Alternatively use DeepSig (open source):
  pip install deepsig
  deepsig -f candidate_proteins/all_small.faa -o signalp_out/deepsig.tsv -k gramn

─────────────────────────────────────────────────────────────────────────────
""")

# ── Step 4: TMHMM / Phobius TM filter instructions ──────────────────────────
print("""
── Next: run TMHMM or Phobius (TM helix prediction) ─────────────────────────
  conda install -c bioconda tmhmm2   # if available
  tmhmm --short < candidate_proteins/all_small.faa > tmhmm_out/all_small.tmhmm

  OR use Phobius web server for small datasets:
  https://phobius.sbc.su.se/ (upload candidate_proteins/all_small.faa)

─────────────────────────────────────────────────────────────────────────────
""")

# ── Step 5: Apply SP + TM filters if output files exist ─────────────────────
from pathlib import Path

sp_file  = Path("signalp_out/prediction_results.txt")    # SignalP6
sp_file2 = Path("signalp_out/deepsig.tsv")               # DeepSig
tm_file  = Path("tmhmm_out/all_small.tmhmm")

with_sp, no_tm = set(), set()

if sp_file.exists():
    import csv
    with open(sp_file) as fh:
        for line in csv.DictReader(fh, delimiter="\t"):
            if line.get("Prediction","").startswith("SP"):
                with_sp.add(line["# ID"])
    print(f"Signal peptide hits (SignalP6): {len(with_sp)}")

elif sp_file2.exists():
    import csv
    with open(sp_file2) as fh:
        for line in csv.DictReader(fh, delimiter="\t"):
            if line.get("prediction","") == "S":   # S = signal peptide in DeepSig
                with_sp.add(line["identifier"])
    print(f"Signal peptide hits (DeepSig): {len(with_sp)}")
else:
    print("No SP filter applied yet — run SignalP6 or DeepSig first")
    with_sp = {r.id for r in all_small}   # pass all through for now

if tm_file.exists():
    with open(tm_file) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 5:
                n_tm = int(parts[4].split("=")[1])
                if n_tm == 0:
                    no_tm.add(parts[0])
    print(f"No-TM proteins: {len(no_tm)}")
else:
    print("No TM filter applied yet — run TMHMM first")
    no_tm = {r.id for r in all_small}   # pass all through for now

filtered = [r for r in all_small if r.id in with_sp and r.id in no_tm]
SeqIO.write(filtered, "candidate_proteins/filtered.faa", "fasta")
print(f"\nProteins passing SP + no-TM: {len(filtered)}")

# ── Step 6: MMseqs2 clustering at 50% identity ───────────────────────────────
print("\nClustering with MMseqs2 (50% identity) ...")
subprocess.run([
    "mmseqs", "easy-cluster",
    "candidate_proteins/filtered.faa",
    "candidate_proteins/nr50",
    "tmp_mmseqs",
    "--min-seq-id", "0.50",
    "-c",           "0.80",
    "--cov-mode",   "1",
    "--threads",    "10",
    "-v",           "1"
])

rep_faa = Path("candidate_proteins/nr50_rep_seq.fasta")
if rep_faa.exists():
    recs = list(SeqIO.parse(rep_faa, "fasta"))
    SeqIO.write(recs, "data/candidate_sponges.faa", "fasta")
    print(f"\nNon-redundant candidate sponge proteins: {len(recs)}")
    print("Saved: data/candidate_sponges.faa")
else:
    # Fallback: just copy filtered without clustering
    Path("data/candidate_sponges.faa").write_text(
        Path("candidate_proteins/filtered.faa").read_text()
    )
    print("MMseqs2 rep file not found — saved filtered.faa as candidate_sponges.faa")
    print(f"Total: {len(filtered)} proteins")
