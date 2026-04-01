#!/usr/bin/env python3
"""
02_download_genomes.py  (v2 - handles GCA_ + CP* accessions)
"""
import subprocess, os, shutil, time
import pandas as pd
from pathlib import Path

os.makedirs("genomes", exist_ok=True)
os.makedirs("data/ncbi_zips", exist_ok=True)

df = pd.read_csv("data/bvbrc_genome_metadata.tsv", sep="\t")

def pick_accession(row):
    for col in ["assembly_accession", "genbank_accessions"]:
        val = str(row.get(col, "")).strip()
        if val and val != "nan":
            return val.split(",")[0].split(";")[0].strip()
    return None

df["accession"] = df.apply(pick_accession, axis=1)
df_valid = df.dropna(subset=["accession"]).drop_duplicates(subset=["accession"])

# Split into GCA_ (datasets CLI) and CP*/other (efetch)
gca_df = df_valid[df_valid["accession"].str.startswith(("GCA_","GCF_"))]
seq_df = df_valid[~df_valid["accession"].str.startswith(("GCA_","GCF_"))]
print(f"GCA/GCF accessions : {len(gca_df)}")
print(f"Sequence accessions: {len(seq_df)}")
print(f"No accession       : {len(df) - len(df_valid)}")

already_done = {p.stem for p in Path("genomes").glob("*.fna")}
print(f"Already downloaded : {len(already_done)}\n")

# ── ROUTE 1: GCA_/GCF_ via NCBI datasets CLI ────────────────────────────────
acc2gid = df_valid.set_index("accession")["genome_id"].to_dict()

gca_todo = [a for a in gca_df["accession"].tolist()
            if acc2gid[a] not in already_done]
print(f"GCA to download: {len(gca_todo)}")

BATCH = 200
failed_gca = []
for i in range(0, len(gca_todo), BATCH):
    batch     = gca_todo[i:i+BATCH]
    bn        = i // BATCH
    acc_file  = f"data/ncbi_zips/gca_batch_{bn:04d}.txt"
    zip_file  = f"data/ncbi_zips/gca_batch_{bn:04d}.zip"
    unzip_dir = f"data/ncbi_zips/gca_unzip_{bn:04d}"

    if not Path(unzip_dir).exists():
        with open(acc_file, "w") as fh:
            fh.write("\n".join(batch))
        cmd = ["datasets","download","genome","accession",
               "--inputfile", acc_file,
               "--include","genome",
               "--no-progressbar",
               "--filename", zip_file]
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            print(f"  Batch {bn} WARN: {r.stderr[:150]}")
            failed_gca.extend(batch)
            continue
        subprocess.run(["unzip","-q","-o", zip_file,"-d", unzip_dir],
                       capture_output=True)

    for fna in Path(unzip_dir).rglob("*.fna"):
        acc_key = fna.parent.name
        gid = acc2gid.get(acc_key)
        if gid:
            dest = Path("genomes") / f"{gid}.fna"
            if not dest.exists():
                shutil.copy(fna, dest)

    print(f"  GCA batch {bn:04d}: done  "
          f"(total so far: {len(list(Path('genomes').glob('*.fna')))})")
    time.sleep(1)

# ── ROUTE 2: CP*/sequence accessions via efetch ──────────────────────────────
seq_todo = [(row["accession"], row["genome_id"])
            for _, row in seq_df.iterrows()
            if row["genome_id"] not in already_done]
print(f"\nSequence accessions to fetch: {len(seq_todo)}")

failed_seq = []
for acc, gid in seq_todo:
    dest = Path("genomes") / f"{gid}.fna"
    if dest.exists():
        continue
    # Use efetch to grab FASTA by nucleotide accession
    cmd = ["efetch", "-db", "nuccore", "-id", acc, "-format", "fasta"]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        if r.returncode == 0 and r.stdout.startswith(">"):
            with open(dest, "w") as fh:
                fh.write(r.stdout)
        else:
            failed_seq.append((acc, gid))
    except subprocess.TimeoutExpired:
        failed_seq.append((acc, gid))
    time.sleep(0.4)   # NCBI rate limit: max 3 req/s without API key

    if len(list(Path("genomes").glob("*.fna"))) % 100 == 0:
        print(f"  ... {len(list(Path('genomes').glob('*.fna')))} genomes downloaded")

# ── Summary ───────────────────────────────────────────────────────────────────
total = len(list(Path("genomes").glob("*.fna")))
print(f"\n{'='*50}")
print(f"Total genomes downloaded : {total}")
print(f"Failed GCA               : {len(failed_gca)}")
print(f"Failed sequence          : {len(failed_seq)}")

if failed_gca or failed_seq:
    with open("data/failed_accessions.txt","w") as fh:
        for a in failed_gca:
            fh.write(f"GCA\t{a}\n")
        for a,g in failed_seq:
            fh.write(f"SEQ\t{a}\t{g}\n")
    print("Failures saved to data/failed_accessions.txt")
