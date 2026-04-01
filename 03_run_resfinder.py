#!/usr/bin/env python3
"""
03_run_resfinder.py  (macOS fixed)
"""
import subprocess, os, csv, sys
import pandas as pd
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed  # Thread not Process

RESFINDER_DB = os.environ.get("RESFINDER_DB", "resfinder_db")
THREADS      = int(os.environ.get("THREADS", 10))

def run_one(fna_path: Path):
    gid     = fna_path.stem
    out_dir = Path("resfinder_results") / gid
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable, "-m", "resfinder",
        "-ifa",  str(fna_path.resolve()),   # ← resolve to absolute path
        "-o",    str(out_dir.resolve()),     # ← resolve to absolute path
        "-db_res", RESFINDER_DB,
        "-acq",
        "-l",  "0.60",
        "-t",  "0.80",
    ]
    # DEBUG: print first command only
    if gid == sorted(Path("genomes").glob("*.fna"))[0].stem:
        print(f"DEBUG CMD: {' '.join(cmd)}", flush=True)
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=180)
    return gid, (r.returncode == 0), r.stderr[:200]

if __name__ == "__main__":
    if not Path(RESFINDER_DB).exists():
        print(f"ERROR: RESFINDER_DB not found at {RESFINDER_DB}")
        sys.exit(1)

    genomes = sorted(Path("genomes").glob("*.fna"))
    already = {p.parent.name for p in
               Path("resfinder_results").rglob("ResFinder_results_tab.txt")}
    todo    = [g for g in genomes if g.stem not in already]

    print(f"Total genomes     : {len(genomes)}")
    print(f"Already done      : {len(already)}")
    print(f"To run            : {len(todo)}")
    print(f"Threads           : {THREADS}\n")

    ok, err = 0, 0
    with ThreadPoolExecutor(max_workers=THREADS) as pool:
        futures = {pool.submit(run_one, g): g for g in todo}
        for i, fut in enumerate(as_completed(futures), 1):
            gid, success, msg = fut.result()
            if success:
                ok += 1
            else:
                err += 1
                print(f"  WARN {gid}: {msg}")
            if i % 100 == 0:
                print(f"  [{i}/{len(todo)}] ok={ok} errors={err}")

    print(f"\nResFinder done. OK={ok}  Errors={err}")

    # ── Summarise ────────────────────────────────────────────────────────────
    BETALACTAM_CLASSES = {
        "Beta-lactam","Carbapenem","Cephalosporin","Penicillin",
        "Monobactam","Cephamycin","Oxacillin","Aminopenicillin"
    }
    rows = []
    for result_file in sorted(Path("resfinder_results").rglob("ResFinder_results_tab.txt")):
        gid   = result_file.parent.name
        genes = []
        with open(result_file) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for line in reader:
                if any(c in line.get("Class","") for c in BETALACTAM_CLASSES):
                    genes.append(line.get("Resistance gene",""))
        rows.append({
            "genome_id":          gid,
            "n_betalactam_genes": len(genes),
            "betalactam_genes":   ";".join(genes) if genes else "NONE"
        })

    df_rf = pd.DataFrame(rows)
    df_rf.to_csv("data/resfinder_summary.tsv", sep="\t", index=False)

    n_none = (df_rf["betalactam_genes"] == "NONE").sum()
    n_with = (df_rf["betalactam_genes"] != "NONE").sum()
    print(f"\nNo known beta-lactam ARGs : {n_none}  ← unexplained candidates")
    print(f"Has known beta-lactam ARGs: {n_with}")
    print(f"\nTop genes found:")
    print(df_rf[df_rf["betalactam_genes"] != "NONE"]["betalactam_genes"]
          .str.split(";").explode().value_counts().head(15).to_string())