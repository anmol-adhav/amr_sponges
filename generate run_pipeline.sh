#!/usr/bin/env bash
# =============================================================================
# run_pipeline.sh — AMR Sponge: Full Pipeline
# Identifies piperacillin/tazobactam sponge candidates in E. coli
#
# Usage:
#   bash run_pipeline.sh [--skip-download] [--skip-fold] [--skip-search]
#
# Requirements:
#   - conda envs: colabfold (Python 3.10 + biopython + pandas + requests)
#   - foldseek binary on PATH
#   - autodock-vina + openbabel on PATH
# =============================================================================

set -euo pipefail

# ── Colour output ─────────────────────────────────────────────────────────────
RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'
BLUE='\033[0;34m'; NC='\033[0m'
log()  { echo -e "${GREEN}[$(date '+%H:%M:%S')] $*${NC}"; }
warn() { echo -e "${YELLOW}[WARN] $*${NC}"; }
err()  { echo -e "${RED}[ERROR] $*${NC}"; exit 1; }
step() { echo -e "\n${BLUE}══════════════════════════════════════${NC}"; \
         echo -e "${BLUE} $*${NC}"; \
         echo -e "${BLUE}══════════════════════════════════════${NC}"; }

# ── Parse flags ───────────────────────────────────────────────────────────────
SKIP_DOWNLOAD=0; SKIP_FOLD=0; SKIP_SEARCH=0
for arg in "$@"; do
  case $arg in
    --skip-download) SKIP_DOWNLOAD=1 ;;
    --skip-fold)     SKIP_FOLD=1 ;;
    --skip-search)   SKIP_SEARCH=1 ;;
  esac
done

# ── Config ────────────────────────────────────────────────────────────────────
CONDA_ENV="colabfold"
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FASTA="$PROJECT_DIR/data/candidates_filtered_A.faa"
ESMFOLD_OUT="$PROJECT_DIR/esmfold_output"
HIGHCONF_DIR="$PROJECT_DIR/esmfold_highconf"
NOHIT_DIR="$PROJECT_DIR/dali_input"
FOLDSEEK_DB="$PROJECT_DIR/foldseek_db/pdb"
AFDB_DB="$PROJECT_DIR/foldseek_db/afdb_swissprot"
RESULTS_DIR="$PROJECT_DIR/results"
DOCKING_DIR="$PROJECT_DIR/docking"
TMP_DIR="$PROJECT_DIR/tmp"
PLDDT_THRESHOLD=0.70
EVALUE_THRESHOLD=0.001
THREADS=8

# ── Preflight checks ──────────────────────────────────────────────────────────
step "Preflight Checks"

[[ -f "$FASTA" ]] || err "Input FASTA not found: $FASTA"
command -v foldseek &>/dev/null || err "foldseek not on PATH. Install from mmseqs.com/foldseek/"
command -v python   &>/dev/null || err "python not found. Activate colabfold conda env first."

# Verify we're in the right conda env
CURRENT_ENV="${CONDA_DEFAULT_ENV:-none}"
if [[ "$CURRENT_ENV" != "$CONDA_ENV" ]]; then
  warn "Expected conda env '$CONDA_ENV', got '$CURRENT_ENV'."
  warn "Run: conda activate $CONDA_ENV && bash run_pipeline.sh"
  exit 1
fi

mkdir -p "$ESMFOLD_OUT" "$HIGHCONF_DIR" "$NOHIT_DIR" \
         "$RESULTS_DIR" "$DOCKING_DIR" "$TMP_DIR" \
         "$PROJECT_DIR/foldseek_db"

log "All checks passed. Project: $PROJECT_DIR"

# ══════════════════════════════════════════════════════════════════════════════
# STEP 1 — Structure Prediction (ESMFold API)
# ══════════════════════════════════════════════════════════════════════════════
step "STEP 1 — Structure Prediction via ESMFold API"

if [[ $SKIP_FOLD -eq 1 ]]; then
  warn "Skipping structure prediction (--skip-fold)"
else
  N_DONE=$(ls "$ESMFOLD_OUT"/*.pdb 2>/dev/null | wc -l | tr -d ' ')
  N_TOTAL=$(grep -c "^>" "$FASTA")
  log "Found $N_DONE / $N_TOTAL structures already predicted."

  if [[ $N_DONE -lt $N_TOTAL ]]; then
    python << 'PYEOF'
import requests, time, os, warnings
from pathlib import Path
from Bio import SeqIO
import urllib3
urllib3.disable_warnings()

FASTA    = os.environ.get("FASTA", "data/candidates_filtered_A.faa")
OUT_DIR  = Path(os.environ.get("ESMFOLD_OUT", "esmfold_output"))
API_URL  = "https://api.esmatlas.com/foldSequence/v1/pdb/"
MAX_LEN  = 400

OUT_DIR.mkdir(exist_ok=True)
recs    = list(SeqIO.parse(FASTA, "fasta"))
already = {p.stem for p in OUT_DIR.glob("*.pdb")}
todo    = [r for r in recs
           if r.id.replace("|","_").replace("/","_").replace(":","_") not in already]

print(f"Total: {len(recs)}  Done: {len(already)}  Todo: {len(todo)}")
failed = []

for idx, rec in enumerate(todo):
    safe_id = rec.id.replace("|","_").replace("/","_").replace(":","_")
    out_pdb = OUT_DIR / f"{safe_id}.pdb"
    seq     = str(rec.seq).replace("*","")[:MAX_LEN]

    for attempt in range(3):
        try:
            r = requests.post(API_URL, data=seq,
                headers={"Content-Type": "application/x-www-form-urlencoded"},
                timeout=60, verify=False)
            if r.status_code == 200 and r.text.startswith("HEADER"):
                out_pdb.write_text(r.text)
                print(f"  [{idx+1}/{len(todo)}] OK {safe_id}")
                break
            else:
                print(f"  WARN {safe_id}: HTTP {r.status_code}")
                time.sleep(10)
        except Exception as e:
            print(f"  ERR {safe_id}: {e}")
            time.sleep(10)
    else:
        failed.append(rec.id)
    time.sleep(3)

print(f"\nDone. Structures: {len(list(OUT_DIR.glob('*.pdb')))}  Failed: {len(failed)}")
PYEOF
  else
    log "All structures already predicted. Skipping."
  fi
fi

# ══════════════════════════════════════════════════════════════════════════════
# STEP 2 — pLDDT Confidence Filtering
# ══════════════════════════════════════════════════════════════════════════════
step "STEP 2 — pLDDT Confidence Filtering (threshold: $PLDDT_THRESHOLD)"

python << PYEOF
from pathlib import Path
import shutil

threshold = float("$PLDDT_THRESHOLD")
esmfold_out  = Path("$ESMFOLD_OUT")
highconf_dir = Path("$HIGHCONF_DIR")
highconf_dir.mkdir(exist_ok=True)

results = []
for pdb in sorted(esmfold_out.glob("*.pdb")):
    b = [float(l[60:66]) for l in pdb.read_text().splitlines() if l.startswith("ATOM")]
    if b:
        results.append((pdb.stem, sum(b)/len(b)))

results.sort(key=lambda x: -x[1])

with open("$RESULTS_DIR/plddt_scores.tsv", "w") as f:
    f.write("protein\tplddt\n")
    for name, score in results:
        f.write(f"{name}\t{score:.4f}\n")

copied = 0
for name, score in results:
    if score >= threshold:
        shutil.copy(esmfold_out / f"{name}.pdb", highconf_dir / f"{name}.pdb")
        copied += 1

high = sum(1 for _, s in results if s >= threshold)
low  = sum(1 for _, s in results if s <  threshold)
print(f"Total: {len(results)}")
print(f"High confidence (pLDDT >= {threshold}): {high}  → {highconf_dir}")
print(f"Low  confidence (pLDDT <  {threshold}): {low}   (excluded)")
PYEOF

# ══════════════════════════════════════════════════════════════════════════════
# STEP 3 — Foldseek Database Download
# ══════════════════════════════════════════════════════════════════════════════
step "STEP 3 — Foldseek Database Setup"

if [[ $SKIP_DOWNLOAD -eq 1 ]]; then
  warn "Skipping database download (--skip-download)"
else
  if [[ ! -f "${FOLDSEEK_DB}.dbtype" ]]; then
    log "Downloading PDB100 database (~2.1 GB)..."
    foldseek databases PDB "$FOLDSEEK_DB" "$TMP_DIR" --threads "$THREADS"
  else
    log "PDB100 database already present."
  fi

  if [[ ! -f "${AFDB_DB}.dbtype" ]]; then
    log "Downloading AlphaFold/Swiss-Prot database (~2 GB)..."
    foldseek databases Alphafold/Swiss-Prot "$AFDB_DB" "$TMP_DIR" --threads "$THREADS"
  else
    log "AlphaFold/Swiss-Prot database already present."
  fi
fi

# ══════════════════════════════════════════════════════════════════════════════
# STEP 4 — Foldseek Structural Search (high-confidence vs PDB)
# ══════════════════════════════════════════════════════════════════════════════
step "STEP 4 — Foldseek Search: High-Confidence vs PDB100"

if [[ $SKIP_SEARCH -eq 1 ]]; then
  warn "Skipping Foldseek search (--skip-search)"
else
  foldseek easy-search \
    "$HIGHCONF_DIR" \
    "$FOLDSEEK_DB" \
    "$RESULTS_DIR/foldseek_hits.tsv" \
    "$TMP_DIR" \
    --exhaustive-search 1 \
    --format-output "query,target,pident,alnlen,evalue,bits,lddt,prob" \
    -e "$EVALUE_THRESHOLD" \
    --threads "$THREADS"

  N_HITS=$(wc -l < "$RESULTS_DIR/foldseek_hits.tsv" | tr -d ' ')
  log "Foldseek complete. Total hits: $N_HITS"
fi

# ══════════════════════════════════════════════════════════════════════════════
# STEP 5 — Annotation, Categorization & β-Lactam Filtering
# ══════════════════════════════════════════════════════════════════════════════
step "STEP 5 — Functional Annotation & β-Lactam Candidate Filtering"

python << 'PYEOF'
import requests, time, json, shutil
import pandas as pd
from pathlib import Path

results_dir  = Path("results")
highconf_dir = Path("esmfold_highconf")
nohit_dir    = Path("dali_input")
nohit_dir.mkdir(exist_ok=True)

# ── Best hit per query ────────────────────────────────────────────────────────
df = pd.read_csv(results_dir / "foldseek_hits.tsv", sep="\t", header=None,
    names=["query","target","pident","alnlen","evalue","bits","lddt","prob"])
best = df.sort_values("bits", ascending=False).groupby("query").first().reset_index()
best["pdb_id"] = best["target"].str.extract(r'^([a-z0-9]+)-')
best.to_csv(results_dir / "foldseek_best_hits.tsv", sep="\t", index=False)

# ── RCSB annotation ───────────────────────────────────────────────────────────
print("Fetching RCSB annotations...")
for pdb in best["pdb_id"].dropna().unique():
    try:
        r = requests.get(f"https://data.rcsb.org/rest/v1/core/entry/{pdb.upper()}", timeout=10)
        if r.status_code == 200:
            title = r.json().get("struct", {}).get("title", "N/A")
            best.loc[best["pdb_id"] == pdb, "annotation"] = title
        time.sleep(0.3)
    except:
        pass
best["annotation"] = best["annotation"].fillna("Unknown")

# ── Functional categorization ─────────────────────────────────────────────────
def categorize(ann):
    ann = str(ann).lower()
    if any(x in ann for x in ["pilus","fim","usher","pilin","gsp","secretion","pseudopilin","eps"]):
        return "Pilus / Secretion System"
    if any(x in ann for x in ["outer membrane","ompw","ompx","tolb","lipoprotein","lpt","mla","trat","pqi"]):
        return "Outer Membrane / Lipid Transport"
    if any(x in ann for x in ["transpeptidase","penicillin","pbp","peptidoglycan","d-ala","endopeptidase","ykud","dpa","hydrolase"]):
        return "Cell Wall / Peptidoglycan"
    if any(x in ann for x in ["superoxide","thioredoxin","copper","merp","cus","dsb","redox","oxidoreductase"]):
        return "Redox / Metal Resistance"
    if any(x in ann for x in ["cyclophilin","periplasmic","glucosidase","kinase","lyase","binding protein"]):
        return "Periplasmic Metabolism / Binding"
    if any(x in ann for x in ["flagell","motb","flil","flga"]):
        return "Flagella / Motility"
    if any(x in ann for x in ["sensor","histidine kinase","response regulator","toxr","cpxp"]):
        return "Signaling / Two-Component System"
    if any(x in ann for x in ["ivy","lysozyme","inhibitor","immunity"]):
        return "Phage / Lysozyme Defense"
    if any(x in ann for x in ["hypothetical","uncharacterized","duf","putative"]):
        return "Hypothetical / Unknown"
    return "Other"

best["category"] = best["annotation"].apply(categorize)
best.sort_values("bits", ascending=False).to_csv(
    results_dir / "foldseek_best_hits_annotated.tsv", sep="\t", index=False)

# ── β-Lactam filter ───────────────────────────────────────────────────────────
blactam_kw = ["transpeptidase","penicillin","pbp","beta-lactam","lactam",
              "dd-peptidase","endopeptidase","cell wall","peptidoglycan",
              "d-ala","carboxypeptidase","ykud","dpa"]
mask = best["annotation"].str.lower().apply(lambda x: any(k in x for k in blactam_kw))
blactam = best[mask].sort_values("bits", ascending=False)
blactam.to_csv(results_dir / "blactam_candidates.tsv", sep="\t", index=False)

# ── Extract no-hit proteins ───────────────────────────────────────────────────
all_hc   = {p.stem for p in highconf_dir.glob("*.pdb")}
hit_ids  = set(best["query"].str.strip())
no_hits  = all_hc - hit_ids
scores   = {}
for line in open(results_dir / "plddt_scores.tsv").readlines()[1:]:
    n, s = line.strip().split("\t"); scores[n] = float(s)

for name in no_hits:
    shutil.copy(highconf_dir / f"{name}.pdb", nohit_dir / f"{name}.pdb")

priority = sorted(no_hits, key=lambda x: -scores.get(x, 0))
with open(results_dir / "dali_priority_list.txt", "w") as f:
    for p in priority:
        f.write(f"{p}\t{scores.get(p,0):.3f}\n")

# ── Summary ───────────────────────────────────────────────────────────────────
print(f"\n{'='*50}")
print(f"ANNOTATION SUMMARY")
print(f"{'='*50}")
print(f"Total high-confidence structures : {len(all_hc)}")
print(f"Structures with Foldseek hit     : {len(hit_ids)}")
print(f"Structures with NO hit (novel?)  : {len(no_hits)}")
print(f"β-Lactam relevant candidates     : {len(blactam)}")
print(f"\nFunctional categories:")
for cat, n in best["category"].value_counts().items():
    print(f"  {n:4d}  {cat}")
print(f"\nTop β-lactam candidates:")
for _, row in blactam.head(5).iterrows():
    print(f"  {row['pdb_id'].upper()}  {row['pident']}%  {row['annotation'][:60]}")
PYEOF

# ══════════════════════════════════════════════════════════════════════════════
# STEP 6 — Foldseek Search: No-Hits vs AlphaFold/Swiss-Prot
# ══════════════════════════════════════════════════════════════════════════════
step "STEP 6 — Foldseek Search: No-Hit Proteins vs AlphaFold/Swiss-Prot"

N_NOHITS=$(ls "$NOHIT_DIR"/*.pdb 2>/dev/null | wc -l | tr -d ' ')
log "Searching $N_NOHITS no-hit structures against AlphaFold/Swiss-Prot..."

if [[ -f "${AFDB_DB}.dbtype" ]]; then
  foldseek easy-search \
    "$NOHIT_DIR" \
    "$AFDB_DB" \
    "$RESULTS_DIR/foldseek_nohits_afdb.tsv" \
    "$TMP_DIR" \
    --exhaustive-search 1 \
    --format-output "query,target,pident,alnlen,evalue,bits,lddt,prob" \
    -e "$EVALUE_THRESHOLD" \
    --threads "$THREADS"
  log "AlphaFold/Swiss-Prot search complete."
else
  warn "AlphaFold/Swiss-Prot DB not found. Skipping (run without --skip-download to fetch it)."
fi

# ══════════════════════════════════════════════════════════════════════════════
# STEP 7 — Ligand Download & Docking Prep
# ══════════════════════════════════════════════════════════════════════════════
step "STEP 7 — Download Ligands (Piperacillin + Tazobactam)"

if ! command -v vina &>/dev/null; then
  warn "AutoDock Vina not found — skipping docking. Install with: conda install -c conda-forge autodock-vina"
else
  if [[ ! -f "$DOCKING_DIR/piperacillin.sdf" ]]; then
    log "Downloading piperacillin (PubChem CID 43672)..."
    curl -s -o "$DOCKING_DIR/piperacillin.sdf" \
      "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/43672/record/SDF/?record_type=3d"
  fi

  if [[ ! -f "$DOCKING_DIR/tazobactam.sdf" ]]; then
    log "Downloading tazobactam (PubChem CID 123630)..."
    curl -s -o "$DOCKING_DIR/tazobactam.sdf" \
      "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/123630/record/SDF/?record_type=3d"
  fi

  # Convert to PDBQT
  for lig in piperacillin tazobactam; do
    if [[ ! -f "$DOCKING_DIR/${lig}.pdbqt" ]]; then
      log "Converting ${lig}.sdf → PDBQT..."
      obabel "$DOCKING_DIR/${lig}.sdf" \
        -O "$DOCKING_DIR/${lig}.pdbqt" \
        --gen3d -p 7.4 2>/dev/null
    fi
  done

  log "Ligands ready in $DOCKING_DIR"
  log "Next: prepare receptor PDBQTs and run vina manually for each candidate."
  log "See METHODS.md §6 for docking parameters."
fi

# ══════════════════════════════════════════════════════════════════════════════
# Done
# ══════════════════════════════════════════════════════════════════════════════
step "Pipeline Complete"
log "Results written to: $RESULTS_DIR"
echo ""
echo "  Key output files:"
echo "    $RESULTS_DIR/plddt_scores.tsv"
echo "    $RESULTS_DIR/foldseek_hits.tsv"
echo "    $RESULTS_DIR/foldseek_best_hits_annotated.tsv"
echo "    $RESULTS_DIR/blactam_candidates.tsv"
echo "    $RESULTS_DIR/dali_priority_list.txt"
echo "    $RESULTS_DIR/foldseek_nohits_afdb.tsv"
echo ""
echo "  Next manual steps:"
echo "    1. Submit dali_input/ top proteins to http://ekhidna2.biocenter.helsinki.fi/dali/"
echo "    2. Run AutoDock Vina on blactam_candidates.tsv top hits"
echo ""
