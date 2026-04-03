# amr_sponges
Structural Identification of Piperacillin/Tazobactam Sponge Candidates in *Escherichia coli*
# AMR-Sponge: Structural Identification of Piperacillin/Tazobactam Sponge Candidates in *Escherichia coli*

[![Python](https://img.shields.io/badge/python-3.10%2B-blue)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![ESMFold](https://img.shields.io/badge/structure-ESMFold-green)](https://esmatlas.com)
[![Foldseek](https://img.shields.io/badge/search-Foldseek-orange)](https://search.foldseek.com)

## Overview

This project identifies *E. coli* proteins that may act as **antibiotic sponges** — periplasmic or outer membrane proteins with structural homology to β-lactam-binding folds that could sequester piperacillin/tazobactam (TZP), titrating the drug away from its primary targets.

Resistance to TZP in *E. coli* is primarily driven by β-lactamase hyperproduction (e.g. *bla*TEM-1 amplification, OXA-1) and IRT β-lactamases. We hypothesize that an additional, underappreciated mechanism involves non-canonical proteins with β-lactam-accessible binding pockets acting as molecular decoys in the periplasm.

---

## Pipeline Overview
E. coli genomes (BV-BRC / NCBI)
│
▼
Candidate protein extraction (TZP-resistant isolates)
│
▼
Structure prediction — ESMFold API (276 proteins)
│
▼
pLDDT confidence filtering (≥ 0.70) → 171 high-confidence structures
│
▼
Structural search — Foldseek vs PDB100
│
├── 120 structures with hits → functional annotation
│ ├── β-lactam relevant hits (5) → AutoDock Vina docking
│ └── Functional categorization (10 categories)
│
└── 51 structures with NO hits → DALI / AlphaFold-SwissProt search
└── Novel fold candidates → priority for experimental 

---

## Key Findings

| Category | n | Top Hit | Identity |
|---|---|---|---|
| Pilus / Secretion System | 32 | FimD (3BWU) | 97.0% |
| Outer Membrane / Lipid Transport | 14 | TolB (6PNV) | 96.8% |
| Periplasmic Metabolism | 14 | Cyclophilin B (1VAI) | 99.3% |
| Redox / Metal Resistance | 10 | Cu,Zn SOD (1ESO) | 100.0% |
| **Cell Wall / Peptidoglycan** | **7** | **DpaA (8IKR)** | **100.0%** |
| No structural hit (novel?) | **51** | — | — |

### Top β-Lactam Sponge Candidates

| Protein | Best PDB Hit | Identity | Function |
|---|---|---|---|
| `NMDK01000077.1_19` | 8IKR DpaA | **100%** | YkuD transpeptidase — β-lactam-accessible active site |
| `CAJSHO010003323.1_1` | 6AZI | **94.8%** | D-Ala-D-Ala endopeptidase crystallized with β-lactam mimic |
| `BGTY01000005.1_212` | 6JFW PA0833 | 55.3% | Enlarged periplasmic PG-binding pocket |

---

## Requirements

```bash
# Core environment
conda create -n amr_sponge python=3.10
conda activate amr_sponge
pip install biopython requests pandas matplotlib

# Structure prediction
conda create -n colabfold python=3.10
conda activate colabfold
pip install colabfold[alphafold]

# Structural search
cd ~
curl -fsSL https://mmseqs.com/foldseek/foldseek-osx-universal.tar.gz -o foldseek.tar.gz
tar -xvf foldseek.tar.gz
export PATH="$HOME/foldseek/bin:$PATH"

# Docking
conda install -c conda-forge autodock-vina openbabel -y
```

---

## Usage

### 1. Structure Prediction (ESMFold API)
```bash
conda activate colabfold
cd amr_sponge_clean
python esmfold_hf.py   # uses ESM Atlas public API, no token required
```

### 2. pLDDT Filtering
```bash
python << 'EOF'
from pathlib import Path
import shutil

results = []
for pdb in sorted(Path('esmfold_output').glob('*.pdb')):
    b = [float(l[60:66]) for l in pdb.read_text().splitlines() if l.startswith('ATOM')]
    if b:
        results.append((pdb.stem, sum(b)/len(b)))

Path('esmfold_highconf').mkdir(exist_ok=True)
for name, score in results:
    if score >= 0.70:
        shutil.copy(f'esmfold_output/{name}.pdb', f'esmfold_highconf/{name}.pdb')
EOF
```

### 3. Foldseek Structural Search
```bash
# Download PDB database
foldseek databases PDB foldseek_db/pdb tmp/ --threads 8

# Run search
foldseek easy-search \
    esmfold_highconf/ \
    foldseek_db/pdb \
    results/foldseek_hits.tsv \
    tmp/ \
    --exhaustive-search 1 \
    --format-output "query,target,pident,alnlen,evalue,bits,lddt,prob" \
    -e 0.001 --threads 8
```

### 4. Functional Annotation & β-Lactam Filtering
```bash
python annotate_hits.py        # fetches titles from RCSB API
python filter_blactam.py       # filters for β-lactam relevant folds
```

### 5. Docking (AutoDock Vina)
```bash
# Download ligands
curl -o docking/piperacillin.sdf \
  "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/43672/record/SDF/?record_type=3d"
curl -o docking/tazobactam.sdf \
  "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/123630/record/SDF/?record_type=3d"

# Run blind docking against top candidates
python run_vina.py
```

---

## Directory Structure
amr_sponge_clean/
├── data/
│ └── candidates_filtered_A.faa # 276 candidate protein sequences
├── esmfold_output/ # All 276 predicted PDB structures
│ └── plddt_scores.tsv # Per-protein confidence scores
├── esmfold_highconf/ # 171 high-confidence structures (pLDDT ≥ 0.70)
├── dali_input/ # 51 no-hit structures for DALI search
├── foldseek_db/ # Local Foldseek databases
├── docking/ # Ligand files and docking results
└── results/
├── foldseek_hits.tsv # All Foldseek hits (9,735 total)
├── foldseek_best_hits.tsv # Best hit per query
├── foldseek_best_hits_annotated.tsv # With RCSB titles
├── blactam_candidates.tsv # β-lactam relevant hits
└── dali_priority_list.txt # No-hit proteins ranked by pLDDT

---

## Background

Piperacillin/tazobactam resistance in *E. coli* is predominantly mediated by β-lactamase hyperproduction via *bla*TEM-1 gene amplification or acquisition of OXA-1/IRT enzymes. [web:108][web:112] This project explores a complementary hypothesis: that periplasmic proteins with β-lactam-accessible transpeptidase or PBP-like folds can act as molecular sponges, sequestering TZP and reducing its effective periplasmic concentration.

The concept of sponge proteins as resistance mediators is well-established in phage-bacteria interactions, where phage-encoded sponge proteins sequester bacterial immune signals to evade defense systems. [web:115] We extend this concept to antibiotic resistance.

---

## Citation

If you use this pipeline, please cite:

- **ESMFold**: Lin et al., *Science* 2023
- **Foldseek**: van Kempen et al., *Nature Biotechnology* 2023
- **AutoDock Vina**: Eberhardt et al., *J. Chem. Inf. Model.* 2021

---

## License

MIT License — see [LICENSE](LICENSE) for details.

## Author

Computational Biology & Bioinformatics  
Copenhagen, Denmark
