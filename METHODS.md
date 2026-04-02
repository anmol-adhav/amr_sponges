# Methods

## 1. Genome Selection and Candidate Extraction

*Escherichia coli* genomes with documented resistance to piperacillin/tazobactam (TZP) were retrieved from the BV-BRC and NCBI databases. Protein sequences were extracted from TZP-resistant isolates (taxon IDs with prefix `562.`) and filtered to produce a non-redundant candidate set. The final input comprised **276 protein sequences** stored in `data/candidates_filtered_A.faa`.

---

## 2. Protein Structure Prediction

### 2.1 ESMFold via ESM Atlas API
Three-dimensional structures were predicted for all 276 candidate sequences using ESMFold v1 via the ESM Atlas public REST API (`https://api.esmatlas.com/foldSequence/v1/pdb/`) with no authentication required. Sequences were truncated to a maximum of 400 residues. Predictions were performed sequentially with a 3-second inter-request delay to respect server rate limits. All predicted structures were saved in PDB format.
Endpoint : https://api.esmatlas.com/foldSequence/v1/pdb/
Method : POST
Headers : Content-Type: application/x-www-form-urlencoded
Max len : 400 residues

### 2.2 Confidence Filtering (pLDDT)
Per-residue predicted local distance difference test (pLDDT) scores were extracted from the B-factor column (columns 61–66) of each predicted PDB file. Mean pLDDT was computed per protein. Structures with mean pLDDT ≥ 0.70 (on a 0–1 scale as returned by the ESM Atlas API) were retained for downstream analysis.

- Total predicted: **276**
- High-confidence (pLDDT ≥ 0.70): **171**
- Low-confidence (excluded): **105**

---

## 3. Structural Homology Search

### 3.1 Foldseek vs PDB100
Structural similarity searches were performed using Foldseek v8 (commit `8dc75c74`) against the PDB100 database (downloaded April 2026, ~2.1 GB). All 171 high-confidence structures were queried in a single batch using the `easy-search` workflow with exhaustive search mode.

```bash
foldseek easy-search \
    esmfold_highconf/ \
    foldseek_db/pdb \
    results/foldseek_hits.tsv \
    tmp/ \
    --exhaustive-search 1 \
    --format-output "query,target,pident,alnlen,evalue,bits,lddt,prob" \
    -e 0.001 \
    --threads 8
```

**Search parameters:**
| Parameter | Value |
|---|---|
| E-value threshold | 0.001 |
| Search mode | Exhaustive (3Di + amino acid) |
| Database | PDB100 (April 2026) |
| Threads | 8 |

A total of **9,735 hits** were returned across 171 queries. The best hit per query (highest bitscore) was selected, yielding **120 proteins with at least one structural match** and **51 with no hit**.

### 3.2 Foldseek vs AlphaFold/Swiss-Prot
The 51 no-hit proteins were re-searched against the AlphaFold/Swiss-Prot database (~2 GB) using identical parameters to improve sensitivity for divergent homologs.

### 3.3 DALI Server (Validation)
The top no-hit candidates by pLDDT score were submitted to the DALI server (http://ekhidna2.biocenter.helsinki.fi/dali/) for cross-validation, searching against the PDB and AlphaFold databases. DALI uses a different structural comparison algorithm (distance matrix alignment) from Foldseek, providing orthogonal sensitivity for distant structural homologs.

---

## 4. Functional Annotation

PDB identifiers from Foldseek hits were annotated using the RCSB PDB REST API:
GET https://data.rcsb.org/rest/v1/core/entry/{PDB_ID}
Field: struct.title

Proteins were classified into 10 functional categories based on keyword matching of RCSB structure titles:

| Category | Keywords |
|---|---|
| Pilus / Secretion System | pilus, fim, gsp, usher, pilin, type 2/3/4 secretion |
| Outer Membrane / Lipid Transport | outer membrane, omp, tolb, mla, lpt, lipoprotein |
| Cell Wall / Peptidoglycan | transpeptidase, pbp, peptidoglycan, d-ala, endopeptidase |
| Redox / Metal Resistance | superoxide, thioredoxin, copper, merp, cus, dsb |
| Periplasmic Metabolism / Binding | cyclophilin, periplasmic, glucosidase, kinase, lyase |
| Flagella / Motility | flagell, motb, flil, flga |
| Signaling / Two-Component | sensor, histidine kinase, response regulator |
| Phage / Lysozyme Defense | ivy, lysozyme, immunity, inhibitor |
| Hypothetical / Unknown | hypothetical, uncharacterized, duf |
| Other | all remaining |

---

## 5. β-Lactam Candidate Prioritization

Foldseek hits were filtered for β-lactam relevance using the following keywords against the RCSB annotation:
transpeptidase, penicillin, pbp, beta-lactam, lactam,
dd-peptidase, endopeptidase, cell wall, peptidoglycan,
d-ala, carboxypeptidase, ykud, dpa

This yielded **5 candidate proteins** with structural homology to experimentally characterized β-lactam-binding folds. The top two candidates for docking were:

1. **`562.22037_NMDK01000077.1_19`** — 100% identity to DpaA (PDB: 8IKR), a YkuD-family transpeptidase with a β-lactam-accessible active site
2. **`562.90314_CAJSHO010003323.1_1`** — 94.8% identity to a D-Ala-D-Ala endopeptidase from *Enterobacter cloacae* (PDB: 6AZI) co-crystallized with a boronic acid β-lactam mimic

---

## 6. Molecular Docking

### 6.1 Ligand Preparation
3D conformers of piperacillin (PubChem CID: 43672) and tazobactam (PubChem CID: 123630) were downloaded from PubChem in SDF format and converted to PDBQT format using OpenBabel:

```bash
obabel piperacillin.sdf -O piperacillin.pdbqt --gen3d
obabel tazobactam.sdf -O tazobactam.pdbqt --gen3d
```

### 6.2 Receptor Preparation
Predicted PDB structures of the top β-lactam candidates were prepared for docking by removing water molecules and adding Gasteiger charges using AutoDockTools/OpenBabel. Hydrogen atoms were added at pH 7.4.

### 6.3 AutoDock Vina Blind Docking
Blind docking (whole-protein search box) was performed using AutoDock Vina to avoid bias toward a predicted binding site:

```bash
vina \
    --receptor receptor.pdbqt \
    --ligand ligand.pdbqt \
    --center_x <cx> --center_y <cy> --center_z <cz> \
    --size_x <sx> --size_y <sy> --size_z <sz> \
    --exhaustiveness 32 \
    --num_modes 9 \
    --out docking/output.pdbqt \
    --log docking/output.log
```

Binding poses were ranked by predicted ΔG (kcal/mol). Poses within the putative active site (within 5 Å of the conserved serine in YkuD-family transpeptidases) were prioritized for visualization.

---

## 7. Software Versions

| Tool | Version | Reference |
|---|---|---|
| Python | 3.10 | |
| ESMFold (API) | v1 | Lin et al., *Science* 2023 |
| Foldseek | 8dc75c74 | van Kempen et al., *Nat. Biotechnol.* 2023 |
| AutoDock Vina | 1.2.x | Eberhardt et al., *J. Chem. Inf. Model.* 2021 |
| OpenBabel | 3.x | O'Boyle et al., *J. Cheminform.* 2011 |
| BioPython | 1.82 | Cock et al., *Bioinformatics* 2009 |
| pandas | 1.5.3 | |
| numpy | 1.24.3 | |
