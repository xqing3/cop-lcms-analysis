# COP LC-MS/MS analysis — code and figures

Analysis code and figure-generation scripts for an untargeted LC-MS/MS study
**identifying cholesterol oxidation products (COPs)**. Cholesterol was heated and the
resulting oxidation products were profiled in positive-ion mode by comparing heated vs.
unheated samples. Six injections (2 replicates × 3 groups) were analyzed: heated
1000 µg/mL, heated 100 µg/mL, and unheated 10 µg/mL.

This repository holds the **code and configuration** needed to reproduce the figures and
the molecular-networking analysis. The underlying raw and processed mass-spectrometry
**data are not stored here** — see [Data availability](#data-availability).

## Pipeline

Thermo `.raw` → **msconvert** → `.mzML` → **MZmine 4.9.14** (`mzmine/COPs_v2.mzbatch`):
feature detection, alignment, ion-identity molecular networking (IIMN), and exports →
**GNPS FBMN** + spectral-library annotation → figures.

The `.mzbatch` is XML and hard-codes absolute input paths to the original `.mzML` files;
adjust those paths before re-running in MZmine 4.9.14.

## Repository layout

```
figures/
  Figure1_notebook.ipynb   Fig 1 — feature counts (A), 3-way Venn (B), PCA panels (C, D)
  Figure3_notebook.ipynb   Fig 3 — molecular-network exploration notebook
  Figure3_network.py       Fig 3 — renders the FBMN/IIMN network using the Cytoscape layout
  rendered/                pre-rendered figure outputs (PNG/PDF)
fbmn/
  COP_style.xml            Cytoscape visual style used for the network figure
  iimn_collapsed.graphml   ion-identity-collapsed network (GNPS IIMN) read by Figure3_network.py
  results.cys              Cytoscape session providing the node layout
mzmine/
  COPs_v2.mzbatch          canonical MZmine 4.9.14 batch (XML)
data/
  output_full_feature_table.csv   master aligned feature table (input to Fig 1)
  output_annotations.csv          spectral-library hits (annotation flags for figures)
  output_metadata.csv             injection metadata
  metadata_updated.tsv            injection → condition / heat-treatment / concentration map
```

## Reproducing the figures

```bash
pip install -r requirements.txt

# Figure 3 (molecular network)
python figures/Figure3_network.py            # writes figures/Fig3_network.{png,pdf}

# Figure 1 (counts, Venn, PCA)
jupyter nbconvert --to notebook --execute figures/Figure1_notebook.ipynb
```

Both scripts resolve their inputs relative to the repository, so they run from a fresh
clone without editing paths. The fixed group palette is teal `#2CA58D` (1000 µg/mL),
amber `#E8A33D` (100 µg/mL), slate `#7E8B99` (unheated).

## Data availability

All raw and processed LC-MS/MS data — converted `.mzML` files, the MS/MS spectral summary
(`.mgf`), the feature quantification table, and sample metadata — are deposited in the
**MassIVE** repository under accession **MSV000099229**
(https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?accession=MSV000099229).

The **Feature-Based Molecular Networking** job is publicly available on **GNPS**:
https://gnps.ucsd.edu/ProteoSAFe/status.jsp?task=b2f1ab175a524a2c811f45fb713f9365

The small derived tables under `data/` are included here only as direct inputs to the
figure scripts.
