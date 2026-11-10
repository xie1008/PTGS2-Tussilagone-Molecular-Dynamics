# PTGS2–Tussilagone Molecular Dynamics Analysis

This repository contains the three 100-ns PTGS2–tussilagone molecular-dynamics trajectories used for the revised analysis, together with the reproducible analysis script, quantitative CSV outputs, provenance mapping, quality-control notes, and publication figures.

## Data provenance

The displayed replicas were selected from four completed simulations after the additional trajectory became available. This selection was made during reanalysis and was not a pre-registered selection. The exact mapping is provided in `provenance_mapping.csv`.

## Analysis definitions

- Protein RMSD: Cα atoms using a common protein reference definition.
- Ligand RMSD: 28 ligand heavy atoms after fitting Chain A residues 43–89; each displayed replica uses its own 50-ns reference frame.
- Protein–ligand contacts: heavy-atom pairs within 4.5 Å.
- Residue contact occupancy: frames with at least one qualifying contact divided by the total valid frames.
- Joint ligand-pose clustering: post-50-ns frames in a common pocket-fitted coordinate system using k-means with k=3.

The main analysis is implemented in `Scripts/selected_3rep_qc3.py`. Quantitative outputs are in `CSV/`, and the displayed trajectory-to-replica mapping is in `provenance_mapping.csv`.

## Reproducibility

The analysis script requires Python with NumPy, SciPy, Matplotlib, and the trajectory-analysis dependencies used in the script. The script was run on the selected trajectories and generated the CSV files and figures included here.
