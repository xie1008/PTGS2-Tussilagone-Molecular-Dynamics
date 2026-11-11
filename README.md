# PTGS2–Tussilagone Molecular Dynamics Analysis

This repository contains the three 100-ns PTGS2–tussilagone molecular-dynamics trajectories used for the revised analysis, together with the reproducible analysis script, quantitative CSV outputs, provenance mapping, quality-control notes, and publication figures.


### Trajectory Provenance and Replica Selection

During data organization, four trajectory files were initially identified. Subsequent verification of the simulation records revealed that the file originally designated as Replica 1 was an inadvertently included trajectory from a previous simulation run rather than a newly generated trajectory from the current set of independent simulations.

Therefore, the previous trajectory (Replica 1) was excluded from the current analysis based on its provenance, independently of its simulation results or binding-stability metrics.

The three newly generated trajectories, originally designated as Replicas 2, 3, and 4, were retained as the three independent replicates for the revised molecular dynamics analysis. For consistency with the manuscript and figures, these trajectories were relabeled as Replicas 1, 2, and 3, respectively.

The correspondence between the original trajectory identifiers and the final replicate labels is documented below.

| Original trajectory identifier | Final manuscript label | Status |
|---|---|---|
| Original Replica 1 | Not included | Previous simulation trajectory inadvertently included during data organization |
| Original Replica 2 | Replica 1 | Included |
| Original Replica 3 | Replica 2 | Included |
| Original Replica 4 | Replica 3 | Included |

No trajectory was excluded on the basis of ligand RMSD, protein RMSD, contact occupancy, binding stability, or other simulation outcomes.

## Analysis definitions

- Protein RMSD: Cα atoms using a common protein reference definition.
- Ligand RMSD: 28 ligand heavy atoms after fitting Chain A residues 43–89; each displayed replica uses its own 50-ns reference frame.
- Protein–ligand contacts: heavy-atom pairs within 4.5 Å.
- Residue contact occupancy: frames with at least one qualifying contact divided by the total valid frames.
- Joint ligand-pose clustering: post-50-ns frames in a common pocket-fitted coordinate system using k-means with k=3.

The main analysis is implemented in `Scripts/selected_3rep_qc3.py`. Quantitative outputs are in `CSV/`, and the displayed trajectory-to-replica mapping is in `provenance_mapping.csv`.

## Reproducibility

The analysis script requires Python with NumPy, SciPy, Matplotlib, and the trajectory-analysis dependencies used in the script. The script was run on the selected trajectories and generated the CSV files and figures included here.
