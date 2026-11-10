# Trajectory and analysis quality-control record

Three displayed trajectories are included in the repository as 100-ns ZIP archives. Their display numbering and original simulation numbering are recorded in `provenance_mapping.csv`.

The selected simulations were run with the same general OpenMM protocol but on different platforms: displayed Replica 1 originated from the Apple M4/OpenCL single-precision run, whereas displayed Replicas 2 and 3 originated from Tesla T4/OpenCL mixed-precision runs.

The analysis used the complete 0–100 ns trajectories. No time interval was removed to improve the trends, and no replicate was excluded from the displayed three-replicate dataset after the stated selection mapping was defined.

The ligand RMSD reference is the frame nearest 50 ns within each displayed trajectory. The ligand selection contains 28 heavy atoms, and protein fitting uses Chain A residues 43–89. The reference-frame self-RMSD check is part of the analysis workflow.

The repository CSV files provide the numerical outputs used for the RMSD, contact, clustering, and representative-pose comparisons. Persistent contacts, within-trajectory stability, and cross-replicate pose convergence are treated as distinct properties and are not interpreted as direct measures of binding affinity.
