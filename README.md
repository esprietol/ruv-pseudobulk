# ruv-pseudobulk

This is a repository for the code used in the paper "Removal of unwanted variation in pseudobulk analysis of single-cell RNA sequencing data and the leveraging of pseudoreplicates".

The datasets used are available at DOI: 10.5281/zenodo.15341158

File description:
- auxf: Auxiliary functions
- DEG_PC_control: Code to find the differentially expressed genes between processing cohorts
- FDRsim_T2noint: Simulations to perform DEA in the control subset, balanced desing with respect to processing cohort, no interaction term is considered.
- FDRsim_T2nointC: Simulations to perform DEA in the control subset, treatment confounded with processing cohort, no interaction term is considered
- FDRsim_PBPS: Simulations to perform DEA in the control subset with PBPS, treatment confounded with processing cohort, no interaction term is considered
- FDRsim_PBPSint: Simulations to perform DEA in the control subset with PBPS, treatment confounded with processing cohort, an interaction term is considered.
- FDRsim_misspecified_ncg: Simulations to perform DEA in the control subset, treatment confounded with processing cohort and negative control genes are misspecified.
- generation PBPS control set: Code to generate the PBPS used in the analysis
- Graphs descriptives control: Code to reproduce the graphs of the paper related to the normalised matrix
- Graphs DEA control: Code to reproduce the FDR graphs presented in the paper
- PBPS: Detailed description of the process to generate pseudobulk pseudosamples
- Trails_allcelltypes: Has the results for the different trails and RUV methods from the celltypes not included in the paper
- Study case: complete analysis, with graphs, of the study case
- LupusDW: Code used to derivate the subsets control and study case



