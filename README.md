# ruv-pseudobulk

This repository contains the code used in the paper *“Removal of unwanted variation in pseudobulk analysis of single-cell RNA sequencing data and the leveraging of pseudoreplicates”*.

The datasets used in the analyses are available [here](https://doi.org/10.5281/zenodo.15341158)

## File Descriptions

- **auxf**: Auxiliary functions used across scripts.  
- **DEG_PC_control**: Code to identify differentially expressed genes between processing cohorts.  
- **FDRsim_T2noint**: Simulations for differential expression analysis (DEA) in the control subset; treatments are balanced with respect to processing cohort, and no interaction term is considered.  
- **FDRsim_T2nointC**: Simulations for DEA in the control subset; treatment is confounded with processing cohort, no interaction term is considered.
- - **FDRsim_intC**: Simulations for DEA in the control subset; treatment is confounded with processing cohort, and an interaction term is considered.  
- **FDRsim_PBPS**: Simulations for DEA in the control subset using PBPS; treatment is confounded with processing cohort, no interaction term is considered.  
- **FDRsim_PBPSint**: Simulations for DEA in the control subset using PBPS; treatment is confounded with processing cohort, and an interaction term is included.  
- **FDRsim_misspecified_ncg**: Simulations for DEA where the treatment is confounded with processing cohort and negative control genes are misspecified.  
- **generation PBPS control set**: Code to generate the PBPS used in the control analyses.  
- **Graphs descriptives control**: Code to reproduce graphs related to the normalized matrices.  
- **Graphs DEA control**: Code to reproduce the FDR-related graphs presented in the paper.  
- **PBPS**: Detailed description of the procedure to generate pseudobulk pseudosamples (PBPS).  
- **Trails_allcelltypes**: Contains results for the three trails and RUV methods from cell types not included in the paper.  
- **Study case**: Complete analysis of the study case, including figures.  
- **LupusDW**: Code used to derive the control and study case subsets.

