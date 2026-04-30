# Proteomic signatures of the *APOE* ε 4 and *APOE* ε 2 genetic variants and Alzheimer's disease
This repository contains code and result files for the analysis of *APOE* genotype-related proteomic signatures in plasma and CSF across multiple cohorts. It includes scripts for data preprocessing, statistical analysis, visualization, and result summarization.

# Files
```bash
.
├── analysis
│   ├── ad_markers.R        # APOE-related proteins and downstream AD phenotypes
│   ├── main_local.R        # Main analysis.
│   ├── main_plots.R        # Plotting scripts (for example, boxplots)
│   ├── sensitivity.R       # Sensitivity analyses
│   ├── setting.R           # Analysis settings for main_local.R
│   ├── snp.R               # SNP summary for coding genes of key mediators
│   └── src.R               # Functions for statistics and visualization
├── generate_ctd.ipynb
├── r2correlation.ipynb     # Direct comparison of measurements within a cohort
├── results_summarize.ipynb # Result summarization and cross-cohort comparison
├── pyproject.toml
├── pysrc                   # Python package for result summarization
│   └── sumAPOE
│       ├── __init__.py
│       ├── loader.py
│       └── utils.py
└── results
    ├── ADNI
    ├── ADNI_MS
    ├── BF2OLINK
    ├── BF2SomaLogic
    ├── GNPC
    ├── ROSMAP
    └── UKBB
```
# Usage
All model results are available in the results folder.

To explore the results in results_summarize.ipynb, the sumAPOE package must first be installed. From the repository root, run:
```bash
cd APOE_Proteomics #or in jupyter notebook
pip install -e . 
```
For analyses and plots in R, creating an RStudio Project for this repository is recommended.

# Reference
Lu, L. et al. Proteomic signatures of the *APOE* ε4 and *APOE* ε2 genetic variants and Alzheimer’s disease. medRxiv 2025.08.04.25332945 (2025) doi:10.1101/2025.08.04.25332945.

# Contact
If you encounter any issues, bugs, or have questions about the code or results, please contact or open an issue:
Lina Lu  
Department of Clinical Sciences Malmö, Lund University  
Email: lina.lu@med.lu.se