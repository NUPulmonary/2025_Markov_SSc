# Profibrotic monocyte-derived alveolar macrophages as a biomarker and therapeutic target in systemic sclerosis-associated interstitial lung disease
https://www.biorxiv.org/content/10.1101/2025.08.07.669006

This repository contains analysis code for this study. Please, direct any questions to the corresponding authors of the study, or feel free to open an issue in this repository.

The code was run with Python 3.9, and we used 4 virtual environments for that: `default, clin, scvi, cellrank`. The first 3 are managed by [pip-tools](https://pypi.org/project/pip-tools/), the last my mamba. The `.in` and `requirements.txt` files describe packages in these environments.

- 00_clinical: (clin), clinical data preprocessing, aggregation, plots & tables
- 01_ssc-bal: (default), scRNAseq samples processing with 10x cellranger
- 02_duke: (default), scRNAseq samples processing
- 03_bal-object: (varies per notebook), scRNAseq BAL object integration, analysis and plots
