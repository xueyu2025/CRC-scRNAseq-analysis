# CRC-scRNAseq-analysis

Reproducible scripts for analyzing public single-cell datasets of colorectal cancer (CRC).  
Focus: metabolic crosstalk between SPP1+ macrophages and epithelial subsets.

## Environment
- **Python 3.10**: scanpy, anndata, pandas, numpy, matplotlib, seaborn
- **R 4.3**: Seurat, AUCell, GSVA/fgsea, NicheNet, inferCNV, tidyverse
Create env:
```bash
conda env create -f environment.yml
conda activate crc-scrna
# or: pip install -r requirements.txt
