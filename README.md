# CRC-scRNAseq-analysis

Reproducible scripts for analyzing public single-cell datasets of colorectal cancer (CRC).  
Focus: metabolic crosstalk between SPP1+ macrophages and epithelial subsets.

## Data
Public datasets only (please cite originals): GSE178341, GSE132256, GSE144735, GSE188711, GSE125553, GSE132257; bulk TCGA-COAD.  
No patient-identifiable data are included.

## Environment
- **Python 3.10**: scanpy, anndata, pandas, numpy, matplotlib, seaborn
- **R 4.3**: Seurat, AUCell, GSVA/fgsea, NicheNet, inferCNV, tidyverse
Create env:
```bash
conda env create -f environment.yml
conda activate crc-scrna
# or: pip install -r requirements.txt
