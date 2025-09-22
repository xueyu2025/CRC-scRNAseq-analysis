import scanpy as sc
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method
import pickle
adata = sc.read_h5ad('EpiImmune-all-rawcount.h5ad')
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata)
from cellphonedb.src.core.methods import cpdb_degs_analysis_method

cpdb_file_path = './v5.0.0/cellphonedb.zip'
meta_file_path = './SubEpi_subimmue_meta.txt'
counts_file_path = adata
degs_file_path = 'SubImmune-SubEpi-subcelltype_diff_genes_TvsN.txt'

cpdb_results = cpdb_degs_analysis_method.call(
    cpdb_file_path = cpdb_file_path,                            # mandatory: CellphoneDB database zip file.
    meta_file_path = meta_file_path,                            # mandatory: tsv file defining barcodes to cell label.
    counts_file_path = counts_file_path,                        # mandatory: normalized count matrix - a path to the counts file, or an in-memory AnnData object
    degs_file_path = degs_file_path,                            # mandatory: tsv file with DEG to account.
    counts_data = 'hgnc_symbol',                                # defines the gene annotation in counts matrix.
    #active_tfs_file_path = active_tf_path,                      # optional: defines cell types and their active TFs.
    #microenvs_file_path = microenvs_file_path,                  # optional (default: None): defines cells per microenvironment.
    score_interactions = True,                                  # optional: whether to score interactions or not. 
    threshold = 0.1,                                            # defines the min % of cells expressing a gene for this to be employed in the analysis.
    result_precision = 3,                                       # Sets the rounding for the mean values in significan_means.
    separator = '|',                                            # Sets the string to employ to separate cells in the results dataframes "cellA|CellB".
    debug = False,                                              # Saves all intermediate tables emplyed during the analysis in pkl format.
    output_path = './output',                                     # Path to save results
    output_suffix = None, #'deg_',                                       # Replaces the timestamp in the output files by a user defined string in the  (default: None)
    threads = 25
    )


with open('cpdb_results-DEG.pkl', 'wb') as f:
   pickle.dump(cpdb_results, f) 
print("Results saved to cpdb_results-DEG.pkl")
import pickle
with open('cpdb_results.pkl', 'rb') as f:
    cpdb_results = pickle.load(f)
print("Results loaded successfully!")

