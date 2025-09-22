#!/usr/bin/env python3
import os
import sys
import subprocess
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)
dbdir = './CRC/00_database/SCENIC'
TF = os.path.join(dbdir, 'hs_hgnc_tfs.txt')
inputloom = os.path.join(script_dir, '..', 'MyeloidCell-Exp-sampled.loom')
out1 = 'adj.sample.tsv'
out2 = 'reg.csv'
out3 = 'sample_SCENIC.loom'
def run(cmd: str):
    print(f"\n➔ Running:\n   {cmd}")
    full = (
        "bash -l -c "
        "\"source ~/mambaforge/etc/profile.d/conda.sh && "
        "conda activate r_env && "
        f"{cmd}\""
    )
    subprocess.run(full, shell=True, check=True)

def main():
    for f in [inputloom, TF]:
        if not os.path.isfile(f):
            print(f"Error: Required file not found: {f}", file=sys.stderr)
            sys.exit(1)

    # —— Step 1: GRN —— 
    grn_cmd = (
        f"pyscenic grn "
        f"--num_workers 20 "
        f"--output {out1} "
        f"--method grnboost2 "
        f"{inputloom} {TF}"
    )
    run(grn_cmd)

    # —— Step 2: Context —— 
    ctx_cmd = (
        f"pyscenic ctx "
        f"{out1} "
        f"{dbdir}/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather "
        f"--annotations_fname {dbdir}/motifs-v9-nr.hgnc-m0.001-o0.0.tbl "
        f"--expression_mtx_fname {inputloom} "
        f"--output {out2} "
        f"--mode dask_multiprocessing "
        f"--num_workers 3 "
        f"--mask_dropouts"
    )
    run(ctx_cmd)

    # —— Step 3: AUCell —— 
    aucell_cmd = (
        f"pyscenic aucell "
        f"{inputloom} "
        f"{out2} "
        f"--output {out3} "
        f"--num_workers 3"
    )
    run(aucell_cmd)

if __name__ == '__main__':
    main()
