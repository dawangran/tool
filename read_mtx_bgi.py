##########################################################
# title:read rna velocity matrix for bgi
# date: 2023-07-20
# author: dawn
##########################################################


import os
import argparse
import scanpy as sc
from scipy.sparse import csr_matrix
from scipy.io import mmread
import numpy as np
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description='Join data path and sample name')
    parser.add_argument('-i', '--input', required=True, help='input path')
    parser.add_argument('-n', '--name', required=True, help='sample name')
    parser.add_argument('-o', '--output', required=True, help='output path')
    args = parser.parse_args()
    return args


def read_bgi_mtx(path,sample_name,out_dir):
    if os.path.exists(path):
        file_list=os.listdir(path)
        for file in file_list:
            file_dir=os.path.join(path,file)
            if file =="barcodes.tsv.gz":
                obs = pd.read_csv(file_dir, header=None, index_col=0, sep="\t")
                obs.index.name = "barcode"
            if file =="features.tsv.gz":
                var = pd.read_csv(file_dir, header=None, index_col=0, sep="\t")
                var.index.name = "gene"
            if file =="spliced.mtx.gz":
                splice_mtx=csr_matrix(mmread(file_dir).T)
            if file =="unspliced.mtx.gz":
                unsplice_mtx=csr_matrix(mmread(file_dir).T)
        mtx = splice_mtx + unsplice_mtx
        adata=sc.AnnData(mtx,obs=obs,var=var)
        adata.obs['sample']=sample_name
        adata.var['gene_ids']=adata.var_names
        adata.layers['spliced']=splice_mtx
        adata.layers['unspliced']=unsplice_mtx
        names=adata.obs_names.tolist()
        samples=adata.obs['sample'].tolist()
        barcodes=[i+"_"+j for i,j in zip(samples,names)]
        adata.obs_names=barcodes
        del adata.obs['sample']
        file_name=sample_name+'.h5ad'
        out_path=os.path.join(out_dir,file_name)
        adata.write(out_path)
        
if __name__=="__main__":
    
    args = parse_args()
    read_bgi_mtx(args.input,args.name,args.output)
    
