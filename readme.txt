#This code aim to provide a program to transform from h5ad to RDS file and keep umap/PCA reduction.
#Three file must in a same folder 
#Test in Python 3.9, R 4.3, Seurat V4.4 (management by conda)
#names: process_sc_data.sh
# function：from .h5ad -> MTX -> Seurat RDS file transform,and PCA / TSNE /UMAP will be keeped
# example: ./process_sc_data.sh -i a.h5ad -o output -p a