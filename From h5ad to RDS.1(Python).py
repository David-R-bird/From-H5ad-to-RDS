import scanpy as sc
import pandas as pd
import numpy as np
from scipy import sparse
import os
import sys
import argparse 

def main():
    # 创建参数解析器
    parser = argparse.ArgumentParser(description='将 .h5ad 文件转换为 MTX 格式')
    parser.add_argument('-i', '--input', required=True, help='输入的 .h5ad 文件路径')
    parser.add_argument('-o', '--output', default='./transform/', help='输出目录')
    parser.add_argument('-p', '--prefix', help='输出文件前缀 (默认为输入文件名前缀)')
    parser.add_argument('-q', '--quiet', action='store_true', help='静默模式')
    
    args = parser.parse_args()  
   

    # 1. 读取数据
    input_file = args.input
    output_dir = args.output
    file_prefix = args.prefix


    # 1. 确保输出目录存在
    os.makedirs(output_dir, exist_ok=True)
    print(f"输出目录: {output_dir} (已存在)")
    print("正在读取 h5ad 文件...")
    # 使用 os.path 模块安全地处理文件名
    file_name = os.path.basename(input_file)  # 获取文件名（去掉路径）
    file_prefix = os.path.splitext(file_name)[0]  # 去掉扩展名，获取前缀
    print(f"输入文件: {input_file}")
    adata = sc.read_h5ad(input_file) # 请确保文件路径正确
    print(f"数据结构: {adata.shape} (细胞数 x 基因数)")

    # 2. 处理并保存细胞注释 (obs)
    print("正在处理并保存细胞注释 (obs)...")
    obs_df = adata.obs.copy()
    # 关键：将所有列强制转换为字符串，避免任何因子/整数/浮点数类型问题
    for col in obs_df.columns:
        obs_df[col] = obs_df[col].astype(str)
    # 保存为制表符分隔的文本文件，R 的 read.delim 读取友好
    obs_df.to_csv(os.path.join(output_dir, f'{file_prefix}_cell_metadata.tsv'), sep='\t', index=True) # index=True 保存细胞ID
    print(f"细胞注释已保存，共 {obs_df.shape[1]} 列。")

    # 3. 处理并保存基因注释 (var)
    print("正在处理并保存基因注释 (var)...")
    var_df = adata.var.copy()
    for col in var_df.columns:
        var_df[col] = var_df[col].astype(str)
    var_df.to_csv(os.path.join(output_dir, f'{file_prefix}_gene_metadata.tsv'), sep='\t', index=True) # index=True 保存基因名
    print(f"基因注释已保存，共 {var_df.shape[1]} 列。")

    # 4. 处理并保存表达矩阵 (X)
    print("正在处理并保存表达矩阵 (X)...")
    # 获取矩阵
    X = adata.X
    # 判断是否为稀疏矩阵
    if sparse.issparse(X):
        print("检测到稀疏矩阵，保存为Matrix Market格式(.mtx)...")
        
        # 1. 保存稀疏矩阵为MTX格式
        from scipy.io import mmwrite
        mtx_path = os.path.join(output_dir, f'{file_prefix}_expression_matrix.mtx')
        # 保存转置，使格式为: 基因×细胞
        mmwrite(mtx_path, X.T)
        
        # 2. 保存基因名（行名）
        genes_path = os.path.join(output_dir, f'{file_prefix}_genes.tsv')
        with open(genes_path, 'w') as f:
            for gene in adata.var_names:
                f.write(f"{gene}\n")
        
        # 3. 保存细胞名（列名）
        cells_path = os.path.join(output_dir, f'{file_prefix}_barcodes.tsv')
        with open(cells_path, 'w') as f:
            for cell in adata.obs_names:
                f.write(f"{cell}\n")
        
        print("已保存MTX格式文件: .mtx, _genes.tsv, _barcodes.tsv")
        
    else:
        print("检测到稠密矩阵，转换为稀疏矩阵并保存为MTX格式...")
        
        # 将稠密矩阵转换为稀疏矩阵
        X_sparse = sparse.csr_matrix(X)
        
        # 1. 保存为MTX格式
        from scipy.io import mmwrite
        mtx_path = os.path.join(output_dir, f'{file_prefix}_expression_matrix.mtx')
        # 保存转置
        mmwrite(mtx_path, X_sparse.T)
        
        # 2. 保存基因名（行名）
        genes_path = os.path.join(output_dir, f'{file_prefix}_genes.tsv')
        with open(genes_path, 'w') as f:
            for gene in adata.var_names:
                f.write(f"{gene}\n")
        
        # 3. 保存细胞名（列名）
        cells_path = os.path.join(output_dir, f'{file_prefix}_barcodes.tsv')
        with open(cells_path, 'w') as f:
            for cell in adata.obs_names:
                f.write(f"{cell}\n")
        
    
        
    print("表达矩阵已保存。")
    # 5. (可选) 保存降维坐标等信息

    if 'X_umap' in adata.obsm_keys():
        print("正在保存 UMAP 坐标...")
        pd.DataFrame(adata.obsm['X_umap'], index=adata.obs_names).to_csv(os.path.join(output_dir, f'{file_prefix}_umap_coordinates.tsv'), sep='\t')
        print("已保存 UMAP 或 PCA 坐标。")
    if 'X_tsne' in adata.obsm_keys():
        print("正在保存 TSNE 坐标...")
        pd.DataFrame(adata.obsm['X_tsne'], index=adata.obs_names).to_csv(os.path.join(output_dir, f'{file_prefix}_tsne_coordinates.tsv'), sep='\t')
        print("已保存 TSNE 坐标。")
    if 'X_pca' in adata.obsm_keys():
        print("正在保存 PCA 坐标...")
        pd.DataFrame(adata.obsm['X_pca'], index=adata.obs_names).to_csv(os.path.join(output_dir, f'{file_prefix}_pca_coordinates.tsv'), sep='\t')
        print("已保存 UMAP 或 PCA 坐标。")
    else:    
        print("未找到 UMAP 或 PCA 坐标，跳过保存。")

    print(f"临时文件已保存至{output_dir}，请使用 R 脚本进行后续处理。")


if __name__ == "__main__":
    sys.exit(main())