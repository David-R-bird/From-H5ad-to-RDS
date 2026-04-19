
#
read_mtx_datasets <- function(file_prefix, output_dir = "output", verbose = TRUE) {
  
  #检查必需的R包
  if (verbose) message("正在检查必要的R包...")
  
  # 检查Seurat包
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("错误: 需要安装Seurat包。请运行: install.packages('Seurat')")
  } else {
    if (verbose) message("✅ Seurat包已安装")
  }
  
  # 检查Matrix包
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("错误: 需要安装Matrix包。请运行: install.packages('Matrix')")
  } else {
    if (verbose) message("✅ Matrix包已安装")
  }
  
  # 加载必需的包
  library(Seurat)
  library(Matrix)
  
  # 输入验证
  if (!is.character(file_prefix)) {
    stop("file_prefix 必须是字符向量")
  }
  
  if (length(file_prefix) == 0) {
    stop("file_prefix 不能为空")
  }
  
  if (!dir.exists(output_dir)) {
    stop("输出目录不存在: ", output_dir)
  }
  
  
  # 1. 读取表达矩阵 (.mtx + 行列名)
  
  mtx_file <- paste0(output_dir, "/", file_prefix, "_expression_matrix.mtx")
  counts_matrix <- readMM(mtx_file)
  
  # 读取基因名（行名）
  gene_file <- paste0(output_dir, "/", file_prefix, "_genes.tsv")
  rownames(counts_matrix) <- readLines(gene_file)
  
  # 读取细胞名（列名）
  cell_file <- paste0(output_dir, "/", file_prefix, "_barcodes.tsv")
  colnames(counts_matrix) <- readLines(cell_file)
  
  #print(counts_matrix[1:3, 1:3])
  
  # 2. 读取细胞注释
  
  cell_meta_file <- paste0(output_dir, "/", file_prefix, "_cell_metadata.tsv")
  cell_metadata <- read.delim(cell_meta_file, row.names = 1, stringsAsFactors = FALSE)
  
  # 3. 读取基因注释
  gene_meta_file <- paste0(output_dir, "/", file_prefix, "_gene_metadata.tsv")
  gene_metadata <- read.delim(gene_meta_file, row.names = 1, stringsAsFactors = FALSE)
  
  #  5. 创建 Seurat 对象
  
  sce <- CreateSeuratObject(
    counts = counts_matrix,
    meta.data = cell_metadata,
    assay = "RNA",
    project = file_prefix
  )
  
  # 6. 添加基因注释
  
  sce[["RNA"]]@meta.features <- gene_metadata
  
  
  # 7. 添加UMAP坐标（如果存在）
  umap_file <- paste0(output_dir, "/", file_prefix, "_umap_coordinates.tsv")
  if (file.exists(umap_file)) {
    umap_coords <- read.delim(umap_file, row.names = 1, stringsAsFactors = FALSE)
    #只能识别"UMAP_1","UMAP_2"
    colnames(umap_coords) <- c("UMAP_1","UMAP_2")
    
    # 确保UMAP坐标中的细胞与Seurat对象中的细胞匹配
    common_cells <- intersect(colnames(sce),rownames(umap_coords))
    is_identical <- identical(rownames(umap_coords),colnames(sce))
    if (length(common_cells) > 0) {
      
      # 重新排序以匹配
      umap_coords <- umap_coords[common_cells, , drop = FALSE]
      
      # 创建DimReduc对象
      sce <- sce[, common_cells] 
      sce[["umap"]] <- CreateDimReducObject(
        embeddings = as.matrix(umap_coords),
        key = "UMAP_",
        assay = "RNA"
      )
      
      if (is_identical) { 
        message(paste0("已添加UMAP,包括所有细胞共",nrow(umap_coords), "个细胞,"
                       ,ncol(umap_coords), "个UMAP维度"))
      }
      else {
        if (length(common_cells) < ncol(sce)){
          message(paste0("已添加UMAP,删除部分不匹配细胞,共",nrow(umap_coords), "个细胞,"
                         ,ncol(umap_coords), "个UMAP维度"))
        }
        else {
          message(paste0("已添加UMAP,包括所有细胞，但顺序已调整,共",
                         nrow(umap_coords), "个细胞,",ncol(umap_coords), " 个UMAP维度"))
        }
      }
    }     else {
      warning("UMAP坐标中的细胞与Seurat对象不匹配，跳过添加UMAP坐标")
    }
  } 
  
  else {
    message("未找到UMAP坐标文件: ", basename(umap_file))
  }

  # 8. 添加Tsne坐标（如果存在）
  tsne_file <- paste0(output_dir, "/", file_prefix, "_tsne_coordinates.tsv")
  if (file.exists(tsne_file)) {
    tsne_coords <- read.delim(tsne_file, row.names = 1, stringsAsFactors = FALSE)
    #只能识别"TSNE_1","TSNE_2"
    colnames(tsne_coords) <- c("TSNE_1","TSNE_2")
    
    # 确保UMAP坐标中的细胞与Seurat对象中的细胞匹配
    common_cells <- intersect(colnames(sce),rownames(tsne_coords))
    is_identical <- identical(rownames(tsne_coords),colnames(sce))
    if (length(common_cells) > 0) {
      
      # 重新排序以匹配
      tsne_coords <- tsne_coords[common_cells, , drop = FALSE]
      
      # 创建DimReduc对象
      sce <- sce[, common_cells] 
      sce[["tsne"]] <- CreateDimReducObject(
        embeddings = as.matrix(tsne_coords),
        key = "TSNE_",
        assay = "RNA"
      )
      
      if (is_identical) { 
        message(paste0("已添加TSNE,包括所有细胞共",nrow(tsne_coords), "个细胞,"
                       ,ncol(tsne_coords), "个TSNE维度"))
      }
      else {
        if (length(common_cells) < ncol(sce)){
          message(paste0("已添加TSNE,删除部分不匹配细胞,共",nrow(tsne_coords), "个细胞,"
                         ,ncol(tsne_coords), "个TSNE维度"))
        }
        else {
          message(paste0("已添加TSNE,包括所有细胞，但顺序已调整,共",
                         nrow(tsne_coords), "个细胞,",ncol(tsne_coords), " 个TSNE维度"))
        }
      }
    }     else {
      warning("TSNE坐标中的细胞与Seurat对象不匹配，跳过添加TSNE坐标")
    }
  } 
  
  else {
    message("未找到TSNE坐标文件: ", basename(tsne_file))
  }
  
  # 9. 添加PCA坐标（如果存在）
  pca_file <- paste0(output_dir, "/", file_prefix, "_pca_coordinates.tsv")
  if (file.exists(pca_file)) {
    pca_coords <- read.delim(pca_file, row.names = 1, stringsAsFactors = FALSE)
    #只能识别"PCA_1","PCA_2"
    #colnames(pca_coords) <- c("PCA_1","PCA_2")
    
    # 确保PCA坐标中的细胞与Seurat对象中的细胞匹配
    common_cells <- intersect(colnames(sce),rownames(pca_coords))
    is_identical <- identical(rownames(pca_coords), colnames(sce))
    if (length(common_cells) > 0) {
      
      # 重新排序以匹配
      pca_coords <- pca_coords[common_cells, , drop = FALSE]
      
      # 创建DimReduc对象
      sce <- sce[, common_cells]
      sce[["pca"]] <- CreateDimReducObject(
        embeddings = as.matrix(pca_coords),
        key = "PC_",
        assay = "RNA"
      )
      
      if (is_identical) { 
        message(paste0("已添加PCA,包括所有细胞共", nrow(pca_coords), "个细胞,"
                       , ncol(pca_coords), "个PCA维度"))
      }
      else {
        if (length(common_cells) < ncol(sce)){
          message(paste0("已添加PCA,,删除部分不匹配细胞,共", nrow(pca_coords), "个细胞,"
                         , ncol(pca_coords), "个PCA维度"))
        }
        else {
          message(paste0("已添加PCA,包括所有细胞，但顺序已调整,共",
                         nrow(pca_coords), "个细胞,", ncol(pca_coords)," 个PCA维度"))
        }
      }
    }     else {
      warning("PCA坐标中的细胞与Seurat对象不匹配，跳过添加PCA坐标")
    }
  } 
  
  else {
    message("未找到PCA坐标文件: ", basename(pca_file))
  }
  
  saveRDS(sce,paste0(output_dir, "/", file_prefix, ".RDS"))
  message("数据集已保存至: ",paste0(output_dir, "/", file_prefix, ".RDS"))
  
}