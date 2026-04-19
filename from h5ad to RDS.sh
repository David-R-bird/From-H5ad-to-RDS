#!/bin/bash
# 脚本名称: process_sc_data.sh
# 功能：自动化串联 .h5ad -> MTX -> Seurat 对象的转换流程
# 使用示例: ./process_sc_data.sh -i sted-ec.h5ad -o output -p sted-ec
# 作者: David
# 建议使用conda管理环境
# 已验证运行环境: Linux/MacOS，安装了 Python 3.9 和 R 4.3，并且安装了必要的包（scanpy, pandas, numpy, scipy, Seurat(V4), Matrix）

# ==================== 1. 参数解析与配置 ====================
# 设置默认值
INPUT_H5AD=""
OUTPUT_DIR="output"
FILE_PREFIX=""
VERBOSE=true
PYTHON_SCRIPT="From h5ad to RDS.1(Python).py"
REMOVE_INTERMEDIATE=true  # 添加这行：默认删除中间文件
R_SCRIPT="From h5ad to RDS.2(R).R"

# 用法说明
usage() {
    echo "用法: $0 -i <input.h5ad> [-o <output_dir>] [-p <file_prefix>] [-v]"
    echo "选项:"
    echo "  -i, --input      输入的 .h5ad 文件路径 (必需)"
    echo "  -o, --output     输出目录 (默认: 'output') (必需)"
    echo "  -p, --prefix     输出文件前缀 (默认: 从输入文件名提取)"
    echo "  -q, --quiet      静默模式，减少输出"
    echo "  -keep, --keep    保留中间文件（默认: 删除）"
    echo "  -h, --help       显示此帮助信息"
    exit 1
}

# 解析命令行参数
while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input)
            INPUT_H5AD="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -p|--prefix)
            FILE_PREFIX="$2"
            shift 2
            ;;
        -q|--quiet)
            VERBOSE=false
            shift
            ;;
        -keep|--keep)  
            REMOVE_INTERMEDIATE=false
            shift
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "未知选项: $1"
            usage
            ;;
    esac
done

# 日志函数
log_info() {
    if [ "$VERBOSE" = true ]; then
        echo "[INFO] $(date '+%Y-%m-%d %H:%M:%S') - $1"
    fi
}

log_error() {
    echo "[ERROR] $(date '+%Y-%m-%d %H:%M:%S') - $1" >&2
}

# ==================== 2. 输入验证 ====================
log_info "开始转换处理流程"

# 检查必需参数
if [ -z "$INPUT_H5AD" ]; then
    log_error "必须指定输入文件 (-i 参数)"
    usage
fi

if [ ! -f "$INPUT_H5AD" ]; then
    log_error "输入文件不存在: $INPUT_H5AD"
    exit 1
fi

# 如果没有指定前缀，从输入文件名提取
if [ -z "$FILE_PREFIX" ]; then
    FILE_PREFIX=$(basename "$INPUT_H5AD" .h5ad)
    log_info "使用自动提取的文件前缀: $FILE_PREFIX"
fi

# 创建输出目录
mkdir -p "$OUTPUT_DIR"
log_info "输出目录: $(realpath "$OUTPUT_DIR")"

# ==================== 3. 环境检查 ====================
log_info "检查运行环境..."

# 检查 Python
if ! command -v python3 &> /dev/null; then
    log_error "未找到 python3，请确保 Python 3 已安装"
    exit 1
fi
log_info "Python 3 已安装: $(python3 --version)"

# 检查 R
if ! command -v Rscript &> /dev/null; then
    log_error "未找到 Rscript，请确保 R 已安装"
    exit 1
fi
log_info "R 已安装: $(Rscript --version 2>&1 | head -1)"

# 检查必要的 Python 包
log_info "检查 Python 依赖..."
if ! python3 -c "import scanpy, pandas, numpy, scipy, os" 2>/dev/null; then
    log_error "缺少必要的 Python 包。请运行: pip install scanpy pandas numpy scipy"
    exit 1
fi
log_info "Python 依赖检查通过"

# 检查必要的 R 包
log_info "检查 R 依赖..."
if ! Rscript -e "if(!require('Seurat')) stop('Seurat 包未安装')" 2>/dev/null; then
    log_error "Seurat 包未安装。请运行: install.packages('Seurat'),(建议使用conda安装R和Seurat，安装命令: conda install -c bioconda r-base r-seurat)"
    exit 1
fi
if ! Rscript -e "if(!require('Matrix')) stop('Matrix 包未安装')" 2>/dev/null; then
    log_error "Matrix 包未安装。请运行: install.packages('Matrix'),(建议使用conda安装R和Seurat，安装命令: conda install -c bioconda r-Matrix)"
    exit 1
fi
log_info "R 依赖检查通过"

# ==================== 4. 第一步：运行 Python 转换脚本 ====================
log_info "开始第一步：运行 Python 转换脚本"

# 检查 Python 脚本是否存在
if [ ! -f "$PYTHON_SCRIPT" ]; then
    log_error "Python 脚本不存在: $PYTHON_SCRIPT"
    log_error "请确保脚本在当前目录"
    exit 1
fi

# 运行 Python 脚本
PYTHON_CMD="python3 \"$PYTHON_SCRIPT\" -i \"$INPUT_H5AD\" -o \"$OUTPUT_DIR\" -p \"$FILE_PREFIX\""
log_info "执行: $PYTHON_CMD"
if ! eval "$PYTHON_CMD"; then
    log_error "Python 脚本执行失败"
    exit 1
fi

log_info "✅ Python 脚本执行完成"

# 验证 Python 脚本的输出文件
log_info "验证 Python 脚本输出文件..."
REQUIRED_PYTHON_OUTPUTS=(
    "${FILE_PREFIX}_expression_matrix.mtx"
    "${FILE_PREFIX}_genes.tsv"
    "${FILE_PREFIX}_barcodes.tsv"
    "${FILE_PREFIX}_cell_metadata.tsv"
    "${FILE_PREFIX}_gene_metadata.tsv"
)

for file in "${REQUIRED_PYTHON_OUTPUTS[@]}"; do
    if [ ! -f "$OUTPUT_DIR/$file" ]; then
        log_error "Python 输出文件缺失: $OUTPUT_DIR/$file"
        exit 1
    fi
    log_info "  ✅ $file"
done

# 检查可选文件
OPTIONAL_FILES=("${FILE_PREFIX}_umap_coordinates.tsv" "${FILE_PREFIX}_pca_coordinates.tsv")
for file in "${OPTIONAL_FILES[@]}"; do
    if [ -f "$OUTPUT_DIR/$file" ]; then
        log_info "  📁 可选文件存在: $file"
    fi
done

# ==================== 5. 第二步：运行 R 处理脚本 ====================
log_info "="
log_info "第二步：运行 R 处理脚本"
log_info "="

# 创建临时的 R 执行脚本
TEMP_R_SCRIPT=$(mktemp /tmp/run_r_function_XXXXXX.R)
# 转换verbose参数为R格式
if [ "$VERBOSE" = true ]; then
  R_VERBOSE="TRUE"
else
  R_VERBOSE="FALSE"
fi

cat > "$TEMP_R_SCRIPT" << EOF
#!/usr/bin/env Rscript
# 临时 R 执行脚本
# 加载 R 函数并执行

cat("=== 开始 R 处理步骤 ===\\n")
cat("时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\\n")
cat("文件前缀:", "$FILE_PREFIX", "\\n")
cat("输出目录:", "$OUTPUT_DIR", "\\n")
cat("\\n")

# 加载 R 函数文件
cat("加载 R 函数文件: $R_SCRIPT\\n")
source("$R_SCRIPT")



# 执行处理
cat("执行 read_mtx_datasets 函数...\\n")
start_time <- Sys.time()

result <- tryCatch({
  read_mtx_datasets(
    file_prefix = "$FILE_PREFIX",
    output_dir = "$OUTPUT_DIR",
    verbose = "$R_VERBOSE"
  )
}, error = function(e) {
  cat("❌ 函数执行失败:", e\$message, "\\n")
  quit(save = "no", status = 1)
})

end_time <- Sys.time()
execution_time <- round(as.numeric(difftime(end_time, start_time, units = "secs")), 1)

cat("\\n")
cat("✅ R 处理步骤完成\\n")
cat("执行时间:", execution_time, "秒\\n")

EOF

# 运行 R 脚本
R_CMD="Rscript \"$TEMP_R_SCRIPT\""
log_info "执行 R 处理步骤..."
if ! eval "$R_CMD"; then
    log_error "R 脚本执行失败"
    rm -f "$TEMP_R_SCRIPT"
    exit 1
fi

# 清理临时文件
rm -f "$TEMP_R_SCRIPT"
log_success "R 脚本执行完成"

# ==================== 6. 清理中间文件 ====================
if [ "$REMOVE_INTERMEDIATE" = true ]; then
    log_info "="
    log_info "清理中间文件"
    log_info "="
    
    # 要清理的中间文件列表
    INTERMEDIATE_FILES=(
        "${FILE_PREFIX}_expression_matrix.mtx"
        "${FILE_PREFIX}_genes.tsv"
        "${FILE_PREFIX}_barcodes.tsv"
        "${FILE_PREFIX}_cell_metadata.tsv"
        "${FILE_PREFIX}_gene_metadata.tsv"
        "${FILE_PREFIX}_umap_coordinates.tsv"
        "${FILE_PREFIX}_tsne_coordinates.tsv"
        "${FILE_PREFIX}_pca_coordinates.tsv"
    )
    
    CLEANED_COUNT=0
    for file in "${INTERMEDIATE_FILES[@]}"; do
        if [ -f "$OUTPUT_DIR/$file" ]; then
            rm -f "$OUTPUT_DIR/$file"
            CLEANED_COUNT=$((CLEANED_COUNT + 1))
            log_info "  清理: $file"
        fi
    done
    
    log_info "已清理 $CLEANED_COUNT 个中间文件"
fi
# ==================== 6. 最终验证 ====================
log_info "进行最终输出验证..."

# 检查最终的 RDS 文件
RDS_FILE="$OUTPUT_DIR/${FILE_PREFIX}.RDS"
if [ -f "$RDS_FILE" ]; then
    FILE_SIZE=$(du -h "$RDS_FILE" | cut -f1)
    log_info "✅ 最终输出文件创建成功: $RDS_FILE ($FILE_SIZE)"
    
    # 可选：验证 RDS 文件是否可以正确加载
    if [ "$VERBOSE" = true ]; then
        log_info "验证 RDS 文件完整性..."
        Rscript -e "
          suppressPackageStartupMessages(library(Seurat))
          tryCatch({
            sce <- readRDS('$RDS_FILE')
            cat('  ✅ 文件可正常加载\\n')
            cat('  ✅ 细胞数:', ncol(sce), '\\n')
            cat('  ✅ 基因数:', nrow(sce), '\\n')
          }, error = function(e) {
            cat('  ❌ 文件加载失败:', e\$message, '\\n')
            quit(save = 'no', status = 1)
          })
        " 2>/dev/null
    fi
else
    log_error "最终输出文件未创建: $RDS_FILE"
    exit 1
fi

# ==================== 7. 生成处理报告 ====================
REPORT_FILE="$OUTPUT_DIR/${FILE_PREFIX}_processing_report.txt"
{
    echo "单细胞数据处理报告"
    echo "========================"
    echo "生成时间: $(date)"
    echo "输入文件: $(realpath "$INPUT_H5AD")"
    echo "输出目录: $(realpath "$OUTPUT_DIR")"
    echo "文件前缀: $FILE_PREFIX"
    echo ""
    echo "生成的文件:"
    echo "-----------"
    find "$OUTPUT_DIR" -name "${FILE_PREFIX}*" -type f | while read file; do
        size=$(du -h "$file" | cut -f1)
        echo "  - $(basename "$file") ($size)"
    done
    echo ""
    echo "处理状态: ✅ 成功完成"
} > "$REPORT_FILE"

log_info "处理报告已保存: $REPORT_FILE"

# ==================== 8. 完成 ====================
log_info "="
log_info "✨ 处理流程全部完成！"
log_info "输入: $INPUT_H5AD"
log_info "输出: $OUTPUT_DIR/${FILE_PREFIX}.RDS"

exit 0