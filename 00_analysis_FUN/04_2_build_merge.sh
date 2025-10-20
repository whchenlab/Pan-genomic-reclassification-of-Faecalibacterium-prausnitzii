#!/bin/bash

# =============================================================================
# HMM分类器构建与合并流程
# =============================================================================

set -e  # 遇到错误立即退出

# =============================================================================
# 步骤1: HMM分类器构建 - 批处理多个基因比对文件
# =============================================================================
echo "构建多个分类单元的HMM模型..."

# 设置目录
INPUT_DIR="split/split_proteins"
OUTPUT_DIR="hmm_classifier"
mkdir -p $OUTPUT_DIR

# 分类单元定义（根据实际情况修改）
FP_STRAINS="F.prausnitzii-ATCC-27768-a F.prausnitzii-ATCC-27768-b F.tardum-CLA-AA-H175 F.intestinale_CLA-AA-H281"
FD_STRAINS="F.duncaniae-A2-165-a F.duncaniae-A2-165-b F.duncaniae-A2-165-c F.duncaniae-JCM-31915"
FL_STRAINS="GCA_000210735 GCA_003287485 GCA_003287505 GCA_003293635"
FO_STRAINS="F.faecis-CLA-JM-H7-B F.taiwanense-HLW78 F.gallinarum_JCM-17207 F.wellingii_HTF-F_a F.wellingii_HTF-F_b F.hattorii_APC922"     

# 创建输出文件夹
mkdir -p "${OUTPUT_DIR}/Fp" "${OUTPUT_DIR}/Fd" "${OUTPUT_DIR}/Fl" "${OUTPUT_DIR}/Other"

# 按分类单元分割单个比对文件
split_alignment() {
    local gene_name=$1
    local group=$2
    local strains=$3
    
    # 创建临时文件
    TMP_FILE=$(mktemp)
    
    # 输入比对文件
    INPUT_ALIGNMENT="${INPUT_DIR}/${gene_name}.aln"
    
    # 检查输入文件是否存在
    if [ ! -f "$INPUT_ALIGNMENT" ]; then
        echo "警告: 比对文件不存在: $INPUT_ALIGNMENT"
        rm "$TMP_FILE"
        return 1
    fi
    
    # 提取属于该类群的序列
    for strain in $strains; do
        grep -A1 "$strain" "$INPUT_ALIGNMENT" | grep -v "^--$" >> "$TMP_FILE"
    done
    
    # 检查是否提取到序列
    if [ ! -s "$TMP_FILE" ]; then
        echo "警告: 未在 $gene_name 中找到 $group 类群的序列"
        rm "$TMP_FILE"
        return 1
    fi
    
    # 保存为单独的比对文件
    OUTPUT_ALIGNMENT="${OUTPUT_DIR}/${group}/${gene_name}_${group}_alignment.aln"
    cp "$TMP_FILE" "$OUTPUT_ALIGNMENT"
    echo "已生成 ${gene_name}_${group} 的比对文件: ${OUTPUT_ALIGNMENT}"
    
    # 清理临时文件
    rm "$TMP_FILE"
    
    # 构建HMM模型
    HMM_FILE="${OUTPUT_DIR}/${group}/${gene_name}_${group}.hmm"
    hmmbuild "$HMM_FILE" "$OUTPUT_ALIGNMENT"
    echo "已生成 ${gene_name}_${group} 的HMM模型: ${HMM_FILE}"
}

# 获取所有基因文件列表
GENE_FILES=$(find "$INPUT_DIR" -name "*.aln" -type f | sed 's#.*/\(.*\)\.aln#\1#')

# 检查是否找到基因文件
if [ -z "$GENE_FILES" ]; then
    echo "错误: 在 $INPUT_DIR 目录下未找到.aln文件！"
    exit 1
fi

# 对每个基因文件执行处理
echo "找到 $(echo "$GENE_FILES" | wc -w) 个基因文件，开始处理..."
for gene in $GENE_FILES; do
    echo -e "\n处理基因: $gene"
    
    # 对每个分类单元执行分割和HMM构建
    split_alignment "$gene" "Fp" "$FP_STRAINS"
    split_alignment "$gene" "Fd" "$FD_STRAINS"
    split_alignment "$gene" "Fl" "$FL_STRAINS"
    split_alignment "$gene" "Other" "$FO_STRAINS"
done

# 为每个分类单元合并并压缩HMM模型
echo -e "\n正在为每个分类单元合并HMM模型..."
for group in Fp Fd Fl Other; do
    GROUP_DIR="${OUTPUT_DIR}/${group}"
    MERGED_HMM="${OUTPUT_DIR}/${group}_models.hmm"
    
    # 检查是否有HMM文件可合并
    if [ $(find "$GROUP_DIR" -name "*.hmm" | wc -l) -eq 0 ]; then
        echo "警告: $group 类群没有HMM模型文件可合并"
        continue
    fi
    
    # 合并HMM模型
    cat "${GROUP_DIR}"/*.hmm > "$MERGED_HMM"
    # 压缩HMM模型
    hmmpress "$MERGED_HMM"
    echo "已生成并压缩 ${group} 类群的合并HMM模型: ${MERGED_HMM}"
done

echo "各分类单元的HMM模型构建完成！"

# =============================================================================
# 步骤2: 合并所有HMM模型为统一分类器
# =============================================================================
echo ""
echo "合并所有HMM模型为统一分类器..."

# 设置目录和文件
HMM_DIR="hmm_classifier"
MERGED_MODEL="${HMM_DIR}/merged.hmm"

# 检查HMM目录是否存在
if [ ! -d "$HMM_DIR" ]; then
    echo "错误: HMM模型目录 '$HMM_DIR' 不存在！"
    exit 1
fi

# 查找并合并HMM模型文件
echo "正在查找HMM模型文件..."
HMM_FILES=$(find "$HMM_DIR" -maxdepth 1 -name "*_models.hmm" -type f)

# 检查是否找到模型文件
if [ -z "$HMM_FILES" ]; then
    echo "错误: 在 '$HMM_DIR' 目录下未找到任何HMM模型文件！"
    exit 1
fi

# 合并HMM模型文件
echo "正在合并HMM模型文件..."
cat $HMM_FILES > "$MERGED_MODEL"

# 检查合并是否成功
if [ ! -f "$MERGED_MODEL" ]; then
    echo "错误: 合并模型文件 '$MERGED_MODEL' 失败！"
    exit 1
fi

# 压缩合并后的模型文件
echo "正在压缩合并后的模型文件..."
hmmpress "$MERGED_MODEL"

# 验证压缩是否成功
if [ ! -f "${MERGED_MODEL}.h3f" ]; then
    echo "错误: 压缩模型文件失败！"
    exit 1
fi

echo "HMM模型合并和压缩完成！"