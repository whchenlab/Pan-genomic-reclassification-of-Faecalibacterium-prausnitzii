#!/bin/bash

# =============================================================================
# 核心基因拆分与翻译流程
# =============================================================================

set -e  # 遇到错误立即退出

# =============================================================================
# 步骤1: 使用Python脚本拆分核心基因对齐文件
# =============================================================================
# 创建Python脚本并执行
python3 << 'EOF'
import re
import os

# Input files
alignment_file = "core_gene_alignment.aln"
header_file = "core_alignment_header.embl"
output_dir = "split_genes"
os.makedirs(output_dir, exist_ok=True)

# Parse header file to get gene positions and labels
features = []
with open(header_file, 'r') as f:
    current_start = None
    current_end = None
    current_label = None
    for line in f:
        line = line.strip()
        pos_match = re.match(r'FT\s+feature\s+(\d+)\.\.(\d+)', line)
        if pos_match:
            if current_start and current_end and current_label:
                features.append((current_start, current_end, current_label))
            current_start = int(pos_match.group(1))
            current_end = int(pos_match.group(2))
            current_label = None
        label_match = re.search(r'/label=([^;\s]+)', line)
        if label_match and current_start and current_end:
            current_label = label_match.group(1)
    if current_start and current_end and current_label:
        features.append((current_start, current_end, current_label))

# Read alignment file and handle multi-line sequences
with open(alignment_file, 'r') as f:
    lines = f.readlines()
    strains_seqs = []
    i = 0
    while i < len(lines):
        if lines[i].startswith('>'):
            strain = lines[i].strip()
            i += 1
            seq = []
            while i < len(lines) and not lines[i].startswith('>'):
                seq.append(lines[i].strip())
                i += 1
            strains_seqs.append((strain, ''.join(seq)))

# Split sequences by gene positions
for start, end, label in features:
    length = end - start + 1
    output_file = os.path.join(output_dir, f"{label}.aln")
    with open(output_file, 'w') as out:
        for strain, seq in strains_seqs:
            if (start - 1 + length) > len(seq):
                print(f"Warning: {label} sequence length mismatch: expected {length}, got {len(seq)}")
                continue
            gene_seq = seq[start-1:start-1+length]
            out.write(f"{strain}\n{gene_seq}\n")
    print(f"Created {output_file} with {len(strains_seqs)} sequences")

print(f"Total split {len(features)} genes")
EOF

# 检查拆分结果
if [ -d "split_genes" ]; then
    gene_count=$(ls split_genes/*.aln 2>/dev/null | wc -l)
    echo "成功拆分出 ${gene_count} 个基因文件"
else
    echo "错误: 基因拆分失败，split_genes 目录未创建"
    exit 1
fi

# =============================================================================
# 步骤2: 将核酸序列翻译为蛋白质序列
# =============================================================================
mkdir -p split_proteins

# 检查transeq命令是否可用
if ! command -v transeq &> /dev/null; then
    echo "错误: 未找到 transeq 命令，请安装 EMBOSS 工具包"
    echo "安装方法: conda install -c bioconda emboss 或 sudo apt-get install emboss"
    exit 1
fi

# 循环处理所有 .aln 文件
processed_count=0
for aln_file in split_genes/*.aln; do
    if [ -f "$aln_file" ]; then
        filename=$(basename "$aln_file")
        out_file="split_proteins/$filename"
        transeq -sequence "$aln_file" -outseq "$out_file" -osformat fasta &>/dev/null
        
        if [ $? -eq 0 ] && [ -f "$out_file" ]; then
            processed_count=$((processed_count + 1))
            if [ $((processed_count % 50)) -eq 0 ]; then
                echo "已翻译 ${processed_count} 个基因..."
            fi
        else
            echo "警告: 翻译失败: $aln_file"
        fi
    fi
done

# 遍历split_proteins中的文件，删除所有星号
cleaned_count=0
for prot_file in split_proteins/*; do
    if [ -f "$prot_file" ]; then
        cp "$prot_file" "${prot_file}.bak"
        sed -i 's/*//g' "$prot_file"
        cleaned_count=$((cleaned_count + 1))
    fi
done