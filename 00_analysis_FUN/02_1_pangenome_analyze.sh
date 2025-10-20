#!/bin/bash

# 第一部分：使用Prokka进行注释
mkdir -p gffs

# 遍历genomes目录下所有.fna文件
for fna_file in genomes/*.fna; do
    # 从文件名提取前缀（移除路径和扩展名）
    prefix=$(basename "$fna_file" .fna)
    
    # 为每个文件创建专属输出目录
    out_dir="prokka_output_${prefix}"
    
    # 运行prokka注释
    prokka --cpus 8 --kingdom Bacteria \
        --outdir "$out_dir" \
        --prefix "$prefix" \
        "$fna_file"
    
    # 将生成的.gff文件移动到gffs目录
    if [ -f "${out_dir}/${prefix}.gff" ]; then
        mv "${out_dir}/${prefix}.gff" gffs/
    else
        echo "Warning: .gff not found for ${prefix}"
    fi
done

echo "All annotations completed. GFF files are in gffs/ directory."

# 第二部分：使用Roary进行泛基因组分析
roary -f roaryresult -p 16 -e -n -v --mafft gffs/*.gff

# 第三部分：使用FastTree构建系统发育树
FastTree -nt -gtr roaryresult/core_gene_alignment.aln > roaryresult/mytree.newick

echo "Pipeline finished. Phylogenetic tree is in roaryresult/mytree.newick"