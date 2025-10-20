#!/bin/bash

# =============================================================================
# 基因分类与eggNOG注释流程
# 1. 从Roary结果中提取core/shell/cloud基因
# 2. 对三类基因分别进行eggNOG注释
# =============================================================================

set -e  # 遇到错误立即退出

# =============================================================================
# 步骤1: 从Roary结果中提取core/shell/cloud基因
# =============================================================================
# 创建Python脚本
cat > extract_genes.py << 'EOF'
import csv
import os

def count_occurrences(row):
    """统计基因在多少个菌株中出现"""
    count = 0
    
    for value in row[14:]:
        if value.strip() != '':
            count += 1
    return count

def find_first_occurrence(row, strains):
    """找到基因第一个出现的菌株及其座位号"""
    
    for i, value in enumerate(row[14:], start=14):
        if value.strip() != '':
            strain_index = i - 14
            return strains[strain_index], value.strip()
    return None, None

def extract_gene_sequence(faa_file, locus_tag):
    """从FAA文件中提取指定座位号的基因序列"""
    if not os.path.exists(faa_file):
        print(f"警告: 文件 {faa_file} 不存在")
        return None
    
    sequence = []
    capture = False
    
    with open(faa_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                
                gene_id = line.split()[0][1:]
                if gene_id == locus_tag:
                    capture = True
                    sequence.append(line.strip())
                else:
                    capture = False
            elif capture:
                sequence.append(line.strip())
    
    if not sequence:
        print(f"警告: 在 {faa_file} 中未找到座位号 {locus_tag} 的基因")
        return None
    
    return '\n'.join(sequence)

def main():
    # 输入文件路径
    csv_file = './roaryresult/gene_presence_absence.csv'
    faa_dir = './faas/'
    
    # 检查输入文件是否存在
    if not os.path.exists(csv_file):
        print(f"错误: 文件 {csv_file} 不存在")
        return
    
    if not os.path.exists(faa_dir):
        print(f"错误: 目录 {faa_dir} 不存在")
        return
    
    # 打开输出文件
    core_file = open('core.faa', 'w')
    shell_file = open('shell.faa', 'w')
    cloud_file = open('cloud.faa', 'w')
    
    try:
        with open(csv_file, 'r') as f:
            reader = csv.reader(f)
            # 读取表头
            header = next(reader)
            
            # 1. 计算总菌株数
            total_strains = len(header) - 14
            print(f"总菌株数: {total_strains}")
            
            # 计算各类型基因的阈值
            core_threshold = total_strains * 0.95
            shell_low_threshold = total_strains * 0.15
            
            print(f"Core基因阈值: ≥ {core_threshold:.1f} 个菌株")
            print(f"Shell基因阈值: {shell_low_threshold:.1f} - {core_threshold:.1f} 个菌株")
            print(f"Cloud基因阈值: < {shell_low_threshold:.1f} 个菌株")
            
            # 获取菌株名称列表
            strains = header[14:]
            
            # 处理每个基因
            gene_count = 0
            for row in reader:
                gene_count += 1
                if gene_count % 100 == 0:
                    print(f"已处理 {gene_count} 个基因...")
                
                # 基因名称在第一列
                gene_name = row[0]
                
                # 统计该基因在多少菌株中出现
                occurrence_count = count_occurrences(row)
                
                # 确定基因类型
                if occurrence_count >= core_threshold:
                    gene_type = 'core'
                    output_file = core_file
                elif occurrence_count >= shell_low_threshold:
                    gene_type = 'shell'
                    output_file = shell_file
                else:
                    gene_type = 'cloud'
                    output_file = cloud_file
                
                # 找到第一个出现该基因的菌株及其座位号
                strain_name, locus_tag = find_first_occurrence(row, strains)
                if not strain_name or not locus_tag:
                    print(f"警告: 基因 {gene_name} 在所有菌株中均未出现，跳过")
                    continue
                
                # 构建FAA文件路径
                faa_file = os.path.join(faa_dir, f"{strain_name}.faa")
                
                # 提取基因序列
                sequence = extract_gene_sequence(faa_file, locus_tag)
                if sequence:
                    # 写入对应的输出文件
                    output_file.write(f"{sequence}\n")
        
        print(f"处理完成，共处理 {gene_count} 个基因")
        print("结果已分别写入 core.faa, shell.faa 和 cloud.faa")
    
    finally:
        # 确保文件被关闭
        core_file.close()
        shell_file.close()
        cloud_file.close()

if __name__ == "__main__":
    main()
EOF

# 运行Python脚本提取基因
python extract_genes.py

# 检查生成的基因文件
echo "检查生成的基因文件:"
for gene_type in core shell cloud; do
    if [ -f "${gene_type}.faa" ]; then
        gene_count=$(grep -c ">" "${gene_type}.faa")
        echo "  ${gene_type}.faa: ${gene_count} 个基因"
    else
        echo "  警告: ${gene_type}.faa 未生成"
    fi
done

# =============================================================================
# 步骤2: 对三类基因分别进行eggNOG注释
# =============================================================================
# eggNOG参数
TAX_SCOPE="2"
DATA_DIR="/mnt/raid6/wuyingjian/biosoft/databases/eggnog-mapper/"
CPU=16

# 检查eggNOG数据库
if [ ! -d "$DATA_DIR" ]; then
    echo "警告: eggNOG数据库目录不存在: $DATA_DIR"
    echo "请检查数据库路径或先下载数据库"
fi

# 对每类基因进行注释
for gene_type in core shell cloud; do
    input_file="${gene_type}.faa"
    output_prefix="${gene_type}_eggnog_annot"
    
    if [ -f "$input_file" ]; then
        echo "正在注释 ${gene_type} 基因..."
        emapper.py -i "$input_file" \
                  --output "$output_prefix" \
                  --tax_scope "$TAX_SCOPE" \
                  --data_dir "$DATA_DIR" \
                  --cpu "$CPU"
        
        if [ $? -eq 0 ]; then
            echo "  ${gene_type} 基因注释完成: ${output_prefix}.emapper.annotations"
        else
            echo "  ${gene_type} 基因注释失败"
        fi
    else
        echo "  跳过 ${gene_type} 基因注释: 输入文件不存在"
    fi
done