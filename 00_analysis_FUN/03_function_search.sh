#!/bin/bash

# 检查BLAST是否安装
if ! command -v blastp &> /dev/null; then
    echo "错误: BLAST未安装。请先安装BLAST+工具包。"
    exit 1
fi

# 创建结果目录
mkdir -p blast_results_fp blast_results_fd blast_results_fl reports logs

# 获取所有目标蛋白质文件
target_proteins=(*.faa)
if [ ${#target_proteins[@]} -eq 0 ]; then
    echo "错误: 当前目录下未找到目标蛋白质文件(*.faa)。"
    exit 1
fi

# BLAST参数
EVALUE="1e-10"        # evalue阈值
IDENTITY="50"        # 最小序列一致性(%)
COVERAGE="60"        # 最小覆盖度(%)
THREADS=8

# 初始化报告文件
report_file="reports/blast_summary_report_$(date +%Y%m%d).txt"
echo "蛋白质BLAST分析报告 - $(date)" > "$report_file"
echo "=========================" >> "$report_file"
echo "参数设置: evalue=$EVALUE, 最小一致性=$IDENTITY%, 最小覆盖度=$COVERAGE%" >> "$report_file"

# 定义待处理文件夹及其对应的结果文件夹
folders=(
    "proteins_fp:blast_results_fp"
    "proteins_fd:blast_results_fd"
    "proteins_fl:blast_results_fl"
)

# 对每个目标蛋白质进行BLAST分析
for target in "${target_proteins[@]}"; do
    target_name=$(basename "$target" .ffn)
    echo -e "\n分析目标蛋白质: $target_name"
    echo -e "\n分析目标蛋白质: $target_name" >> "$report_file"
    
    # 遍历每个文件夹
    for folder_pair in "${folders[@]}"; do
        IFS=':' read -r protein_folder result_folder <<< "$folder_pair"
        
        # 获取该文件夹中的所有细菌蛋白质组文件
        bacteria_proteins=("$protein_folder"/*.faa)
        if [ ${#bacteria_proteins[@]} -eq 0 ]; then
            echo "警告: $protein_folder目录下未找到细菌蛋白质组文件(*.faa)。"
            continue
        fi
        
        echo -e "\n  处理文件夹: $protein_folder"
        echo -e "\n  处理文件夹: $protein_folder" >> "$report_file"
        
        found_count=0
        not_found_count=0
        
        for protein in "${bacteria_proteins[@]}"; do
            protein_name=$(basename "$protein" .faa)
            log_file="logs/${target_name}_vs_${protein_name}_${protein_folder//\//_}.log"
            echo -n "    正在比对 $protein_name..."
            
            # 创建蛋白质数据库
            makeblastdb -in "$protein" -dbtype prot -out "temp_db" > "$log_file" 2>&1
            
            # 运行BLASTP
            blastp -query "$target" -db "temp_db" \
                   -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \
                   -evalue "$EVALUE" \
                   -out "temp_blast.txt" \
                   -num_threads "$THREADS" > /dev/null 2>&1
            
            if [ -s "temp_blast.txt" ]; then
                # 提取查询序列长度
                qlen=$(awk 'NR==1 {print $13}' temp_blast.txt)
                
                # 处理qlen为0或空的情况
                if [[ -z "$qlen" || "$qlen" -eq 0 ]]; then
                    echo " 无效(查询长度为0)"
                    ((not_found_count++))
                    cp "temp_blast.txt" "${result_folder}/${target_name}_vs_${protein_name}_invalid.txt"
                    continue
                fi
                
                # 计算覆盖度
                awk -v qlen="$qlen" -v cov="$COVERAGE" -v id="$IDENTITY" '
                    BEGIN {max_end=0; min_start=1000000000}
                    $3 >= id {
                        if ($7 < min_start) min_start=$7  # qstart
                        if ($8 > max_end) max_end=$8      # qend
                    }
                    END {
                        if (min_start <= max_end) {
                            coverage = (max_end - min_start + 1) * 100 / qlen
                            if (coverage >= cov) print coverage
                        }
                    }' temp_blast.txt > coverage.txt
                
                if [ -s coverage.txt ]; then
                    echo " 存在"
                    ((found_count++))
                    cp "temp_blast.txt" "${result_folder}/${target_name}_vs_${protein_name}.txt"
                    echo "覆盖度: $(cat coverage.txt)" >> "${result_folder}/${target_name}_vs_${protein_name}.txt"
                else
                    echo " 不存在(低覆盖度)"
                    ((not_found_count++))
                    cp "temp_blast.txt" "${result_folder}/${target_name}_vs_${protein_name}_low.txt"
                fi
            else
                echo " 不存在(无显著比对)"
                ((not_found_count++))
            fi
            
            # 保存日志
            echo "比对 $target_name vs $protein_name 完成" >> "$log_file"
            echo "结果: $(if [ -s "temp_blast.txt" ]; then echo "有比对"; else echo "无比对"; fi)" >> "$log_file"
            
            # 清理临时文件
            rm -f temp_db.* temp_blast.txt coverage.txt
        done
        
        echo "    总结: $target_name 在 $found_count 株菌中存在，在 $not_found_count 株菌中不存在"
        echo "    总结: $target_name 在 $found_count 株菌中存在，在 $not_found_count 株菌中不存在" >> "$report_file"
    done
done

# 生成总体总结
echo -e "\n\n总体总结" >> "$report_file"
echo "=========================" >> "$report_file"
echo "分析日期: $(date)" >> "$report_file"
echo "目标基因总数: ${#target_proteins[@]}" >> "$report_file"
echo "分析的细菌菌株总数:" >> "$report_file"
for folder_pair in "${folders[@]}"; do
    IFS=':' read -r protein_folder _ <<< "$folder_pair"
    bacteria_proteins=("$protein_folder"/*.ffn)
    echo "  - $protein_folder: ${#bacteria_proteins[@]}" >> "$report_file"
done

echo -e "\nBLAST分析完成！"
echo "详细比对结果保存在以下目录下:"
for folder_pair in "${folders[@]}"; do
    IFS=':' read -r _ result_folder <<< "$folder_pair"
    echo "  - $result_folder"
done
echo "总结报告保存在: $report_file"