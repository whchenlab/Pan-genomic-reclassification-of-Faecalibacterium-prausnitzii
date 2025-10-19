#!/bin/bash

# BUSCO参数
LINEAGE="bacteria_odb10"
MODE="genome"
CPU=4

# 初始化报告文件
echo -e "Status\tBin_ID\tLevel\tComplete\tFragmented\tDuplicated\tReason" > pass_report_busco.txt
echo -e "Status\tBin_ID\tLevel\tComplete\tFragmented\tDuplicated\tReason" > fail_report_busco.txt

process_busco() {
    local fna_file="$1"
    local level="$2"
    local bin_id=$(basename "$fna_file" .fna)
    
    # 运行BUSCO
    local output_dir="busco_results/${bin_id}_${level}"
    busco -i "$fna_file" -l "$LINEAGE" -o "$output_dir" -m "$MODE" --cpu "$CPU" &>/dev/null
    
    # 检查BUSCO结果
    local summary_file="${output_dir}/run_${LINEAGE}/short_summary.txt"
    if [ ! -f "$summary_file" ]; then
        echo -e "Fail\t$bin_id\t$level\t-\t-\t-\tBUSCO运行失败" >> fail_report_busco.txt
        return 1
    fi

    # 提取指标
    local line=$(grep "C:" "$summary_file")
    local complete=$(echo "$line" | sed -E 's/.*C:([0-9]+\.?[0-9]*)%.*/\1/')
    local fragmented=$(echo "$line" | sed -E 's/.*F:([0-9]+\.?[0-9]*)%.*/\1/')
    local duplicated=$(echo "$line" | sed -E 's/.*D:([0-9]+\.?[0-9]*)%.*/\1/')

    # 构建原因
    local reason=""
    (( $(echo "$complete < 90" | bc -l) )) && reason+="Complete <90%; "
    (( $(echo "$fragmented > 10" | bc -l) )) && reason+="Fragmented >10%; "
    (( $(echo "$duplicated > 5" | bc -l) )) && reason+="Duplicated >5%; "

    # 写入报告
    if [ -z "$reason" ]; then
        echo -e "Pass\t$bin_id\t$level\t$complete\t$fragmented\t$duplicated\t-" >> pass_report_busco.txt
    else
        echo -e "Fail\t$bin_id\t$level\t$complete\t$fragmented\t$duplicated\t${reason%%; }" >> fail_report_busco.txt
    fi
}

# 主循环
for level in genomes; do
    if [ -d "./$level" ]; then
        for fna in "./$level"/*.fna; do
            [ -f "$fna" ] && process_busco "$fna" "$level"
        done
    fi
done

echo "分析完成，结果已写入报告文件。"