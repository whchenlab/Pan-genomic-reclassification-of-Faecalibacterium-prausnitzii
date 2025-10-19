#!/bin/bash

# checkm参数
level="genomes"
comp_thresh=95                 # 完整性阈值
contamination_threshold=5      # 污染阈值
strain_threshold=10            # 菌株异质性阈值
threads=8                      # 线程数

# 递归深度问题
export PYTHONRECURSIONLIMIT=10000

# 初始化报告文件
report_pass="pass_report_checkm.txt"
report_fail="fail_report_checkm.txt"
for report in "$report_pass" "$report_fail"; do
    if [ ! -f "$report" ]; then
        echo -e "Status\tBin_ID\tLevel\tCompleteness\tContamination\tStrain_heterogeneity\tReason" > "$report"
    fi
done

# 核心处理函数
process_bin() {
    local bin_id="$1"
    local level="$2"
    local checkm_output_dir="checkm_results/$level"
    
    # 实时检查每个bin的CheckM结果
    if grep -q "^${bin_id}" "$checkm_output_dir/quality_report.tsv"; then
        # 从现有报告提取数据
        local line=$(grep "^${bin_id}" "$checkm_output_dir/quality_report.tsv")
        IFS=$'\t' read -ra cols <<< "$line"
        
        completeness="${cols[5]}"
        contamination="${cols[6]}"
        strain_heterogeneity="${cols[7]}"
    else
        # 如果未找到则标记为未处理
        echo -e "FAIL\t${bin_id}\t${level}\t0.00\t0.00\t0.00\t未找到CheckM结果" >> "$report_fail"
        return 1
    fi

    # 质量控制
    reason=""
    status="PASS"
    
    (( $(echo "$completeness < $comp_thresh" | bc -l) )) && { reason+="完整性不足; "; status="FAIL"; }
    (( $(echo "$contamination > $contamination_threshold" | bc -l) )) && { reason+="污染过高; "; status="FAIL"; }
    (( $(echo "$strain_heterogeneity > $strain_threshold" | bc -l) )) && { reason+="菌株异质性过高; "; status="FAIL"; }

    # 实时写入报告
    record=$(printf "%s\t%s\t%.2f\t%.2f\t%.2f" "$bin_id" "$level" "$completeness" "$contamination" "$strain_heterogeneity")
    if [ "$status" == "PASS" ]; then
        echo -e "PASS\t$record\t-" >> "$report_pass"
    else
        echo -e "FAIL\t$record\t${reason%; }" >> "$report_fail"
    fi
}

# 主流程
main() {
    # 0. 检查输入目录
    if [ ! -d "$level" ]; then
        echo "错误：输入目录 $level 不存在！"
        exit 1
    fi

    # 1. 运行CheckM lineage_wf
    checkm_output_dir="checkm_results/$level"
    mkdir -p "$checkm_output_dir"
    
    if [ ! -f "$checkm_output_dir/lineage.ms" ]; then
        echo "启动CheckM lineage_wf..."
        if ! checkm lineage_wf -t $threads -x fna "$level" "$checkm_output_dir"; then
            echo "CheckM lineage_wf运行失败！请检查错误日志。"
            exit 1
        fi
    else
        echo "检测到已存在lineage.ms，跳过CheckM lineage_wf"
    fi

    # 2. 生成质量报告
    quality_report="$checkm_output_dir/quality_report.tsv"
    if [ ! -f "$quality_report" ]; then
        echo "生成CheckM质量报告..."
        if ! checkm qa "$checkm_output_dir/lineage.ms" "$checkm_output_dir" \
            -o 2 --tab_table -f "$quality_report"; then
            echo "CheckM qa运行失败！请检查错误日志。"
            exit 1
        fi
    else
        echo "检测到已存在质量报告，直接解析"
    fi

    # 3. 实时处理每个bin（并行加速）
    echo "开始解析质量报告..."
    bin_ids=($(tail -n +2 "$quality_report" | cut -f1))
    total=${#bin_ids[@]}
    
    for ((i=0; i<$total; i++)); do
        bin_id="${bin_ids[$i]}"
        printf "处理进度: %d/%d (%.1f%%)\r" $((i+1)) $total $((100*(i+1)/total))
        process_bin "$bin_id" "$level"
    done
    echo  # 换行

    # 4. 统计结果
    pass_count=$(wc -l < "$report_pass")
    fail_count=$(wc -l < "$report_fail")
    echo "===== 最终统计 ====="
    echo "合格基因组数: $((pass_count-1))"  # 减去表头
    echo "不合格基因组数: $((fail_count-1))"
}

# 执行主流程
main