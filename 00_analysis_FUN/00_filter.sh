#!/bin/bash

# =============================================================================
# 基因组质量筛选流程
# 1. BUSCO质量评估
# 2. CheckM质量评估  
# 3. FastANI基因组相似性分析
# 4. GTDB-Tk分类学鉴定
# =============================================================================

set -e  # 遇到任何错误立即退出

# 设置目录
INPUT_DIR="genomes"
OUTPUT_DIR="analysis_results"
mkdir -p "$OUTPUT_DIR"

# =============================================================================
# 步骤1: BUSCO质量评估
# =============================================================================
# BUSCO参数
LINEAGE="bacteria_odb10"
MODE="genome"
CPU=4

# 初始化报告文件
echo -e "Status\tBin_ID\tLevel\tComplete\tFragmented\tDuplicated\tReason" > "$OUTPUT_DIR/pass_report_busco.txt"
echo -e "Status\tBin_ID\tLevel\tComplete\tFragmented\tDuplicated\tReason" > "$OUTPUT_DIR/fail_report_busco.txt"

process_busco() {
    local fna_file="$1"
    local level="$2"
    local bin_id=$(basename "$fna_file" .fna)
    
    # 运行BUSCO
    local output_dir="$OUTPUT_DIR/busco_results/${bin_id}_${level}"
    mkdir -p "$(dirname "$output_dir")"
    
    echo "运行BUSCO: $bin_id"
    busco -i "$fna_file" -l "$LINEAGE" -o "$output_dir" -m "$MODE" --cpu "$CPU" &>/dev/null
    
    # 检查BUSCO结果
    local summary_file="${output_dir}/run_${LINEAGE}/short_summary.txt"
    if [ ! -f "$summary_file" ]; then
        echo -e "Fail\t$bin_id\t$level\t-\t-\t-\tBUSCO运行失败" >> "$OUTPUT_DIR/fail_report_busco.txt"
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
        echo -e "Pass\t$bin_id\t$level\t$complete\t$fragmented\t$duplicated\t-" >> "$OUTPUT_DIR/pass_report_busco.txt"
    else
        echo -e "Fail\t$bin_id\t$level\t$complete\t$fragmented\t$duplicated\t${reason%%; }" >> "$OUTPUT_DIR/fail_report_busco.txt"
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
# =============================================================================
# 步骤2: CheckM质量评估
# =============================================================================
# checkm参数
level="genomes"
comp_thresh=95                 # 完整性阈值
contamination_threshold=5      # 污染阈值
strain_threshold=10            # 菌株异质性阈值
threads=8

# 递归深度问题
export PYTHONRECURSIONLIMIT=10000

# 初始化报告文件
report_pass="$OUTPUT_DIR/pass_report_checkm.txt"
report_fail="$OUTPUT_DIR/fail_report_checkm.txt"
for report in "$report_pass" "$report_fail"; do
    if [ ! -f "$report" ]; then
        echo -e "Status\tBin_ID\tLevel\tCompleteness\tContamination\tStrain_heterogeneity\tReason" > "$report"
    fi
done

# 核心处理函数
process_bin() {
    local bin_id="$1"
    local level="$2"
    local checkm_output_dir="$OUTPUT_DIR/checkm_results/$level"
    
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
# 0. 检查输入目录
if [ ! -d "$level" ]; then
    echo "错误：输入目录 $level 不存在！"
    exit 1
fi

# 1. 运行CheckM lineage_wf
checkm_output_dir="$OUTPUT_DIR/checkm_results/$level"
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

# 3. 实时处理每个bin
echo "开始解析质量报告..."
bin_ids=($(tail -n +2 "$quality_report" | cut -f1))
total=${#bin_ids[@]}

for ((i=0; i<$total; i++)); do
    bin_id="${bin_ids[$i]}"
    printf "处理进度: %d/%d (%.1f%%)\r" $((i+1)) $total $((100*(i+1)/total))
    process_bin "$bin_id" "$level"
done
echo

# 4. 统计结果
pass_count=$(wc -l < "$report_pass")
fail_count=$(wc -l < "$report_fail")
echo "CheckM分析完成!"
echo "合格基因组数: $((pass_count-1))"
echo "不合格基因组数: $((fail_count-1))"

# =============================================================================
# 步骤3: FastANI基因组相似性分析
# =============================================================================
# 创建基因组列表文件
genomes_list="$OUTPUT_DIR/genomes.txt"
find "$INPUT_DIR" -name "*.fna" -o -name "*.fasta" > "$genomes_list"

genome_count=$(wc -l < "$genomes_list")
if [ "$genome_count" -eq 0 ]; then
    echo "错误：在 $INPUT_DIR 中未找到任何基因组文件"
    exit 1
fi

echo "找到 $genome_count 个基因组文件"

# 运行FastANI
echo "运行FastANI..."
fastANI --ql "$genomes_list" --rl "$genomes_list" -o "$OUTPUT_DIR/ani_matrix.txt" -t 10

if [ ! -f "$OUTPUT_DIR/ani_matrix.txt" ]; then
    echo "FastANI运行失败，未生成输出文件"
    exit 1
fi
# =============================================================================
# 步骤4: GTDB-Tk分类学鉴定
# =============================================================================
gtdbtk_output="$OUTPUT_DIR/gtdbtk_results"
mkdir -p "$gtdbtk_output"

# 检查输入文件
file_count=$(find "$INPUT_DIR" -name "*.fna" -o -name "*.fasta" | wc -l)
if [ "$file_count" -eq 0 ]; then
    echo "错误：在 $INPUT_DIR 中未找到任何基因组文件"
    exit 1
fi

# 运行GTDB-Tk
echo "运行GTDB-Tk classify_wf..."
gtdbtk classify_wf \
    --genome_dir "$INPUT_DIR" \
    --extension '.fna' \
    --out_dir "$gtdbtk_output" \
    --cpus 4 \
    --min_perc_aa 10 \
    --min_af 0.65 \
    --skip_ani_screen