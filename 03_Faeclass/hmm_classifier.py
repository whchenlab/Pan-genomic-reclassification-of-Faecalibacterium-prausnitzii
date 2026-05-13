#!/usr/bin/env python3
"""
HMM分类器运行脚本
用于对 test 文件夹中的所有 .faa 文件执行 HMM 搜索并保存结果，然后汇总分类结果。
"""

import os
import sys
import subprocess
import argparse
import glob
from collections import defaultdict

# ==================== 参数解析 ====================
def parse_args():
    parser = argparse.ArgumentParser(description="运行 HMM 搜索并汇总分类结果")
    parser.add_argument("TEST_DIR", help="包含 .faa 文件的测试序列目录（必填）")
    parser.add_argument("--HMM_DIR", default="hmm_classifier",
                        help="HMM 模型所在目录（默认：hmm_classifier）")
    parser.add_argument("--OUTPUT_DIR", default="hmm_results",
                        help="结果输出目录（默认：hmm_results）")
    return parser.parse_args()

# ==================== 定义物种映射 ====================
SPECIES_PATTERNS = {
    "_Fp_alignment$": "F.prausnitzii",
    "_Fd_alignment$": "F.duncaniae",
    "_Fl_alignment$": "F.longum",
    "_Ff_alignment$": "F.faecis",
    "_Fw_alignment$": "F.wellingii",
    "_Fb_alignment$": "F.butyricigenerans",
    "_Fh_alignment$": "F.hattorii",
    "_Fla_alignment$": "F.langellae",
    "_Fun1_alignment$": "F.centralis",
    "_Fun2_alignment$": "F.indigenum",
    "_Fun3_alignment$": "F.bellum",
    "_Fun4_alignment$": "F.protidificans",
}

VALID_SPECIES = list(SPECIES_PATTERNS.values())

# ==================== 辅助函数 ====================
def check_hmm_press(model_path):
    """检查 HMM 模型是否已压缩（hmmpress）"""
    for ext in [".h3f", ".h3i", ".h3m", ".h3p"]:
        if not os.path.exists(model_path + ext):
            return False
    return True

def run_hmmsearch(model_file, query_faa, output_tbl, output_dom):
    """运行 hmmsearch，返回是否成功"""
    cmd = [
        "hmmsearch",
        "--tblout", output_tbl,
        "--domtblout", output_dom,
        model_file,
        query_faa
    ]
    try:
        # 重定向标准输出到 /dev/null（静默模式）
        with open(os.devnull, 'w') as devnull:
            subprocess.run(cmd, check=True, stdout=devnull, stderr=subprocess.PIPE)
        return True
    except subprocess.CalledProcessError as e:
        print(f"错误: hmmsearch 执行失败，返回码 {e.returncode}", file=sys.stderr)
        return False

def parse_hmm_results(tbl_file):
    """
    解析单个 hmmsearch 的 --tblout 结果文件
    返回一个字典：{物种名: (total_score, count)}
    规则：每个 query 只保留最高分，然后按物种汇总
    """
    # 存储每个 query 的最高分及对应的物种
    query_max_score = {}
    query_species = {}

    with open(tbl_file, 'r') as f:
        for line in f:
            if line.startswith('#') or line.strip() == '':
                continue
            parts = line.split()
            if len(parts) < 7:
                continue
            # 第三列: query name, 第六列: full sequence score
            query = parts[2]
            try:
                score = float(parts[5])   # 注意：索引从0开始，第五列是score
            except ValueError:
                continue
            if score <= 0:
                continue

            # 根据 query 名称匹配物种
            species = None
            for pattern, sp in SPECIES_PATTERNS.items():
                import re
                if re.search(pattern, query):
                    species = sp
                    break
            if species is None:
                continue

            # 保留每个 query 的最高分
            if query not in query_max_score or score > query_max_score[query]:
                query_max_score[query] = score
                query_species[query] = species

    # 按物种汇总
    species_total = defaultdict(float)
    species_count = defaultdict(int)
    for query, score in query_max_score.items():
        sp = query_species[query]
        species_total[sp] += score
        species_count[sp] += 1

    return species_total, species_count

def determine_species(species_total):
    """
    根据各物种的总分，选择总分最高的物种，并根据平均分阈值决定最终分类
    返回 (max_species, max_total, mean_score)
    """
    if not species_total:
        return "Outgroup", 0.0, 0.0

    # 找出总分最高的物种
    max_species = max(species_total, key=lambda sp: species_total[sp])
    max_total = species_total[max_species]
    # 注意：需要知道该物种的 count，但由于 species_total 只有总分，我们需要从原始汇总数据获取 count
    # 这个函数会在外部调用时传入 count 字典，所以这里只做占位，实际在调用处计算平均分
    # 为了方便，改为传入 (total, count) 的字典
    return max_species, max_total

def process_summary(details_dir, output_dir):
    """处理 details 目录下的所有 *_hmm_results.txt 文件，生成 summary.txt"""
    summary_file = os.path.join(output_dir, "summary.txt")
    with open(summary_file, 'w') as out_f:
        out_f.write("Genome\tMScore\tMean_MScore\tSpecies\n")

        # 查找所有结果文件
        pattern = os.path.join(details_dir, "*_hmm_results.txt")
        result_files = glob.glob(pattern)
        if not result_files:
            print(f"警告: 在 {details_dir} 中没有找到任何 *_hmm_results.txt 文件", file=sys.stderr)
            return

        for tbl_file in result_files:
            # 提取基因组名称（与 bash 脚本逻辑一致）
            basename = os.path.basename(tbl_file)
            # 去掉 _hmm_results.txt 后缀
            genome_name = basename.replace("_hmm_results.txt", "")

            # 解析结果
            species_total, species_count = parse_hmm_results(tbl_file)

            if not species_total:
                # 没有任何有效匹配
                out_f.write(f"{genome_name}\t0.00\t0.00\tOutgroup\n")
                continue

            # 找出总分最高的物种
            max_species = max(species_total, key=lambda sp: species_total[sp])
            max_total = species_total[max_species]
            mean_score = max_total / species_count[max_species] if species_count[max_species] > 0 else 0.0

            # 根据阈值决定最终分类
            if mean_score < 120:
                final_species = "Outgroup"
            elif mean_score < 139:
                final_species = "F.else"
            else:
                final_species = max_species

            out_f.write(f"{genome_name}\t{max_total:.2f}\t{mean_score:.2f}\t{final_species}\n")
            print(f"处理完成: {tbl_file} -> 添加到汇总")

    print(f"最终汇总文件已生成：{summary_file}")
    # 统计记录数（除去表头）
    with open(summary_file, 'r') as f:
        line_count = sum(1 for _ in f) - 1
    print(f"总记录数: {line_count} 条")

# ==================== 主函数 ====================
def main():
    args = parse_args()
    test_dir = args.TEST_DIR
    hmm_dir = args.HMM_DIR
    output_dir = args.OUTPUT_DIR
    details_dir = os.path.join(output_dir, "details")

    # 1. 检查输入目录
    if not os.path.isdir(test_dir):
        print(f"错误: 测试序列目录 '{test_dir}' 不存在！", file=sys.stderr)
        sys.exit(1)

    # 2. 检查合并的 HMM 模型
    merged_model = os.path.join(hmm_dir, "merged.hmm")
    if not os.path.isfile(merged_model):
        print(f"错误: HMM 模型文件 '{merged_model}' 不存在！", file=sys.stderr)
        sys.exit(1)

    if not check_hmm_press(merged_model):
        print(f"错误: HMM 模型文件 '{merged_model}' 未压缩！请先运行 hmmpress。", file=sys.stderr)
        sys.exit(1)

    # 3. 创建输出目录
    os.makedirs(details_dir, exist_ok=True)

    # 4. 查找所有 .faa 文件
    faa_files = glob.glob(os.path.join(test_dir, "*.faa"))
    if not faa_files:
        print(f"错误: 在 '{test_dir}' 目录下未找到任何 .faa 文件！", file=sys.stderr)
        sys.exit(1)

    total_files = len(faa_files)
    print(f"找到 {total_files} 个 .faa 文件，开始处理...")

    # 5. 对每个文件运行 hmmsearch
    for idx, faa_file in enumerate(faa_files, start=1):
        print(f"\n[{idx}/{total_files}] 处理文件: {faa_file}")
        base_name = os.path.splitext(os.path.basename(faa_file))[0]
        results_file = os.path.join(details_dir, f"{base_name}_hmm_results.txt")
        dom_results_file = os.path.join(details_dir, f"{base_name}_hmm_results_dom.txt")

        print("  正在运行 HMM 搜索...")
        success = run_hmmsearch(merged_model, faa_file, results_file, dom_results_file)
        if not success:
            print(f"  错误: 对 {faa_file} 的 HMM 搜索失败！", file=sys.stderr)
            continue

        print(f"  结果已保存至: {results_file}")
        print(f"  结构域详情已保存至: {dom_results_file}")

    print("\nHMM 搜索完成！")
    print(f"结果保存在: {output_dir}/")
    print(f"详细结果保存在: {details_dir}/")

    # 6. 生成汇总文件
    process_summary(details_dir, output_dir)

if __name__ == "__main__":
    main()