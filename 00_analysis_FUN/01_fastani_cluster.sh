#!/bin/bash

# =============================================================================
# ANI聚类分析流程
# =============================================================================

set -e  # 遇到错误立即退出

# =============================================================================
# 函数：执行单个文件夹的ANI分析和聚类
# =============================================================================
run_ani_clustering() {
    local folder="$1"
    local folder_name="$2"
    
    echo "处理文件夹: $folder (标识: $folder_name)"
    cd "$folder"
    
    # 步骤1: fastANI分析
    ls -1 > genomes.txt  
    fastANI --ql genomes.txt --rl genomes.txt -o ani_matrix.txt -t 10    
    # 步骤2: 聚类分析
    
    # 配置参数
    ANI_MATRIX="ani_matrix.txt"
    FASTA_DIR="."
    OUTPUT_FILE="../${folder_name}_cluster_results.txt"
    CLUSTER_DIR="../${folder_name}_clusters"
    THREADS=$(nproc)
    
    # 预处理ANI矩阵
    filtered_pairs=$(mktemp)
    awk -v threshold=95 '
        NF >=3 && $3 >= threshold {print $1,$2}
    ' "$ANI_MATRIX" > "$filtered_pairs"
    
    # 构建菌株哈希表
    declare -A strain_map
    while IFS= read -r -d $'\0' file; do
        strain=$(basename "$file")
        strain_map["$strain"]=1
    done < <(find "$FASTA_DIR" -maxdepth 1 -name "*.fna" -print0)
    
    # 并查集操作
    clusters_raw=$(mktemp)
    awk -v strains="${!strain_map[*]}" '
    BEGIN {
        split(strains, arr, " ")
        for (i in arr) exists[arr[i]] = 1
        
        for (s in exists) {
            parent[s] = s
            rank[s] = 1
        }
    }
    
    function find(s) {
        while (parent[s] != s) {
            parent[s] = parent[parent[s]]
            s = parent[s]
        }
        return s
    }
    
    function union(s1, s2) {
        root1 = find(s1)
        root2 = find(s2)
        if (root1 == root2) return
        
        if (rank[root1] < rank[root2]) {
            parent[root1] = root2
            rank[root2] += rank[root1]
        } else {
            parent[root2] = root1
            rank[root1] += rank[root2]
        }
    }
    
    NR==FNR {  # 处理过滤后的配对
        if ($1 in exists && $2 in exists) union($1, $2)
        next
    }
    
    END {  # 生成聚类结果
        delete clusters
        for (s in exists) {
            root = find(s)
            clusters[root] = clusters[root] ? clusters[root]","s : s
        }
        
        for (root in clusters) {
            split(clusters[root], members, ",")
            printf "%d\t%s\t%s\n", length(members), root, clusters[root]
        }
    }
    ' "$filtered_pairs" | sort -k1,1nr -k2,2 > "$clusters_raw"
    
    # 写入聚类目录
    mkdir -p "$CLUSTER_DIR"
    rm -f "$OUTPUT_FILE"
    type_num=1
    while IFS=$'\t' read -r count key members; do
        # 写入结果文件
        sorted_members=$(tr ',' '\n' <<< "$members" | sort | tr '\n' ',' | sed 's/,$//')
        echo -e "Type${type_num}\t${sorted_members}" >> "$OUTPUT_FILE"
        
        # 移动文件到独立的聚类目录
        type_subdir="${CLUSTER_DIR}/Type${type_num}"
        mkdir -p "$type_subdir"
        echo "$members" | tr ',' '\n' | xargs -P $THREADS -I {} sh -c "
        src=\"{}\"
        if [ -f \"\$src\" ]; then
            cp -f \"\$src\" \"$type_subdir/\"
        fi
    "
        
        ((type_num++))
    done < "$clusters_raw"
    
    # 清理临时文件
    rm -f "$filtered_pairs" "$clusters_raw"
    
    echo "[$folder_name] 处理完成!"
    echo "[$folder_name] 聚类结果: $OUTPUT_FILE"
    echo "[$folder_name] 聚类目录: $CLUSTER_DIR"
    
    # 返回上级目录
    cd ..
}

# =============================================================================
# 主流程
# =============================================================================

# 检查必要的文件夹是否存在
if [ ! -d "genome_genus" ] || [ ! -d "genome_fp" ]; then
    echo "错误: 需要 genome_genus 和 genome_fp 文件夹"
    echo "请确保这两个文件夹存在于当前目录"
    exit 1
fi

# 检查FastANI是否安装
if ! command -v fastANI &> /dev/null; then
    echo "错误: 未找到 fastANI，请先安装"
    exit 1
fi

# 执行genome_genus的分析
run_ani_clustering "genome_genus" "genome_genus"

# 执行genome_fp的分析  
run_ani_clustering "genome_fp" "genome_fp"