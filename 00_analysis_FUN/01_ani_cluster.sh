#!/bin/bash

# 配置参数
ANI_MATRIX="ani_matrix_genus.txt"  # ANI矩阵文件
FASTA_DIR="genomes_genus"            # 原始基因组文件夹
OUTPUT_FILE="cluster_results_95.txt"  # 聚类结果文件
THREADS=$(nproc)                   # 自动获取CPU核心数

# 预处理ANI矩阵
echo "[1/4] 预处理ANI矩阵..."
filtered_pairs=$(mktemp)
awk -v threshold=95 '
    NF >=3 && $3 >= threshold {print $1,$2}
' "$ANI_MATRIX" > "$filtered_pairs"

# 构建菌株哈希表
echo "[2/4] 构建菌株索引..."
declare -A strain_map
while IFS= read -r -d $'\0' file; do
    strain=$(basename "$file")
    strain_map["$strain"]=1
done < <(find "$FASTA_DIR" -maxdepth 1 -name "*.fna" -print0)

# 并查集操作
echo "[3/4] 执行聚类分析..."
clusters_raw=$(mktemp)
awk -v strains="${!strain_map[*]}" '
BEGIN {
    split(strains, arr, " ")
    for (i in arr) exists[arr[i]] = 1
    
    # 并查集实现
    for (s in exists) {
        parent[s] = s
        rank[s] = 1
    }
}

function find(s) {
    while (parent[s] != s) {
        parent[s] = parent[parent[s]]  # 路径压缩
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

# 移动文件
echo "[4/4] 组织结果文件..."
rm -f "$OUTPUT_FILE"
type_num=1
while IFS=$'\t' read -r count key members; do
    # 写入结果文件
    sorted_members=$(tr ',' '\n' <<< "$members" | sort | tr '\n' ',' | sed 's/,$//')
    echo -e "Type${type_num}\t${sorted_members}" >> "$OUTPUT_FILE"
    
    # 并行移动文件
    type_dir="Type${type_num}"
    mkdir -p "$type_dir"
    echo "$members" | tr ',' '\n' | xargs -P $THREADS -I {} sh -c "
    src=\"$FASTA_DIR/{}\"
    if [ -f \"\$src\" ]; then
        cp -f \"\$src\" \"$type_dir/\"
    fi
"
    
    ((type_num++))
done < "$clusters_raw"

# 清理临时文件
rm -f "$filtered_pairs" "$clusters_raw"

echo "处理完成！结果保存至 $OUTPUT_FILE"