#!/bin/bash

# Run FastANI
genomes_list="genomes.txt"
find "$INPUT_DIR" -name "*.fna" -o -name "*.fasta" > "$genomes_list"
fastANI --ql "$genomes_list" --rl "$genomes_list" -o "ani_matrix.txt" -t 10

# Extract 95-ANI pairs
ANI_MATRIX="ani_matrix.txt"
OUTPUT_FILE="cluster_results.txt"

filtered_pairs=$(mktemp)
awk -v threshold=95 '
    NR>1 && $3 >= threshold {
        gsub(/.*\//, "", $1)
        gsub(/.*\//, "", $2)
        print $1, $2
    }
' "$ANI_MATRIX" > "$filtered_pairs"

strain_list_file=$(mktemp)
awk 'NR>1 {gsub(/.*\//,"",$1); gsub(/.*\//,"",$2); print $1; print $2}' "$ANI_MATRIX" | sort -u > "$strain_list_file"

# Union-find clustering
clusters_raw=$(mktemp)
awk '
FNR==NR {
    strain_map[$0] = 1
    next
}
{
    strain1 = $1
    strain2 = $2
    if (strain1 in strain_map && strain2 in strain_map) {
        if (!(strain1 in parent)) { parent[strain1] = strain1; rank[strain1] = 1 }
        if (!(strain2 in parent)) { parent[strain2] = strain2; rank[strain2] = 1 }
        union(strain1, strain2)
    }
}
function find(s) {
    if (!(s in parent)) { parent[s] = s; rank[s] = 1; return s }
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
END {
    for (s in strain_map) {
        if (!(s in parent)) {
            parent[s] = s
            rank[s] = 1
        }
    }
    delete cluster_members
    for (s in strain_map) {
        root = find(s)
        if (cluster_members[root]) {
            cluster_members[root] = cluster_members[root] "," s
        } else {
            cluster_members[root] = s
        }
    }
    for (root in cluster_members) {
        split(cluster_members[root], members, ",")
        count = length(members)
        printf "%d\t%s\t%s\n", count, root, cluster_members[root]
    }
}
' "$strain_list_file" "$filtered_pairs" | sort -k1,1nr -k2,2 > "$clusters_raw"

# Write species cluster file
rm -f "$OUTPUT_FILE"
type_num=1
if [ -s "$clusters_raw" ]; then
    while IFS=$'\t' read -r count key members; do
        sorted_members=$(echo "$members" | tr ',' '\n' | sort | tr '\n' ',' | sed 's/,$//')
        echo -e "Type${type_num}\t${sorted_members}" >> "$OUTPUT_FILE"
        ((type_num++))
    done < "$clusters_raw"
fi
rm -f "$filtered_pairs" "$clusters_raw" "$strain_list_file"