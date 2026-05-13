#!/bin/bash

# Set input and output
INPUT_DIR="genomes"
OUTPUT_DIR="analysis_results"
mkdir -p "$OUTPUT_DIR"

# =============================================================================
# Step1: BUSCO
# =============================================================================
LINEAGE="bacteria_odb10"
MODE="genome"
CPU=4

# Set report file
echo -e "Status\tBin_ID\tLevel\tComplete\tFragmented\tDuplicated\tReason" > "$OUTPUT_DIR/pass_report_busco.txt"
echo -e "Status\tBin_ID\tLevel\tComplete\tFragmented\tDuplicated\tReason" > "$OUTPUT_DIR/fail_report_busco.txt"

# Main function
process_busco() {
    local fna_file="$1"
    local level="$2"
    local bin_id=$(basename "$fna_file" .fna)
    
    # Run BUSCO
    local output_dir="$OUTPUT_DIR/busco_results/${bin_id}_${level}"
    mkdir -p "$(dirname "$output_dir")"
    busco -i "$fna_file" -l "$LINEAGE" -o "$output_dir" -m "$MODE" --cpu "$CPU" &>/dev/null
    
    # Extract results
    local summary_file="${output_dir}/run_${LINEAGE}/short_summary.txt"
    local line=$(grep "C:" "$summary_file")
    local complete=$(echo "$line" | sed -E 's/.*C:([0-9]+\.?[0-9]*)%.*/\1/')
    local fragmented=$(echo "$line" | sed -E 's/.*F:([0-9]+\.?[0-9]*)%.*/\1/')
    local duplicated=$(echo "$line" | sed -E 's/.*D:([0-9]+\.?[0-9]*)%.*/\1/')

    # Explain reason
    local reason=""
    (( $(echo "$complete < 50" | bc -l) )) && reason+="Complete <50%; "
    (( $(echo "$fragmented > 10" | bc -l) )) && reason+="Fragmented >10%; "
    (( $(echo "$duplicated > 5" | bc -l) )) && reason+="Duplicated >5%; "

    # Write report
    if [ -z "$reason" ]; then
        echo -e "Pass\t$bin_id\t$level\t$complete\t$fragmented\t$duplicated\t-" >> "$OUTPUT_DIR/pass_report_busco.txt"
    else
        echo -e "Fail\t$bin_id\t$level\t$complete\t$fragmented\t$duplicated\t${reason%%; }" >> "$OUTPUT_DIR/fail_report_busco.txt"
    fi
}

# Main loop
for level in genomes; do
    if [ -d "./$level" ]; then
        for fna in "./$level"/*.fna; do
            [ -f "$fna" ] && process_busco "$fna" "$level"
        done
    fi
done
# =============================================================================
# Step2: CheckM
# =============================================================================
level="genomes"
comp_thresh=70
contamination_threshold=5
strain_threshold=100
threads=16
export PYTHONRECURSIONLIMIT=10000

# Set report file
echo -e "Status\tBin_ID\tLevel\tCompleteness\tContamination\tStrain_heterogeneity\tReason" > "$OUTPUT_DIR/pass_report_checkm.txt"
echo -e "Status\tBin_ID\tLevel\tCompleteness\tContamination\tStrain_heterogeneity\tReason" > "$OUTPUT_DIR/fail_report_checkm.txt"

# Main function
process_checkm() {
    local bin_id="$1"
    local level="$2"
    local checkm_output_dir="$OUTPUT_DIR/checkm_results/$level"
    
    # Extract results
    if grep -q "^${bin_id}" "$checkm_output_dir/quality_report.tsv"; then
        local line=$(grep "^${bin_id}" "$checkm_output_dir/quality_report.tsv")
        IFS=$'\t' read -ra cols <<< "$line"
        completeness="${cols[5]}"
        contamination="${cols[6]}"
        strain_heterogeneity="${cols[7]}"
    else
        echo -e "FAIL\t${bin_id}\t${level}\t0.00\t0.00\t0.00\tfind no result" >> "$report_fail"
        return 1
    fi

    # Explain reason
    reason=""
    status="PASS"
    (( $(echo "$completeness < $comp_thresh" | bc -l) )) && { reason+="completeness<70; "; status="FAIL"; }
    (( $(echo "$contamination > $contamination_threshold" | bc -l) )) && { reason+="contamination>5; "; status="FAIL"; }
    (( $(echo "$strain_heterogeneity > $strain_threshold" | bc -l) )) && { reason+="heterogeneity>100; "; status="FAIL"; }

    # Write report
    record=$(printf "%s\t%s\t%.2f\t%.2f\t%.2f" "$bin_id" "$level" "$completeness" "$contamination" "$strain_heterogeneity")
    if [ "$status" == "PASS" ]; then
        echo -e "PASS\t$record\t-" >> "$report_pass"
    else
        echo -e "FAIL\t$record\t${reason%; }" >> "$report_fail"
    fi
}

# Run CheckM
checkm_output_dir="$OUTPUT_DIR/checkm_results/$level"
mkdir -p "$checkm_output_dir"
quality_report="$checkm_output_dir/quality_report.tsv"

if [ ! -f "$quality_report" ]; then
    if ! checkm qa "$checkm_output_dir/lineage.ms" "$checkm_output_dir" \
        -o 2 --tab_table -f "$quality_report"; then
        echo "Fail to run CheckM"
        exit 1
    fi
fi

# Run main function
bin_ids=($(tail -n +2 "$quality_report" | cut -f1))
total=${#bin_ids[@]}

for ((i=0; i<$total; i++)); do
    bin_id="${bin_ids[$i]}"
    process_checkm "$bin_id" "$level"
done
echo

# =============================================================================
# Step3: GUNC
# =============================================================================
gunc_output="$OUTPUT_DIR/gunc_results"
mkdir -p "$gunc_output"
gunc_db="gunc_db_progenomes2.1.dmnd"

# Run GUNC
gunc run --inpud_dir "$INPUT_DIR" --out_dir "$gunc_output" --threads 32 --db_file "$gunc_db"

# =============================================================================
# Step4: GTDB-Tk
# =============================================================================
gtdbtk_output="$OUTPUT_DIR/gtdbtk_results"
mkdir -p "$gtdbtk_output"

# Run GTDB-Tk
echo "Run GTDB-Tk classify_wf..."
gtdbtk classify_wf \
  --genome_dir "$INPUT_DIR" \
  -x fna \
  --out_dir "$gtdbtk_output" \
  --cpus 32 \

# =============================================================================
# Step5: FastANI
# =============================================================================
# Run FastANI
echo "Run FastANI..."
genomes_list="$OUTPUT_DIR/genomes.txt"
find "$INPUT_DIR" -name "*.fna" -o -name "*.fasta" > "$genomes_list"
fastANI --ql "$genomes_list" --rl "$genomes_list" -o "$OUTPUT_DIR/ani_matrix.txt" -t 32 --minFraction 0.6

# Extract redundant pairs
ANI_MATRIX="$OUTPUT_DIR/ani_matrix.txt"
OUTPUT_FILE="$OUTPUT_DIR/redundant_cluster.txt"

filtered_pairs=$(mktemp)
awk -v threshold=99.9 '
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

# Write redundant cluster file
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