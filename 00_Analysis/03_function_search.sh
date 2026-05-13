#!/bin/bash

mkdir -p blast_results logs
folders=("faas:blast_results")
target_proteins=(*.faa)
if [ ${#target_proteins[@]} -eq 0 ]; then
    echo "Error: no function coding sequence(*.faa) is found in current dir"
    exit 1
fi

# Set Blastp threshold
EVALUE="10"
IDENTITY="80"
COVERAGE="70"
THREADS=16

# Init summary file
report_file="blast_summary_$(date +%Y%m%d).txt"
echo "Function coding sequences - $(date)" > "$report_file"
echo "=========================" >> "$report_file"
echo "Thresholds: evalue=$EVALUE, identity=$IDENTITY%, coverage=$COVERAGE%" >> "$report_file"

# Run Blastp
for target in "${target_proteins[@]}"; do
    target_name=$(basename "$target" .ffn)
    echo -e "\nAnalyzing target gene: $target_name" >> "$report_file"
    
    # Search in genus/species cluster
    for folder_pair in "${folders[@]}"; do
        IFS=':' read -r protein_folder result_folder <<< "$folder_pair"
        bacteria_proteins=("$protein_folder"/*.faa)
        if [ ${#bacteria_proteins[@]} -eq 0 ]; then
            echo "Warning: no genome(.faa) file is found in $protein_folder"
            continue
        fi   
        echo -e "\n  Processing cluster: $protein_folder" >> "$report_file"
        
        found_count=0
        not_found_count=0
        
        for protein in "${bacteria_proteins[@]}"; do
            protein_name=$(basename "$protein" .faa)
            log_file="logs/${target_name}_vs_${protein_name}_${protein_folder//\//_}.log"
            echo -n "    Starting blast $protein_name..."
            
            # Construct genome databases
            makeblastdb -in "$protein" -dbtype prot -out "temp_db" > "$log_file" 2>&1
            
            # Run Blastp
            blastp -query "$target" -db "temp_db" \
                   -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \
                   -evalue "$EVALUE" \
                   -out "temp_blast.txt" > /dev/null 2>&1
            
            if [ -s "temp_blast.txt" ]; then
                # Extract query length
                qlen=$(awk 'NR==1 {print $13}' temp_blast.txt)
                if [[ -z "$qlen" || "$qlen" -eq 0 ]]; then
                    echo " invalid 0 length"
                    ((not_found_count++))
                    cp "temp_blast.txt" "${result_folder}/${target_name}_vs_${protein_name}_invalid.txt"
                    continue
                fi
                
                # Calculate coverage
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
                    echo " Present"
                    ((found_count++))
                    cp "temp_blast.txt" "${result_folder}/${target_name}_vs_${protein_name}.txt"
                    echo "Coverage: $(cat coverage.txt)" >> "${result_folder}/${target_name}_vs_${protein_name}.txt"
                else
                    echo " Absent(Low coverage)"
                    ((not_found_count++))
                    cp "temp_blast.txt" "${result_folder}/${target_name}_vs_${protein_name}_low.txt"
                fi
            else
                echo " Absent(No hitting)"
                ((not_found_count++))
            fi
            
            echo "Blastp $target_name vs $protein_name finish" >> "$log_file"
            echo "Results: $(if [ -s "temp_blast.txt" ]; then echo "matching"; else echo "not matching"; fi)" >> "$log_file"
            rm -f temp_db.* temp_blast.txt coverage.txt
        done
        
        echo "    Summarize: $target_name is present in $found_count strains, is absent in $not_found_count strains" >> "$report_file"
    done
done