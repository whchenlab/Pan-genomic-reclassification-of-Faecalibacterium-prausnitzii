#!/bin/bash

# Set input core genes
INPUT_DIR="split/split_proteins"

# Define species clusters
FP_STRAINS="F.prausnitzii-ATCC-27768-a F.prausnitzii-ATCC-27768-b F.tardum-CLA-AA-H175 F.intestinale_CLA-AA-H281"
FD_STRAINS="F.duncaniae-A2-165-a F.duncaniae-A2-165-b F.duncaniae-A2-165-c F.duncaniae-JCM-31915"
FL_STRAINS="GCA_000210735 GCA_003287485 GCA_003287505 GCA_003293635"
FF_STRAINS="F.faecis-CLA-JM-H7-B F.taiwanense-HLW78 GCA_000166035 GCA_002549755 GCA_002549895"     
FW_STRAINS="F.wellingii-HTF-F-a Fw_GCA_018381455 Fw_GCA_018381555 Fw_GCA_037478275" 
FH_STRAINS="GCA_008671395 GCA_965147195 GCA_002549975 GCA_001406355"
FB_STRAINS="GCA_020687345 GCA_040094315 GCA_048528205 GCA_900758465"
FLA_STRAINS="GCA_902470615 GCA_032584235"
FUN1_STRAINS="GCA_041326985 GCA_048402285 GCA_902490285 GCA_934541075"
FUN2_STRAINS="GCA_017377625 GCA_002550035 GCA_003478405 GCA_019969425"
FUN3_STRAINS="GCA_003287415 GCA_900539065 GCA_003287495 GCA_045655365"
FUN4_STRAINS="GCA_048396855 GCA_018374075 GCA_048526205 GCA_037248125"

# Set output classifier
OUTPUT_DIR="hmm_classifier"
mkdir -p $OUTPUT_DIR
mkdir -p "${OUTPUT_DIR}/Fp" "${OUTPUT_DIR}/Fd" "${OUTPUT_DIR}/Fl" "${OUTPUT_DIR}/Ff" "${OUTPUT_DIR}/Fw" "${OUTPUT_DIR}/Fh" "${OUTPUT_DIR}/Fb" "${OUTPUT_DIR}/Fla" "${OUTPUT_DIR}/Fun1" "${OUTPUT_DIR}/Fun2" "${OUTPUT_DIR}/Fun3" "${OUTPUT_DIR}/Fun4"

# Function: split core-gene by species-cluster
split_alignment() {
    local gene_name=$1
    local group=$2
    local strains=$3
    
    TMP_FILE=$(mktemp)
    INPUT_ALIGNMENT="${INPUT_DIR}/${gene_name}.aln"
    
    # Extract genes in each species cluster
    for strain in $strains; do
        grep -A1 "$strain" "$INPUT_ALIGNMENT" | grep -v "^--$" >> "$TMP_FILE"
    done
    
    # Save gene_species file
    OUTPUT_ALIGNMENT="${OUTPUT_DIR}/${group}/${gene_name}_${group}_alignment.aln"
    cp "$TMP_FILE" "$OUTPUT_ALIGNMENT"
    echo "Save ${gene_name}_${group} alignment file: ${OUTPUT_ALIGNMENT}"
    rm "$TMP_FILE"
    
    # Construct hmm model
    HMM_FILE="${OUTPUT_DIR}/${group}/${gene_name}_${group}.hmm"
    hmmbuild "$HMM_FILE" "$OUTPUT_ALIGNMENT"
    echo "Save ${gene_name}_${group} HMM model: ${HMM_FILE}"
}

# Process core genes
GENE_FILES=$(find "$INPUT_DIR" -name "*.aln" -type f | sed 's#.*/\(.*\)\.aln#\1#')

echo "Find $(echo "$GENE_FILES" | wc -w) core genes, start processing..."
for gene in $GENE_FILES; do
    echo -e "\nProcess: $gene"
    # Split core-gene by species-cluster
    split_alignment "$gene" "Fp" "$FP_STRAINS"
    split_alignment "$gene" "Fd" "$FD_STRAINS"
    split_alignment "$gene" "Fl" "$FL_STRAINS"
    split_alignment "$gene" "Ff" "$FF_STRAINS"
    split_alignment "$gene" "Fw" "$FW_STRAINS"
    split_alignment "$gene" "Fh" "$FH_STRAINS"
    split_alignment "$gene" "Fb" "$FB_STRAINS"
    split_alignment "$gene" "Fla" "$FLA_STRAINS"
    split_alignment "$gene" "Fun1" "$FUN1_STRAINS"
    split_alignment "$gene" "Fun2" "$FUN2_STRAINS"
    split_alignment "$gene" "Fun3" "$FUN3_STRAINS"
    split_alignment "$gene" "Fun4" "$FUN4_STRAINS"
done

# Merge model for each species cluster
for group in Fp Fd Fl Ff Fw Fh Fb Fla Fun1 Fun2 Fun3 Fun4; do
    GROUP_DIR="${OUTPUT_DIR}/${group}"
    MERGED_HMM="${OUTPUT_DIR}/${group}_models.hmm"
    
    cat "${GROUP_DIR}"/*.hmm > "$MERGED_HMM"
    hmmpress "$MERGED_HMM"
    
    echo "Construct merged HMM model for ${group} species: ${MERGED_HMM}"
done

# Merge all species models
MERGED_ALL_MODEL="${OUTPUT_DIR}/merged.hmm"
HMM_FILES=$(find "$OUTPUT_DIR" -maxdepth 1 -name "*_models.hmm" -type f)

cat $HMM_FILES > "$MERGED_ALL_MODEL"
hmmpress "$MERGED_ALL_MODEL"

echo "HMM model for all species is saved at: ${MERGED_ALL_MODEL}"