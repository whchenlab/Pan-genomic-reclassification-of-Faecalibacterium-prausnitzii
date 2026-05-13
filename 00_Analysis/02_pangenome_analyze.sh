#!/bin/bash

# Prokka annotation
mkdir -p gffs faas

find filtered_genomes/ -maxdepth 1 -type f \( -name "*.fna" -o -name "*.fa" -o -name "*.fasta" \) | while read -r fna_file; do
    #Set individual output dir
    prefix=$(basename "$fna_file" | sed 's/\.[^.]*$//')
    out_dir="prokka_output_${prefix}"
    #Run prokka
    prokka --cpus 8 --kingdom Bacteria \
        --outdir "$out_dir" \
        --prefix "$prefix" \
        "$fna_file" || { echo "Prokka failed for $prefix"; continue; }
    #Copy annotation files
    cp "${out_dir}/${prefix}.gff" gffs/ 2>/dev/null || echo "Warning: .gff not found for ${prefix}"
    cp "${out_dir}/${prefix}.faa" faas/ 2>/dev/null || echo "Warning: .faa not found for ${prefix}"
done

# Roary pangenome analysis
roary -f roaryresult -p 64 -i 50 -cd 90 -e -n -v --mafft gffs/*gff

# Panaroo pangenome analysis
panaroo -i gffs/*.gff -o panarooresults --clean-mode strict -a core --merge_paralogs --core_threshold 0.90 -t 64

# Extract pangene sequences
python3 02_1_extract_pangene.py
mkdir pangene
mv core.faa shell.faa cloud.faa pangene

# Eggnog-mapper annotation
for faa_file in pangene/*.faa; do
    base_name=$(basename "${faa_file}" .faa)
    emapper.py -i "${faa_file}" --data_dir ./eggnog-mapper-data/ \
      -m diamond --tax_scope bacteria --cpu 16 \
      -o "eggnog_results/${base_name}_eggnog"
done