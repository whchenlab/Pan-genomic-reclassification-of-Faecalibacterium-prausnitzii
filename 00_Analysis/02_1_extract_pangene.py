import csv
import os

def find_all_occurrences(row, strains):
    """Find present coding sequences by strain"""
    occurrences = []
    # Start at column 4
    for i, value in enumerate(row[3:], start=3):
        if value.strip() != '':
            strain_index = i - 3
            occurrences.append((strains[strain_index], value.strip()))
    return occurrences

def extract_gene_sequence_from_strain(faa_file, locus_tag):
    """Extract coding sequences from faa files"""
    if not os.path.exists(faa_file):
        return None, f"file {faa_file} is absent!"
    
    sequence = []
    capture = False
    
    with open(faa_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                # Check annotation for target sequence
                gene_id = line.split()[0][1:]
                if gene_id == locus_tag:
                    capture = True
                    sequence.append(line.strip())
                else:
                    capture = False
            elif capture:
                sequence.append(line.strip())
    
    if not sequence:
        return None, f"sequence {locus_tag} is absent in file {faa_file}"
    
    return '\n'.join(sequence), None

def extract_gene_sequence_auto(faa_dir, occurrences):
    """Find present coding sequences by strains"""
    if not occurrences:
        return None, "target sequence is absent in all strains"
    
    # Try all strains until target sequence is found
    for strain_name, locus_tag in occurrences:
        faa_file = os.path.join(faa_dir, f"{strain_name}.faa")
        sequence, error = extract_gene_sequence_from_strain(faa_file, locus_tag)
        
        if sequence:
            return sequence, None
    
    return None, f"target sequence is absent in {len(occurrences)} strains"

def main():
    csv_file = './panarooresults/gene_presence_absence.csv'
    faa_dir = './faas/'
    
    if not os.path.exists(csv_file) or not os.path.exists(faa_dir):
        return
    
    with open('core.faa', 'w') as core_file, \
         open('shell.faa', 'w') as shell_file, \
         open('cloud.faa', 'w') as cloud_file:
        
        with open(csv_file, 'r') as f:
            # sum strain number and threshold
            reader = csv.reader(f)
            header = next(reader)
            total_strains = len(header) - 3
            core_threshold = total_strains * 0.90
            shell_low_threshold = total_strains * 0.15
            strains = header[3:]
            
            for row in reader:
                gene_name = row[0]
                occurrences = find_all_occurrences(row, strains)
                occurrence_count = len(occurrences)
                # Determin pan-gene type
                if occurrence_count >= core_threshold:
                    output_file = core_file
                elif occurrence_count >= shell_low_threshold:
                    output_file = shell_file
                else:
                    output_file = cloud_file
                
                if not occurrences:
                    continue
                
                sequence, _ = extract_gene_sequence_auto(faa_dir, occurrences)
                if sequence:
                    output_file.write(f"{sequence}\n")

if __name__ == "__main__":
    main()