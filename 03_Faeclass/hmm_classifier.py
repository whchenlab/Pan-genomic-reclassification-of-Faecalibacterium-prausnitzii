#!/usr/bin/env python3
"""
Script to run HMM classifier: perform HMM search on all .faa files in the test directory,
save the results, and generate a classification summary.
"""

import os
import sys
import subprocess
import argparse
import glob
from collections import defaultdict

# ==================== Parse arguments ====================
def parse_args():
    parser = argparse.ArgumentParser(description="Run HMM search and summarize classification results")
    parser.add_argument("TEST_DIR", help="Directory containing test sequence .faa files (required)")
    parser.add_argument("--HMM_DIR", default="hmm_classifier",
                        help="Directory containing HMM models (default: hmm_classifier)")
    parser.add_argument("--OUTPUT_DIR", default="hmm_results",
                        help="Output directory for results (default: hmm_results)")
    return parser.parse_args()

# ==================== Species mapping ====================
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

# ==================== Helper functions ====================
def check_hmm_press(model_path):
    """Check if the HMM model has been compressed with hmmpress"""
    for ext in [".h3f", ".h3i", ".h3m", ".h3p"]:
        if not os.path.exists(model_path + ext):
            return False
    return True

def run_hmmsearch(model_file, query_faa, output_tbl, output_dom):
    """Run hmmsearch, return True if successful"""
    cmd = [
        "hmmsearch",
        "--tblout", output_tbl,
        "--domtblout", output_dom,
        model_file,
        query_faa
    ]
    try:
        # Redirect stdout to /dev/null (silent mode)
        with open(os.devnull, 'w') as devnull:
            subprocess.run(cmd, check=True, stdout=devnull, stderr=subprocess.PIPE)
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error: hmmsearch execution failed, return code {e.returncode}", file=sys.stderr)
        return False

def parse_hmm_results(tbl_file):
    """
    Parse a single hmmsearch --tblout result file
    Returns a dictionary: {species: (total_score, count)}
    Rule: keep the highest score for each query, then aggregate by species
    """
    # Store the highest score and species for each query
    query_max_score = {}
    query_species = {}

    with open(tbl_file, 'r') as f:
        for line in f:
            if line.startswith('#') or line.strip() == '':
                continue
            parts = line.split()
            if len(parts) < 7:
                continue
            # column 3: query name, column 6: full sequence score
            query = parts[2]
            try:
                score = float(parts[5])
            except ValueError:
                continue
            if score <= 0:
                continue

            # Match species based on query name
            species = None
            for pattern, sp in SPECIES_PATTERNS.items():
                import re
                if re.search(pattern, query):
                    species = sp
                    break
            if species is None:
                continue

            # Keep the highest score for each query
            if query not in query_max_score or score > query_max_score[query]:
                query_max_score[query] = score
                query_species[query] = species

    # Aggregate by species
    species_total = defaultdict(float)
    species_count = defaultdict(int)
    for query, score in query_max_score.items():
        sp = query_species[query]
        species_total[sp] += score
        species_count[sp] += 1

    return species_total, species_count

def determine_species(species_total):
    """
    Select the species with the highest total score, and decide final classification based on mean score threshold.
    Returns (max_species, max_total, mean_score)
    """
    if not species_total:
        return "Outgroup", 0.0, 0.0

    # Find species with highest total score
    max_species = max(species_total, key=lambda sp: species_total[sp])
    max_total = species_total[max_species]
    return max_species, max_total

def process_summary(details_dir, output_dir):
    """Process all *_hmm_results.txt files in the details directory and generate summary.txt"""
    summary_file = os.path.join(output_dir, "summary.txt")
    with open(summary_file, 'w') as out_f:
        out_f.write("Genome\tMScore\tMean_MScore\tSpecies\n")

        # Find all result files
        pattern = os.path.join(details_dir, "*_hmm_results.txt")
        result_files = glob.glob(pattern)
        if not result_files:
            print(f"Warning: No *_hmm_results.txt files found in {details_dir}", file=sys.stderr)
            return

        for tbl_file in result_files:
            # Extract genome name (consistent with bash script logic)
            basename = os.path.basename(tbl_file)
            # Remove _hmm_results.txt suffix
            genome_name = basename.replace("_hmm_results.txt", "")

            # Parse results
            species_total, species_count = parse_hmm_results(tbl_file)

            if not species_total:
                # No valid matches
                out_f.write(f"{genome_name}\t0.00\t0.00\tOutgroup\n")
                continue

            # Find species with highest total score
            max_species = max(species_total, key=lambda sp: species_total[sp])
            max_total = species_total[max_species]
            mean_score = max_total / species_count[max_species] if species_count[max_species] > 0 else 0.0

            # Decide final classification based on thresholds
            if mean_score < 120:
                final_species = "Outgroup"
            elif mean_score < 139:
                final_species = "F.else"
            else:
                final_species = max_species

            out_f.write(f"{genome_name}\t{max_total:.2f}\t{mean_score:.2f}\t{final_species}\n")
            print(f"Processed: {tbl_file} -> added to summary")

    print(f"Final summary file generated: {summary_file}")
    # Count records (excluding header)
    with open(summary_file, 'r') as f:
        line_count = sum(1 for _ in f) - 1
    print(f"Total records: {line_count}")

# ==================== Main function ====================
def main():
    args = parse_args()
    test_dir = args.TEST_DIR
    hmm_dir = args.HMM_DIR
    output_dir = args.OUTPUT_DIR
    details_dir = os.path.join(output_dir, "details")

    # 1. Check input directory
    if not os.path.isdir(test_dir):
        print(f"Error: Test sequence directory '{test_dir}' does not exist!", file=sys.stderr)
        sys.exit(1)

    # 2. Check merged HMM model
    merged_model = os.path.join(hmm_dir, "merged.hmm")
    if not os.path.isfile(merged_model):
        print(f"Error: HMM model file '{merged_model}' does not exist!", file=sys.stderr)
        sys.exit(1)

    if not check_hmm_press(merged_model):
        print(f"Error: HMM model file '{merged_model}' is not compressed! Please run hmmpress.", file=sys.stderr)
        sys.exit(1)

    # 3. Create output directory
    os.makedirs(details_dir, exist_ok=True)

    # 4. Find all .faa files
    faa_files = glob.glob(os.path.join(test_dir, "*.faa"))
    if not faa_files:
        print(f"Error: No .faa files found in '{test_dir}'!", file=sys.stderr)
        sys.exit(1)

    total_files = len(faa_files)
    print(f"Found {total_files} .faa files, starting processing...")

    # 5. Run hmmsearch for each file
    for idx, faa_file in enumerate(faa_files, start=1):
        print(f"\n[{idx}/{total_files}] Processing file: {faa_file}")
        base_name = os.path.splitext(os.path.basename(faa_file))[0]
        results_file = os.path.join(details_dir, f"{base_name}_hmm_results.txt")
        dom_results_file = os.path.join(details_dir, f"{base_name}_hmm_results_dom.txt")

        print("  Running HMM search...")
        success = run_hmmsearch(merged_model, faa_file, results_file, dom_results_file)
        if not success:
            print(f"  Error: HMM search failed for {faa_file}!", file=sys.stderr)
            continue

        print(f"  Results saved to: {results_file}")
        print(f"  Domain details saved to: {dom_results_file}")

    print("\nHMM search completed!")
    print(f"Results stored in: {output_dir}/")
    print(f"Detailed results stored in: {details_dir}/")

    # 6. Generate summary file
    process_summary(details_dir, output_dir)

if __name__ == "__main__":
    main()