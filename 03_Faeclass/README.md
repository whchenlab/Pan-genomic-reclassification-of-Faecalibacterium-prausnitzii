# FaeClass:  Faecalibacterium Species Distinguish

FaeClass, a core gene-based classifier, enables robust discrimination among closely related Faecalibacterium species, including *F.longum, F.prausnitzii, F.duncaniae, F.taiwanense, F.indigenum, F.centralis, F.hattorii, F.butyricigenerans, F.bellum, F.protidificans, F.wellingii, F.langellae*, and other *Faecalibacterium* species.

## Table of Contents

* FaeClass: Faecalibacterium Species Distinguish

  * Table of Contents

  * Software requirement

  * Installation

  * Usage

  * Output Explanation

  * Citation

## Software requirement

* Python 3.6+

* HMMER 3.1+

## Installation

To install FaeClass, download the model directory (hmm_classifier) and the execution script (hmm_classifier.py) from the GitHub repository.

* **hmm_classifier/**

  * **Fl, Fp, Fd, ...**: Species-specific directories.

    * **.aln**: Core gene alignment files for each species group.

    * **.hmm**: Individual single-gene HMM models for each species group.

  * **_models.hmm**: Single-species HMM databases.

  * **_models.hmm.h3***: Compressed index files for each single-species HMM database.

  * **merged.hmm**: Complete HMM database combining all species groups.

  * **merged.hmm.h3***: Compressed index files for the complete merged database.

* **hmm_classifier.py**: Main script implementing the classification functionality.

## Usage

* **Basic Usage**

  python3 hmm_classifier.py /path/to/genome_dir

* **Advanced Usage with Custom Parameters**

  python3 hmm_classifier.py /path/to/genome_dir --HMM_DIR my_models --OUTPUT_DIR my_results

* **Parameter Description**

  * /path/to/genome_dir: Input directory containing sequence files for classification.

  * --HMM_DIR: (Optional) Custom directory containing HMM database files.

  * --OUTPUT_DIR: (Optional) Custom directory for storing classification results.

* **Input Requirements**

  * Files should be standard protein sequence (FASTA format) with the ".faa" file extension, either from public databases or derived from de novo genome annotation pipelines.

## Output Explanation

All results will be saved in the hmm_results directory (customizable via the --output-dir parameter).

* **Directory Structure**

  hmm_results/&#x20;

  ├── details/                           # Detailed Results Directory&#x20;

  │   ├── [Sample1]_hmm_results.txt      # Raw HMM Search Results&#x20;

  │   ├── [Sample1]_hmm_results_dom.txt  # HMM Domain-Level Details&#x20;

  │   ├── [Sample2]_hmm_results.txt&#x20;

  │   └── ...&#x20;

  └── summary.txt                       # Final Summary File

* **Running Example**

  ```markup
  cd FaeClass
  python3 hmm_classifier.py test
  ```

  * **[Sample]_hmm_results.txt**

    · Complete sequence alignment results from HMM search, containing global alignment statistics for each target sequence against HMM models.

    |         Column        |       Description       |
    | :-------------------: | :---------------------: |
    |      target name      |   Target sequence name  |
    |       query name      |      HMM model name     |
    |        E-value        |  Full sequence E-value  |
    |         score         | Full sequence bit score |
    | description of target |  Functional description |

  ![](README_md_files/68dc7200-bcb7-11f0-9a45-433d18e6f797.jpeg?v=1&type=image)

  * **[Sample]_hmm_results_dom.txt**

    · Domain-level alignment results from HMM search, providing detailed alignment information for each individual domain.

    |       Column      |          Description         |
    | :---------------: | :--------------------------: |
    |      c-Evalue     |      Conditional E-value     |
    |      i-Evalue     |      Independent E-value     |
    | hmm coord from/to |      HMM model positions     |
    | ali coord from/to | Alignment sequence positions |

  * **summary.txt**

    · Summary table of classification result for all samples, including final species identification and score information.

    |    Column   |                     Description                     |
    | :---------: | :-------------------------------------------------: |
    |    Genome   |                  Genome sample name                 |
    |    Msocre   |       Total score of the best-matching species      |
    | Mean_Mscore | Average score per gene in the best-matching species |
    |   Species   |              The best-matching species              |

  ```markup
  Genome	MScore	Mean_MScore	Species
  AF27-11BH	73639.40	150.90	F.faecis
  AF31-14AC	72918.20	148.81	F.prausnitzii
  AF32-8AC	72848.00	149.28	F.longum
  AF52-21	71691.70	150.61	F.butyricigenerans
  AHMP21	73619.20	150.24	F.indigenum
  AM33-14AC	73594.60	151.12	F.duncaniae
  AM37-13AC	72974.50	148.93	F.prausnitzii
  AM43-5AT	73538.50	150.69	F.faecis
  APC918_95b	73099.40	149.18	F.prausnitzii
  APC922_41-1	71243.90	150.94	F.hattorii
  APC923_51-1	72651.00	148.57	F.prausnitzii
  APC923_61-1	72881.40	151.21	F.bellum
  APC924_119	73056.00	149.09	F.prausnitzii
  APC942_8-14-2	72566.00	150.55	F.bellum
  ATCC_27768	73342.50	149.68	F.prausnitzii
  ```

  ## Citation

  If you use FaeClass, please cite:

  Li et al. [期刊]. Pan-genomic reclassification of Faecalibacterium prausnitzii sensu lato reveals F. longum as a dominant, functionally distinct, and health-associated gut anaerobe. Reference click here

