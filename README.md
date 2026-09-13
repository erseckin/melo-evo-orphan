# melo-evo-orphan

Scripts and command lines used to identify orphan proteins in *Meloidogyne* species and investigate their possible evolutionary origins.

This repository accompanies the study on the evolution and characteristics of orphan genes in root-knot nematodes of the genus *Meloidogyne*.

## Overview

Orphan genes are genes without detectable homologs outside a defined taxonomic group. In this project, we identify genus-specific orphan proteins in *Meloidogyne* and combine comparative proteomics, synteny-based analyses, expression evidence, protein feature extraction, and machine learning to characterize their origin and properties.

The main goals of the pipeline are to:

1. identify *Meloidogyne* orphan proteins using comparative proteome analyses;
2. filter and validate orphan candidates using expression evidence;
3. infer probable mechanisms of origin, including *de novo* emergence and high divergence from ancestral genes;
4. generate feature tables describing orphan and non-orphan proteins;
5. train machine learning classifiers to distinguish orphan from non-orphan proteins and *de novo* from highly diverged orphan proteins.

## Workflow summary

The analysis is organized around the following main steps.

### 1. Proteome collection and quality control

Predicted proteomes were collected for several *Meloidogyne* species and outgroup nematodes. Proteome quality was assessed using tools such as BUSCO and OMArk before downstream comparative analyses.

### 2. Comparative proteome analysis

Orthology and homology relationships were inferred using comparative proteome tools, including OrthoFinder and SonicParanoid. These analyses were used to identify *Meloidogyne* proteins without detectable homologs outside the genus.

### 3. Orphan protein identification

Candidate orphan proteins were defined as proteins present in *Meloidogyne* but lacking detectable homologs in non-*Meloidogyne* nematodes and other compared taxa. Additional filtering was applied to retain high-confidence candidates and reduce the effect of annotation artifacts.

### 4. Expression support

Transcriptomic and other omics data were used to identify expressed orphan proteins. This step helps distinguish biologically supported orphan genes from unsupported gene predictions.

### 5. Evolutionary origin inference

Synteny and comparative genomic context were used to infer probable mechanisms of orphan gene origin. Candidate orphan genes were classified into categories:

- putative *de novo* genes;
- highly diverged genes;
- ambiguous origin.

### 6. Feature table generation

Protein-level features were calculated and merged into final feature tables. These features include, among others:

- sequence length;
- amino acid composition;
- physicochemical properties;
- predicted localization;
- signal peptide prediction;
- transcriptional evidence;
- proteomic or ribosome profiling evidence when available;
- orphan status and inferred origin.

### 7. Machine learning classifiers

Random Forest classifiers were trained to test whether sequence-derived and predicted protein features can distinguish:

1. orphan proteins from non-orphan proteins;
2. candidate *de novo* orphan proteins from highly diverged orphan proteins.

These classifiers provide a fast strategy for prioritizing candidate orphan proteins and predicting the likely origin of orphan genes in newly available *Meloidogyne* proteomes.

## Folders

### `cmd/`

This folder contains command lines used to run external programs and reproduce the main computational analyses. These commands are provided for transparency and reproducibility.

Typical analyses include:

- proteome quality control;
- orthology inference;
- comparative proteome analyses;
- synteny-related analyses;
- downstream formatting and filtering steps.

### `treat/`

This folder contains scripts used to process and reformat outputs from external tools. These scripts are used to clean intermediate results, merge tables, add annotations, and generate final datasets used in downstream analyses.

### `classifier/`

This folder contains scripts related to machine learning analyses, including the training and evaluation of Random Forest classifiers for orphan gene identification and origin inference.

## Requirements

The pipeline uses a combination of external bioinformatics tools and Python/R scripts. Depending on the analysis step, the following tools may be required:

- Python 3
- R
- OrthoFinder
- SonicParanoid
- BUSCO
- OMArk
- MCScanX / synteny-related tools
- standard Unix command-line tools

Python packages may include:

- pandas
- numpy
- scikit-learn
- matplotlib
- seaborn
- biopython

Additional dependencies may be required by specific scripts.

Moreover, additional data needed to test some steps can be found on: https://entrepot.recherche.data.gouv.fr/dataverse/orphan-rkn

## Usage

This repository is organized as a collection of analysis scripts and command lines rather than as a single automated software package. The general usage is:

1. prepare input proteomes and annotation files;
2. run the command lines in `cmd/` for quality control and comparative analyses;
3. process intermediate results using scripts in `treat/`;
4. generate final feature tables;
5. run machine learning analyses from `classifier/`.

For each script, check the input paths, output paths.

## Input data

The main input files are predicted proteomes, genome annotations, orthology inference outputs, synteny outputs, and omics-derived evidence tables.

Large input datasets are not necessarily included in this repository. Users should adapt the scripts to their own local paths and file organization.

## Output data

The pipeline produces several types of output files, including:

- lists of candidate orphan proteins;
- orphan orthogroups;
- classifications of orphan origin;
- processed synteny tables;
- protein feature tables;
- machine learning training and prediction outputs.

## Citation

If you use this repository or adapt parts of the code, please cite the associated manuscript:

```text
Seçkin, E. et al. Orphan genes shape the genome and parasitic arsenal of root-knot nematodes. bioRxiv, 2026. doi: 10.64898/2025.12.19.695360. 
```

## Contact

For questions about the repository or the analysis, please contact:

```text
Ercan Seçkin ercan.seckin@inrae.fr
```


Please check the repository license before reuse. If no license is provided, the code is shared for transparency but reuse rights are not explicitly granted.
