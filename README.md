# ASV Count Table Workflow: Generating Insights from Amplicon Sequencing

### *By Anika Nag*

---------------------

### Repository Contents- 
## Directories

  + `raw_fastq`: Contains raw FASTQ files used as input for the workflow. 
  
  + `results`: Generated from the pipeline and contains all the outputs, including intermediate and final results. 
        
  + `work`: Generated from the pipeline and stores intermediate files created during Nextflow execution. 
  
  
## Files

  + `ASVs_main.nf`: The main Nextflow script for the ASV workflow, defining the processes and workflow.
  
  + `nextflow.config`: Configuration file for Nextflow that sets up parameters and Docker containers and settings.
  
  + `README.md`: This file, providing an overview of the ASV count table workflow and usage instructions.
  
  + `.gitignore`: Specifies files and directories to exclude from version control. 

---------------------

## Project

This Nextflow workflow processes raw paired-ed RASTQ files and generates an Amplicon Sequenxe Variant (ASV) count matrix. In addition, it also generates many intermediate QC'd files. The workflow automates steps such as quality control, merging, filtering, dereplication, chimera removal and creation of ASV count matrix. This output is ideal for various downstream microbial ecology analyses, such as ASV-species diversity matrix. 

--------------------

## Overview

Microbial communities play an essential role in diverse environments, including human health, soil health, and marine ecosystems. Characterizing these communities often requires high-resolution data about microbial diversity and abundance. Amplicon Sequence Variants (ASVs) provide a more precise representation of microbial community composition than traditional clustering approaches (e.g., OTUs).

This workflow implements steps to generate an ASV count matrix using tools such as `fastp` and `vsearch`. It processes raw FASTQ files into a format ready for ecological and statistical analysis.

### Workflow Summary

This workflow includes the following steps:

1. **Quality Control**: Reads are trimmed to remove low-quality sequences and adapters(`fastp`).
2. **Read Merging**: Paired-end reads are merged into single reads (`vsearch`).
3. **Filtering**: Reads are filtered based on length and quality (`vsearch`).
4. **Concatenation**: Filtered reads are concatenated into a single file for dereplication.
5. **Dereplication**: Unique sequences are identified and their abundances are recorded.
6. **ASV Denoising**: Noise is reduced and ASVs are identified. 
7. **Chimera Removal**: Chimeric sequences are identified and removed.
8. **Singleton Removal**: Low-abundance sequences (singletons) are removed.
9. **Count Matrix Generation**: Reads are mapped back to ASVs to create the count matrix.

Here is visualization of the workflow:

![Workflow Illustration](BIOF501_Project/VSearch_501)

---------------------

## Dependencies

Installing and running this workflow requires the following tools:
- `Nextflow`
- `Docker` (for containerized execution)


The key software dependencies include:
- `fastp` (for quality control)
- `vsearch` (for merging, filtering, dereplication, and chimera detection)

--------------------

## Installation

### Clone the repository:

```bash
git clone git@github.com:anika-n-01/BIOF501_Project.git
cd BIOF501_Project
```
### Usage

To run the pipeline on the FASTQ files in raw_fastq directory, execute the following command:

```bash
nextflow run ASVs_main.nf -with-report report.html -with-timeline timeline.html -with-trace trace.txt -with-dag dag.png
```


