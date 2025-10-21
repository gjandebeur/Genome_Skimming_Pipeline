# Genome_Skimming_Pipeline
This repository contains the workflow for basecalling, alignment, and modkit analysis of genetic data. Pipeline was built for RNA but can be altered for DNA with relative ease.
Repository contains the scripts to install and run software to identify modification counts across genomic samples.


## Whats Included
  -README.md 
  
  -installation.sh 
  
  -pipeline.sh 
  
  -RNA_ModkitCounts.R 

## Requirements

Scripts are designed to be run on OSCER (OU's supercomputer) using SLURM job scheduler.

Conda is used in this repository as a way to load in dependencies for each software.

Apptainer is used to fulfill version mismatches for GLIBC2.18, a version OSCER has not been updated to run on yet. The apptainer also allows running of each software using its necessary dependencies that are preloaded (see installation script).


# Quick Start
  After loading .sh scripts onto OSCER in respective directories. 
  First run the installation script using "sbatch installation.sh"

  -This will begin by setting up the apptainer needed downstream. The apptainer will be produced in the scratch directory for storage purposes. 
  
  *scratch/ is deleted weekly if the filenames aren't changed or files aren't executed*

### Installation of softwares
  
  The next steps install 3 different softwares: Dorado (basecaller), minimap2 (aligner), modkit (modification detection)

All installations will occur within the specific directory the script was run in, these can be adjusted for organization purposes.


## Running the Pipeline

To run the pipeline, first alter the #USER INPUTS# section of the pipeline.sh script.
This contains the input pod5 directory of raw reads. the sample name, paths to each installation discussed above, and importantly, the parameters

Parameters are currently set for RNA sequencing. To change to DNA alter; 

  -MODEL (https://github.com/nanoporetech/dorado) 

  -MODIFICATIONS (change to DNA mods 5mC/5hmC...)  
  
  -ALIGNMENT_TYPE (change 'rna' to 'dna')

Importantly, for statistical purposes, change:
    
  -MIN-QSCORE =10. Currently set at 10 which is top 90%, but 20 is standard for DNA mods which removes any below 99%. This is the initial QC basecall filterings 
  -For the R script, all modifications need to be switched to DNA for their respective codes.

Pipeline is set up to load necessary dependencies and source to the users conda environment, and to run fully through the modification detection. 

**Important outcome files** 

basecalled.fastq (.bam) - both dorado outputs 

aligned.sam - minimap2 output 

alignment_sorted - SAMtools sorted aligned file 

**stats.txt - SAMtools flagstat output (total reads, aligned reads, alignment %)**

## After Pipeline Runs

To extract the stats out of the pileup output from modkit, run the rna_modificationcount.R Rscript. Modify the codes and any other key features to focus on DNA rather than RNA. 
**I recommend removing samples without a minimum 10 n_valid_cov and keeping modkit thresholds at 0.99 for data accuracy**
