#!/bin/bash
#PBS -l walltime=16:00:00
#PBS -l select=1:ncpus=8:mem=64gb
#PBS -N ridge_regression

# Change to working directory
cd /rds/general/project/hda_24-25/live/TDS/Group09/cluster_protein_analysis/multiple_ridge_reg/ridge_reg_scripts

# Load your conda environment
source /rds/general/user/kja24/home/anaconda3/bin/activate r413

# Run the script
Rscript ridge_reg.r