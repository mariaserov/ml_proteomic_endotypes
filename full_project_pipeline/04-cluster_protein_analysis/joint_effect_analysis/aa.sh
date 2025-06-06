#!/bin/bash
#PBS -l walltime=16:00:00
#PBS -l select=1:ncpus=12:mem=120gb
#PBS -N joint_effect_aa
#PBS -o /rds/general/project/hda_24-25/live/TDS/Group09/cluster_protein_analysis/joint_effect_analysis/output/aa.log
#PBS -e /rds/general/project/hda_24-25/live/TDS/Group09/cluster_protein_analysis/joint_effect_analysis/output/aa.err

# Change to working directory
cd /rds/general/project/hda_24-25/live/TDS/Group09/cluster_protein_analysis/joint_effect_analysis

# Load your conda environment
source /rds/general/user/kja24/home/anaconda3/bin/activate r413

# Run the script
Rscript joint_effect_aa.r