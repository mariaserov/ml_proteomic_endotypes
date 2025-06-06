#!/bin/bash
#PBS -l walltime=16:00:00
#PBS -l select=1:ncpus=8:mem=64gb
#PBS -N preprocessing

cd /rds/general/project/hda_24-25/live/TDS/Group09/cluster_protein_analysis/joint_effect_analysis

source /rds/general/user/kja24/home/anaconda3/bin/activate r413

Rscript preprocessing.r