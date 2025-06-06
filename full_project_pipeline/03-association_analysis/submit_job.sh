#!/bin/bash
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=1:mem=100gb
#PBS -N outcome_association_job

# Change directory to the location where the job was submitted 
cd /rds/general/project/hda_24-25/live/TDS/Group09/outcome_definition/Scripts

# Load the required Anaconda module and activate your environment
module load anaconda3/personal
source activate r413

# Define the path to the folder containing the R scripts
def_path="/rds/general/user/jdc124/projects/hda_24-25/live/TDS/Group09/full_pipeline/4_outcome_association"

# Loop over each .R file in the folder (non-recursively) and run it
for rfile in "$def_path"/*.R; do
    if [ -f "$rfile" ]; then
        echo "Running $rfile"
        Rscript "$rfile"
    fi
done

