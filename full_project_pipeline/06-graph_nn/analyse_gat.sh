#PBS -l walltime=5:00:00
#PBS -l select=1:ncpus=12:mem=100gb
#PBS -N GAT_Analysis_IS
#PBS -j oe
#PBS -o /rds/general/project/hda_24-25/live/TDS/Group09/graph_nn/gat_output_is/job_output_analysis.log

# Create output directory if it doesn't exist
mkdir -p /rds/general/project/hda_24-25/live/TDS/Group09/graph_nn/gat_output_is

cd /rds/general/project/hda_24-25/live/TDS/Group09/graph_nn

# Add memory monitoring
echo "Available memory before run:"
free -h

source /rds/general/user/kja24/home/anaconda3/bin/activate r413

# Run the Python script (updated path)
python scripts/analyse_gat.py

echo "Available memory after run:"
free -h