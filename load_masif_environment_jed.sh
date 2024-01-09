export MSMS_BIN=/work/upcorreia/bin/apbs-pdb2pqr/bin/msms
export REDUCE_HET_DICT=/work/upcorreia/bin/reduce/reduce_wwPDB_het_dict.txt
export PATH=$PATH:/work/upcorreia/bin/reduce/
module load gcc

# APBS variables
APBS_BIN=/work/upcorreia/bin/apbs-pdb2pqr/bin/apbs
MULTIVALUE_BIN=/work/upcorreia/bin/apbs-pdb2pqr/tools/bin/multivalue
PDB2PQR_BIN=/work/upcorreia/bin/apbs-pdb2pqr/pdb2pqr/pdb2pqr.py
export APBS_BIN
export MULTIVALUE_BIN
export PDB2PQR_BIN
export LD_LIBRARY_PATH=/home/gainza/lpdi_fs/programs/apbs/APBS-1.5-linux64/lib/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/home/shxiao/anaconda3/envs/masif/lib/:$LD_LIBRARY_PATH

# masif_seed_search_root=/scratch/shxiao/masif_seed/masif_seed_search
# masif_root=$masif_seed_search_root/..
# masif_target_root=$masif_seed_search_root/data/masif_targets/
# export masif_db_root=/scratch/shxiao/masif_seed/masif/
# masif_source=$masif_root/source/
# masif_data=$masif_root/data/
# export masif_root
# export masif_target_root

export PATH="$HOME/.local/bin:$PATH"
export PYTHONPATH=$PYTHONPATH:$masif_source:`pwd`

echo Setting up virtual python environment...
source /home/shxiao/anaconda3/etc/profile.d/conda.sh
conda deactivate
conda activate masif
echo masif environment succesfully loaded!