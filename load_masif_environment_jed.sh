export MSMS_BIN=/work/upcorreia/users/gainza/scripts/msms/msms
export REDUCE_HET_DICT=/home/shxiao/.local/reduce_wwPDB_het_dict.txt
export PATH=$PATH:/home/shxiao/.local/reduce/
module load gcc

# APBS variables
APBS_BIN=/home/gainza/lpdi_fs/programs/apbs/APBS-1.5-linux64/bin/apbs
MULTIVALUE_BIN=/home/gainza/lpdi_fs/programs/apbs/APBS-1.5-linux64/share/apbs/tools/bin/multivalue
PDB2PQR_BIN=/home/gainza/lpdi_fs/programs/apbs/pdb2pqr-linux-bin64-2.1.1/pdb2pqr
export APBS_BIN
export MULTIVALUE_BIN
export PDB2PQR_BIN
export LD_LIBRARY_PATH=/home/gainza/lpdi_fs/programs/apbs/APBS-1.5-linux64/lib/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/home/shxiao/anaconda3/envs/masif/lib/:$LD_LIBRARY_PATH

masif_seed_search_root=/scratch/shxiao/masif_seed/masif_seed_search
masif_root=$masif_seed_search_root/..
masif_target_root=$masif_seed_search_root/data/masif_targets/
export masif_db_root=/scratch/shxiao/masif_seed/masif/
masif_source=$masif_root/source/
masif_data=$masif_root/data/
export masif_root
export masif_target_root

export PATH="$HOME/.local/bin:$PATH"
export PYTHONPATH=$PYTHONPATH:$masif_source:`pwd`

echo Setting up virtual python environment...
source /home/shxiao/anaconda3/etc/profile.d/conda.sh
conda deactivate
conda activate masif
echo masif environment succesfully loaded!