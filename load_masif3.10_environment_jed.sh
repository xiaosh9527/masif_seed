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

export PATH="$HOME/.local/bin:$PATH"

echo Setting up virtual python environment...
source /home/shxiao/anaconda3/etc/profile.d/conda.sh
conda deactivate
conda activate masif_3.10
echo masif environment succesfully loaded!