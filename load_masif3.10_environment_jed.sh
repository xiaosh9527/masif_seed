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

export PATH="$HOME/.local/bin:$PATH"
export PYTHONPATH=$PYTHONPATH:$masif_source:`pwd`

echo Setting up virtual python environment...
source /home/shxiao/anaconda3/etc/profile.d/conda.sh
conda deactivate
conda activate masif_3.10
echo masif environment succesfully loaded!