conf_file=$1
ligands_folder=$2
mkdir -p results logs

for LIG in $(ls $ligands_folder/*.pdbqt); do
  BASENAME=$(basename "$LIG" .pdbqt)
  vina \
    --ligand "$LIG" \
    --config $conf_file \
    --out results/"${BASENAME}"_docked.pdbqt > logs/"${BASENAME}".log
done
