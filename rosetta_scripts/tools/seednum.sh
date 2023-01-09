###
# This file is adding a unique seed number ('_S[ID]') to each input PDB containing a seed in complex with the target
###

i=1
for file in $(ls ./input/*.pdb)
  do  base=$(basename $file .pdb)
      mv $file $base\_S$i.pdb
      i=`expr $i + 1`
      echo $i
  done
