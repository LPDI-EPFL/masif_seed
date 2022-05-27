repo_root=$(git rev-parse --show-toplevel)
masif_root=$repo_root/$masif_seed_search_root/masif/
masif_seed_search_root=$repo_root/masif_seed_search/
masif_target_root=$masif_seed_search_root/data/masif_targets_seed_benchmark/
export masif_db_root=/home/gainza/lpdi_fs/masif_design_paper/masif/
masif_source=$masif_root/source/
masif_data=$masif_root/data/
export masif_root
export masif_target_root
export PYTHONPATH=$PYTHONPATH:$masif_source:`pwd`
if [[ $# -eq 0 ]] ; then
    echo 'You must pass as an argument the PDB id of the target and its chain (e.g. 3R2X_AB)'
    exit 0
fi
pdbid=$(echo $1 | cut -d"_" -f1)
chain1=$(echo $1 | cut -d"_" -f2)
chain2=$(echo $1 | cut -d"_" -f3)
> decoy_list_$pdbid\.tmp
#cp ../benchmark_lists/peptide_decoys_random_1000.txt ./decoy_list_$pdbid\.tmp
echo $pdbid\_$chain2 >> decoy_list_$pdbid\.tmp
python -W ignore $masif_seed_search_root/source/masif_seed_search_nn.py params_gt_site $pdbid\_$chain1 decoy_list_$pdbid\.tmp
