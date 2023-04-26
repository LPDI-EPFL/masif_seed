masif_root=$(git rev-parse --show-toplevel)
masif_source=$masif_root/masif/source/
masif_data=$masif_root/masif/data/
# Compute geometric invariant descriptors, as implemented in Yin et al. PNAS 2009
export PYTHONPATH=$PYTHONPATH:$masif_source:$masif_data/masif_ppi_search/
python $masif_source/gif_descriptors/compute_gif_descriptors.py lists/ransac_benchmark_list.txt
