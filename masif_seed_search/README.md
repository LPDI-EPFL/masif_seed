# masif_seed_search
This directory contains the code and scripts to perform seed searching with masif-seed. It relies on scripts/code on the masif base as well.

# Procedure to run masif seed search. 

## Preliminaries

I strongly recommend you _always_ use tmux when using MaSIF. 

Create a new session under the name masif:

```
tmux new -t masif
```

## Preparing your structure for masif seed search on a GPU

You can do this protocol on a CPU but it is going to be slower than a GPU for some parts (specifically the part of computing the descriptors). 

Go into the actual run directory.
```
# cd masif_seed_search/data/masif_targets/
```

You must now select a chain (or multiple chains) and a PDB id to target. You can also target a specific PDB file. 

Let's suppose you want to target the binding site of protein 4QVF, chain A (this is BclxL, which we know binds a peptide). You can run this protocol as follows: 

```
# ./run_target_protocol.sh 4QVF_A
``` 

This process should generate a site prediction for MaSIF-site and the descriptors for masif-search. If you wish to visualize the predictions for MaSIF-site you can do so using the masif pymol plugin. The predictions are located in the subdirectory: 

```
output/all_feat_3l/pred_surfaces/4QVF_A.ply
```

## Performing an actual search.

This steps can be performed on a CPU or a GPU - the difference is smaller between the two for this process and CPU is just fine.  

```
cd targets/4QVF_A/
```

```
./run.sh 4QVF_A
```

## Configuring parameters

There are some parameters that may improve your search in 'params_peptides.py':

The main criteria for speed is the 'descriptor distance cutoff'. This cutoff determines which fingerprints are further considered and which are completely discarded. The lower the value, the faster the search (and the higher the number of false negatives):
```
# Score cutoffs -- these are empirical values, if they are too loose, then you get a lot of results.
# Descriptor distance cutoff for the patch. All scores below this value are accepted for further processing.
params['desc_dist_cutoff'] = 1.7 # Recommended values: [1.5-2.0] (lower is stricter)
```

Another important cutoff value is the interface cutoff. All patches are scored due to their interface propensity. One can assume that during the search we want to find peptides with a high interface cutoff. You can lower this value to increase the number of candidates:

```
# Interface cutoff value, all patches with a score below this value are discarded.
params['iface_cutoff'] = 0.75 # Recommended values: [0.75-0.95] range (higher is stricter)
```

Finally, an important value is the neural network score cutoff. This neural network cutoff is the second stage ranking. The lower the value, the more candidates are output. In practical terms, a value above 0.8 is good:

```
# Neural network score cutoff - Discard anything below this score
params['nn_score_cutoff'] = 0.8 # Recommended values: [0.8-0.95] (higher is stricter)
```

The number of sites to target in the protein (generally one should be fine? ):

```
# Number of sites to target in the protein
params['num_sites'] = 1
```

## Clustering helical peptides. 

The distribution provides a few scripts to cluster helical results in a specific site. To do this I recommend reducing strictness of the desc_dist_cutoff and nn_score_cutoff parameters to increase the number of results, for example to the following values: 
```
params['desc_dist_cutoff'] = 2.5
params['nn_score_cutoff'] = 0.90
```

The following files should be on your ```masif_seed_search/data/masif_targets/targets/template``` template directory: 
```
align_all_to_all_fixed_length.slurm 
align_all_to_all_setup.sh
MDS_HELICES_RMSD.ipynb
```

Copy them to your run directory and then run in the following order: 

```
./align_all_to_all_setup.sh
sbatch align_all_to_all_fixed_length.slurm 
```

Once the batch job finishes, the data can be plotted in the MDS_HELICES_RMSD.ipynb jupyter notebook. 

## Features in development

The goal is to improve the quality, speed, and usability of masif-seed-search. The following features are in process or in advanced stages of development: 

+ Improved electrostatics (currently these have low influence) 




