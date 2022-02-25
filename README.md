# Masif-seed - An enhanced representation of protein structures that enables de novo design of protein interactions

This repository contains code to design de novo binders based on surface fingerprints. The code was used to perform the experiments in: [citation]

## Table of Contents: 

- [Description](#description)
- [System and hardware requirements](#system-and-hardware-requirements)
- [Running through a docker container](#running-through-a-docker-container)
- [Method overview](#Method-overview)
     * [MaSIF data preparation](#MaSIF-data-preparation)
- [Step-by-step example](#code-overview)
- [Reproducing the benchmarks](#reproducing-the-benchmark)
- [Running the code on PD-L1, PD-1, and RBD](#running-the-code-on-pd-l1,-pd-1,-and-RBD)
- [PyMOL plugin](#PyMOL-plugin)
- [Docker container](#Docker-container)
- [License](#License)
- [Reference](#Reference)

## Description

## System and hardware requirements

MaSIF-seed has been tested on Linux (Red Hat Enterprise Linux Server release 7.4, with a Intel(R) Xeon(R) CPU E5-2650 v2 @ 2.60GHz 
processesor and 16GB of memory allotment). We find some of the code extremely slow on the newer M1 processors.
To reproduce the experiments in the paper, the entire datasets for all proteins consume several terabytes. 

Currently, MaSIF takes a few seconds to preprocess every protein. We find the main bottleneck to be the APBS computation for surface charges,
which can likely be optimize. Nevertheless, we recommend a distributed cluster to 
preprocess the data for large datasets of proteins. If retraining, we strongly recommend using a GPU to 
train or evaluate the trained models as it can be up to 100 times faster than a CPU; for inference, we find a CPU is enough.

## Running through a docker container

Since Masif-seed relies on a few external programs (msms, APBS) and libraries (pyMesh, tensorflow, scipy, open3D), 
we strongly recommend you use the Dockerfile and Docker container. 



## Preliminaries

Log into deneb1 or 2 (remember your choice)

I strongly recommend you _always_ use tmux on the server when using MaSIF. 
Add the following line to your ~/.bashrc file: 

```
module load tmux 
```

Then create a new session under the name masif:

```
tmux new -t masif
```

You have now created a long lived session in deneb.

## Downloading masif.

Attach your tmux session if you it is not attached yet: 

```
# tmux attach -t masif
```

Load the masif environment:
``` 
# source /work/upcorreia/bin/load_masif_environment.sh
```

and create a directory where you will run masif seed search on. 
```
# mkdir masif_runs
# cd masif_runs
```

clone the masif repository: 

``` 
# git clone https://github.com/LPDI-EPFL/masif
Cloning into 'masif'...
...

```

The masif-seed-search methods are located on a different Github repository, to protect it while we publish this paper. Clone the masif seed search:

```
# cd masif
# git clone https://github.com/pablogainza/masif_seed_search/
```

Your directory structure within the parent masif repository should look like this: 

```
# ls
citation.bib  comparison  data  docker_tutorial.md  img  LICENSE  masif_seed_search  README.md  requirements.txt  source
```


## Preparing your structure for masif seed search on a GPU

You can do this protocol on a CPU but it is going to be slower than a GPU for some parts (specifically the part of computing the descriptors). 


Request an interactive GPU session. 
```
Sinteract -p gpu -q gpu -g gpu:1 -t 40:0:0 -m 32000
```

Load the MaSIF environment for GPUs: 
```
source /work/upcorreia/bin/load_masif_environment_gpu.sh
```

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

This steps can be performed on a CPU or a GPU - the difference is smaller between the two for this process. 

First request a cpu or gpu. You may need it for a while, perhaps 24 hours.

```
Sinteract -p gpu -q gpu -g gpu:1 -t 24:0:0 -m 32000
```

If you can't get a GPU, during weekdays you can access our exclusive partition:
```
Sinteract -p gpu -r lpdi-gpu -q gpu -g gpu:1 -t 8:0:0 -m 32000
```
Go into the run subdirectory that was created: 

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



## License

MaSIF is released under an [Apache v2.0 license](LICENSE).

## Reference
If you use this code, please use the bibtex entry in [citation.bib](citation.bib)

