# Masif-seed - An enhanced representation of protein structures that enables de novo design of protein interactions

This repository contains code to design de novo binders based on surface fingerprints. The code was used to perform the experiments in: [citation]

## Table of Contents: 

- [Description](#description)
- [Method overview](#Method-overview)
- [System and hardware requirements](#system-and-hardware-requirements)
- [Running through a docker container](#running-through-a-docker-container)
     * [MaSIF data preparation](#MaSIF-data-preparation)
- [Step-by-step example](#code-overview)
- [Reproducing the benchmarks](#reproducing-the-benchmark)
- [Running the code on PD-L1, PD-1, and RBD](#running-the-code-on-pd-l1,-pd-1,-and-RBD)
- [PyMOL plugin](#PyMOL-plugin)
- [Docker container](#Docker-container)
- [License](#License)
- [Reference](#Reference)

## Description

## Method overview

## System and hardware requirements

MaSIF-seed has been tested on Linux. We find some of the code extremely slow on the newer M1 processors.
To reproduce the experiments in the paper, the entire datasets for all proteins consume several terabytes. 

Currently, MaSIF takes a few seconds to preprocess every protein. We find the main bottleneck to be the APBS computation for surface charges,
which can likely be optimize. Nevertheless, we recommend a distributed cluster to 
preprocess the data for large datasets of proteins. A GPU is strongly recommended, especially for computing masif-search and masif-site 
inference, as it can be 60 times faster!

## Running through a docker container

Since Masif-seed relies on a few external programs (msms, APBS) and libraries (pyMesh, tensorflow, scipy, open3D), 
we strongly recommend you use the Dockerfile and Docker container. 

```
git clone http://https://github.com/LPDI-EPFL/masif_seed/
cd masif_seed
docker build . -t masif_seed 
docker run -it -v $PWD:$PWD masif_seed
```

## Step-by-step example

We will test masif seed using one example consisting of a single helix (BH3) and a receptor (Bcl-xL):

```
cd masif/data/masif_peptides/
```

First split the BIM BH3, as crystallized in PDB id: 4QVF, chain B into helices: 

```
./data_extract_helix_one.sh 4QVF_B
```

This example contains exactly one helix, `4QVF000_B`. Now precompute the features for this 
helix, including the geodesic coordinates: 

```
./data_precompute_patches_one.sh 4QVF000_B
```

Finally, compute the Masif-site prediction and the Masif-search descriptors. 

```
./predict_site.sh 4QVF000_B
./compute_descriptors.sh 4QVF000_B
```

Once the site predictions and descriptors on `4QVF000_B` have been computed, we 
can focus on the target. 

```
cd ../../../
cd masif_seed_search/data/masif_targets/
```

The features, Masif-site and Masif-search descriptors must be computed as well for the target, 
as well as a surface with per-vertex coloring. 

```
./run_target_protocol.sh 4QVF_A
cd targets
cp -r template/ bclxl/
cd bclxl
```

Finally, run the script to match Bclxl to all precomputed peptides: 

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

