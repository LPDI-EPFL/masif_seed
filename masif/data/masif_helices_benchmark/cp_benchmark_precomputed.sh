#!/bin/bash
while read p;
do
    cp data_preparation_source/01-benchmark_pdbs/$p\.pdb data_preparation/01-benchmark_pdbs/
    cp data_preparation_source/01-benchmark_surfaces/$p\.ply data_preparation/01-benchmark_surfaces/ 
done < lists/peptide_decoys_random_1000.txt
