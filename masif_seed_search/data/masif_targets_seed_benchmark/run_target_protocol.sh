#!/bin/bash
./data_prepare_one.sh $1
./predict_site.sh $1
./color_site.sh $1
./compute_descriptors.sh $1

