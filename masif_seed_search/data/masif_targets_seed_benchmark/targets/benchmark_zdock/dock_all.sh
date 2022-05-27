#!/bin/sh
while read ppi_pair;
do
	echo "Starting $ppi_pair"
	./dock_one.sh $ppi_pair &
        sleep 15
done < ../benchmark_lists/merged_receptor_helix_list.txt
