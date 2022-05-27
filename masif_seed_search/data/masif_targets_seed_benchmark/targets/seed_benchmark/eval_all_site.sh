#!/bin/bash
while read pdbid_chain;
do
    pdbid=$(echo $pdbid_chain | cut -d"_" -f1)
    recchain=$(echo $pdbid_chain | cut -d"_" -f2)
    cd out_peptides_site/$pdbid\_$recchain
    python3 ../../eval_one.py $pdbid_chain
    cd ../../
    continue
    num_results=$(find . -name "*.score" | grep $pdbid | wc -l)
    echo Num_results: $num_results
    if (( num_results > 0));
    then
#result=$(cat $(find . -name "*.score") |grep -v "clashing_heavy: 5" | grep -v "clashing_heavy: 4"| grep -v "clashing_heavy: 3" |sort -k 7 -n | tac |nl | grep $pdbid| head -1)
        result=$(cat $(find . -name "*.score") |sort -k 7 -n | tac |nl | grep $pdbid| head -1)
        helix_id=$(echo $result | cut -d"," -f1| cut -d" " -f3 | sed -e 's/ //')
        pointid=$(echo $result | cut -d"," -f2 | cut -d":" -f2 | sed -e 's/ //')
        rmsd=$(python3 ../../measure_rmsd.py $pdbid_chain $helix_id $helix_id\_$pointid)
        echo "Success: $pdbid_chain: $result iRMSD: $rmsd"
    else 
        echo "Failed for $pdbid_chain"
    fi
    cd ../../
done < lists/merged_receptor_helix_list.txt

