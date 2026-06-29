#!/bin/bash 
  
rm iptm_scores.txt plddt_scores.txt

proteins=$(ls -d */)
#proteins="haemoglobin"

for protein in $proteins
do
    cd $protein
    combinations=$(ls -d */)

    for combination in $combinations
    do
        # check if it's run yet 
        did_it_run=$(ls ${combination}pdb/ | wc -l | tail -n 1)
        if [ "$did_it_run" -eq 0 ]; then
            echo "skipping" 
        else
            cd $combination
            cd pdb
            pwd
            jsons=$(ls *scores*.json)

            for json in $jsons
            do
                current_folder=$(pwd)
                iptm=$(grep iptm $json | awk '{print $NF}')
                echo "$current_folder $iptm" >> ../../../iptm2_scores.txt
            done

            pdbs=$(ls *unrelaxed*.pdb)

            for pdb in $pdbs
            do
                current_folder=$(pwd)
                plddt=$(awk '{ total += $11; count++ } END { print total/count }' ${pdb})
            # iptm=$(grep iptm $json | awk '{print $NF}')
                echo "$current_folder $plddt" >> ../../../plddt2_scores.txt
            done
            cd ../
            cd ../
        fi
    done
    cd ../
done
                                                                                                                                                                                                                                                                                                                                                                               