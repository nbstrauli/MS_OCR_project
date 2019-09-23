#!/bin/bash

patient_id=$(/home/hwu/dotfiles/app/jq .patient config.json)

## to obtain data from mysql database and perform cleaning, QC control and finally write to file
stage="111111111111"
echo $stage
Rscript get_data_tcr.R > /tmp/downstream.log

stage="22222222222"
echo $stage
Rscript edge_list_from_distance.R >> /tmp/downstream.log

stage="333333333"
echo $stage
Rscript write_all_clusters.R 
Rscript write_cd4_and_cd8.R  

### generate intermediate connctome
stage="444444444"
echo $stage
Rscript connectome_script_all_clsuter_and_write_number_of_clusters_and_write_overlapping_among_subsets.R  >> /tmp/downstream.log

stage="555555555"
echo $stage
## get simplfied connectome
Rscript connectome_script_CSF_clusters_plus_CSF_to_PB.R >> /tmp/downstream.log

stage="6666666666"
echo $stage
## assess overlap among the subset
Rscript write_subset_overlapping.R >> /tmp/downstream.log


