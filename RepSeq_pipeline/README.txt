The code herein comprises the bioinformatic pipeline for processing immune repertoire (RepSeq) data. The author of this code is Hao Wu. The code is divided into two parts: 'upstream' and 'downstream'. Below gives the order of commands to run the full pipeline.

## Hao's order for running both the Upstream and Downstream Pipeline ##

python ./upstream/upstream.py
python ./upstream/last_change_date.py
python ./upstream/seq_count.py >./upstream/seq_count.csv
python ./upstream/last_change_date.py

# update the runs_and_sample.txt in ~/downstream/result
# the copy log have to be exists, all the downstream are dependent on this
python ./upstream/copy_mixcr_result.py 

python ./upstream/barcode_tables.py 

# shm
# TODO: plot, make shm indenpendt
python ./upstream/shm_changed.py  && python ./upstream/copy_shm_result.py 
# dir hash
python ./upstream/dir_hash.py 

#be carefull about the time to the pipeline runs
# this one have to be in the end, or there will be bizzard behavour
cd ./downstream/ &&  python  run_bcr.py && python run_tcr.py &&  python update_symbol_link.py
