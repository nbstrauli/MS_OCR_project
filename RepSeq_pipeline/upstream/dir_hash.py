from __future__ import print_function
# Author: Hao Wu <echowuhao@gmail.com>
# for each dir, generate the first 6 char of roothash
import os, glob
import os,sys
import json
import logging
from os.path import join, dirname, basename
import time 
import hashlib
import shutil

ANALYSIS_ROOT ='/home/hwu/analysis'
MIXCR_CSV_ROOT ='/home/hwu/downstream/result'



if __name__ == '__main__':
    dir_hash_fh = open(join(MIXCR_CSV_ROOT, 'status', 'data', 'dir_hash.csv'), 'w+')
    for root, dirs, files in os.walk(ANALYSIS_ROOT):
        # when top is FastaCreat_out*
         # we are in the dir of fasta file
        mixcr_result_dir = basename(root)
        if mixcr_result_dir.startswith('mixcr_trb_extracted') or mixcr_result_dir.startswith('mixcr_igh_extracted'): 
            csv_txt_pattern = os.path.join(root, "*.csv")
            # calculate hash to diffeireate different dirs with same csv file name
            # i.e. same sample are sequence in different time
            root_hash = hashlib.sha224(root.encode('utf-8')).hexdigest()[0:6]
            print("%s,%s" % (root_hash, root), file= dir_hash_fh)

