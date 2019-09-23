#
# Author: Hao Wu <echowuhao@gmail.com>
#
# create barcode file in the dir
# python barcode.py /home/hwu/analysis/109-Stimulation_61515a_Proliferating_T_and_B_cells_dupli_1_8-18-2015_149_119
# will create a barcode.json file in the dir /home/hwu/analysis/109-Stimulation_61515a_Proliferating_T_and_B_cells_dupli_1_8-18-2015_149_119
#
import os, glob
import os,sys
import json
import logging
import time 
import datetime as dt
import shutil
from os.path import join, dirname, basename

from hcvb_lab import make_dir_if_not_exist, hash_of_dir 

ANALYSIS_ROOT ='/home/hwu/analysis'
MIXCR_CSV_ROOT ='/home/hwu/downstream/result'

def copy_csv_file(source_dir_hash , source_dir, csv_file):
    file_name  = basename(csv_file)
    file_name_splited = file_name.split('_')
    patient_dir = os.path.join(MIXCR_CSV_ROOT, file_name_splited[0] )
    patient_visit_dir  = os.path.join(MIXCR_CSV_ROOT, '_'.join(file_name_splited[0:2]))
    for d in [patient_dir, patient_visit_dir]: 
        new_file_name = file_name.replace('.shm', '_' +source_dir_hash + '.shm')
        new_dir_name = join(d, 'shm')
        make_dir_if_not_exist(new_dir_name)
        new_file_path  = join(new_dir_name, new_file_name)
        if not os.path.isfile(new_file_path):
            print(new_file_name)
            shutil.copyfile(join(source_dir,csv_file),new_file_path)

if __name__ == '__main__':
    counter = 0
    new_dirs = []
    for root, dirs, files in os.walk(ANALYSIS_ROOT):
        # when top is FastaCreat_out*
         # we are in the dir of fasta file
        mixcr_result_dir = basename(root)
        if mixcr_result_dir.startswith('shm'): 
            csv_txt_pattern = os.path.join(root, "*.shm")
            # calculate hash to diffeireate different dirs with same csv file name
            # i.e. same sample are sequence in different time
            root_hash = hash_of_dir(root)
            names = glob.glob(csv_txt_pattern)
            for csv_file in names:
                counter = counter  + 1;
                copy_csv_file(root_hash,root, csv_file)


