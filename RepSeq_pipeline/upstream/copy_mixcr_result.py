#
# Author: Hao Wu <echowuhao@gmail.com>
#
# create barcode file in the dir
# python barcode.py /home/hwu/analysis/109-Stimulation_61515a_Proliferating_T_and_B_cells_dupli_1_8-18-2015_149_119
# will create a barcode.json file in the dir /home/hwu/analysis/109-Stimulation_61515a_Proliferating_T_and_B_cells_dupli_1_8-18-2015_149_119
#
import os, glob,sys
import json
import logging
from os.path import join, dirname, basename
import time 
import datetime as dt
import hashlib
import shutil
from hcvb_lab import make_dir_if_not_exist, hash_of_dir, is_bcr

ANALYSIS_ROOT ='/home/hwu/analysis'
LOG_ROOT= "/home/hwu/analysis_logs/"
MIXCR_CSV_ROOT ='/home/hwu/downstream/result'



def copy_csv_file(source_dir_hash , source_dir, csv_file):
    file_name  = basename(csv_file)
    file_name_splited = file_name.split('_')
    patient_dir = os.path.join(MIXCR_CSV_ROOT, file_name_splited[0] )
    patient_visit_dir  = os.path.join(MIXCR_CSV_ROOT, '_'.join(file_name_splited[0:2]))

    # find if a sequence is bcr or tcr
    bcr_or_tcr= ''
    if is_bcr(file_name):
        bcr_or_tcr = 'bcr'
    else:
        bcr_or_tcr = 'tcr'

    for d in [patient_dir, patient_visit_dir]: 
        new_file_name = file_name.replace('.csv', '_' +source_dir_hash + '.csv')
        new_dir_name = join(d, bcr_or_tcr)
        make_dir_if_not_exist(new_dir_name)
        new_file_path  = join(new_dir_name, new_file_name)
        if not os.path.isfile(new_file_path):
            logging.info("copying %s" ,csv_file)
            shutil.copyfile(join(source_dir,csv_file),new_file_path)




if __name__ == '__main__':
    LOG_FILENAME =  os.path.join(LOG_ROOT, "copy_log", "%s.log" % (dt.datetime.today().strftime("%Y-%m-%d")))
    logging.basicConfig(filename=LOG_FILENAME,
                        level=logging.DEBUG,
                        filemode='a',
                        format='%(asctime)s %(message)s')

    logging.getLogger().addHandler(logging.StreamHandler())
    logging.info('%s start', __file__)


    counter = 0
    # dir to escape
    escape_dirs = ['00000']
    for root, dirs, files in os.walk(ANALYSIS_ROOT):
        # when top is FastaCreat_out*
         # we are in the dir of fasta file
        mixcr_result_dir = basename(root)
        analysis_dir = basename(dirname(root))
        if mixcr_result_dir.startswith('mixcr_trb_extracted') or mixcr_result_dir.startswith('mixcr_igh_extracted') and analysis_dir not in escape_dirs:
            csv_txt_pattern = os.path.join(root, "*.csv")
            # calculate hash to diffeireate different dirs with same csv file name
            # i.e. same sample are sequence in different time
            root_hash = hash_of_dir(root)
            names = glob.glob(csv_txt_pattern)
            for csv_file in names:
                counter = counter  + 1;
                copy_csv_file(root_hash,root, csv_file)

    logging.info( "number of total csv file is %s" , counter)

