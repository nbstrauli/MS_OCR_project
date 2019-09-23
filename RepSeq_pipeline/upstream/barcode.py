#
# Author: Hao Wu <echowuhao@gmail.com>
#
# create barcode file in the dir
# python barcode.py /home/hwu/analysis/109-Stimulation_61515a_Proliferating_T_and_B_cells_dupli_1_8-18-2015_149_119
# will create a barcode.json file in the dir /home/hwu/analysis/109-Stimulation_61515a_Proliferating_T_and_B_cells_dupli_1_8-18-2015_149_119
#
import os, glob
import json
import re, time
import os, sys
import json
import argparse
from shutil import copyfile


from os.path import join, dirname, basename
from pprint import pprint
from multiprocessing import Pool

ANALYSIS_ROOT ='/home/hwu/analysis'
RAW_DATA_ROOT ='/home/hwu/raw_data'


def analysis_dir_to_raw_data_dir(analysis_dir):
    """we striped the Auto_user part in post_sync.py
    so add it back
    """
    raw_data_dir= os.path.basename(analysis_dir)
    if os.path.basename(analysis_dir)[0].isdigit():
        raw_data_dir =  'Auto_user_SN2-' + raw_data_dir
    return raw_data_dir

def create_barcode_json(work_dir):    
    base_dir = os.path.basename(work_dir)
    barcode_path= join('/home/hwu/barcode', base_dir, 'barcode.json')
    if os.path.exists(barcode_path):
        copyfile(barcode_path, os.path.join(work_dir, 'barcode.json'))
        return

    expmeta_json_file = 'ion_params_00.json'
    expmeta_path = os.path.join(RAW_DATA_ROOT, analysis_dir_to_raw_data_dir(work_dir), expmeta_json_file)
    #read expmeta json
    barcode_dict = {}
    with open(expmeta_path) as f:
          barcode_dict = json.load(f)
    barcode_info = barcode_dict['barcodeInfo']

    barcode_json = {}

    for barcode in barcode_info:
        sample_name = barcode_info[barcode]['sample'] 
        if sample_name != 'none':
            barcode_json[barcode[-3:]] = sample_name.replace(' ','_')

    with open(os.path.join(work_dir, 'barcode.json'), 'w') as f:
          print("creating barcode.json file")
          f.write(json.dumps(barcode_json, sort_keys=True,indent = 4))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("work_dir")
    args = parser.parse_args()    
    print(args.work_dir) 
    create_barcode_json(args.work_dir)
