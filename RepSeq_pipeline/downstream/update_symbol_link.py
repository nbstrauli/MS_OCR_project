import json
import os,sys,glob
import argparse
import os,sys
from os.path import join, dirname, basename
import time 
import hcvb_lab
from shutil import copyfile
import csv
import re
from datetime import datetime

MIXCR_CSV_ROOT ='/home/hwu/downstream/result'
folder_pattern = re.compile("^[0-9]+\w*$")


def date_from_filename(filename_with_date):
    date_str = filename_with_date[10:]
    date_object = datetime.strptime(date_str , '%b_%d_%Y')
    return date_object


def find_new_dir_bcr(dir_path):
    mixcr_dirs = [ folder_name for folder_name in os.listdir(dir_path) if folder_name.startswith('mixcr_bcr')
            and folder_name != 'mixcr_bcr_final'
            ]
    if mixcr_dirs:
        sorted_dir = sorted(mixcr_dirs, key = date_from_filename)
        link_target = join(dir_path, 'mixcr_bcr_final')
        print(link_target)
        if os.path.exists(link_target):
            os.unlink(link_target)
        os.symlink(join(dir_path, sorted_dir[-1]),link_target)

def find_new_dir_tcr(dir_path):
    mixcr_dirs = [ folder_name for folder_name in os.listdir(dir_path) if folder_name.startswith('mixcr_tcr')
            and folder_name != 'mixcr_tcr_final'
            ]
    if mixcr_dirs:
        sorted_dir = sorted(mixcr_dirs, key = date_from_filename)
        link_target = join(dir_path, 'mixcr_tcr_final')
        print(link_target)
        if os.path.exists(link_target):
            os.unlink(link_target)
        os.symlink(join(dir_path, sorted_dir[-1]),link_target)

for folder_name in os.listdir(MIXCR_CSV_ROOT):
    full_dir_path = join(MIXCR_CSV_ROOT,folder_name)
    if os.path.isdir(full_dir_path) and  folder_pattern.match(folder_name):
        print(folder_name)
        find_new_dir_bcr(full_dir_path)
        find_new_dir_tcr(full_dir_path)




