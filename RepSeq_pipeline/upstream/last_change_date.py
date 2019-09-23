#!/usr/bin/python
from __future__ import print_function
import os,sys,csv
import time 
import re
import glob
import argparse

ANALYSIS_ROOT ='/home/hwu/analysis'
MIXCR_CSV_ROOT = '/home/hwu/downstream/result'
import os.path, time

def get_last_modified_time(target_file):
    return "%s" % time.ctime(os.path.getmtime(target_file))


if __name__ == '__main__':
    #run_mixcr_for_dir(work_dir)
    fh = open(os.path.join(MIXCR_CSV_ROOT, 'status', 'data', 'last_changed.csv'), 'w+')
    print("%s,%s" % ('sample_name', 'date'), file=fh)
    for root, dirs, files in os.walk(ANALYSIS_ROOT):
        if os.path.basename(root) == 'mixcr_igh_extracted': 
            run_dir = os.path.dirname(root)
            file_and_date = [(os.path.splitext(f)[0],get_last_modified_time(os.path.join(root,f)))for f in files]
            for fd in file_and_date:
                print("%s,%s" % fd, file=fh)
