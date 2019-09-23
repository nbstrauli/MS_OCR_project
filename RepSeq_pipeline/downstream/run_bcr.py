import json
import os,sys
import argparse
import os,sys
from os.path import join, dirname, basename
import time 
import datetime as dt
import hcvb_lab
import subprocess

ANALYSIS_ROOT ='/home/hwu/analysis'
LOG_ROOT= "/home/hwu/analysis_logs/"

MIXCR_CSV_ROOT ='/home/hwu/downstream/result'


def get_revision():
    """get the pipepline code head"""
    binary_head = subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD'])
    return binary_head.decode('ascii').strip()

### 
def run_for_id(patient_id, start_time):
    template_json = {'patient': patient_id, 'length_cutoff' : 8, 'dbname' : 'bcr'}
    specific_config = [(count, distance, method) 
        for count in range(1, 3) 
        for distance in range(0, 3) 
        for method in ['lv', 'hamming']]
    for config_tuple in specific_config:
        config_json = template_json
        count_cut, distance_cut, distance_method = config_tuple
        config_json['count_cutoff'] = count_cut
        config_json['distance_cutoff'] = distance_cut
        config_json['method'] = distance_method
        config_json['date'] = start_time 
        config_json['pipeline_version'] = get_revision()
        with open('config.json', 'w') as f:
            f.write(json.dumps(config_json, sort_keys=True,indent = 4))
        os.system('sh /home/hwu/downstream/mixcr_bcr.sh')

if __name__ == '__main__':
     copy_log_of_today = os.path.join(LOG_ROOT, "copy_log", "%s.log" % (dt.datetime.today().strftime("%Y-%m-%d")))
     changed_ids = hcvb_lab.get_changed_bcr_ids(copy_log_of_today)
     start_time = dt.datetime.now().strftime('%b_%d_%Y')
     for patient_id in changed_ids:
        print(patient_id)
        print("concate data")
        os.system("Rscript /home/hwu/downstream/prepross/concat_csv.R %s" % (patient_id))
        print("generating plot")
        os.system("Rscript /home/hwu/downstream/prepross/plot_vj.R %s" % (patient_id))
        run_for_id(patient_id, start_time)

