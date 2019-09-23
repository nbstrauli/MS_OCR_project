import json
import os,sys
import argparse
import os,sys
from os.path import join, dirname, basename
import time 
import datetime as dt
import hcvb_lab

ANALYSIS_ROOT ='/home/hwu/analysis'
LOG_ROOT= "/home/hwu/analysis_logs/"

MIXCR_CSV_ROOT ='/home/hwu/downstream/result'

def run_for_id(patient_id, start_time):
    template_json = {'patient': patient_id, 'length_cutoff' : 8, 'dbname' : 'tcr'}
    specific_config = [(5, 0 , 'hamming'),(10, 0 , 'hamming'), (1, 0 , 'hamming')]
    for config_tuple in specific_config:
        config_json = template_json
        count_cut, distance_cut, distance_method = config_tuple
        config_json['count_cutoff'] = count_cut
        config_json['distance_cutoff'] = distance_cut
        config_json['method'] = distance_method
        config_json['date'] = start_time 
        print(config_json)
        with open('config.json', 'w') as f:
            f.write(json.dumps(config_json, sort_keys=True,indent = 4))
        os.system('sh /home/hwu/downstream/mixcr_tcr.sh')

if __name__ == '__main__':
     copy_log_of_today = os.path.join(LOG_ROOT, "copy_log", "%s.log" % (dt.datetime.today().strftime("%Y-%m-%d")))
     changed_ids = hcvb_lab.get_changed_tcr_ids(copy_log_of_today)
     start_time = dt.datetime.now().strftime('%b_%d_%Y')
     for patient_id in changed_ids:
        print(patient_id)
        print("concate data")
        os.system("Rscript /home/hwu/downstream/prepross/concat_csv.R %s" % (patient_id))
        print("generating plot")
        os.system("Rscript /home/hwu/downstream/prepross/plot_vj.R %s" % (patient_id))
        run_for_id(patient_id, start_time)

