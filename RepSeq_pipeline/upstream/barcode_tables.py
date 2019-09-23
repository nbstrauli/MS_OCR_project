# directory - samplename
#
import glob, os
import csv, json
import hcvb_lab
import pandas as pd
import numpy as np

from hcvb_lab import get_changed_run, copy_log_of_today

ANALYSIS_PATH = '/home/hwu/analysis/'
RESULT_PATH = '/home/hwu/downstream/result/'

def reads_of_fastq_file(file_path):
    ''' try raw_fastq fold, if not find, try 
    raw_fasta fodler'''
    if os.path.exists(file_path):
        num_lines = sum(1 for line in open(file_path))
        return(num_lines/4)
    elif os.path.exists(file_path.replace('raw_fastq', 'raw_fasta')):
        num_lines = sum(1 for line in open(file_path))
        return(num_lines/2)
    else:
        print("checked both fastq and fasta folder, path is {} not exists".format(file_path))
        return 0

def shanon_diversity(file_path):
    if not os.path.exists(file_path):
        print(file_path)
        print("!!!!!!!!!! file not exits !!!!!!!!!!!!!!!!!!")
        return (0,0,0)

    df = pd.read_csv(file_path)
    total_cdr3 = df['count'].sum()
    grouped = df.groupby(['vgene', 'jgene', 'cdr3'])
    df_g= grouped['count'].agg({'count' : 'sum'})
    df_g['p'] = df_g['count']/ df_g['count'].sum()
    df_g['lnpp'] = np.log(df_g['p']) * df_g['p']
    diveristy = 0 -df_g['lnpp'].sum()
    uniq_cdr3 = len(df_g)
    return (uniq_cdr3,total_cdr3,diveristy)


outputFile = open(os.path.join(RESULT_PATH, 'status/data/runs_and_samples.csv'), 'w')
outputWriter = csv.writer(outputFile)

outputWriter.writerow(['barcode', 'sample', 'raw_read','uniq_cdr3_count', 'total_cdr3_count', 'diveristy','directory'])

def recreate_barcode():
    for folder in glob.glob(ANALYSIS_PATH + '*'):
        barcode_file = os.path.join(folder, 'barcode.json')
        with open(barcode_file) as json_data:
            sample_to_barcode = json.load(json_data)
            for k,v  in sample_to_barcode.items():
                reads_num = reads_of_fastq_file(os.path.join(folder, 'raw_fastq', v + '.fastq'))
                if hcvb_lab.is_bcr(v):# if bcr, so go to bcr lib
                    vdj_file_name = os.path.join(folder, 'mixcr_igh_extracted', v + '.csv') 
                else: # tcr
                    vdj_file_name = os.path.join(folder, 'mixcr_trb_extracted', v + '.csv') 
                (uniq_cdr3, total_cdr3, diveristy)  = shanon_diversity(vdj_file_name)
                outputWriter.writerow([k,v, reads_num, uniq_cdr3, total_cdr3, diveristy,  os.path.basename(folder)])


if __name__ == '__main__':
    run_dir = get_changed_run(copy_log_of_today())
    if True or run_dir:
        recreate_barcode()
    else:
        print("no new run today {}".format(copy_log_of_today()))
