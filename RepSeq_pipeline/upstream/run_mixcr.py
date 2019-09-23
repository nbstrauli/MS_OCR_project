import json,os,glob
from datetime import datetime, timedelta
#from airflow.operators import PythonOperator
#from airflow.models import DAG
from hcvb_lab import is_bcr
import time
import argparse

def run_mixcr(input_file_name, chain):
    mixcr = '/home/hwu/app/mixcr-2.1.3/mixcr'
    mixcr_dir = "mixcr_" + chain.lower()
    output_file_name = input_file_name.replace("raw_fastq",mixcr_dir).replace('.fastq', '.txt')
    align_file_name = output_file_name.replace('.txt', '.vdjca')
    assemle_file_name = output_file_name.replace('.txt', '.clns')
    align_report_file_name = output_file_name.replace('.txt', '.align')
    assemble_report_file_name = output_file_name.replace('.txt', '.assemble')
    align_command = "%s align --chains %s  -r %s %s %s;" % (mixcr, chain ,align_report_file_name , input_file_name,align_file_name)
    assemle_command = "%s assemble -r %s %s %s ;" % (mixcr,assemble_report_file_name, align_file_name,assemle_file_name)
    export_command = " %s exportClones -vHit -jHit -dHit -count -aaFeature CDR3 %s %s;" % (mixcr,assemle_file_name ,output_file_name) 
    command = align_command + assemle_command + export_command
    print(input_file_name)
    print(command)
    os.system(command)

def run_mixcr_for_dir(work_dir):
    # Make sure the output dir is created
    print("mixcr start")
    start = time.time()
    fastq = os.path.join(work_dir, "raw_fastq/*.fastq")
    names = glob.glob(fastq)
    for name in names:
        barcode_name = os.path.basename(name)
        if is_bcr(barcode_name):
            run_mixcr(name, 'IGH')
        else:
            run_mixcr(name, 'TRB')

    print("mixcr done")
    print("total time used %f" % (time.time() - start))
	
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("work_dir")
    args = parser.parse_args()    
    print(args.work_dir) 
    work_dir = args.work_dir
    run_mixcr_for_dir(work_dir)


