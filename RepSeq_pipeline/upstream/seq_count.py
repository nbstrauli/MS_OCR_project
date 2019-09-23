import os,sys,glob
import logging

ANALYSIS_ROOT ='/home/hwu/analysis'

def seq_count_for_dir(work_dir):
    fastq = os.path.join(work_dir, "raw_fastq/*.fastq")
    names = glob.glob(fastq)
    for name in names:
        num_lines = sum(1 for line in open(name))
        print(','.join([os.path.basename(name),str(num_lines)]))

	
print("sample_name, seq_count")
for run_dir in os.listdir(ANALYSIS_ROOT):
    seq_count_for_dir(os.path.join(ANALYSIS_ROOT,run_dir))


        
