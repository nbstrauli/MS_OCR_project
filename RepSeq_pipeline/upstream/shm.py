#!/usr/bin/python
import os,sys,csv
import time 
import re
import glob
import argparse

def find_barcode(mixcr_file):
    return os.path.basename(mixcr_file).rstrip(".txt")


def export_shm(export_file_name, shm_file_name):

    input_fh = open(export_file_name, 'r')
    mixcr_reader = csv.reader(input_fh, delimiter='\t')
    # This skips the first row of the CSV file.
    # csvreader.next() also works in Python 2.
    #
    # in rare case, file size is zero
    barcode = find_barcode(export_file_name)

    if os.path.getsize(export_file_name):
        next(mixcr_reader)
        output_fh = open(shm_file_name, 'wb')
        writer = csv.writer(output_fh)

        header_info = ['shm', 'barcode']
        writer.writerow(header_info)
        regex = r"\d+"
        for row in mixcr_reader:
            J_align  =row[6]
            J_mute = J_align.split('|')[5]
            J_mute_count = 0
            if J_mute:
                J_mute_count =len(re.findall(regex, J_mute))
            writer.writerow([J_mute_count, barcode])
        output_fh.close()
    else:
        print("!!! NO valid alignments in sample %s" % barcode)
    input_fh.close()



def run_mixcr_bcr(aligned_file_name, export_file_name):
    mixcr = '/home/hwu/app/mixcr-2.1.3/mixcr'
    command = " %s exportAlignments %s %s;" % (mixcr,aligned_file_name, export_file_name) 
    os.system(command)

def run_mixcr_for_dir(work_dir):
    # Make sure the output dir is created
    print("export align start")
    start = time.time()
    mixcr_igh = os.path.join(work_dir, "mixcr_igh/*.vdjca")
    names =  glob.glob(mixcr_igh)
    for name in names:
        export_file_name = name.replace("mixcr_igh",
            "mixcr_igh_export_align").replace('.vdjca', '.txt')
        run_mixcr_bcr(name, export_file_name)
        shm_file_name = export_file_name.replace('mixcr_igh_export_align', 'shm').replace('.txt', '.shm')
        export_shm(export_file_name, shm_file_name)

    print("export align done")
    print("total time used %f" % (time.time() - start))
	

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("work_dir")
    args = parser.parse_args()    
    print(args.work_dir) 
    work_dir = args.work_dir
    run_mixcr_for_dir(work_dir)

