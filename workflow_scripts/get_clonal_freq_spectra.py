import os
import sys
sys.path.insert(0, '../python_classes')
from immune_time_series_class import immune_time_series

def get_freq_spectra(input_dirpath, output_dirpath, truncate_to, bcr_or_tcr):
    """
    Uses 'immune_time_series_class.py' to write the frequency spectra of each of the samples to disk.
    """

    ####### parameters #######
    min_num_reads = 1000
    #these are only used if use_counts_or_freqs='freqs'
    #this gives the width (span) of each bin in the freq spectrum
    #This is currently set to 1/2355 because the sample with the
    #lowest number of reads in OPERA/ORATORIO BCR data is 2355
    #lowest number of reads in OPERA/ORATORIO TCR data is 18876
    if bcr_or_tcr == 'bcr':
        freq_bin_width = 1/2354.
    elif bcr_or_tcr == 'tcr':
        freq_bin_width = 1/18875.
    ####### parameters #######

    if input_dirpath[-1] != '/':
        input_dirpath += '/'
    if output_dirpath[-1] != '/':
        output_dirpath += '/'
    if not os.path.exists(output_dirpath):
        os.makedirs(output_dirpath)
    for i in os.listdir(input_dirpath):
        if i[0] == '.' or i[:6] == 'README':
            continue
        output_dirpath_isotype = '%s%s/' % (output_dirpath, i)
        if not os.path.exists(output_dirpath_isotype):
            os.makedirs(output_dirpath_isotype)
        for j in os.listdir(input_dirpath + i):
            if j[0] == '.' or j[:6]== 'README':
                continue
            input_dirpath_patient = '%s%s/%s/' % (input_dirpath, i, j)
            print input_dirpath_patient
            output_dirpath_patient = '%s%s/' % (output_dirpath_isotype, j)
            if not os.path.exists(output_dirpath_patient):
                os.makedirs(output_dirpath_patient)
            samples = immune_time_series(input_dirpath_patient, 'count', 'V_gene', None, 'J_gene', 'CDR3')
            freq_spectra, bin_edges_list = samples.get_clonal_freq_spectra(define='v_j_cdr3', use_counts_or_freqs='freqs', truncate_to=truncate_to, write_output=output_dirpath_patient, min_num_reads=min_num_reads, freq_bin_width=freq_bin_width)
    return

if __name__ == '__main__':
    get_freq_spectra(input_dirpath='/Users/nstrauli/data/ms/clone_seqs_fasta_format/bcr/opera', output_dirpath='/Users/nstrauli/data/ms/clonal_diversity/clonal_freq_spectra/bcr/opera', truncate_to=None, bcr_or_tcr='bcr')
    get_freq_spectra(input_dirpath='/Users/nstrauli/data/ms/clone_seqs_fasta_format/tcr/opera', output_dirpath='/Users/nstrauli/data/ms/clonal_diversity/clonal_freq_spectra/tcr/opera', truncate_to=None, bcr_or_tcr='tcr')
    get_freq_spectra(input_dirpath='/Users/nstrauli/data/ms/clone_seqs_fasta_format/bcr/oratorio', output_dirpath='/Users/nstrauli/data/ms/clonal_diversity/clonal_freq_spectra/bcr/oratorio', truncate_to=None, bcr_or_tcr='bcr')
    get_freq_spectra(input_dirpath='/Users/nstrauli/data/ms/clone_seqs_fasta_format/tcr/oratorio', output_dirpath='/Users/nstrauli/data/ms/clonal_diversity/clonal_freq_spectra/tcr/oratorio', truncate_to=None, bcr_or_tcr='tcr')
    get_freq_spectra(input_dirpath='/Users/nstrauli/data/ms/clone_seqs_fasta_format_clustered_data/bcr/opera', output_dirpath='/Users/nstrauli/data/ms/clonal_diversity/clonal_freq_spectra_clustered_data/bcr/opera', truncate_to=None, bcr_or_tcr='bcr')
    get_freq_spectra(input_dirpath='/Users/nstrauli/data/ms/clone_seqs_fasta_format_clustered_data/tcr/opera', output_dirpath='/Users/nstrauli/data/ms/clonal_diversity/clonal_freq_spectra_clustered_data/tcr/opera', truncate_to=None, bcr_or_tcr='tcr')
