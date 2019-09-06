import sys
import os
sys.path.insert(0, '../python_classes/')
from immune_time_series_class import immune_time_series
import subprocess

def get_data(input_data_dirpath, gene_class, include_counts, normalize_by_total=True, cell_count_data=None):
	"""
	cell_count_data - If None (default), then will not normalize by the estimated cell counts for the given data. If this is defined then it should be the path to the shannon diversity file that contains the approx cell counts for the OPERA data.
	normalize_by_total - If true then the counts for each gene are normalized by the total number of reads in the sample. This variable is disregarded if cell_count_data is defined.
	"""

	#if using cell_count_data then get that data from the file before proceeding
	if cell_count_data:
		cell_count_dic = {}
		filein = open(cell_count_data, "r")
		filein.readline()
		for i in filein:
			line = i[:-1].split('\t')
			pat_id = line[0]
			try:
				cell_count_dic[pat_id]
			except KeyError:
				cell_count_dic[pat_id] = {}
			tpoint = round(float(line[2]))
			cell_count = line[8]
			cell_count_dic[pat_id][tpoint] = cell_count
		filein.close()

	gene_dic = {}
	all_gene_names = set()
	for i in os.listdir(input_data_dirpath):
		if i[:6] == 'README' or i[0] == '.':
			continue
		patient = i

		####### skip some patients ############
		#patient 46313 has incomplete data
		if patient == '46313':
			continue
		#patient 39913 was treated with Gilenya
		if patient == '39913':
			continue
		#######################################

		gene_dic[patient] = {'treatments':[], 'tpoints':[], 'totals':[], 'cell_counts':[], 'gene_data':{}}
		time_ordered_files = sorted([[float(j[:-6].split("_")[1]), j] for j in os.listdir(input_data_dirpath + i) if j[0] != '.'])
		time_ordered_files = [j[1] for j in time_ordered_files]

		#need to remove the time-points that don't have cell counts
		if cell_count_data:
			time_ordered_files_new = []
			for j in time_ordered_files:
				samp_info = j[:-6].split('_')
				cell_count = cell_count_dic[patient][round(float(samp_info[1]))]
				if cell_count != 'NA':
					time_ordered_files_new.append(j)
			time_ordered_files = time_ordered_files_new[:]

		for tpoint_rank, j in enumerate(time_ordered_files):
			samp_info = j[:-6].split('_')
			if cell_count_data:
				cell_count = cell_count_dic[patient][round(float(samp_info[1]))]
				gene_dic[patient]['cell_counts'].append(cell_count)
			gene_dic[patient]['treatments'].append(samp_info[0])
			gene_dic[patient]['tpoints'].append(samp_info[1])
			gene_dic[patient]['totals'].append(0)
			filein = open('%s%s/%s' % (input_data_dirpath, patient, j), "r")
			for k in filein:
				if k[0] == '>':
					clone_data = k[1:-1].split('|')
					if gene_class == 'v_gene':
						gene_name = clone_data[1].split('=')[1]
					elif gene_class == 'j_gene':
						gene_name = clone_data[2].split('=')[1]
					all_gene_names.update([gene_name])
					if include_counts == True:
						count = int(clone_data[4].split('=')[1])
					elif include_counts == False:
						count = 1
					gene_dic[patient]['totals'][-1] += count
					try:
						gene_dic[patient]['gene_data'][gene_name][tpoint_rank] += count
					except KeyError:
						gene_dic[patient]['gene_data'][gene_name] = [0 for l in xrange(len(time_ordered_files))]
						gene_dic[patient]['gene_data'][gene_name][tpoint_rank] = count
			filein.close()

	#add in 0's for the genes that were not observed in a given patient
	for i in gene_dic:
		patient = i
		pat_set = set(gene_dic[patient]['gene_data'])
		for j in all_gene_names-pat_set:#for each gene that is missing from the patient, but was observed in others
			gene_dic[patient]['gene_data'][j] = [0 for k in xrange(len(gene_dic[patient]['tpoints']))]
	if cell_count_data:
		for i in gene_dic:
			for j in gene_dic[i]['gene_data']:
				for k in xrange(len(gene_dic[i]['tpoints'])):
					#this should give the number of B cells that had a given gene
					gene_dic[i]['gene_data'][j][k] = (gene_dic[i]['gene_data'][j][k] / float(gene_dic[i]['totals'][k])) * float(gene_dic[i]['cell_counts'][k])
	elif normalize_by_total == True:
		for i in gene_dic:
			for j in gene_dic[i]['gene_data']:
				for k in xrange(len(gene_dic[i]['tpoints'])):
					gene_dic[i]['gene_data'][j][k] = gene_dic[i]['gene_data'][j][k] / float(gene_dic[i]['totals'][k])
	return gene_dic, all_gene_names

def get_data_simulate(input_data_dirpath, gene_class, include_counts, cell_count_data, normalize_by_total=True):
	"""Runs simulations, and get data from each one."""
	
	### parameters ###
	element_to_simulate = 'clone'
	trials = 100
	# trials = 1
	update_freqs_with_serial_samp = True
	clone_def = 'v_j_cdr3'
	use_comp_cluster = True
	# use_comp_cluster = False
	#temp_dirpath = '/hernandez/mandrill/users/nstrauli/data/ms/temp_stuff'
	temp_dirpath = '/Users/nstrauli/data/ms/temp_stuff'
	times_per_generation = [0.001, 0.01, 0.05, 0.1, 0.5, 1., 2., 5., 10., None]
	#times_per_generation = [None,0.5,1.,2.,5.,10.]
	# times_per_generation = [0.1]
	##################

	#get cell count data (i.e. population size for each time-point)
	cell_count_dic = {}
	filein = open(cell_count_data, "r")
	filein.readline()
	for i in filein:
		line = i[:-1].split('\t')
		pat_id = line[0]
		try:
			cell_count_dic[pat_id]
		except KeyError:
			cell_count_dic[pat_id] = {}
		tpoint = round(float(line[2]))
		cell_count = line[8]
		cell_count_dic[pat_id][tpoint] = cell_count
	filein.close()

	gene_dics = {}
	for i in times_per_generation:
		gene_dics[str(i)] = [{} for j in xrange(trials)]#this will be a list of gene_dic's
	all_gene_names = set()
	#for each patient
	for i in os.listdir(input_data_dirpath):
		if i[:6] == 'README' or i[0] == '.':
			continue
		patient = i

		####### skip some patients ############
		#patient 46313 has incomplete data
		if patient == '46313':
			continue
		#patient 39913 was treated with Gilenya
		if patient == '39913':
			continue
		#######################################

		pop_counts = [cell_count_dic[patient][j] for j in sorted(cell_count_dic[patient])]
		os.chdir('../abr_hiv_coevo')
		samples = immune_time_series(input_data_dirpath + i, count_attribute_name='count', vgene_name='V_gene', dgene_name=None, jgene_name='J_gene', cdr3_name='CDR3', ignor_allele_info=True)

		#need to remove the time-points that don't have population count data
		tpoints_to_remove = [j for j in xrange(len(pop_counts)) if pop_counts[j]=='NA']
		samples.remove_time_points(tpoints_to_remove)
		time_ordered_files = sorted([[float(j[:-6].split("_")[1]), j] for j in os.listdir(input_data_dirpath + i) if j[0] != '.'])
		time_ordered_files = [j[1] for j in time_ordered_files]
		time_ordered_files_new = [time_ordered_files[j] for j in xrange(len(time_ordered_files)) if pop_counts[j]!='NA']
		time_ordered_files = time_ordered_files_new[:]
		pop_counts_new = [j for j in pop_counts if j!='NA']
		pop_counts = pop_counts_new[:]

		#get info for each patient
		cell_counts = []
		treatments = []
		tpoints = []
		for j in time_ordered_files:
			samp_info = j[:-6].split('_')
			cell_count = int(cell_count_dic[patient][round(float(samp_info[1]))])
			cell_counts.append(cell_count)
			treatments.append(samp_info[0])
			tpoints.append(samp_info[1])

		print 'simulating data for:', input_data_dirpath + i
		element_names, sim_element_counts_record = samples.simulate_time_points(element_to_simulate=element_to_simulate, trials=trials, serial_samp_pop_counts=cell_counts, update_freqs_with_serial_samp=update_freqs_with_serial_samp, output_dirpath=None, clone_def=clone_def, use_comp_cluster=use_comp_cluster, temp_dirpath=temp_dirpath, time_per_generation=times_per_generation, store_full_data_in_ram=False)
		os.chdir('../ms')

		#foreach generation time parameter value
		for gen_time_index, gen_time in enumerate(gene_dics):
			#foreach independent trial
			for j in xrange(len(sim_element_counts_record[gen_time_index])):
				gene_dics[gen_time][j][patient] = {'treatments':treatments, 'tpoints':tpoints, 'cell_counts':cell_counts, 'gene_data':{}}
				sim_data_filepath = sim_element_counts_record[gen_time_index][j]
				filein = open(sim_data_filepath, "r")
				one_sim = [[int(l) for l in k[:-1].split(',')] for k in filein]
				#for each time-point
				for k in xrange(len(one_sim)):
					total = 0.
					#for each clone
					for l in xrange(len(one_sim[k])):
						if gene_class == 'v_gene':
							gene_name = element_names[l].split('_')[0]
						elif gene_class == 'j_gene':
							gene_name = element_names[l].split('_')[1]
						all_gene_names.update([gene_name])
						if include_counts:
							count = one_sim[k][l]
						else:
							if one_sim[k][l] == 0:
								count = 0
							else:
								count = 1
						total += count
						try:
							gene_dics[gen_time][j][patient]['gene_data'][gene_name][k] += count
						except KeyError:
							gene_dics[gen_time][j][patient]['gene_data'][gene_name] = [0 for m in tpoints]
							gene_dics[gen_time][j][patient]['gene_data'][gene_name][k] += count
					#normalize by total, if desired
					if normalize_by_total:
						#for each gene name
						for l in gene_dics[gen_time][j][patient]['gene_data']:
							gene_dics[gen_time][j][patient]['gene_data'][l][k] /= total
				subprocess.call(['rm', sim_data_filepath])

	#add in 0's for the genes that were not observed in a given patient
	for gen_time in gene_dics:
		for i in xrange(len(gene_dics[gen_time])):
			for j in gene_dics[gen_time][i]:
				patient = j
				pat_set = set(gene_dics[gen_time][i][patient]['gene_data'])
				for k in all_gene_names-pat_set:#for each gene that is missing from the patient, but was observed in others
					gene_dics[gen_time][i][patient]['gene_data'][k] = [0 for l in xrange(len(gene_dics[gen_time][i][patient]['tpoints']))]

	return gene_dics, all_gene_names

def write_data(gene_dic, all_gene_names, output_dirpath, cell_type, study, isotype, gene_class, count_status, sort_data=True):
	#write data for regression:
	out_dir_regression = '%sregression/%s/%s/%s/%s/%s/' % (output_dirpath, count_status, study, cell_type, isotype, gene_class)
	if not os.path.exists(out_dir_regression):
		os.makedirs(out_dir_regression)
	for i in all_gene_names:
		gene_name = i
		output_filepath = out_dir_regression + gene_name
		fileout = open(output_filepath, "w")
		fileout.write('patient_id\tweeks_since_start\ttreatment\tgene_level\n')
		for j in gene_dic:
			patient = j
			for tpoint_rank, k in enumerate(gene_dic[patient]['gene_data'][gene_name]):
				fileout.write('%s\t%s\t%s\t%s\n' % (patient, gene_dic[patient]['tpoints'][tpoint_rank], gene_dic[patient]['treatments'][tpoint_rank], k))
		fileout.close()
	#write data for time-series plots:
	out_dir_tseries = '%stime_series_viz/%s/%s/%s/%s/%s/' % (output_dirpath, count_status, study, cell_type, isotype, gene_class)
	if not os.path.exists(out_dir_tseries):
		os.makedirs(out_dir_tseries)
	for i in gene_dic:
		patient = i
		output_filepath = out_dir_tseries + patient
		fileout = open(output_filepath, "w")
		fileout.write('\t%s\n' % '\t'.join([str(j) for j in gene_dic[patient]['tpoints']]))
		if sort_data == True:
			sorted_gene_data = []
			for j in gene_dic[patient]['gene_data']:
				total = sum(gene_dic[patient]['gene_data'][j])
				sorted_gene_data.append([total, j] + gene_dic[patient]['gene_data'][j])
			for j in sorted(sorted_gene_data):
				fileout.write('%s\t%s\n' % (j[1], '\t'.join([str(k) for k in j[2:]])))
		else:
			for j in gene_dic[patient]['gene_data']:
				fileout.write('%s\t%s\n' % (j, '\t'.join([str(k) for k in gene_dic[patient]['gene_data'][j]])))
		fileout.close()
	return

def write_data_sim(gene_dics, all_gene_names, output_dirpath, cell_type, study, isotype, gene_class, count_status, sort_data=True):
	
	#write data formatted for the regression analysis
	for gen_time in gene_dics:
		for trial in xrange(len(gene_dics[gen_time])):
			out_dir_regression = '%sregression/%s/%s/%s/%s/%s/%s/%s/' % (output_dirpath, count_status, study, cell_type, isotype, gene_class, gen_time, trial)
			if not os.path.exists(out_dir_regression):
				os.makedirs(out_dir_regression)
			gene_dic = gene_dics[gen_time][trial]
			for i in all_gene_names:
				gene_name = i
				output_filepath = out_dir_regression + gene_name
				fileout = open(output_filepath, "w")
				fileout.write('patient_id\tweeks_since_start\ttreatment\tgene_level\n')
				for j in gene_dic:
					patient = j
					for tpoint_rank, k in enumerate(gene_dic[patient]['gene_data'][gene_name]):
						fileout.write('%s\t%s\t%s\t%s\n' % (patient, gene_dic[patient]['tpoints'][tpoint_rank], gene_dic[patient]['treatments'][tpoint_rank], k))
				fileout.close()
	
	#write data formatted for time-series visualization
	for gen_time in gene_dics:#for each generation time parameter value
		#for each patient
		for i in gene_dics[gen_time][0]:
			patient = i
			out_dir_tseries = '%stime_series_viz/%s/%s/%s/%s/%s/%s/%s/' % (output_dirpath, count_status, study, cell_type, isotype, gene_class, patient, gen_time)
			if not os.path.exists(out_dir_tseries):
				os.makedirs(out_dir_tseries)
			#for each sim trial
			for trial in xrange(len(gene_dics[gen_time])):
				gene_dic = gene_dics[gen_time][trial]
				output_filepath = out_dir_tseries + str(trial)
				fileout = open(output_filepath, "w")
				fileout.write('\t%s\n' % '\t'.join([str(j) for j in gene_dic[patient]['tpoints']]))
				if sort_data == True:
					sorted_gene_data = []
					for j in gene_dic[patient]['gene_data']:
						total = sum(gene_dic[patient]['gene_data'][j])
						sorted_gene_data.append([total, j] + gene_dic[patient]['gene_data'][j])
					for j in sorted(sorted_gene_data):
						fileout.write('%s\t%s\n' % (j[1], '\t'.join([str(k) for k in j[2:]])))
				else:
					for j in gene_dic[patient]['gene_data']:
						fileout.write('%s\t%s\n' % (j, '\t'.join([str(k) for k in gene_dic[patient]['gene_data'][j]])))
				fileout.close()

	return

def get_data_foreach_dataset(clone_seqs_dirpath, cell_count_data_dirpath, output_dirpath_observed, output_dirpath_simulation, simulate=True):
	"""'cell_count_data_dirpath' - This should be the 'shannon_diversity' directory. The data files in this directory contain cell count information (i.e. number B/T cells)."""
	if clone_seqs_dirpath[-1] != '/':
		clone_seqs_dirpath += '/'
	if cell_count_data_dirpath[-1] != '/':
		cell_count_data_dirpath += '/'
	if output_dirpath_observed[-1] != '/':
		output_dirpath_observed += '/'
	if not os.path.exists(output_dirpath_observed):
		os.makedirs(output_dirpath_observed)
	if output_dirpath_simulation[-1] != '/':
		output_dirpath_simulation += '/'
	if not os.path.exists(output_dirpath_simulation):
		os.makedirs(output_dirpath_simulation)
	for i in os.listdir(clone_seqs_dirpath):
		if i[:6] == 'README' or i[0] == '.':
			continue
		cell_type = i
		for j in os.listdir('%s%s' % (clone_seqs_dirpath, cell_type)):
			if j[:6] == 'README' or j[0] == '.':
				continue
			study = j
			for k in os.listdir('%s%s/%s' % (clone_seqs_dirpath, cell_type, study)):
				if k[:6] == 'README' or k[0] == '.':
					continue
				isotype = k
				print cell_type, study, isotype

				# if cell_type == 'bcr' or study == 'oratorio' or isotype == 'IgG':
				# 	continue

				#get the cell count file path right quick
				input_data_dirpath = '%s%s/%s/%s/' % (clone_seqs_dirpath, cell_type, study, isotype)
				#get data including clonal counts
				gene_dic, all_gene_names = get_data(input_data_dirpath=input_data_dirpath, gene_class='v_gene', include_counts=True)
				write_data(gene_dic, all_gene_names, output_dirpath_observed, cell_type, study, isotype, 'v_gene', 'using_clone_counts')
				gene_dic, all_gene_names = get_data(input_data_dirpath=input_data_dirpath, gene_class='j_gene', include_counts=True)
				write_data(gene_dic, all_gene_names, output_dirpath_observed, cell_type, study, isotype, 'j_gene', 'using_clone_counts')
				#get data excluding clone counts (i.e. clones are presence/absense)
				gene_dic, all_gene_names = get_data(input_data_dirpath=input_data_dirpath, gene_class='v_gene', include_counts=False)
				write_data(gene_dic, all_gene_names, output_dirpath_observed, cell_type, study, isotype, 'v_gene', 'unique_clones_only')
				gene_dic, all_gene_names = get_data(input_data_dirpath=input_data_dirpath, gene_class='j_gene', include_counts=False)
				write_data(gene_dic, all_gene_names, output_dirpath_observed, cell_type, study, isotype, 'j_gene', 'unique_clones_only')
				#if study == OPERA get the cell count data file path (don't currently have cell count data for ORATORIO)
				if study == 'opera':
					isotype_dic = {'IgG':'igg', 'IgM':'igm', 'CD4':'cd4', 'CD8':'cd8'}
					cell_count_filepath = '%s%s/%s.txt' % (cell_count_data_dirpath, cell_type, isotype_dic[isotype])
					gene_dic, all_gene_names = get_data(input_data_dirpath=input_data_dirpath, gene_class='v_gene', include_counts=True, cell_count_data=cell_count_filepath)
					write_data(gene_dic, all_gene_names, output_dirpath_observed, cell_type, study, isotype, 'v_gene', 'using_cell_counts')
					gene_dic, all_gene_names = get_data(input_data_dirpath=input_data_dirpath, gene_class='j_gene', include_counts=True, cell_count_data=cell_count_filepath)
					write_data(gene_dic, all_gene_names, output_dirpath_observed, cell_type, study, isotype, 'j_gene', 'using_cell_counts')
					if simulate:
						#simulate data including clonal counts
						gene_dics, all_gene_names = get_data_simulate(input_data_dirpath=input_data_dirpath, gene_class='v_gene', include_counts=True, cell_count_data=cell_count_filepath)
						write_data_sim(gene_dics=gene_dics, all_gene_names=all_gene_names, output_dirpath=output_dirpath_simulation, cell_type=cell_type, study=study, isotype=isotype, gene_class='v_gene', count_status='using_clone_counts')
						gene_dics, all_gene_names = get_data_simulate(input_data_dirpath=input_data_dirpath, gene_class='j_gene', include_counts=True, cell_count_data=cell_count_filepath)
						write_data_sim(gene_dics=gene_dics, all_gene_names=all_gene_names, output_dirpath=output_dirpath_simulation, cell_type=cell_type, study=study, isotype=isotype, gene_class='j_gene', count_status='using_clone_counts')
						#simulate data excluding clone counts (i.e. clones are presence/absense)
						gene_dics, all_gene_names = get_data_simulate(input_data_dirpath=input_data_dirpath, gene_class='v_gene', include_counts=False, cell_count_data=cell_count_filepath)
						write_data_sim(gene_dics=gene_dics, all_gene_names=all_gene_names, output_dirpath=output_dirpath_simulation, cell_type=cell_type, study=study, isotype=isotype, gene_class='v_gene', count_status='unique_clones_only')
						gene_dics, all_gene_names = get_data_simulate(input_data_dirpath=input_data_dirpath, gene_class='j_gene', include_counts=False, cell_count_data=cell_count_filepath)
						write_data_sim(gene_dics=gene_dics, all_gene_names=all_gene_names, output_dirpath=output_dirpath_simulation, cell_type=cell_type, study=study, isotype=isotype, gene_class='j_gene', count_status='unique_clones_only')
	return

if __name__ == '__main__':
	get_data_foreach_dataset(clone_seqs_dirpath='/hernandez/mandrill/users/nstrauli/data/ms/clone_seqs_fasta_format_clustered_data', cell_count_data_dirpath='/hernandez/mandrill/users/nstrauli/data/ms/clonal_diversity/shannon_diversity', output_dirpath_observed='/hernandez/mandrill/users/nstrauli/data/ms/gene_usage_clustered_data/observed_values', output_dirpath_simulation='/hernandez/mandrill/users/nstrauli/data/ms/gene_usage_clustered_data/simulations', simulate=True)
