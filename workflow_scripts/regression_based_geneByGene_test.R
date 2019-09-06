run_model = function(input_data){
	t = read.table(input_data, header=TRUE)
	t$patient_id = as.factor(t$patient_id)
	t$treatment = relevel(t$treatment, ref='Untreated')#make 'Untreated' the reference level
	lin_model = lm(gene_level ~ weeks_since_start + treatment + patient_id + patient_id:weeks_since_start, data=t)
	coefficients = coef(summary(lin_model))
	Ocrelizumab_pval = coefficients[3,4]
	Ocrelizumab_effect = coefficients[3,1]
	ocre_results = list(Ocrelizumab_pval, Ocrelizumab_effect)
	return (ocre_results)
}

plot_trajs = function(filein, gene_name, p_val){
	library(ggplot2)
	t = read.table(filein, header=TRUE)
	t$patient_id = as.factor(t$patient_id)
	p = ggplot(data=t, aes(x=weeks_since_start, y=gene_level))
	p = p + geom_line(aes(group=patient_id))
	p = p + geom_point(aes(colour=treatment), size=3)
	p = p + ggtitle(paste(gene_name, '; p value = ', p_val, sep='')) + xlab('Time Since Start (weeks)') + ylab('Proportion of Unique Clones')
	pdf(paste(filein, '_plot.pdf', sep=''))
	plot (p)
	dev.off()
	return
}

get_pvals = function(input_dirpath, make_plots=TRUE){

	####################
	sig_threshold = 0.05
	####################

	p_values = c()
	gene_names = c()
	effect_sizes = c()
	for (i in list.files(input_dirpath)){
		len_i = nchar(i)
		if ((substring(i, 1, 1) == '.') | (substring(i, len_i-3) == '.pdf') | (i == 'regression_results.txt')){
			next
		}
		gene_name = i
		gene_names = c(gene_names, gene_name)
		input_filepath = paste(input_dirpath, i, sep='')
		results = run_model(input_filepath)
		p_val = results[[1]]
		p_values = c(p_values, p_val)
		beta = results[[2]]
		effect_sizes = c(effect_sizes, beta)
	}
	p_values = p.adjust(p_values, method='BH')
	df = data.frame(gene_names, effect_sizes, p_values)
	df = df[order(df$p_values),]
	output_filepath = paste(input_dirpath, 'regression_results.txt', sep='')
	write.table(df, file=output_filepath, row.names=FALSE, quote=FALSE, sep='\t')
	if (make_plots == TRUE){
		for (i in 1:length(df$p_values)){
			if (is.nan(df$p_values[i])){
				next
			}
			if (df$p_values[i] <= sig_threshold){
				filein = paste(input_dirpath, df$gene_names[i], sep='')
				plot_trajs(filein, df$gene_names[i], df$p_values[i])
			}
		}
	}
}

cycle_thru_all_files = function(input_dirpath, simulation_data){
	#input_dirpath should be path to regression formated gene usage data
	#simulation_data is a boolean. If TRUE then it tells the script to cycle though file structure formatted for simulation data (i.e. there are 2 extra directory levels)
	for (i in list.files(input_dirpath)){
		len_i = nchar(i)
		if ((substring(i, 1, 1) == '.') | (substring(i, len_i-3) == '.pdf') | (substring(i, 1, 6) == 'README')){
			next
		}
		for (j in list.files(paste(input_dirpath, i, sep=''))){
			len_j = nchar(j)
			if ((substring(j, 1, 1) == '.') | (substring(j, len_j-3) == '.pdf') | (substring(j, 1, 6) == 'README')){
				next
			}
			for (k in list.files(paste(input_dirpath, i, '/', j, sep=''))){
				len_k = nchar(k)
				if ((substring(k, 1, 1) == '.') | (substring(k, len_k-3) == '.pdf') | (substring(k, 1, 6) == 'README')){
					next
				}
				for (l in list.files(paste(input_dirpath, i, '/', j, '/', k, sep=''))){
					len_l = nchar(l)
					if ((substring(l, 1, 1) == '.') | (substring(l, len_l-3) == '.pdf') | (substring(l, 1, 6) == 'README')){
						next
					}
					for (m in list.files(paste(input_dirpath, i, '/', j, '/', k, '/', l, sep=''))){
						len_m = nchar(m)
						if ((substring(m, 1, 1) == '.') | (substring(m, len_m-3) == '.pdf') | (substring(m, 1, 6) == 'README')){
							next
						}
						if (simulation_data == TRUE){
							for (n in list.files(paste(input_dirpath, i, '/', j, '/', k, '/', l, '/', m, sep=''))){
								len_n = nchar(n)
								if ((substring(n, 1, 1) == '.') | (substring(n, len_n-3) == '.pdf') | (substring(n, 1, 6) == 'README')){
									next
								}
								for (o in list.files(paste(input_dirpath, i, '/', j, '/', k, '/', l, '/', m, '/', n, sep=''))){
									len_o = nchar(o)
									if ((substring(o, 1, 1) == '.') | (substring(o, len_o-3) == '.pdf') | (substring(o, 1, 6) == 'README')){
										next
									}
									simulation_dirpath = paste(input_dirpath, i, '/', j, '/', k, '/', l, '/', m, '/', n, '/', o, '/', sep='')
									print (simulation_dirpath)
									get_pvals(simulation_dirpath, make_plots=FALSE)
								}
							}
						}
						else if (simulation_data == FALSE){
							geneClass_dirpath = paste(input_dirpath, i, '/', j, '/', k, '/', l, '/', m, '/', sep='')
							print (geneClass_dirpath)
							get_pvals(geneClass_dirpath)
						}
					}
				}
			}
		}
	}
}

args = commandArgs(trailingOnly=TRUE)
if (args[2] == 'TRUE'){
	simulation_data = TRUE
} else if (args[2] == 'FALSE'){
	simulation_data = FALSE
}
cycle_thru_all_files(args[1], simulation_data)
