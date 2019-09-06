run_models = function(input_data_path){
	t = read.table(input_data_path, header=TRUE)
	attach(t)
	t$patient_id = as.factor(t$patient_id)
	t$weeks_since_start_bin = as.factor(t$weeks_since_start_bin)
	t$treatment = relevel(t$treatment, ref='Untreated')#make 'Untreated' the reference level
	sorted_tpoint_levels = sort((as.numeric(levels(t$weeks_since_start_bin))))
	t$weeks_since_start_bin = factor(t$weeks_since_start_bin, levels=sorted_tpoint_levels, ordered=TRUE)#make it an ordered factor
	tpoint_rank = c()#make a new col that has the rank of each time-point (as a numeric vector)
	for (i in t$weeks_since_start_bin){
		if (i == 0){tpoint_rank = c(tpoint_rank, 0)}
		else if (i == 12){tpoint_rank = c(tpoint_rank, 2)}
		else if (i == 36){tpoint_rank = c(tpoint_rank, 3)}
		else if (i == 60){tpoint_rank = c(tpoint_rank, 4)}
		else if (i == 84){tpoint_rank = c(tpoint_rank, 5)}
		else if (i == 96){tpoint_rank = c(tpoint_rank, 6)}
		else if (i == 100){tpoint_rank = c(tpoint_rank, 7)}
		else if (i == 112){tpoint_rank = c(tpoint_rank, 8)}
		else if (i == 122){tpoint_rank = c(tpoint_rank, 9)}
		else if (i == 146){tpoint_rank = c(tpoint_rank, 10)}
	}
	t = data.frame(t, tpoint_rank)
	#remove patient 46313
	t = subset(t, patient_id!=46313)

	## exploratory plots of data ##
	library(ggplot2)

	#effect of treatment on diversity
	p = ggplot(data=t, aes(x=treatment, y=shannon_diversity))
	p = p + geom_point(position=position_jitter(w=0.1, h=0)) + ggtitle('Effect of Treatment on Diversity') + xlab('Treatment Drug') + ylab('Shannon Diversity')
	pdf(paste(input_data_path, '_plotTreatmentVsDiversity.pdf', sep=''))
	plot(p)
	dev.off()
	#grouping Untreated and Rebif treatment categories
	p = ggplot(data=t, aes(x=ocrelizumab_or_not, y=shannon_diversity))
	p = p + geom_point(position=position_jitter(w=0.1, h=0)) + ggtitle('Effect of Ocrelizumab on Diversity') + xlab('Treatment Drug') + ylab('Shannon Diversity')
	pdf(paste(input_data_path, '_plotOcrelizumabVsDiversity.pdf', sep=''))
	plot(p)
	dev.off()
	#effect of time on diversity
	p = ggplot(data=t, aes(x=weeks_since_start, y=shannon_diversity))
	p = p + geom_point() + facet_grid(~ocrelizumab_or_not)
	p = p + ggtitle('Effect of Time on Diversity') + xlab('Weeks Since Start') + ylab('Shannon Diversity')
	pdf(paste(input_data_path, '_plotTimeVsDiversity.pdf', sep=''))
	plot(p)
	dev.off()

	## linear models of data ##

	#simplest linear model
	lin_model = lm(shannon_diversity ~ weeks_since_start + treatment, data=t, na.action=na.exclude)
	print ('#############################')
	print ('Simple lin model')
	a = anova(lin_model)
	print (a)
	print (summary(lin_model))
	print ('#############################')
	output_dir = paste(input_data_path, '_lin_model/', sep='')
	dir.create(output_dir, showWarnings=FALSE)
	plot.model_fit(lin_model, t, output_dir, model_fit_grouping='treatment')

	#linear model with interaction b/t treatment and time
	lin_model_inter = lm(shannon_diversity ~ weeks_since_start + treatment + weeks_since_start:treatment, data=t, na.action=na.exclude)
	output_dir = paste(input_data_path, '_lin_model_inter/', sep='')
	dir.create(output_dir, showWarnings=FALSE)
	plot.model_fit(lin_model_inter, t, output_dir, model_fit_grouping='treatment')

	#model testing basic lin model with lin model with interaction variable
	library(lmtest)
	a = lrtest(lin_model, lin_model_inter)
	print ('#############################')
	print ('Lin model vs lin model with time:treament interaction')
	print(a)
	aic = AIC(lin_model, lin_model_inter)
	print (aic)
	print ('#############################')

	#linear model with patient indicator variable
	lin_model_patient = lm(shannon_diversity ~ weeks_since_start + treatment + patient_id, data=t)
	print ('#############################')
	print ('lin model with patient covariate')
	a = anova(lin_model_patient)
	print (a)
	print (summary(lin_model_patient))
	print ('#############################')
	output_dir = paste(input_data_path, '_lin_model_patient/', sep='')
	dir.create(output_dir, showWarnings=FALSE)
	plot.model_fit(lin_model_patient, t, output_dir)

	#model testing basic lin model with one with patient indicator
	a = lrtest(lin_model, lin_model_patient)
	print ('#############################')
	print ('Lin model vs. lin model with patient covariate')
	print (a)
	aic = AIC(lin_model, lin_model_patient)
	print (aic)
	print ('#############################')

	#linear model with patient indicator and patient:time interaction
	lin_model_patientInter = lm(shannon_diversity ~ weeks_since_start + treatment + patient_id + patient_id:weeks_since_start, data=t)
	print ('#############################')
	print ('lin model with patient:time interaction')
	#using likelihood ratio test to test each categorical variable's sig
	null_mod = lm(shannon_diversity ~ weeks_since_start + patient_id + patient_id:weeks_since_start, data=t)
	a = lrtest(null_mod, lin_model_patientInter)
	print (a)
	null_mod = lm(shannon_diversity ~ weeks_since_start + treatment + patient_id:weeks_since_start, data=t)
	a = lrtest(null_mod, lin_model_patientInter)
	print (a)
	null_mod = lm(shannon_diversity ~ weeks_since_start + treatment + patient_id, data=t)
	a = lrtest(null_mod, lin_model_patientInter)
	print (a)

	a = anova(lin_model_patientInter)
	print (a)
	print (summary(lin_model_patientInter))
	print ('#############################')
	output_dir = paste(input_data_path, '_lin_model_patientInter/', sep='')
	dir.create(output_dir, showWarnings=FALSE)
	coefs_to_get = c('weeks_since_start', 'treatmentRebif', 'treatmentOcrelizumab', 'patient_id16912', 'patient_id20712', 'patient_id27812', 'patient_id30312', 'patient_id30712', 'patient_id32312', 'patient_id33012')
	effect_and_sig = get_effect_sig(lin_model_patientInter, coefs_to_get)
	factors_to_test = c('treatment', 'patient_id', 'patient_id:weeks_since_start')
	like_rat_test_pvals = like_rat_test_factors(lin_model_patientInter, factors_to_test, t)
	table_output_filepath = paste(output_dir, 'coef_effect_and_sig_table.txt')
	write_output_table(coefs_to_get, effect_and_sig[[1]], effect_and_sig[[2]], factors_to_test, like_rat_test_pvals, table_output_filepath)
	plot.model_fit(lin_model_patientInter, t, output_dir)

	#model testing lin model with patient info against the same but with patient:time interaction
	a = lrtest(lin_model_patient, lin_model_patientInter)
	print ('#############################')
	print ('Lin model with patient covar. vs the same thing but with patient:time interaction')
	print (a)
	aic = AIC(lin_model_patient, lin_model_patientInter)
	print (aic)
	print ('#############################')

	#linear model with heteroscedastic errors over time
	library(nlme)
	lin_model_heteroError = gls(shannon_diversity ~ weeks_since_start + treatment + patient_id + patient_id:weeks_since_start, weights=varIdent(form=~1|patient_id), data=t, na.action=na.exclude)
	output_dir = paste(input_data_path, '_lin_model_heteroError/', sep='')
	dir.create(output_dir, showWarnings=FALSE)
	plot.model_fit(lin_model_heteroError, t, output_dir)

	#testing basic lin model with one with hetero. error
	a = lrtest(lin_model_patientInter, lin_model_heteroError)
	print ('#############################')
	print ('lin model with patient:time inter. vs same but with hetero error')
	print (a)
	aic = AIC(lin_model_patientInter, lin_model_heteroError)
	print (aic)
	print ('#############################')

	#linear model with correlated errors over time
	t$weeks_since_start_bin = as.numeric(as.character(t$weeks_since_start_bin))#need to change this to a numeric vector to make the correlation structure (corCAR1) work
	lin_model_corrError = gls(shannon_diversity ~ weeks_since_start + treatment + patient_id + patient_id:weeks_since_start, correlation=corCAR1(form = ~ weeks_since_start_bin|patient_id), data=t, na.action=na.exclude, control=glsControl(singular.ok=TRUE))
	print ('#############################')
	print ('lin model with patinet:time inter. and corr. errors over time')
	a = anova(lin_model_corrError)
	print (a)
	print (summary(lin_model_corrError))
	print ('#############################')
	#print (corMatrix(lin_model_heteroError$modelStruct$corStruct))#this prints the correlation matrices for each individual
	output_dir = paste(input_data_path, '_lin_model_corrError/', sep='')
	dir.create(output_dir, showWarnings=FALSE)
	plot.model_fit(lin_model_corrError, t, output_dir)

	#model testing lin model to lin model with cor. errors
	a = lrtest(lin_model_patientInter, lin_model_corrError)
	print ('#############################')
	print ('lin model with patient:time inter. vs same but with corr. error')
	print (a)
	aic = AIC(lin_model_patientInter, lin_model_corrError)
	print (aic)
	print ('#############################')

	## linear mixed models of data ##

	#linear mixed model with homoscedastic errors
	t$weeks_since_start_bin = as.factor(t$weeks_since_start_bin)
	sorted_tpoint_levels = sort((as.numeric(levels(t$weeks_since_start_bin))))
	t$weeks_since_start_bin = factor(t$weeks_since_start_bin, levels=sorted_tpoint_levels, ordered=TRUE)#make it an ordered factor
	mixed_lin_model = lme(shannon_diversity ~ weeks_since_start + treatment, random=~1|patient_id, data=t, na.action=na.exclude)
	output_dir = paste(input_data_path, '_mixed_lin_model/', sep='')
	dir.create(output_dir, showWarnings=FALSE)
	plot.model_fit(mixed_lin_model, t, output_dir)

	#model testing lin model with patient:time interaction against mixed effects model where patient is a random effect
	a = lrtest(lin_model_patientInter, mixed_lin_model)
	print ('#############################')
	print ('with patient:time inter. vs mixed effects model where patient is random covariate')
	print (a)
	aic = AIC(lin_model_patientInter, mixed_lin_model)
	print (aic)
	print ('#############################')
}

plot.model_fit = function(model, data, output_dir, model_fit_grouping='patient_id'){
	#plot residuals
	df = data.frame(fitted_values=fitted(model), residuals=resid(model, type='pearson'))
	p = ggplot(data=df, aes(x=fitted_values, y=residuals))
	p = p + geom_point() + xlab('Fitted Values') + ylab('Standardized Residuals')
	pdf(paste(output_dir, 'plotLinModelResids.pdf', sep=''))
	plot(p)
	dev.off()
	#plot residuals with time
	data_trunc = data[complete.cases(data),]
	df = data.frame(residuals=resid(model, type='pearson'), weeks_since_start=data_trunc$weeks_since_start_bin)
	df$weeks_since_start = as.factor(df$weeks_since_start)
	p = ggplot(data=df, aes(x=weeks_since_start, y=residuals))
	p = p + geom_point(position=position_jitter(w=0.1, h=0))
	p = p + ggtitle('Heteroscedasticity With Time') + xlab('jitter(Weeks Since Start)') + ylab('Standardized Residuals')
	pdf(paste(output_dir, 'plotLinModelResidsWithTime.pdf', sep=''))
	plot(p)
	dev.off()
	#plot residuals by patient
	df = data.frame(residuals=resid(model, type='pearson'), patient=data_trunc$patient_id)
	df$patient = as.factor(df$patient)
	p = ggplot(data=df, aes(x=patient, y=residuals))
	p = p + geom_point(position=position_jitter(w=0.1, h=0))
	p = p + ggtitle('Heteroscedasticity by Patient') + xlab('jitter(Patient ID)') + ylab('Standardized Residuals')
	pdf(paste(output_dir, 'plotLinModelResidsByPatient.pdf', sep=''))
	plot(p)
	dev.off()
	#QQ plot
	pdf(paste(output_dir, 'plotQQnorm.pdf', sep=''))
	qqnorm(resid(model, type='pearson'), lwd=2)
	dev.off()
	#adjacent time-point residual correlation plot
	patient_id = data_trunc$patient_id[1]
	resid_for_pat = c()
	resids_prior_tpoint = c()
	resids_next_tpoint = c()
	for (i in 1:length(data_trunc$shannon_diversity)){
		pat_id_2 = data_trunc$patient_id[i]
		if (pat_id_2 == patient_id){
			resid_for_pat = c(resid_for_pat, resid(model, type='pearson')[i])
		}
		else if (pat_id_2 != patient_id){
			resids_prior_tpoint = c(resids_prior_tpoint, resid_for_pat[1:length(resid_for_pat)-1])
			resids_next_tpoint = c(resids_next_tpoint, resid_for_pat[2:length(resid_for_pat)])
			patient_id = pat_id_2
			resid_for_pat = c(resid(model, type='pearson')[i])
		}
	}
	resids_prior_tpoint = c(resids_prior_tpoint, resid_for_pat[1:length(resid_for_pat)-1])
	resids_next_tpoint = c(resids_next_tpoint, resid_for_pat[2:length(resid_for_pat)])
	c = cor.test(resids_prior_tpoint, resids_next_tpoint, type='pearson')
	df = data.frame(x=resids_prior_tpoint, y=resids_next_tpoint)
	p = ggplot(data=df, aes(x=x, y=y))
	p = p + geom_point() + xlab('Standardized Residuals') + ylab('Standardized Residuals in Next Time-Point') + ggtitle(paste('Pearson Correlation = ', c$estimate, '; p = ', c$p.value, sep=''))
	pdf(paste(output_dir, 'plotResidCorr.pdf', sep=''))
	plot(p)
	dev.off()
	#plot the fitted values
	df = data.frame(Patients=data_trunc$patient_id, fitted_values=fitted(model), time_point=data_trunc$weeks_since_start, time_point_bin=data_trunc$weeks_since_start_bin, observed_values=data_trunc$shannon_diversity, treatment=data_trunc$treatment, treatment_regimen=data_trunc$treatment_regimen)
	p = ggplot(data=df, aes(x=time_point_bin, y=observed_values))
	p = p + geom_point(aes(colour=treatment), size=3)
	if (model_fit_grouping == 'patient_id'){
		p = p + geom_line(aes(group=Patients, y=fitted_values, colour=treatment_regimen))
	}
	else if (model_fit_grouping == 'treatment'){
		p = p + geom_line(aes(group=treatment, y=fitted_values, colour=treatment_regimen))
	}
	p = p + geom_point(aes(y=fitted_values), colour='black', size=3)
	p = p + ggtitle('Fitted Values') + xlab('Weeks Since Start') + ylab('Shannon Diversity')
	pdf(paste(output_dir, 'plotFittedModel.pdf', sep=''))
	plot(p)
	dev.off()
}

get_effect_sig = function(lin_model, coefs_to_get){
	p_vals = c()
	effects = c()
	coefs = summary(lin_model)$coefficients
	for (i in coefs_to_get){
		row_index = match(i, rownames(coefs))
		p_vals = c(p_vals, coefs[row_index, 4])
		effects = c(effects, coefs[row_index, 1])
	}
	return (list(effects, p_vals))
}

like_rat_test_factors = function(lin_model, factors_to_test, dat){
	library(lmtest)
	p_vals = c()
	for (i in factors_to_test){
		small_form = formula(paste('. ~ . -', i, sep=''))
		small_model = update(lin_model, small_form, data=dat)
		a = lrtest(small_model, lin_model)
		p_vals = c(p_vals, a$'Pr(>Chisq)'[2])
	}
	return(p_vals)
}

write_output_table = function(coefs_to_get, coef_effects, coef_pvals, factors_to_test, like_rat_test_pvals, table_output_filepath){
	effects = c(coef_effects[1], '-', coef_effects[2:3], '-', coef_effects[4:10], '-')
	Chi_square = c('', 'x', '', '', 'x', rep('', 7), 'x')
	t = c('x', '', 'x', 'x', '', rep('x', 7), '')
	p = c(coef_pvals[1], like_rat_test_pvals[1], coef_pvals[2], coef_pvals[3], like_rat_test_pvals[2], coef_pvals[4:10], like_rat_test_pvals[3])
	df = data.frame(effects, Chi_square, t, p)
	rownames(df) = c('time', 'treatment', 'rebif', 'ocrelizumab', 'patient', '16912', '20712', '27812', '30312', '30712', '32312', '33012', 'patinet:time')
	write.table(df, table_output_filepath, col.names=NA, row.names=TRUE, quote=FALSE, sep='\t')
}

args = commandArgs(trailingOnly=TRUE)
run_models(args[1])
