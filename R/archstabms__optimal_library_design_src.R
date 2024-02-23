
#' archstabms__optimal_library_design_src
#'
#' Library design for Src optimising WT-like abundance and activity.
#'
#' @param mochi_outpath path to ddPCA data (required)
#' @param position_offset residue position offset (default:0)
#' @param max_num_mutants maximum number of mutations (required)
#' @param explore_mutations vector of mutation ids to explore or "all_starting" (default:"all_starting")
#' @param explore_positions vector of positions to explore or "all_starting" (default:"all_starting")
#' @param RT constant (default:0.001987*(273+30))
#' @param sample_size number of random samples (default:1000) 
#' @param optimal_perc optimal percentage inverse minimum geometric mean growth rate deviation from WT (default:0.7)
#' @param min_num_doubles minimum number of doubles per single for candidate mutations (default:20)
#' @param colour_scheme colour scheme file (required)
#' @param outpath output path for plots and saved objects (required)
#'
#' @return Nothing
#' @export
#' @import data.table
archstabms__optimal_library_design_src <- function(
	mochi_outpath,
	position_offset = 0,
	max_num_mutants, 
	explore_mutations = "all_starting", 
	explore_positions = "all_starting",
	RT = 0.001987*(273+30), 
	sample_size = 1000,
	optimal_perc = 0.7,
	min_num_doubles = 20,
	colour_scheme,
	outpath
	){

  #Domain name
  domain_name <- rev(unlist(strsplit(basename(outpath), "_")))[1]

	#Display status
	message(paste("\n\n*******", paste("running stage: archstabms__optimal_library_design_src for", domain_name), "*******\n\n"))

	#Create output directory
	archstabms__create_dir(archstabms_dir = outpath)

	#Fitness
	all_variants <- fread(file.path(mochi_outpath, "task_1", "predictions", "predicted_phenotypes_all.txt"))

  #Fitness - abundance
  gra_lm <- all_variants[phenotype==1,lm(growthrate~fitness)]
  wta_gr <- all_variants[phenotype==1 & WT==T,growthrate] #WT growth rate

  #Fitness - activity
  grb_lm <- all_variants[phenotype==6,lm(growthrate~fitness)]
  wtb_gr <- all_variants[phenotype==6 & WT==T,growthrate] #WT growth rate

	#WT nucleotide sequence
	nt_wtseq <- fread(file.path(mochi_outpath, "WT.txt"), header = F)[,V1]

  #Growth rate thresholds
  modpar_a <- fread(file.path(mochi_outpath, "task_1", "weights", "linears_weights_Abundance_b1.txt"))
  modpar_b <- fread(file.path(mochi_outpath, "task_1", "weights", "linears_weights_Activity_b1.txt"))
  gr_thresholds <- list(
    "Abundance" = as.numeric(predict(gra_lm, newdata = data.frame(fitness = 0.5*modpar_a[fold==1,kernel]+modpar_a[fold==1,bias]))),
    "Activity" = as.numeric(predict(grb_lm, newdata = data.frame(fitness = 0.5*modpar_b[fold==1,kernel]+modpar_b[fold==1,bias]))))
		
	#Initialise variables
	d_ddg_pred_sumsq <- 0
	selected_mutants_strong <- c()
	comb_dt_strong <- data.table()
	selected_mutants <- c()
	folding_median_values <- data.table(number_of_mutants = 2:max_num_mutants)
	activity_median_values <- data.table(number_of_mutants = 2:max_num_mutants)

	#Determine candidate mutations (confident activity and folding ddGs, reachable by single nt substitutions, at least 20 doubles with single)
	#WT nucleotide sequence
	nt_wtseq <- fread(file.path(mochi_outpath, "WT.txt"), header = F)[,V1]
	#Model results
	ddg_dt <- merge(
		fread(file.path(mochi_outpath, "task_1", "weights", "weights_Folding.txt"))[,.(id_ref, Pos, Pos_ref, f_ddg_pred = .SD[[1]], f_ddg_pred_conf = .SD[[2]]<1),,.SDcols = c("mean_kcal/mol", "ci95_kcal/mol")],
		fread(file.path(mochi_outpath, "task_1", "weights", "weights_Activity.txt"))[,.(id_ref, Pos, Pos_ref, b_ddg_pred = .SD[[1]], b_ddg_pred_conf = .SD[[2]]<1),,.SDcols = c("mean_kcal/mol", "ci95_kcal/mol")], by = c("id_ref", "Pos", "Pos_ref"))
	#Add dg
	ddg_dt[, f_dg_pred := f_ddg_pred]
	ddg_dt[id_ref != 'WT', f_dg_pred := f_ddg_pred + ddg_dt[id_ref == 'WT',f_ddg_pred]]
	ddg_dt[, b_dg_pred := b_ddg_pred]
	ddg_dt[id_ref != 'WT', b_dg_pred := b_ddg_pred + ddg_dt[id_ref == 'WT',b_ddg_pred]]

	#WT amino acid sequence
	aa_wtseq_split <- unlist(strsplit(all_variants[WT==T,aa_seq][1], ""))
	#Table of singles in doubles
	pred_dt <- all_variants[Nham_aa==2 & phenotype %in% 1:5]
	pred_dt[, Pos_ref1 := which(unlist(strsplit(aa_seq, ""))!=aa_wtseq_split)[1],aa_seq]
	pred_dt[, Pos_ref2 := which(unlist(strsplit(aa_seq, ""))!=aa_wtseq_split)[2],aa_seq]
	pred_dt[, WT_AA1 := aa_wtseq_split[Pos_ref1],aa_seq]
	pred_dt[, WT_AA2 := aa_wtseq_split[Pos_ref2],aa_seq]
	pred_dt[, Mut1 := substr(aa_seq, Pos_ref1, Pos_ref1),aa_seq]
	pred_dt[, Mut2 := substr(aa_seq, Pos_ref2, Pos_ref2),aa_seq]
	pred_dt[, id_ref1 := paste0(WT_AA1, Pos_ref1, Mut1),aa_seq]
	pred_dt[, id_ref2 := paste0(WT_AA2, Pos_ref2, Mut2),aa_seq]
	num_tab <- table(unlist(pred_dt[,.(id_ref1, id_ref2)]))

	#Number of doubles per single
	ddg_dt[id_ref!="WT", num_doubles := 0]
	ddg_dt[id_ref!="WT" & id_ref %in% names(num_tab), num_doubles := num_tab[id_ref]]

	#Explore positions
	if(explore_positions[1]=='all_starting'){
		explore_positions <- unique(ddg_dt[id!="WT",Pos_ref])
	}

	#Reachable mutants
	reachable_singles <- archstabms__reachable_single_AA_mutants(nt_wtseq, position_offset)

	#Candidate windows
	window_size <- 22
	window_results <- list()
	for(i in 1:(length(aa_wtseq_split)-window_size+1)){
		candidate_muts <- ddg_dt[f_ddg_pred_conf==T & b_ddg_pred_conf==T & Pos_ref %in% i:(i+window_size-1) & id_ref %in% reachable_singles[["reachable_muts"]] & num_doubles>=min_num_doubles,]
		window_results[[i]] <- data.table(
			"start_pos" = i,
			"total_muts" = candidate_muts[,.N],
			"total_pos" = candidate_muts[,length(unique(Pos))])
	}
	window_results_dt <- rbindlist(window_results)
	save(window_results_dt, file = file.path(outpath, "candidate_windows.Rdata"))

	#Plot Candidate windows
	pos_offset <- 264
	plot_dt <- copy(window_results_dt)
	plot_dt[, start_pos_ref := start_pos+pos_offset]
	d <- ggplot2::ggplot(plot_dt, ggplot2::aes(start_pos_ref, total_muts)) +
		ggplot2::geom_line(linetype = 1) + 
		ggplot2::geom_line(ggplot2::aes(start_pos_ref, total_pos), color = colour_scheme[["shade 0"]][[1]]) + 
		ggplot2::theme_classic() +
		ggplot2::geom_point(size = 1) +
		ggplot2::geom_vline(xintercept = explore_positions[1]+pos_offset, linetype = 1) +
		ggplot2::geom_hline(yintercept = 15, linetype = 1) +
		ggplot2::geom_hline(yintercept = 22, linetype = 2) +
		ggplot2::ylab("Number of candidate mutations") +
		ggplot2::xlab("Window start residue")
	ggplot2::ggsave(file.path(outpath, "candidate_windows.pdf"), d, width = 7, height = 3, useDingbats=FALSE)

	#Candidate mutations
	ddgs_dt <- ddg_dt[f_ddg_pred_conf==T & b_ddg_pred_conf==T & Pos_ref %in% explore_positions & id_ref %in% reachable_singles[["reachable_muts"]] & num_doubles>=min_num_doubles,]
		
	#Set starting mutations
	mutations_id <- c()
	if(explore_mutations=="all_starting"){
		mutations_id <- ddgs_dt[,id_ref]
	}else{
		mutations_id <- explore_mutations
	}

	#Iterate over all starting mutations
	set.seed(123456)
	for(mutation in mutations_id){
		print(paste0("mutation ", which(mutations_id==mutation), "/",length(mutations_id)))
		selected_id <- c(mutation)
		median_values_activity <- c()
		median_values_folding <- c()

		#Folding and activity free energy combinatorial data.tables
		#WT
		combf_dt <- data.table(WT = rep(ddg_dt[id_ref=="WT",f_dg_pred], sample_size))
		combb_dt <- data.table(WT = rep(ddg_dt[id_ref=="WT",b_dg_pred], sample_size))
		#First mutant
		combf_dt[, paste0(mutation) := sample(0:1, sample_size, replace = TRUE)]
		combb_dt[, paste0(mutation) := combf_dt[,.SD[[1]]*ddgs_dt[id_ref==selected_id[1],b_ddg_pred],,.SDcols = mutation]]
		combf_dt[, paste0(mutation) := .SD[[1]]*ddgs_dt[id_ref==selected_id[1],f_ddg_pred],,.SDcols = mutation]

		#Greedily add new mutations up to max_num_mutants
		for(n_mut in 2:max_num_mutants){
			#Already selected mutations
			selected_dt<-ddgs_dt[id_ref %in% selected_id]
			#Free energies
			sum_ddg_folding<-sum(selected_dt[,f_ddg_pred])+ddg_dt[id_ref=="WT",f_dg_pred]
			sum_ddg_activity<-sum(selected_dt[,b_ddg_pred])+ddg_dt[id_ref=="WT",b_dg_pred]

			#Remaining mutations in other positions
			new_possible_mut<-ddgs_dt[!Pos %in% selected_dt[,Pos]]

			#No mutations left to select - break
			if(new_possible_mut[,.N]==0){
				break
			}
			
			#Compute abundance and activity growth rates of all current order +1 mutants
			#Free energies
			new_possible_mut[, f_dg_pred_comb := sum_ddg_folding+f_ddg_pred]
			new_possible_mut[, b_dg_pred_comb := sum_ddg_activity+b_ddg_pred]
			#Fitness
			new_possible_mut[, fitness_folding := archstabms__predict_fitness(mochi_outpath, f_dg_pred_comb, 0, RT = RT, phenotype_names = c("Abundance", "Activity"), dataset_name = "b1")[["fitness_folding"]]] #fitness_binding
			new_possible_mut[, fitness_activity := archstabms__predict_fitness(mochi_outpath, f_dg_pred_comb, b_dg_pred_comb, RT = RT, phenotype_names = c("Abundance", "Activity"), dataset_name = "b1")[["fitness_binding"]]] #fitness_binding
			#Growth rates deviation from WT
			new_possible_mut[, growthrate_folding := abs(wta_gr - gra_lm[["coefficients"]][1] + gra_lm[["coefficients"]][2]*fitness_folding)]
			new_possible_mut[, growthrate_activity := abs(wtb_gr - grb_lm[["coefficients"]][1] + grb_lm[["coefficients"]][2]*fitness_activity)]

			#Compute geometric means of folding and activity growth rates
			new_possible_mut[, geom_mean := (growthrate_folding*growthrate_activity)**(1/2)]

			#Select optimal mutation (minimum deviation from WT)
			selected_id<-c(selected_id, new_possible_mut[order(geom_mean)][1, id_ref])

			#Compute median abundance and activity growth rates of entire library
			#Randomly sample additional mutant
			combf_dt[, paste0(selected_id[n_mut]) := sample(rep(c(0, 1), each = sample_size/2))]
			combb_dt[, paste0(selected_id[n_mut]) := combf_dt[,.SD[[1]]*ddgs_dt[id_ref == selected_id[n_mut],b_ddg_pred],,.SDcols = selected_id[n_mut]]]
			combf_dt[, paste0(selected_id[n_mut]) := .SD[[1]]*ddgs_dt[id_ref == selected_id[n_mut],f_ddg_pred],,.SDcols = selected_id[n_mut]]
			#Compute abundance and activity growth rates
			gr_dt <- data.table()
			#Free energies
			gr_dt[, f_dg_pred := rowSums(combf_dt)]
			gr_dt[, b_dg_pred := rowSums(combb_dt)]
			#Fitness
			temp_fitness <- archstabms__predict_fitness(mochi_outpath, gr_dt[, f_dg_pred], gr_dt[, b_dg_pred], RT = RT, phenotype_names = c("Abundance", "Activity"), dataset_name = "b1")
			gr_dt[, a_fitness := temp_fitness[["fitness_folding"]]]
			gr_dt[, b_fitness := temp_fitness[["fitness_binding"]]]
			#Growth rates
			gr_dt[, a_growthrate := gra_lm[["coefficients"]][1] + gra_lm[["coefficients"]][2]*a_fitness]
			gr_dt[, b_growthrate := grb_lm[["coefficients"]][1] + grb_lm[["coefficients"]][2]*b_fitness]

			#Append the absolute deviation of the medians of the new libraries from the WT
			median_values_folding <- c(median_values_folding, abs(wta_gr - gr_dt[,median(a_growthrate)]))
			median_values_activity <- c(median_values_activity, abs(wtb_gr - gr_dt[,median(b_growthrate)]))
		}
		#Store
		folding_median_values[1:length(median_values_folding), paste0(mutation) := median_values_folding]
		activity_median_values[1:length(median_values_activity), paste0(mutation) := median_values_activity]
	}
	
	if(explore_mutations=="all_starting"){
		#Prepare table of all computed values and store
		total <- merge(
			melt(folding_median_values, measure.vars = mutations_id, variable.name = "mutation", value.name = "folding"),
			melt(activity_median_values, measure.vars = mutations_id, variable.name = "mutation", value.name = "folding"), 
			by=c("number_of_mutants","mutation"))
		names(total) <- c("mutant_order","mutation", "Abundance", "Activity")
		#Calculate geometric mean of deviation from WT
		total[, geom_mean := (Abundance*Activity)**(1/2)]
		save(total, file = file.path(outpath, "activity_folding_means_all_computed.Rdata"))
		
		#Select best geometric mean library for each order mutant and store
		final_table <- total[, .SD[which.min(geom_mean)], by = mutant_order]
		save(final_table, file = file.path(outpath, "activity_folding_BestGeomMean_bestEachStep.Rdata"))
		write.table(final_table, file = file.path(outpath, "activity_folding_BestGeomMean_bestEachStep.csv"), sep = "\t", row.names = F)
		
		#Select optimal library
		print("Optimal highest order mutant combination:")
		opt_table <- final_table[mutant_order==max_num_mutants][.N]
		print(opt_table)
		optimal_num_mutants <- opt_table[,mutant_order]
		optimal_first_mut <- opt_table[,mutation]
		
		#Plot - mutant order vs median growth rates progression
		table_plot <- melt(
			final_table, 
			measure.vars = c("Abundance", "Activity", "geom_mean"),
			variable.name = "assay", 
			value.name = "median_predicted_growth_rate")
		write.table(table_plot, file = file.path(outpath, "plot_table.csv"), sep = "\t", row.names = F)
		cc <- c(colour_scheme[["shade 0"]][[3]], colour_scheme[["shade 0"]][[1]], 'black')
		names(cc) <- c("Abundance", "Activity", "geom_mean")
		d <- ggplot2::ggplot(table_plot, ggplot2::aes(mutant_order, median_predicted_growth_rate, color = assay)) +
			ggplot2::geom_line(linetype = 1) + 
			ggplot2::theme_classic() +
			ggplot2::geom_point(size = 1) +
			ggplot2::geom_vline(xintercept = optimal_num_mutants, linetype = 2) +
			ggplot2::geom_hline(yintercept = final_table[,min(geom_mean)/optimal_perc], linetype = 2) +
			ggplot2::ylab("Predicted median growth rate dev.") +
			ggplot2::xlab("AA substitution order") +
			ggplot2::scale_colour_manual(values = cc)
		ggplot2::ggsave(file.path(outpath, "bestGeomMeanGR_perOrderMutants.pdf"), d, width = 7, height = 3, useDingbats=FALSE)
				
		#Store medians for that starting position, for each mutant order library
		return(list(
			optimal_first_mut = optimal_first_mut, 
			optimal_num_mutants = optimal_num_mutants, 
			folding_medians = folding_median_values, 
			activity_medians = activity_median_values))

	}else{
		#Specified starting mutation
		print("set of selected mutations:")
		print(selected_id)
		write.table(selected_id, file = file.path(outpath, "selected_ids.csv"), sep = "\t", row.names = F)

		#Mutant oligo
		nt_mutseq <- archstabms__get_codons(nt_wtseq) 
		for(i in selected_id){
			nt_mutseq[as.integer(substr(i, 2, nchar(i)-1))-position_offset] <- reachable_singles[['reachable_muts_codons']][[i]][1]
		}
		nt_mutseq <- paste0(nt_mutseq, collapse = "")
		#Oigo sequence with ambiguity codes
		oligo_seq <- paste0(sapply(sapply(apply(cbind(unlist(strsplit(toupper(nt_wtseq), "")), unlist(strsplit(toupper(nt_mutseq), ""))), 1, unique), paste0, collapse = ""), archstabms__get_ambiguity_code), collapse = "")
		write.table(oligo_seq, file = file.path(outpath, "oligo_seq.txt"), sep = "\n", row.names = F, col.names = F, quote = F)

	  #Random variant data.table
	  set.seed(123456)

	  folding_energies_dict <- as.list(ddg_dt[selected_id,f_ddg_pred])
	  names(folding_energies_dict) <- selected_id
	  activity_energies_dict <- as.list(ddg_dt[selected_id,b_ddg_pred])
	  names(activity_energies_dict) <- selected_id

	  #Set variant ids
	  rand_dt <- data.table(varid = 1:sample_size)
	  #Select random mutation at each step
	  rand_dt <- rand_dt[, mutlist := paste0(sample(selected_id, max_num_mutants), collapse = ","),varid]
	  #Energies
	  rand_dt[, mutcoeff_0 := ddg_dt[id_ref=="WT",f_dg_pred]]
	  rand_dt[, mutcoefb_0 := ddg_dt[id_ref=="WT",b_dg_pred]]
	  for(i in 1:max_num_mutants){
	    rand_dt[, paste0("mutcoeff_", i) := unlist(folding_energies_dict[as.character(sapply(strsplit(mutlist, ","), '[', i))])]
	    rand_dt[, paste0("mutcoefb_", i) := unlist(activity_energies_dict[as.character(sapply(strsplit(mutlist, ","), '[', i))])]
	  }
	  #Predict growth rate
	  for(i in 1:max_num_mutants){
	    rand_dt[, paste0("mutenergyf_", i) := apply(.SD, 1, sum),,.SDcols = paste0("mutcoeff_", 0:i)]
	    rand_dt[, paste0("mutenergyb_", i) := apply(.SD, 1, sum),,.SDcols = paste0("mutcoefb_", 0:i)]
	    rand_dt[, paste0("mutfitnessa_", i) := archstabms__predict_fitness(mochi_outpath = mochi_outpath, folding_energy = .SD[[1]], binding_energy = NULL, RT = RT, phenotype_names = c("Abundance", "Activity"), dataset_name = "b1")[["fitness_folding"]],,.SDcols = paste0("mutenergyf_", i)]
	    rand_dt[, paste0("mutfitnessb_", i) := archstabms__predict_fitness(mochi_outpath = mochi_outpath, folding_energy = .SD[[1]], binding_energy = .SD[[2]], RT = RT, phenotype_names = c("Abundance", "Activity"), dataset_name = "b1")[["fitness_binding"]],,.SDcols = paste0(c("mutenergyf_", "mutenergyb_"), i)]
	    rand_dt[, paste0("mutgra_", i) := gra_lm[["coefficients"]][1] + gra_lm[["coefficients"]][2] * .SD[[1]],,.SDcols = paste0("mutfitnessa_", i)]
	    rand_dt[, paste0("mutgrb_", i) := grb_lm[["coefficients"]][1] + grb_lm[["coefficients"]][2] * .SD[[1]],,.SDcols = paste0("mutfitnessb_", i)]
	  }

	  #Plot - abundance
	  plot_dt <- melt(rand_dt[,.SD,,.SDcols = paste0("mutgra_", 1:max_num_mutants)], measure.vars = paste0("mutgra_", 1:max_num_mutants))
	  plot_dt[, mut_order := factor(sapply(strsplit(as.character(variable), "_"), '[', 2), levels = as.character(c(1:max_num_mutants)))]
	  #Percentage of folded variants
	  text_dt <- plot_dt[,.(label = round(sum(value>gr_thresholds[["Abundance"]])/.N*100, 0), value = 0.17),mut_order]
	  cc <- scales::seq_gradient_pal(colour_scheme[["shade 0"]][3], colour_scheme[["shade 0"]][1], "Lab")(seq(0,1,length.out=plot_dt[, max(as.integer(as.character(mut_order)))]))
	  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(mut_order, value, fill = mut_order, color = mut_order)) +
	    ggplot2::geom_violin(scale = 'width') +
	    # ggplot2::geom_boxplot() +
	    ggplot2::scale_colour_manual(values = cc) +
	    ggplot2::scale_fill_manual(values = cc) +
	    ggplot2::geom_hline(yintercept = wta_gr) +
	    ggplot2::geom_hline(yintercept = gr_thresholds[["Abundance"]], linetype = 2) +
	    ggplot2::theme_classic() +
	    ggplot2::geom_text(data = plot_dt[,.(label = paste("n = ", .N, sep=""), mut_order),][1], ggplot2::aes(label=label, x=-Inf, y=-Inf, hjust = 0, vjust = 0), color = "black") +
	    ggplot2::geom_text(data = text_dt, ggplot2::aes(label=label)) +
	    ggplot2::xlab("AA substitution order") +
	    ggplot2::ylab("Predicted growth rate")
	  ggplot2::ggsave(file.path(outpath, paste0("growthrate_violins_", domain_name, "_abundance_Nham_aa.pdf")), d, width = 10, height = 3, useDingbats=FALSE)

	  #Plot - activity
	  plot_dt <- melt(rand_dt[,.SD,,.SDcols = paste0("mutgrb_", 1:max_num_mutants)], measure.vars = paste0("mutgrb_", 1:max_num_mutants))
	  plot_dt[, mut_order := factor(sapply(strsplit(as.character(variable), "_"), '[', 2), levels = as.character(c(1:max_num_mutants)))]
	  #Percentage of folded variants
	  text_dt <- plot_dt[,.(label = round(sum(value>gr_thresholds[["Activity"]])/.N*100, 0), value = 0.17),mut_order]
	  cc <- scales::seq_gradient_pal(colour_scheme[["shade 0"]][3], colour_scheme[["shade 0"]][1], "Lab")(seq(0,1,length.out=plot_dt[, max(as.integer(as.character(mut_order)))]))
	  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(mut_order, value, fill = mut_order, color = mut_order)) +
	    ggplot2::geom_violin(scale = 'width') +
	    # ggplot2::geom_boxplot() +
	    ggplot2::scale_colour_manual(values = cc) +
	    ggplot2::scale_fill_manual(values = cc) +
	    ggplot2::geom_hline(yintercept = wtb_gr) +
	    ggplot2::geom_hline(yintercept = gr_thresholds[["Activity"]], linetype = 2) +
	    ggplot2::theme_classic() +
	    ggplot2::geom_text(data = plot_dt[,.(label = paste("n = ", .N, sep=""), mut_order),][1], ggplot2::aes(label=label, x=-Inf, y=-Inf, hjust = 0, vjust = 0), color = "black") +
	    ggplot2::geom_text(data = text_dt, ggplot2::aes(label=label)) +
	    ggplot2::xlab("AA substitution order") +
	    ggplot2::ylab("Predicted growth rate")
	  ggplot2::ggsave(file.path(outpath, paste0("growthrate_violins_", domain_name, "_activity_Nham_aa.pdf")), d, width = 10, height = 3, useDingbats=FALSE)

		return(list(
			ids = selected_id, 
			folding_medians = median_values_folding, 
			activity_medians = median_values_activity))
	}
}
