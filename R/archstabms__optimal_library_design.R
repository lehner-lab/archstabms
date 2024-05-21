
#' archstabms__optimal_library_design
#'
#' Library design.
#'
#' @param ddpca_outpath path to ddPCA data (required)
#' @param position_offset residue position offset (default:0)
#' @param fitness_abundance_file path to abundance fitness data (required)
#' @param fitness_binding_file path to binding fitness data (required)
#' @param max_num_mutants maximum number of mutations (required)
#' @param explore_mutations vector of mutation ids to explore or "all_starting" (default:"all_starting")
#' @param explore_positions vector of positions to explore or "all_starting" (default:"all_starting")
#' @param RT constant (default:0.001987*(273+30))
#' @param sample_size number of random samples (default:1000) 
#' @param optimal_perc optimal percentage maximum geometric mean growth rate (default:0.7)
#' @param min_num_doubles minimum number of doubles per single for candidate mutations (default:20)
#' @param colour_scheme colour scheme file (required)
#' @param outpath output path for plots and saved objects (required)
#'
#' @return Nothing
#' @export
#' @import data.table
archstabms__optimal_library_design <- function(
	ddpca_outpath,
	position_offset = 0,
	fitness_abundance_file,
	fitness_binding_file,
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
	message(paste("\n\n*******", paste("running stage: archstabms__optimal_library_design for", domain_name), "*******\n\n"))

	#Create output directory
	archstabms__create_dir(archstabms_dir = outpath)

  #Fitness - abundance
  load(fitness_abundance_file)
  gra_lm <- all_variants[,lm(growthrate~fitness)]
  wta_gr <- all_variants[WT==T,growthrate] #WT growth rate

  #Fitness - binding
  load(fitness_binding_file)
  grb_lm <- all_variants[,lm(growthrate~fitness)]
  wtb_gr <- all_variants[WT==T,growthrate] #WT growth rate

	#WT nucleotide sequence
	nt_wtseq <- all_variants[WT==T,nt_seq]

  #Growth rate thresholds
  modpar <- fread(file.path(ddpca_outpath, "model_parameters_0.txt"), header = F)[,V1]
  modpar_list <- as.list(as.numeric(modpar[seq(2, length(modpar), 2)]))
  names(modpar_list) <- modpar[seq(1, length(modpar), 2)]
  gr_thresholds <- list(
    "Abundance" = as.numeric(predict(gra_lm, newdata = data.frame(fitness = 0.5*modpar_list[["folding_linear_kernel"]]+modpar_list[["folding_linear_bias"]]))),
    "Binding" = as.numeric(predict(grb_lm, newdata = data.frame(fitness = 0.5*modpar_list[["binding_linear_kernel"]]+modpar_list[["binding_linear_bias"]]))))
		
	#Initialise variables
	d_ddg_pred_sumsq <- 0
	selected_mutants_strong <- c()
	comb_dt_strong <- data.table()
	selected_mutants <- c()
	folding_median_values <- data.table(number_of_mutants = 2:max_num_mutants)
	binding_median_values <- data.table(number_of_mutants = 2:max_num_mutants)

	#Determine candidate mutations (confident binding and folding ddGs, reachable by single nt substitutions, at least 20 doubles with single)
	#WT nucleotide sequence
	load(fitness_abundance_file)
	nt_wtseq <- all_variants[WT==T,nt_seq]
	#Model results
	pred_dt <- fread(file.path(ddpca_outpath, "../model_results.txt"))
	ddg_dt <- fread(file.path(ddpca_outpath, "../dg_singles.txt"))
	pred_double1_dt <- copy(pred_dt[mut_order==2 & dataset_binding==0])
	pred_double1_dt[, id_ref_double := sapply(strsplit(id_ref, ","), "[", 1)]
	pred_double2_dt <- copy(pred_dt[mut_order==2 & dataset_binding==0])
	pred_double2_dt[, id_ref_double := sapply(strsplit(id_ref, ","), "[", 2)]
	pred_double_dt <- rbind(pred_double1_dt, pred_double2_dt)
	pred_double_dt[, num_doubles := .N,id_ref_double,]
	pred_double_dt <- pred_double_dt[!duplicated(id_ref_double),.(id_ref = id_ref_double, num_doubles)]
	ddg_dt <- merge(ddg_dt, pred_double_dt, by = "id_ref", all.x = T) #Merge with ddGs

	#Explore positions
	if(explore_positions[1]=='all_starting'){
		explore_positions <- unique(ddg_dt[id!="-0-",Pos_ref])
	}

	#Reachable mutants
	reachable_singles <- archstabms__reachable_single_AA_mutants(nt_wtseq, position_offset)

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
		median_values_binding <- c()
		median_values_folding <- c()

		#Folding and binding free energy combinatorial data.tables
		#WT
		combf_dt <- data.table(WT = rep(ddg_dt[id_ref=="-0-",f_dg_pred], sample_size))
		combb_dt <- data.table(WT = rep(ddg_dt[id_ref=="-0-",b_dg_pred], sample_size))
		#First mutant
		combf_dt[, paste0(mutation) := sample(0:1, sample_size, replace = TRUE)]
		combb_dt[, paste0(mutation) := combf_dt[,.SD[[1]]*ddgs_dt[id_ref==selected_id[1],b_ddg_pred],,.SDcols = mutation]]
		combf_dt[, paste0(mutation) := .SD[[1]]*ddgs_dt[id_ref==selected_id[1],f_ddg_pred],,.SDcols = mutation]

		#Greedily add new mutations up to max_num_mutants
		for(n_mut in 2:max_num_mutants){
			#Already selected mutations
			selected_dt<-ddgs_dt[id_ref %in% selected_id]
			#Free energies
			sum_ddg_folding<-sum(selected_dt[,f_ddg_pred])+ddg_dt[id_ref=="-0-",f_dg_pred]
			sum_ddg_binding<-sum(selected_dt[,b_ddg_pred])+ddg_dt[id_ref=="-0-",b_dg_pred]

			#Remaining mutations in other positions
			new_possible_mut<-ddgs_dt[!Pos %in% selected_dt[,Pos]]

			#No mutations left to select - break
			if(new_possible_mut[,.N]==0){
				break
			}
			
			#Compute abundance and binding growth rates of all current order +1 mutants
			#Free energies
			new_possible_mut[, f_dg_pred_comb := sum_ddg_folding+f_ddg_pred]
			new_possible_mut[, b_dg_pred_comb := sum_ddg_binding+b_ddg_pred]
			#Fitness
			new_possible_mut[, fitness_folding := archstabms__predict_fitness_ddPCA(ddpca_outpath, f_dg_pred_comb, 0, RT = RT)[["fitness_folding"]]] #fitness_binding
			new_possible_mut[, fitness_binding := archstabms__predict_fitness_ddPCA(ddpca_outpath, f_dg_pred_comb, b_dg_pred_comb, RT = RT)[["fitness_binding"]]] #fitness_binding
			#Growth rates
			new_possible_mut[, growthrate_folding := gra_lm[["coefficients"]][1] + gra_lm[["coefficients"]][2]*fitness_folding]
			new_possible_mut[, growthrate_binding := grb_lm[["coefficients"]][1] + grb_lm[["coefficients"]][2]*fitness_binding]

			#Compute geometric means of folding and binding growth rates
			new_possible_mut[, geom_mean := (growthrate_folding*growthrate_binding)**(1/2)]

			#Select optimal mutation
			selected_id<-c(selected_id, new_possible_mut[order(-geom_mean)][1, id_ref])

			#Compute median abundance and binding growth rates of entire library
			#Randomly sample additional mutant
			combf_dt[, paste0(selected_id[n_mut]) := sample(rep(c(0, 1), each = sample_size/2))]
			combb_dt[, paste0(selected_id[n_mut]) := combf_dt[,.SD[[1]]*ddgs_dt[id_ref == selected_id[n_mut],b_ddg_pred],,.SDcols = selected_id[n_mut]]]
			combf_dt[, paste0(selected_id[n_mut]) := .SD[[1]]*ddgs_dt[id_ref == selected_id[n_mut],f_ddg_pred],,.SDcols = selected_id[n_mut]]
			#Compute abundance and binding growth rates
			gr_dt <- data.table()
			#Free energies
			gr_dt[, f_dg_pred := rowSums(combf_dt)]
			gr_dt[, b_dg_pred := rowSums(combb_dt)]
			#Fitness
			temp_fitness <- archstabms__predict_fitness_ddPCA(ddpca_outpath, gr_dt[, f_dg_pred], gr_dt[, b_dg_pred], RT = RT)
			gr_dt[, a_fitness := temp_fitness[["fitness_folding"]]]
			gr_dt[, b_fitness := temp_fitness[["fitness_binding"]]]
			#Growth rates
			gr_dt[, a_growthrate := gra_lm[["coefficients"]][1] + gra_lm[["coefficients"]][2]*a_fitness]
			gr_dt[, b_growthrate := grb_lm[["coefficients"]][1] + grb_lm[["coefficients"]][2]*b_fitness]

			#Append the medians of the new libraries
			median_values_folding <- c(median_values_folding, gr_dt[,median(a_growthrate)])
			median_values_binding <- c(median_values_binding, gr_dt[,median(b_growthrate)])
		}
		#Store
		folding_median_values[1:length(median_values_folding), paste0(mutation) := median_values_folding]
		binding_median_values[1:length(median_values_binding), paste0(mutation) := median_values_binding]
	}
	
	if(explore_mutations=="all_starting"){
		#Prepare table of all computed values and store
		total <- merge(
			melt(folding_median_values, measure.vars = mutations_id, variable.name = "mutation", value.name = "folding"),
			melt(binding_median_values, measure.vars = mutations_id, variable.name = "mutation", value.name = "binding"), 
			by=c("number_of_mutants","mutation"))
		names(total) <- c("mutant_order","mutation", "Abundance", "Binding")
		#Set minimum growth rate to zero
		total[Abundance<0, Abundance := 0]
		total[Binding<0, Binding := 0]
		#Calculate geometric mean
		total[, geom_mean := (Abundance*Binding)**(1/2)]
		save(total, file = file.path(outpath, "binding_folding_means_all_computed.Rdata"))
		
		#Select best geometric mean library for each order mutant and store
		final_table <- total[, .SD[which.max(geom_mean)], by = mutant_order]
		save(final_table, file = file.path(outpath, "binding_folding_BestGeomMean_bestEachStep.Rdata"))
		write.table(final_table, file = file.path(outpath, "binding_folding_BestGeomMean_bestEachStep.csv"), sep = "\t", row.names = F)
		
		#Select optimal library
		print("Optimal higher order mutant combination:")
		opt_table <- final_table[geom_mean > (optimal_perc*max(geom_mean))][.N]
		print(opt_table)
		optimal_num_mutants <- opt_table[,mutant_order]
		optimal_first_mut <- opt_table[,mutation]
		
		#Plot - mutant order vs median growth rates progression
		table_plot <- melt(
			final_table, 
			measure.vars = c("Abundance", "Binding", "geom_mean"),
			variable.name = "assay", 
			value.name = "median_predicted_growth_rate")
		write.table(table_plot, file = file.path(outpath, "plot_table.csv"), sep = "\t", row.names = F)
		cc <- c(colour_scheme[["shade 0"]][[3]], colour_scheme[["shade 0"]][[1]], 'black')
		names(cc) <- c("Abundance", "Binding", "geom_mean")
		d <- ggplot2::ggplot(table_plot, ggplot2::aes(mutant_order, median_predicted_growth_rate, color = assay)) +
			ggplot2::geom_line(linetype = 1) + 
			ggplot2::theme_classic() +
			ggplot2::geom_point(size = 1) +
			ggplot2::geom_vline(xintercept = optimal_num_mutants, linetype = 2) +
			ggplot2::geom_hline(yintercept = final_table[,optimal_perc*max(geom_mean)], linetype = 2) +
			ggplot2::ylab("Predicted median growth rate") +
			ggplot2::xlab("AA substitution order") +
			ggplot2::scale_colour_manual(values = cc)
		ggplot2::ggsave(file.path(outpath, "bestGeomMeanGR_perOrderMutants.pdf"), d, width = 7, height = 3, useDingbats=FALSE)
				
		#Store medians for that starting position, for each mutant order library
		return(list(
			optimal_first_mut = optimal_first_mut, 
			optimal_num_mutants = optimal_num_mutants, 
			folding_medians = folding_median_values, 
			binding_medians = binding_median_values))

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
	  binding_energies_dict <- as.list(ddg_dt[selected_id,b_ddg_pred])
	  names(binding_energies_dict) <- selected_id

	  #Set variant ids
	  rand_dt <- data.table(varid = 1:sample_size)
	  #Select random mutation at each step
	  rand_dt <- rand_dt[, mutlist := paste0(sample(selected_id, max_num_mutants), collapse = ","),varid]
	  #Energies
	  rand_dt[, mutcoeff_0 := ddg_dt[id_ref=="-0-",f_dg_pred]]
	  rand_dt[, mutcoefb_0 := ddg_dt[id_ref=="-0-",b_dg_pred]]
	  for(i in 1:max_num_mutants){
	    rand_dt[, paste0("mutcoeff_", i) := unlist(folding_energies_dict[as.character(sapply(strsplit(mutlist, ","), '[', i))])]
	    rand_dt[, paste0("mutcoefb_", i) := unlist(binding_energies_dict[as.character(sapply(strsplit(mutlist, ","), '[', i))])]
	  }
	  #Predict growth rate
	  for(i in 1:max_num_mutants){
	    rand_dt[, paste0("mutenergyf_", i) := apply(.SD, 1, sum),,.SDcols = paste0("mutcoeff_", 0:i)]
	    rand_dt[, paste0("mutenergyb_", i) := apply(.SD, 1, sum),,.SDcols = paste0("mutcoefb_", 0:i)]
	    rand_dt[, paste0("mutfitnessa_", i) := archstabms__predict_fitness_ddPCA(mochi_outpath = ddpca_outpath, folding_energy = .SD[[1]], binding_energy = NULL, RT = RT)[["fitness_folding"]],,.SDcols = paste0("mutenergyf_", i)]
	    rand_dt[, paste0("mutfitnessb_", i) := archstabms__predict_fitness_ddPCA(mochi_outpath = ddpca_outpath, folding_energy = .SD[[1]], binding_energy = .SD[[2]], RT = RT)[["fitness_binding"]],,.SDcols = paste0(c("mutenergyf_", "mutenergyb_"), i)]
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

	  #Plot - binding
	  plot_dt <- melt(rand_dt[,.SD,,.SDcols = paste0("mutgrb_", 1:max_num_mutants)], measure.vars = paste0("mutgrb_", 1:max_num_mutants))
	  plot_dt[, mut_order := factor(sapply(strsplit(as.character(variable), "_"), '[', 2), levels = as.character(c(1:max_num_mutants)))]
	  #Percentage of folded variants
	  text_dt <- plot_dt[,.(label = round(sum(value>gr_thresholds[["Binding"]])/.N*100, 0), value = 0.17),mut_order]
	  cc <- scales::seq_gradient_pal(colour_scheme[["shade 0"]][3], colour_scheme[["shade 0"]][1], "Lab")(seq(0,1,length.out=plot_dt[, max(as.integer(as.character(mut_order)))]))
	  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(mut_order, value, fill = mut_order, color = mut_order)) +
	    ggplot2::geom_violin(scale = 'width') +
	    # ggplot2::geom_boxplot() +
	    ggplot2::scale_colour_manual(values = cc) +
	    ggplot2::scale_fill_manual(values = cc) +
	    ggplot2::geom_hline(yintercept = wtb_gr) +
	    ggplot2::geom_hline(yintercept = gr_thresholds[["Binding"]], linetype = 2) +
	    ggplot2::theme_classic() +
	    ggplot2::geom_text(data = plot_dt[,.(label = paste("n = ", .N, sep=""), mut_order),][1], ggplot2::aes(label=label, x=-Inf, y=-Inf, hjust = 0, vjust = 0), color = "black") +
	    ggplot2::geom_text(data = text_dt, ggplot2::aes(label=label)) +
	    ggplot2::xlab("AA substitution order") +
	    ggplot2::ylab("Predicted growth rate")
	  ggplot2::ggsave(file.path(outpath, paste0("growthrate_violins_", domain_name, "_binding_Nham_aa.pdf")), d, width = 10, height = 3, useDingbats=FALSE)

		return(list(
			ids = selected_id, 
			folding_medians = median_values_folding, 
			binding_medians = median_values_binding))
	}
}
