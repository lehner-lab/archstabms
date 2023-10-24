
#' archstabms__window_library_design
#'
#' Saturation combinatorial mutagenesis library design.
#'
#' @param fitness_abundance_file path to abundance fitness data (required)
#' @param pdb_file path to PDB file (required)
#' @param pdb_chain_query query chain id (default:A)
#' @param pdb_chain_target target chain id (default:B)
#' @param final_window_size final window size (default:22)
#' @param final_window_pos final window start position (default:10)
#' @param colour_scheme colour scheme file (required)
#' @param outpath output path for plots and saved objects (required)
#'
#' @return Nothing
#' @export
#' @import data.table
archstabms__window_library_design <- function(
  fitness_abundance_file,
	pdb_file,
	pdb_chain_query = "A",
	pdb_chain_target = "B",
	final_window_size = 22,
	final_window_pos = 10,	
	colour_scheme,
	outpath
	){

  #Domain name
  domain_name <- rev(unlist(strsplit(basename(outpath), "_")))[1]

	#Display status
	message(paste("\n\n*******", paste("running stage: archstabms__window_library_design for", domain_name), "*******\n\n"))

	#Create output directory
	archstabms__create_dir(archstabms_dir = outpath)

  #Fitness - abundance
  load(fitness_abundance_file)

	#File names
	pdb_ligand_file <- pdb_file
	pdb_RSASA_file <- gsub(".pdb$", "_RSASA.pdb", pdb_file)
	pdb_depth_file <- gsub(".pdb$", "_residue_depth.pdb", pdb_file)
	pdb_ss_file <- pdb_file
	pdb_cbeta_file <- gsub(".pdb$", "_RSASA.pdb", pdb_file)
	pdb_angle_file <- gsub(".pdb$", "_RSASA.pdb", pdb_file)

	#Calculate minimum ligand distances
	ligdist_dt <- archstabms__pairwise_distances_from_PDB(
		input_file = pdb_ligand_file,
		chain_query = pdb_chain_query,
		chain_target = pdb_chain_target)
	ligdist_dt[, Pos_ref := as.character(Pos)]

	#Calculate pairwise distances - scHAmin
	dist_dt_scHA <- archstabms__distance_matrix_from_PDB(
		input_file = pdb_ligand_file,
		metric = 'scHA')
	dist_dt_calpha <- archstabms__distance_matrix_from_PDB(
		input_file = pdb_ligand_file,
		metric = 'calpha')

	#Merge distances and contacts
	dist_dt <- merge(dist_dt_scHA, dist_dt_calpha, by = c("Pos1", "Pos2"))
	dist_dt[, Pos_ref := paste0(Pos1, "_", Pos2)]

	#Calculate pairwise distances - backbone
	dist_dt[, backbone := abs(Pos1-Pos2)]

	#Get relative SASA
	sasa_dt <- archstabms__temperature_factor_from_PDB(
		input_file = pdb_RSASA_file,
		chain = pdb_chain_query)
	sasa_dt <- sasa_dt[,.(Pos_ref = as.character(Pos), RSASA = bfactor)]

	#Get residue depth
	depth_dt <- archstabms__temperature_factor_from_PDB(
		input_file = pdb_depth_file,
		chain = pdb_chain_query)
	depth_dt <- depth_dt[,.(Pos_ref = as.character(Pos), depth = bfactor)]

	#Get secondary structure
	ss_dt <- archstabms__secondary_structure_from_PDB(
		input_file = pdb_ss_file,
		chain = pdb_chain_query)
	ss_dt[, Pos_ref := as.character(Pos)]

	#Get residue depth
	angle_dt <- archstabms__relative_angle_from_PDB(
		input_file = pdb_angle_file,
		chain = pdb_chain_query)
	angle_dt <- angle_dt[,.(Pos_ref = as.character(Pos), relative_angle = bfactor)]

	#Get ligand-relative annotations
	bi_residues <- ligdist_dt[scHAmin_ligand<5,Pos]
	ss_residues <- as.integer(unique(unlist(dist_dt_scHA[scHAmin<5][Pos1 %in% bi_residues | Pos2 %in% bi_residues, .(Pos1, Pos2)])))
	ba_residues <- unique(c(bi_residues-1, bi_residues+1))
	anno_dt <- copy(ligdist_dt[,.(Pos, Pos_ref)])
	anno_dt[, anno_ligand := 'distal']
	anno_dt[Pos %in% ba_residues, anno_ligand := 'adjacent']
	anno_dt[Pos %in% ss_residues, anno_ligand := 'second_shell']
	anno_dt[Pos %in% bi_residues, anno_ligand := 'binding_interface']

	#Residue features
	res_dt <- merge(ligdist_dt, sasa_dt[,.SD,,.SDcols = names(sasa_dt)[names(sasa_dt)!="Pos"]], by = "Pos_ref", all.x = T)
	res_dt <- merge(res_dt, depth_dt[,.SD,,.SDcols = names(depth_dt)[names(depth_dt)!="Pos"]], by = "Pos_ref", all.x = T)
	res_dt <- merge(res_dt, ss_dt[,.SD,,.SDcols = names(ss_dt)[names(ss_dt)!="Pos"]], by = "Pos_ref", all.x = T)
	res_dt <- merge(res_dt, angle_dt[,.SD,,.SDcols = names(angle_dt)[names(angle_dt)!="Pos"]], by = "Pos_ref", all.x = T)
	res_dt <- merge(res_dt, anno_dt[,.SD,,.SDcols = names(anno_dt)[names(anno_dt)!="Pos"]], by = "Pos_ref", all.x = T)

	#Define position class
	res_dt[RSASA<0.25, Pos_class := "core"]
	res_dt[RSASA>=0.25, Pos_class := "surface"]
	res_dt[scHAmin_ligand<5, Pos_class := "binding_interface"]

	#Surface beta strand non-binding interface residues
	res_subset <- res_dt[Pos_class != 'binding_interface'][order(as.integer(Pos_ref)),Pos_ref]

	### How many positions within all possible residue windows of length N (20-22) meet single AA substitution criteria?
	### 1. no/intermediate fitness effect (within 1/3rd fitness distribution IQR of WT)
	### 2. reacheable by single nucleotide substitution
	### 3. not in binding interface
	###########################

	reg_result_list <- list()
	iqr_prop <- 1/3
	window_sizes <- 20:22

	#Load singles
	singles_dt <- singles
	wtseq <- all_variants[WT==T,nt_seq]

	#Loop over all window sizes
	for(mutlenAA in window_sizes){
		reg_result_list[[as.character(mutlenAA)]] <- list()
		#Loop over all starting positions
		for(start_pos in 1:(nchar(wtseq)/3-(mutlenAA-1))){
			mutseq <- paste0(archstabms__get_codons(wtseq)[start_pos:(start_pos+mutlenAA-1)], collapse = "")
			mutpos <- start_pos:(start_pos+mutlenAA-1)
			mutpos <- mutpos[mutpos %in% res_subset]
			#Single AA mutants reacheable with single nucleotide substitutions
			reacheable_info <- archstabms__reachable_single_AA_mutants(mutseq, start_pos-1)
			singles_dt[, mut_code := paste0(WT_AA, Pos, Mut)]
			singles_dt[, Nham_nt1 := mut_code %in% reacheable_info[["reachable_muts"]]]
			#x*IQR range
			iqrx <- singles_dt[Mut!="*",iqr_prop * (quantile(fitness, 0.75) - quantile(fitness, 0.25))]
			#Number of unique positions with candidate mutations
			reg_result_list[[as.character(mutlenAA)]][[as.character(start_pos)]] <- length(unique(singles_dt[Pos %in% mutpos & (fitness+1.96*sigma) < iqrx & (fitness-1.96*sigma) > -iqrx & Nham_nt1, Pos]))
		}
		reg_result_list[[as.character(mutlenAA)]] <- data.table(
			region_length = mutlenAA, 
			residue_start = names(reg_result_list[[as.character(mutlenAA)]]), 
			candidate_positions = unlist(reg_result_list[[as.character(mutlenAA)]]))
	}

	#Plot
	plot_dt <- rbindlist(reg_result_list)
	plot_dt[, residue_start := factor(residue_start, levels = unique(residue_start[order(as.integer(residue_start))]))] 
	cc <- c('lightgrey', 'darkgrey', 'black')
	names(cc) <- plot_dt[,unique(region_length)]
	d <- ggplot2::ggplot(plot_dt,ggplot2::aes(x = residue_start, y = candidate_positions, fill = as.factor(region_length))) +
	  ggplot2::geom_col(position = ggplot2::position_dodge()) +
	  ggplot2::xlab("Starting residue position") +
	  ggplot2::ylab("Num. candidate positions") +
	  ggplot2::labs(fill = "Variable region\nlength (residues)") +
		ggplot2::geom_hline(ggplot2::aes(yintercept = 15), linetype = 2) +
		ggplot2::scale_fill_manual(values = cc) +
	  ggplot2::theme_classic()
	ggplot2::ggsave(file.path(outpath, paste0("num_candidate_mutant_pos_window.pdf")), d, width = 8, height = 3, useDingbats=FALSE)

	### Final mutations
	###########################

	mutseq <- paste0(archstabms__get_codons(wtseq)[final_window_pos:(final_window_pos+final_window_size-1)], collapse = "")
	mutpos <- final_window_pos:(final_window_pos+final_window_size-1)
	mutpos <- mutpos[mutpos %in% res_subset]
	#Single AA mutants reacheable with single nucleotide substitutions
	reacheable_info <- archstabms__reachable_single_AA_mutants(mutseq, final_window_pos-1)
	singles_dt[, mut_code := paste0(WT_AA, Pos, Mut)]
	singles_dt[, Nham_nt1 := mut_code %in% reacheable_info[["reachable_muts"]]]
	#x*IQR range
	iqrx <- singles_dt[Mut!="*",iqr_prop * (quantile(fitness, 0.75) - quantile(fitness, 0.25))]
	#Candidate mutations
	cand_mut <-singles_dt[Pos %in% mutpos & (fitness+1.96*sigma) < iqrx & (fitness-1.96*sigma) > -iqrx & Nham_nt1,]
	#Select random mutation per candidate position
	set.seed(1)
	write.table(cand_mut[,sample(mut_code, 1),Pos][order(Pos)][,V1], file = file.path(outpath, "final_singles.txt"), sep = "\n", quote = F, row.names = F, col.names = F)

}
