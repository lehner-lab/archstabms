
#' archstabms__saturation_library_design
#'
#' Saturation combinatorial mutagenesis library design.
#'
#' @param pdb_file path to PDB file (required)
#' @param pdb_chain_query query chain id (default:A)
#' @param pdb_chain_target target chain id (default:B)
#' @param colour_scheme colour scheme file (required)
#' @param outpath output path for plots and saved objects (required)
#'
#' @return Nothing
#' @export
#' @import data.table
archstabms__saturation_library_design <- function(
	pdb_file,
	pdb_chain_query = "A",
	pdb_chain_target = "B",
	colour_scheme,
	outpath
	){

  #Domain name
  domain_name <- rev(unlist(strsplit(basename(outpath), "_")))[1]

	#Display status
	message(paste("\n\n*******", paste("running stage: archstabms__saturation_library_design for", domain_name), "*******\n\n"))

	#Create output directory
	archstabms__create_dir(archstabms_dir = outpath)

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
	res_subset <- res_dt[Pos_class == 'surface' & SS == "sheet"][order(as.integer(Pos_ref)),Pos_ref]

	#Relevant pairs
	dist_dt[, Pos_ref1 := sapply(strsplit(Pos_ref, "_"), '[', 1)]
	dist_dt[, Pos_ref2 := sapply(strsplit(Pos_ref, "_"), '[', 2)]
	res_pair_subset <- dist_dt[Pos_ref1 %in% res_subset & Pos_ref2 %in% res_subset]

	#Contact map
	plot_mat <- matrix(NA, nrow = max(as.integer(res_subset)), ncol = max(as.integer(res_subset)))
	for(i in 1:res_pair_subset[,.N]){
		plot_mat[res_pair_subset[i,as.integer(Pos_ref1)], res_pair_subset[i,as.integer(Pos_ref2)]] <- as.integer(res_pair_subset[i,scHAmin]<5)
		plot_mat[res_pair_subset[i,as.integer(Pos_ref2)], res_pair_subset[i,as.integer(Pos_ref1)]] <- as.integer(res_pair_subset[i,scHAmin]<5)
	}
	rownames(plot_mat) <- paste0("Pos", 1:max(as.integer(res_subset)))
	colnames(plot_mat) <- paste0("Pos", 1:max(as.integer(res_subset)))
	plot_mat <- plot_mat[apply(!is.na(plot_mat), 1, sum, na.rm = T)!=0, apply(!is.na(plot_mat), 2, sum, na.rm = T)!=0]
	#Plot contact map
	archstabms__tile_heatmap_wrapper(plot_mat, file.path(outpath, "contact_heatmap.pdf"), xlab = 'Residue', ylab = "Residue",
		cluster = 'both', width = 4, height = 3.5, xaxis_angle = 90, xaxis_size = NULL, xaxis_hjust = 1)

}
