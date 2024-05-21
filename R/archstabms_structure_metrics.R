
#' archstabms_structure_metrics
#'
#' Add 3D structure metrics.
#'
#' @param input_file path to MoCHI thermo model fit results (required)
#' @param outpath output path for plots and saved objects (required)
#' @param pdb_file path to PDB file (required)
#' @param getcontacts_WT_file path to getcontacts WT file (required)
#' @param getcontacts_counts_file path to getcontacts counts file (required)
#' @param pdb_chain_query query chain id (default:A)
#' @param pdb_chain_target target chain id (default:B)
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return Nothing
#' @export
#' @import data.table
archstabms_structure_metrics <- function(
	input_file,
	outpath,
	pdb_file,
	getcontacts_WT_file,
	getcontacts_counts_file,
	pdb_chain_query = "A",
	pdb_chain_target = "B",
	execute = TRUE
	){

	#Return if analysis not executed
  if(!execute | !file.exists(input_file)){
		return()
	}

	#Display status
	message(paste("\n\n*******", paste("running stage: archstabms_structure_metrics for", basename(pdb_file)), "*******\n\n"))

	#Create output directory
	archstabms__create_dir(archstabms_dir = outpath)

	#File names
	pdb_ligand_file <- pdb_file
	pdb_RSASA_file <- gsub(".pdb$", "_RSASA.pdb", pdb_file)
	# pdb_depth_file <- gsub(".pdb$", "_residue_depth.pdb", pdb_file)
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

	#WT contacts
	gc_dt <- fread(getcontacts_WT_file, skip = 2)
	names(gc_dt) <- c("z", "contact", "p1", "p2")
	gc_dt[, Pos1temp := as.integer(sapply(strsplit(p1, ":"), '[', 3))]
	gc_dt[, Pos2temp := as.integer(sapply(strsplit(p2, ":"), '[', 3))]
	gc_dt[, Pos1 := apply(.SD, 1, min),,.SDcols = c("Pos1temp", "Pos2temp")]
	gc_dt[, Pos2 := apply(.SD, 1, max),,.SDcols = c("Pos1temp", "Pos2temp")]
	gcc_dt <- gc_dt[,.(count = .N),.(Pos1, Pos2, contact)]
	gcc_dt <- dcast(gcc_dt, Pos1 + Pos2 ~ contact, value.var = "count", fill = 0)
	gcc_dt[, contacts := apply(.SD, 1, sum),,.SDcols = names(gcc_dt)[-c(1,2)]]
	names(gcc_dt)[-c(1,2)] <- paste0("WT_", names(gcc_dt)[-c(1,2)])

	#Merge distances and contacts
	dist_dt <- merge(dist_dt_scHA, dist_dt_calpha, by = c("Pos1", "Pos2"))
	dist_dt <- merge(dist_dt, gcc_dt, by = c("Pos1", "Pos2"), all.x = T)
	nm1 <- names(gcc_dt)[-c(1,2)]
	setDT(dist_dt)[, (nm1) := lapply(.SD, nafill, fill = 0), .SDcols = nm1]
	dist_dt[, Pos_ref := paste0(Pos1, "_", Pos2)]

	#Calculate pairwise distances - backbone
	dist_dt[, backbone := abs(Pos1-Pos2)]

	#Get relative SASA
	sasa_dt <- archstabms__temperature_factor_from_PDB(
		input_file = pdb_RSASA_file,
		chain = pdb_chain_query)
	sasa_dt <- sasa_dt[,.(Pos_ref = as.character(Pos), RSASA = bfactor)]

	# #Get residue depth
	# depth_dt <- archstabms__temperature_factor_from_PDB(
	# 	input_file = pdb_depth_file,
	# 	chain = pdb_chain_query)
	# depth_dt <- depth_dt[,.(Pos_ref = as.character(Pos), depth = bfactor)]

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

	#Load free energies - order 0
	dg_dt <- fread(input_file)
	dg_dt[, Pos_ref := as.character(Pos_ref)]
	dg_dt_list <- list()
	dg_dt_list[['0']] <- dg_dt[coef_order==0]

	#Load free energies - order 1
	dg_dt_list[['1']] <- dg_dt[coef_order==1]

	#Merge with free energies
	dg_dt_list[['1']] <- merge(dg_dt_list[['1']], ligdist_dt[,.SD,,.SDcols = names(ligdist_dt)[names(ligdist_dt)!="Pos"]], by = "Pos_ref", all.x = T)
	dg_dt_list[['1']] <- merge(dg_dt_list[['1']], sasa_dt[,.SD,,.SDcols = names(sasa_dt)[names(sasa_dt)!="Pos"]], by = "Pos_ref", all.x = T)
	# dg_dt_list[['1']] <- merge(dg_dt_list[['1']], depth_dt[,.SD,,.SDcols = names(depth_dt)[names(depth_dt)!="Pos"]], by = "Pos_ref", all.x = T)
	dg_dt_list[['1']] <- merge(dg_dt_list[['1']], ss_dt[,.SD,,.SDcols = names(ss_dt)[names(ss_dt)!="Pos"]], by = "Pos_ref", all.x = T)
	dg_dt_list[['1']] <- merge(dg_dt_list[['1']], angle_dt[,.SD,,.SDcols = names(angle_dt)[names(angle_dt)!="Pos"]], by = "Pos_ref", all.x = T)
	dg_dt_list[['1']] <- merge(dg_dt_list[['1']], anno_dt[,.SD,,.SDcols = names(anno_dt)[names(anno_dt)!="Pos"]], by = "Pos_ref", all.x = T)

	#Define position class
	dg_dt_list[['1']][RSASA<0.25, Pos_class := "core"]
	dg_dt_list[['1']][RSASA>=0.25, Pos_class := "surface"]
	dg_dt_list[['1']][scHAmin_ligand<5, Pos_class := "binding_interface"]
  #Position class dict
  pclass_dict <- as.list(dg_dt_list[['1']][,Pos_class])
  names(pclass_dict) <- as.character(dg_dt_list[['1']][,Pos_ref])
  #Secondary structure dict
  ss_dict <- as.list(dg_dt_list[['1']][,SS])
  names(ss_dict) <- as.character(dg_dt_list[['1']][,Pos_ref])
  #Ligand distance dict
  ld_dict <- as.list(dg_dt_list[['1']][,scHAmin_ligand])
  names(ld_dict) <- as.character(dg_dt_list[['1']][,Pos_ref])

	#Load free energies - order 2
	dg_dt_list[['2']] <- dg_dt[coef_order==2]

	if(dg_dt_list[['2']][,.N]!=0){
		#Define position class
  	dg_dt_list[['2']][, Pos1_class := pclass_dict[[sapply(strsplit(Pos_ref, "_"), '[', 1)]],Pos_ref]
  	dg_dt_list[['2']][, Pos2_class := pclass_dict[[sapply(strsplit(Pos_ref, "_"), '[', 2)]],Pos_ref]
		#Define secondary structure
  	dg_dt_list[['2']][, SS1 := ss_dict[[sapply(strsplit(Pos_ref, "_"), '[', 1)]],Pos_ref]
  	dg_dt_list[['2']][, SS2 := ss_dict[[sapply(strsplit(Pos_ref, "_"), '[', 2)]],Pos_ref]
		#Define ligand distance
  	dg_dt_list[['2']][, scHAmin_ligand1 := ld_dict[[sapply(strsplit(Pos_ref, "_"), '[', 1)]],Pos_ref]
  	dg_dt_list[['2']][, scHAmin_ligand2 := ld_dict[[sapply(strsplit(Pos_ref, "_"), '[', 2)]],Pos_ref]
		#Merge with free energies
		dg_dt_list[['2']] <- merge(dg_dt_list[['2']], dist_dt[,.SD,,.SDcols = names(dist_dt)[!names(dist_dt) %in% c("Pos1", "Pos2")]], by = c("Pos_ref"), all.x = T)
	}

	#Save dGs and ddGs
	write.table(rbindlist(dg_dt_list, fill = T),
		file = file.path(outpath, "model_coefficients.txt"), 
		quote = F, sep = "\t", row.names = F)
}
