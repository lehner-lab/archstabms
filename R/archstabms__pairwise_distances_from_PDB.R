
#' archstabms__pairwise_distances_from_PDB
#'
#' Calculate pairwise (side-chain) heavy atom distances or minimum thereof.
#'
#' @param input_file path to PDB file (required)
#' @param chain_query query chain id (default:A)
#' @param chain_target target chain id (default:B)
#' @param min_dist minimum distance (default:TRUE)
#'
#' @return data.table with minimum inter-chain (side-chain) heavy atom distances or list of data.tables with distances (if min_dist==FALSE)
#' @export
#' @import data.table
archstabms__pairwise_distances_from_PDB <- function(
	input_file,
	chain_query = "A",
	chain_target = "B",
	min_dist = TRUE
	){
	
	#load PDB structure
	sink(file = "/dev/null")
	pdb <- bio3d::read.pdb(input_file, rm.alt = TRUE)
	sink()

	### Atom selections
	###########################

	#Protein atoms
	sele_protein <- bio3d::atom.select(pdb, "protein", verbose=FALSE)
	#Hydrogen atoms
	sele_H <-bio3d::atom.select(pdb, "h", verbose=FALSE)
	#Water atoms
	sele_water <- bio3d::atom.select(pdb, "water", verbose=FALSE)
	#Side chain atoms
	sele_sc <- bio3d::atom.select(pdb, "sidechain", verbose=FALSE)
	#C-alpha atoms
	sele_ca <- bio3d::atom.select(pdb, "calpha", verbose=FALSE)
	#Glycine c-alpha atoms
	sele_glyca <- bio3d::atom.select(pdb, resid = "GLY", string = "calpha", verbose=FALSE)

	### Combine atom selections
	###########################

	#Heavy atoms
	sele_HA <- bio3d::combine.select(sele_protein, sele_H, sele_water, operator = "-", verbose=FALSE)

	#Side chain heavy atoms + c-alpha for glycine
	sele_prot_sc <- bio3d::combine.select(sele_protein, sele_sc, operator = "AND", verbose=FALSE)
	sele_prot_sc_glyca <- bio3d::combine.select(sele_prot_sc, sele_glyca, operator = "OR", verbose=FALSE)
	sele_scHA <- bio3d::combine.select(sele_prot_sc_glyca, sele_H, sele_water, operator = "-", verbose=FALSE)

	#List
	sele_list <- list(
		"HA" = sele_HA,
		"scHA" = sele_scHA,
		"calpha" = sele_ca)

	### Calculate pairwise distances or reductions thereof
	###########################

	result_list <- list()
	result_dt <- data.table()
	for(metric in names(sele_list)){
		#Distance matrix
		pdb_sub <- bio3d::trim.pdb(pdb, sele_list[[metric]])
		dist_mat <- bio3d::dm.xyz(pdb_sub$xyz, grpby=apply(pdb_sub$atom[,c("resno", "chain")], 1, paste, collapse = "_"), scut=0, mask.lower = FALSE)
		resno_sub <- unique(pdb_sub$atom[,c("resno", "chain")])

		#Check if reduction specified
		if(!min_dist){
			#Distance matrix
			result_list[[metric]] <- dist_mat[resno_sub[,"chain"]==chain_query,resno_sub[,"chain"]==chain_target]
		}else{
			#Reduce distance matrix
			ligand_dist <- dist_mat[resno_sub[,"chain"]==chain_query,resno_sub[,"chain"]==chain_target]
			#Absolute residue number
			ligand_dist_dt <- data.table(Pos = resno_sub[resno_sub[,"chain"]==chain_query,"resno"])
			#Reduce distance
			ligand_dist_dt[, paste0(metric, "min_ligand") := apply(ligand_dist, 1, min)]
			if(nrow(result_dt)==0){
				result_dt <- ligand_dist_dt
			}else{
				result_dt <- merge(result_dt, ligand_dist_dt, by = "Pos", all = T)
			}
		}
	}

	#Return
	if(!min_dist){
		return(result_list)
	}else{
		return(result_dt)
	}

}

