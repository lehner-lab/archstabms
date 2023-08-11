
#' archstabms__distance_matrix_from_PDB
#'
#' Residue distance matrix from PDB.
#'
#' @param input_file path to PDB file (required)
#' @param residue_subset subset of residues (default:all)
#' @param query_chain query chain (default:A)
#' @param metric distance metric (default:'scHA')
#'
#' @return data.table with minimum residue distances
#' @export
#' @import data.table
archstabms__distance_matrix_from_PDB <- function(
	input_file,
	residue_subset="all",
	query_chain="A",
	metric="scHA"
	){

	#load PDB structure
	sink(file = "/dev/null")
	pdb <- bio3d::read.pdb(input_file, rm.alt = TRUE)
	sink()

	#All residue numbers
	pdb_dt <- as.data.table(pdb$atom)
	resnos <- pdb_dt[chain==query_chain & type=="ATOM"][order(resno)][,unique(resno)]

	#Residue subset
	if(residue_subset[1]=="all"){
			residue_subset <- resnos
	}

	#Distance matrix
	dist_mat <- archstabms__pairwise_distances_from_PDB(
		input_file, 
		chain_query = query_chain, 
		chain_target = query_chain,
		min_dist = FALSE)[[metric]]
	#Subset
	dist_mat <- dist_mat[which(resnos %in% residue_subset),which(resnos %in% residue_subset)]
	#Melt
	dist_dt <- data.table(reshape2::melt(dist_mat))
	colnames(dist_dt) <- c("Pos1", "Pos2", paste0(metric, "min"))

	return(dist_dt[Pos1<Pos2])
}
