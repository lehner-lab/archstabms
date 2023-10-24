
#' archstabms_coupling_heatmap
#'
#' Coupling scatterplots.
#'
#' @param input_file path to MoCHI thermo model fit results (required)
#' @param trait trait name (default:"Folding")
#' @param outpath output path for plots and saved objects (required)
#'
#' @return Nothing
#' @export
#' @import data.table
archstabms_coupling_heatmap <- function(
	input_file,
  trait="Folding",
	outpath
	){

  #Return if input_file doesn't exist
  if(!file.exists(input_file)){
    return()
  }

	#Load free energies
	dg_dt_singles <- fread(input_file)[coef_order==1]

  #Return if trait doesn't exist
  if(!trait %in% dg_dt_singles[,trait_name]){
    return()
  }

	#Display status
	message(paste("\n\n*******", paste("running stage: archstabms_coupling_heatmap"), "*******\n\n"))

	#Create output directory
	archstabms__create_dir(archstabms_dir = outpath)

	#Load free energies
	dg_dt_singles <- dg_dt_singles[trait_name==trait,]
	ref_ids <- dg_dt_singles[,id_ref]
	names(ref_ids) <- dg_dt_singles[,Pos_ref]
	dg_dt <- fread(input_file)[coef_order==2 & trait_name==trait][!duplicated(Pos_ref)]

	#Plot
	temp_dt <- copy(dg_dt)
	temp_dt[, Pos_ref1 := as.integer(sapply(strsplit(Pos_ref, "_"), '[', 1))]
	temp_dt[, Pos_ref2 := as.integer(sapply(strsplit(Pos_ref, "_"), '[', 2))]
	heat_mat <- matrix(NA, nrow = temp_dt[,max(Pos_ref2)], ncol = temp_dt[,max(Pos_ref2)])
	heat_txt <- matrix("", nrow = temp_dt[,max(Pos_ref2)], ncol = temp_dt[,max(Pos_ref2)])
	for(i in 1:temp_dt[,.N]){
	  heat_mat[temp_dt[i,Pos_ref2],temp_dt[i,Pos_ref1]] <- temp_dt[i,.SD[[1]],,.SDcols = 'mean_kcal/mol']
	  heat_mat[temp_dt[i,Pos_ref1],temp_dt[i,Pos_ref2]] <- temp_dt[i,.SD[[1]],,.SDcols = 'mean_kcal/mol']
	  heat_txt[temp_dt[i,Pos_ref2],temp_dt[i,Pos_ref1]] <- temp_dt[i,c('', '*')[as.integer(scHAmin<8)+1]]
	  heat_txt[temp_dt[i,Pos_ref1],temp_dt[i,Pos_ref2]] <- temp_dt[i,c('', '*')[as.integer(scHAmin<8)+1]]
	}
	#Remove unnecessary rows and columns
	pos_window <- temp_dt[,min(Pos_ref1)]:temp_dt[,max(Pos_ref2)]
	heat_mat <- heat_mat[pos_window,pos_window]
	heat_txt <- heat_txt[pos_window,pos_window]
	rownames(heat_mat) <- pos_window
	rownames(heat_mat)[rownames(heat_mat) %in% names(ref_ids)] <- ref_ids[rownames(heat_mat)[rownames(heat_mat) %in% names(ref_ids)]]
	colnames(heat_mat) <- rownames(heat_mat)
	archstabms__tile_heatmap_wrapper(
	  input_matrix = heat_mat,
	  input_matrix_text = heat_txt,
	  text_size = 5,
	  output_file = file.path(outpath, "heatmap_order2.pdf"),
	  cluster = 'none',
	  width = 7, height = 6,
	  xaxis_angle=90, 
	  xaxis_hjust=1, 
	  xaxis_vjust=0.5, 
	  xaxis_size=7, 
	  yaxis_size=7,
	  na_colour="lightgrey",
	  xlab = "",
	  ylab = "")
	archstabms__tile_heatmap_wrapper(
	  input_matrix = abs(heat_mat),
	  input_matrix_text = heat_txt,
	  text_size = 5,
	  output_file = file.path(outpath, "heatmap_order2_abs.pdf"),
	  cluster = 'none',
	  width = 7, height = 6,
	  xaxis_angle=90, 
	  xaxis_hjust=1, 
	  xaxis_vjust=0.5, 
	  xaxis_size=7, 
	  yaxis_size=7,
	  colour_high='red',
	  na_colour="lightgrey",
	  xlab = "",
	  ylab = "")
	archstabms__tile_heatmap_wrapper(
	  input_matrix = heat_mat,
	  input_matrix_text = heat_txt,
	  text_size = 5,
	  output_file = file.path(outpath, "heatmap_order2_clip.pdf"),
	  cluster = 'none',
	  width = 7, height = 6,
	  xaxis_angle=90, 
	  xaxis_hjust=1, 
	  xaxis_vjust=0.5, 
	  xaxis_size=7, 
	  yaxis_size=7,
	  na_colour="lightgrey",
	  xlab = "",
	  ylab = "",
	  colour_clip=0.3,
	  colour_limits = c(-0.3, 0.3))
	archstabms__tile_heatmap_wrapper(
	  input_matrix = abs(heat_mat),
	  input_matrix_text = heat_txt,
	  text_size = 5,
	  output_file = file.path(outpath, "heatmap_order2_abs_clip.pdf"),
	  cluster = 'none',
	  width = 7, height = 6,
	  xaxis_angle=90, 
	  xaxis_hjust=1, 
	  xaxis_vjust=0.5, 
	  xaxis_size=7, 
	  yaxis_size=7,
	  colour_high='red',
	  na_colour="lightgrey",
	  xlab = "",
	  ylab = "",
	  colour_clip=0.3,
	  colour_limits = c(-0.3, 0.3))
}
