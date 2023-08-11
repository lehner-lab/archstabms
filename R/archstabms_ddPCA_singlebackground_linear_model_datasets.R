
#' archstabms_ddPCA_singlebackground_linear_model_datasets
#'
#' Biophysical model from shallow ddPCA data for multiple datasets.
#'
#' @param dataset_names character vector of dataset names (required)
#' @param seq_position_offset list of sequence position offsets (required)
#' @param base_dir Base directory (required)
#' @param output_dir Output directory (required)
#' @param stagenum stage number (required)
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return Nothing
#' @export
archstabms_ddPCA_singlebackground_linear_model_datasets <- function(
  dataset_names,
  seq_position_offset,
  base_dir,
  output_dir,
  stagenum,
  execute = TRUE
  ){

  #Return if analysis not executed
  if(!execute){
    return()
  }

  #Set up reticulate
  require("reticulate")
  use_condaenv(condaenv = "pymochi", conda = "~/conda/bin/conda")
  source_python(system.file("py", "load_mochi_data.py", package = "archstabms"))

  for(i in dataset_names){
    #order 1
    archstabms_ddPCA_singlebackground_linear_model(
      mochi_outpath = file.path(base_dir, "Data", "mochi", paste0(i, "o1")),
      ddpca_outpath = file.path(base_dir, "Data", "fitness", "ddPCA", "JD_GRB2_NM2_stabilityPCA_dimsum128_fitness_replicates.RData"),
      position_offset = seq_position_offset[[i]],
      load_mochi_data = load_mochi_data,
      outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_ddPCA_singlebackground_linear_model_", i, "o1"), stagenum=stagenum, base_dir=output_dir))
  }
}