
#' archstabms_coupling_heatmap_datasets
#'
#' Coupling scatterplots for multiple datasets.
#'
#' @param dataset_names character vector of dataset names (required)
#' @param base_dir Base directory (required)
#' @param output_dir Output directory (required)
#' @param stagenum stage number (required)
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return Nothing
#' @export
archstabms_coupling_heatmap_datasets <- function(
  dataset_names,
  base_dir,
  output_dir,
  stagenum,
  execute = TRUE
  ){

  #Return if analysis not executed
  if(!execute){
    return()
  }

  for(i in dataset_names){
    #order 2
    archstabms_coupling_heatmap(
      input_file = file.path(base_dir, paste0("003", "_archstabms_structure_metrics_", i, "o2"), "model_coefficients.txt"),
      outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_coupling_heatmap_", i, "o2"), stagenum=stagenum, base_dir=output_dir))
    #order 2 - binding
    archstabms_coupling_heatmap(
      input_file = file.path(base_dir, paste0("003", "_archstabms_structure_metrics_", i, "o2"), "model_coefficients.txt"),
      trait = "Binding",
      outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_coupling_heatmap_", i, "o2b"), stagenum=stagenum, base_dir=output_dir))
  }
}
