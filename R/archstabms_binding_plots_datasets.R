
#' archstabms_binding_plots_datasets
#'
#' Coupling scatterplots for multiple datasets.
#'
#' @param dataset_names character vector of dataset names (required)
#' @param base_dir Base directory (required)
#' @param output_dir Output directory (required)
#' @param stagenum stage number (required)
#' @param colour_scheme colour scheme file (required)
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return Nothing
#' @export
archstabms_binding_plots_datasets <- function(
  dataset_names,
  base_dir,
  output_dir,
  stagenum,
  colour_scheme,
  execute = TRUE
  ){

  #Return if analysis not executed
  if(!execute){
    return()
  }

  for(i in dataset_names){
    #1st order coefficient
    archstabms_binding_plots(
      dataset_name = i,
      coef_order = 1,
      base_dir = base_dir,
      colour_scheme = colour_scheme,
      outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_binding_plots_", i, "o1"), stagenum=stagenum, base_dir=output_dir))
    #2nd order coefficient
    archstabms_binding_plots(
      dataset_name = i,
      coef_order = 2,
      base_dir = base_dir,
      colour_scheme = colour_scheme,
      outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_binding_plots_", i, "o2"), stagenum=stagenum, base_dir=output_dir))
  }
}
