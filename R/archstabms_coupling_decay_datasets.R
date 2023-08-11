
#' archstabms_coupling_decay_datasets
#'
#' Coupling predictive decay plots for multiple datasets.
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
archstabms_coupling_decay_datasets <- function(
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
    #order 1,2,3 - up to 13th order mutants
    archstabms_coupling_decay(
      input_dir_list = list(
        "order1" = file.path(base_dir, "Data", "mochi", paste0(i, "o1-morder13sparse")),
        "order2" = file.path(base_dir, "Data", "mochi", paste0(i, "o2-morder13sparse")),
        "order3" = file.path(base_dir, "Data", "mochi", paste0(i, "o3-morder13sparse"))),
      min_mut_order = 14,
      test_variants = TRUE,
      outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_coupling_decay_sparse_", i), stagenum=stagenum, base_dir=output_dir),
      colour_scheme = colour_scheme)
    #order 1,2,3 - equal number of random variants
    archstabms_coupling_decay(
      input_dir_list = list(
        "order1" = file.path(base_dir, "Data", "mochi", paste0(i, "o1-morder13sparserand")),
        "order2" = file.path(base_dir, "Data", "mochi", paste0(i, "o2-morder13sparserand")),
        "order3" = file.path(base_dir, "Data", "mochi", paste0(i, "o3-morder13sparserand"))),
      min_mut_order = 14,
      test_variants = TRUE,
      outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_coupling_decay_sparse_random_", i), stagenum=stagenum, base_dir=output_dir),
      colour_scheme = colour_scheme)

    # #order 1,2,3 - up to 13th order mutants
    # archstabms_coupling_decay(
    #   input_dir_list = list(
    #     "order1" = file.path(base_dir, "Data", "mochi", paste0(i, "o1-morder13")),
    #     "order2" = file.path(base_dir, "Data", "mochi", paste0(i, "o2-morder13")),
    #     "order3" = file.path(base_dir, "Data", "mochi", paste0(i, "o3-morder13"))),
    #   min_mut_order = 14,
    #   test_variants = FALSE,
    #   outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_coupling_decay_train_", i), stagenum=stagenum, base_dir=output_dir),
    #   colour_scheme = colour_scheme)
    #order 1,2,3 - equal number of random variants
    archstabms_coupling_decay(
      input_dir_list = list(
        "order1" = file.path(base_dir, "Data", "mochi", paste0(i, "o1-morder13sparserand")),
        "order2" = file.path(base_dir, "Data", "mochi", paste0(i, "o2-morder13sparserand")),
        "order3" = file.path(base_dir, "Data", "mochi", paste0(i, "o3-morder13sparserand"))),
      min_mut_order = 10,
      test_variants = FALSE,
      outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_coupling_decay_sparse_random_train_", i), stagenum=stagenum, base_dir=output_dir),
      colour_scheme = colour_scheme)

  }
}
