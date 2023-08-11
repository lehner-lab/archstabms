
#' archstabms_structure_metrics_datasets
#'
#' Add 3D structure metrics of multiple datasets.
#'
#' @param dataset_names character vector of dataset names (required)
#' @param pdb_file path to PDB file (required)
#' @param base_dir Base directory (required)
#' @param output_dir Output directory (required)
#' @param stagenum stage number (required)
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return Nothing
#' @export
archstabms_structure_metrics_datasets <- function(
  dataset_names,
  pdb_file,
  base_dir,
  output_dir,
  stagenum,
  execute = TRUE
  ){

  for(i in dataset_names){
    #order 1
    archstabms_structure_metrics(
      input_file = file.path(base_dir, paste0("001", paste0("_archstabms_thermo_model_results_", i, "o1")), "model_coefficients.txt"),
      outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_structure_metrics_", i, "o1"), stagenum=stagenum, base_dir=output_dir),
      pdb_file = pdb_file,
      getcontacts_WT_file = file.path(base_dir, "Data", "getcontacts", "2vwfH_A54P_noligand.tsv"),
      getcontacts_counts_file = file.path(base_dir, "Data", "getcontacts", "getcontacts_counts.txt"),
      execute = execute)
    #order 1 -linear
    archstabms_structure_metrics(
      input_file = file.path(base_dir, paste0("001", paste0("_archstabms_thermo_model_results_", i, "o1-linear")), "model_coefficients.txt"),
      outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_structure_metrics_", i, "o1-linear"), stagenum=stagenum, base_dir=output_dir),
      pdb_file = pdb_file,
      getcontacts_WT_file = file.path(base_dir, "Data", "getcontacts", "2vwfH_A54P_noligand.tsv"),
      getcontacts_counts_file = file.path(base_dir, "Data", "getcontacts", "getcontacts_counts.txt"),
      execute = execute)
    #order 2
    archstabms_structure_metrics(
      input_file = file.path(base_dir, paste0("001", paste0("_archstabms_thermo_model_results_", i, "o2")), "model_coefficients.txt"),
      outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_structure_metrics_", i, "o2"), stagenum=stagenum, base_dir=output_dir),
      pdb_file = pdb_file,
      getcontacts_WT_file = file.path(base_dir, "Data", "getcontacts", "2vwfH_A54P_noligand.tsv"),
      getcontacts_counts_file = file.path(base_dir, "Data", "getcontacts", "getcontacts_counts.txt"),
      execute = execute)
    #order 2 -linear
    archstabms_structure_metrics(
      input_file = file.path(base_dir, paste0("001", paste0("_archstabms_thermo_model_results_", i, "o2-linear")), "model_coefficients.txt"),
      outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_structure_metrics_", i, "o2-linear"), stagenum=stagenum, base_dir=output_dir),
      pdb_file = pdb_file,
      getcontacts_WT_file = file.path(base_dir, "Data", "getcontacts", "2vwfH_A54P_noligand.tsv"),
      getcontacts_counts_file = file.path(base_dir, "Data", "getcontacts", "getcontacts_counts.txt"),
      execute = execute)
    # #order 3
    # archstabms_structure_metrics(
    #   input_file = file.path(base_dir, paste0("001", paste0("_archstabms_thermo_model_results_", i, "o3")), "model_coefficients.txt"),
    #   outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_structure_metrics_", i, "o3"), stagenum=stagenum, base_dir=output_dir),
    #   pdb_file = pdb_file,
    #   execute = execute)
    # #order 3 -linear
    # archstabms_structure_metrics(
    #   input_file = file.path(base_dir, paste0("001", paste0("_archstabms_thermo_model_results_", i, "o3-linear")), "model_coefficients.txt"),
    #   outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_structure_metrics_", i, "o3-linear"), stagenum=stagenum, base_dir=output_dir),
    #   pdb_file = pdb_file,
    #   execute = execute)

    # #order 1 - ensemble
    # archstabms_structure_metrics(
    #   input_file = file.path(base_dir, paste0("001", paste0("_archstabms_thermo_model_results_", i, "o1-ensemble")), "model_coefficients.txt"),
    #   outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_structure_metrics_", i, "o1-ensemble"), stagenum=stagenum, base_dir=output_dir),
    #   pdb_file = pdb_file,
    #   execute = execute)
    # #order 1 -linear - ensemble
    # archstabms_structure_metrics(
    #   input_file = file.path(base_dir, paste0("001", paste0("_archstabms_thermo_model_results_", i, "o1-linear-ensemble")), "model_coefficients.txt"),
    #   outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_structure_metrics_", i, "o1-linear-ensemble"), stagenum=stagenum, base_dir=output_dir),
    #   pdb_file = pdb_file,
    #   execute = execute)
    # #order 2 - ensemble
    # archstabms_structure_metrics(
    #   input_file = file.path(base_dir, paste0("001", paste0("_archstabms_thermo_model_results_", i, "o2-ensemble")), "model_coefficients.txt"),
    #   outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_structure_metrics_", i, "o2-ensemble"), stagenum=stagenum, base_dir=output_dir),
    #   pdb_file = pdb_file,
    #   execute = execute)
    # #order 2 -linear - ensemble
    # archstabms_structure_metrics(
    #   input_file = file.path(base_dir, paste0("001", paste0("_archstabms_thermo_model_results_", i, "o2-linear-ensemble")), "model_coefficients.txt"),
    #   outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_structure_metrics_", i, "o2-linear-ensemble"), stagenum=stagenum, base_dir=output_dir),
    #   pdb_file = pdb_file,
    #   execute = execute)
    # #order 3 - ensemble
    # archstabms_structure_metrics(
    #   input_file = file.path(base_dir, paste0("001", paste0("_archstabms_thermo_model_results_", i, "o3-ensemble")), "model_coefficients.txt"),
    #   outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_structure_metrics_", i, "o3-ensemble"), stagenum=stagenum, base_dir=output_dir),
    #   pdb_file = pdb_file,
    #   execute = execute)
    # #order 3 -linear - ensemble
    # archstabms_structure_metrics(
    #   input_file = file.path(base_dir, paste0("001", paste0("_archstabms_thermo_model_results_", i, "o3-linear-ensemble")), "model_coefficients.txt"),
    #   outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_structure_metrics_", i, "o3-linear-ensemble"), stagenum=stagenum, base_dir=output_dir),
    #   pdb_file = pdb_file,
    #   execute = execute)

    # #order 3 - sparse
    # archstabms_structure_metrics(
    #   input_file = file.path(base_dir, paste0("001", paste0("_archstabms_thermo_model_results_", i, "o3-sparse")), "model_coefficients.txt"),
    #   outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_structure_metrics_", i, "o3-sparse"), stagenum=stagenum, base_dir=output_dir),
    #   pdb_file = pdb_file,
    #   execute = execute)
    # #order 3 -linear - sparse
    # archstabms_structure_metrics(
    #   input_file = file.path(base_dir, paste0("001", paste0("_archstabms_thermo_model_results_", i, "o3-sparse-linear")), "model_coefficients.txt"),
    #   outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_structure_metrics_", i, "o3-sparse-linear"), stagenum=stagenum, base_dir=output_dir),
    #   pdb_file = pdb_file,
    #   execute = execute)
    # #order 3 - ensemble - sparse
    # archstabms_structure_metrics(
    #   input_file = file.path(base_dir, paste0("001", paste0("_archstabms_thermo_model_results_", i, "o3-sparse-ensemble")), "model_coefficients.txt"),
    #   outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_structure_metrics_", i, "o3-sparse-ensemble"), stagenum=stagenum, base_dir=output_dir),
    #   pdb_file = pdb_file,
    #   execute = execute)
    # #order 3 -linear - ensemble - sparse
    # archstabms_structure_metrics(
    #   input_file = file.path(base_dir, paste0("001", paste0("_archstabms_thermo_model_results_", i, "o3-sparse-linear-ensemble")), "model_coefficients.txt"),
    #   outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_structure_metrics_", i, "o3-sparse-linear-ensemble"), stagenum=stagenum, base_dir=output_dir),
    #   pdb_file = pdb_file,
    #   execute = execute)
  }
}