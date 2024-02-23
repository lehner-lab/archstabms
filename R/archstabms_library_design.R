
#' archstabms_library_design
#'
#' Library design.
#'
#' @param base_dir Base directory (required)
#' @param output_dir Output directory (required)
#' @param stagenum stage number (required)
#' @param colour_scheme colour scheme file (required)
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return Nothing
#' @export
#' @import data.table
archstabms_library_design <- function(
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

  #Set up reticulate
  require("reticulate")
  use_condaenv(condaenv = "pymochi", conda = "~/conda/bin/conda")
  source_python(system.file("py", "load_mochi_data.py", package = "archstabms"))

  ### GRB2-SH3
  ###########################

  #Random library design
  archstabms__random_library_design(
   	ddpca_outpath = file.path(base_dir, "Data", "mochi", "ddPCA", "mochi__fit_tmodel_3state_sparse_dimsum128"),
    fitness_abundance_file = file.path(base_dir, "Data", "fitness", "ddPCA", "JD_GRB2_NM2_stabilityPCA_dimsum128_fitness_replicates.RData"),
    fitness_binding_file = file.path(base_dir, "Data", "fitness", "ddPCA", "JD_GRB2_epPCA_bindingPCA_dimsum128_fitness_replicates.RData"),
	  max_num_mutants = 10, 
	  load_mochi_data = load_mochi_data,
	  colour_scheme = colour_scheme,
    outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_library_design_random_GRB2"), stagenum=stagenum, base_dir=output_dir))

  #Explore all first mutations up to maximum order
  exp_results <- archstabms__optimal_library_design(
    ddpca_outpath = file.path(base_dir, "Data", "mochi", "ddPCA", "mochi__fit_tmodel_3state_sparse_dimsum128"),
    fitness_abundance_file = file.path(base_dir, "Data", "fitness", "ddPCA", "JD_GRB2_NM2_stabilityPCA_dimsum128_fitness_replicates.RData"),
    fitness_binding_file = file.path(base_dir, "Data", "fitness", "ddPCA", "JD_GRB2_epPCA_bindingPCA_dimsum128_fitness_replicates.RData"),
  	max_num_mutants = 56,
    colour_scheme = colour_scheme,
    outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_library_design_explore_GRB2"), stagenum=stagenum, base_dir=output_dir))

  #Optimised library design
  opt_results <- archstabms__optimal_library_design(
    ddpca_outpath = file.path(base_dir, "Data", "mochi", "ddPCA", "mochi__fit_tmodel_3state_sparse_dimsum128"),
    fitness_abundance_file = file.path(base_dir, "Data", "fitness", "ddPCA", "JD_GRB2_NM2_stabilityPCA_dimsum128_fitness_replicates.RData"),
    fitness_binding_file = file.path(base_dir, "Data", "fitness", "ddPCA", "JD_GRB2_epPCA_bindingPCA_dimsum128_fitness_replicates.RData"),
    max_num_mutants = exp_results[['optimal_num_mutants']], 
    # max_num_mutants = 34, 
    explore_mutations = exp_results[['optimal_first_mut']],
    # explore_mutations = "Y51H",
    sample_size = 10000,
    colour_scheme = colour_scheme,
    outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_library_design_optimise_GRB2"), stagenum=stagenum, base_dir=output_dir))

  #Saturation combinatorial library design
  archstabms__saturation_library_design(
    pdb_file = file.path(base_dir, "Data", "pdb", "2vwf.pdb"),
    colour_scheme = colour_scheme,
    outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_library_design_saturation_GRB2"), stagenum=stagenum, base_dir=output_dir))

  #Window combinatorial library design
  archstabms__window_library_design(
    fitness_abundance_file = file.path(base_dir, "Data", "fitness", "ddPCA", "JD_GRB2_NM2_stabilityPCA_dimsum128_fitness_replicates.RData"),
    pdb_file = file.path(base_dir, "Data", "pdb", "2vwf.pdb"),
    colour_scheme = colour_scheme,
    outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_library_design_window_GRB2"), stagenum=stagenum, base_dir=output_dir))

  ### PSD95-PDZ3
  ###########################

  #Random library design
  archstabms__random_library_design(
    ddpca_outpath = file.path(base_dir, "Data", "mochi", "ddPCA_PDZ", "mochi__fit_tmodel_3state_sparse_dimsum128"),
    fitness_abundance_file = file.path(base_dir, "Data", "fitness", "ddPCA_PDZ", "JD_PDZ_NM2_stabilityPCA_dimsum128_filtered_fitness_replicates.RData"),
    fitness_binding_file = file.path(base_dir, "Data", "fitness", "ddPCA_PDZ", "JD_PDZ_NM2_bindingPCA_dimsum128_filtered_fitness_replicates.RData"),
    max_num_mutants = 10, 
    load_mochi_data = load_mochi_data,
    colour_scheme = colour_scheme,
    outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_library_design_random_PSD95"), stagenum=stagenum, base_dir=output_dir))

  #Explore all first mutations up to maximum order
  exp_results <- archstabms__optimal_library_design(
    ddpca_outpath = file.path(base_dir, "Data", "mochi", "ddPCA_PDZ", "mochi__fit_tmodel_3state_sparse_dimsum128"),
    position_offset = 310,
    fitness_abundance_file = file.path(base_dir, "Data", "fitness", "ddPCA_PDZ", "JD_PDZ_NM2_stabilityPCA_dimsum128_filtered_fitness_replicates.RData"),
    fitness_binding_file = file.path(base_dir, "Data", "fitness", "ddPCA_PDZ", "JD_PDZ_NM2_bindingPCA_dimsum128_filtered_fitness_replicates.RData"),
    max_num_mutants = 77,
    explore_positions = c(311:348, 357:394),
    min_num_doubles = 1,
    colour_scheme = colour_scheme,
    outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_library_design_explore_PSD95"), stagenum=stagenum, base_dir=output_dir))

  #Optimised library design
  opt_results <- archstabms__optimal_library_design(
    ddpca_outpath = file.path(base_dir, "Data", "mochi", "ddPCA_PDZ", "mochi__fit_tmodel_3state_sparse_dimsum128"),
    position_offset = 310,
    fitness_abundance_file = file.path(base_dir, "Data", "fitness", "ddPCA_PDZ", "JD_PDZ_NM2_stabilityPCA_dimsum128_filtered_fitness_replicates.RData"),
    fitness_binding_file = file.path(base_dir, "Data", "fitness", "ddPCA_PDZ", "JD_PDZ_NM2_bindingPCA_dimsum128_filtered_fitness_replicates.RData"),
    max_num_mutants = exp_results[['optimal_num_mutants']], 
    # max_num_mutants = 60, 
    explore_mutations = exp_results[['optimal_first_mut']],
    # explore_mutations = "P394A",
    explore_positions = c(311:348, 357:394),
    sample_size = 10000,
    min_num_doubles = 1,
    colour_scheme = colour_scheme,
    outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_library_design_optimise_PSD95"), stagenum=stagenum, base_dir=output_dir))

  ### SRC 1
  ###########################

  #Explore all first mutations up to maximum order
  exp_results <- archstabms__optimal_library_design_src(
    mochi_outpath = file.path(base_dir, "Data", "mochi", "ddPCA_SRC"),
    position_offset = 0,
    max_num_mutants = 15,
    explore_positions = 64:85,
    min_num_doubles = 7,
    colour_scheme = colour_scheme,
    outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_library_design_explore_SRC_64-85"), stagenum=stagenum, base_dir=output_dir))

  #Optimised library design
  opt_results <- archstabms__optimal_library_design_src(
    mochi_outpath = file.path(base_dir, "Data", "mochi", "ddPCA_SRC"),
    position_offset = 0,
    max_num_mutants = 15, 
    explore_mutations = exp_results[['optimal_first_mut']],
    explore_positions = 64:85,
    min_num_doubles = 7,
    sample_size = 10000,
    colour_scheme = colour_scheme,
    outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_library_design_optimise_SRC_64-85"), stagenum=stagenum, base_dir=output_dir))

  ### SRC 2
  ###########################

  #Explore all first mutations up to maximum order
  exp_results <- archstabms__optimal_library_design_src(
    mochi_outpath = file.path(base_dir, "Data", "mochi", "ddPCA_SRC"),
    position_offset = 0,
    max_num_mutants = 15,
    explore_positions = 183:204,
    min_num_doubles = 7,
    colour_scheme = colour_scheme,
    outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_library_design_explore_SRC_183-204"), stagenum=stagenum, base_dir=output_dir))

  #Optimised library design
  opt_results <- archstabms__optimal_library_design_src(
    mochi_outpath = file.path(base_dir, "Data", "mochi", "ddPCA_SRC"),
    position_offset = 0,
    max_num_mutants = 15, 
    explore_mutations = exp_results[['optimal_first_mut']],
    explore_positions = 183:204,
    min_num_doubles = 7,
    sample_size = 10000,
    colour_scheme = colour_scheme,
    outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_library_design_optimise_SRC_183-204"), stagenum=stagenum, base_dir=output_dir))

}
