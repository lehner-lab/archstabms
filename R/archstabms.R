
#' archstabms
#'
#' Main analysis script.
#'
#' @param startStage Start at a specified analysis stage (default:1)
#' @param stopStage Stop at a specified analysis stage (default:15)
#' @param base_dir Base directory for all input files (default:NB private CRG server path; change accordingly)
#' @param output_dir Output directory for all output files (default:same as base_dir)
#'
#' @return Nothing
#' @export
archstabms <- function(
	startStage = 1,
	stopStage = 6,
	base_dir = "/users/project/prj004631/afaure/DMS/Results/archstabms_proj",
	output_dir = ""
	){

	# startStage=1
	# stopStage=6
	# base_dir = "/users/project/prj004631/afaure/DMS/Results/archstabms_proj"
	# output_dir = ""

	colour_scheme <- list(
		"shade 0" = list(
			"#F4270C",
			"#F4AD0C",
			"#1B38A6",
			"#09B636"),
		"shade 1" = list(
			"#FFB0A5",
			"#FFE4A5",
			"#9DACE3",
			"#97E9AD"),
		"shade 2" = list(
			"#FF6A56",
			"#FFCB56",
			"#4C63B7",
			"#43C766"),
		"shade 3" = list(
			"#A31300",
			"#A37200",
			"#0C226F",
			"#007A20"),
		"shade 4" = list(
			"#410800",
			"#412D00",
			"#020B2C",
			"#00300D"),
		"shade 5" = list(
			"#CCCCCC",
			"#FF991F",
			"#5CB8FF",
			"#B22222"
		))

	#First and last analysis stages
	first_stage <- startStage
	last_stage <- stopStage

	#Default output directory
	if(output_dir==""){
		output_dir <- base_dir
	}

	### Fit thermo models
	###########################

	# stagenum <- 0
	# archstabms_fit_thermo_model(
	#   base_dir = base_dir,
	#   execute = (first_stage <= stagenum & last_stage >= stagenum))

	### Fit linear models
	###########################

	stagenum <- 0
	# archstabms_fit_linear_model_datasets(
	# 	dataset_names = c(
	# 		"CM1",
	# 		"CM2",
	# 		"CM3",
	# 		"CM4",
	# 		"CM5",
	# 		"CM1345",
	# 		"CM6"),
	# 	seq_position_offset = list(
	# 		"CM1" = 9,
	# 		"CM2" = 25,
	# 		"CM3" = 0,
	# 		"CM4" = 0,
	# 		"CM5" = 2,
	# 		"CM1345" = 0,
	# 		"CM6" = 0),
	#   base_dir = base_dir,
	# 	output_dir = output_dir,
	# 	stagenum = stagenum,
	#   execute = (first_stage <= stagenum & last_stage >= stagenum))

	# archstabms_ddPCA_singlebackground_linear_model_datasets(
	# 	dataset_names = c(
	# 		"CM1",
	# 		"CM2",
	# 		"CM3",
	# 		"CM4",
	# 		"CM5",
	# 		"CM1345",
	# 		"CM6"),
	# 	seq_position_offset = list(
	# 		"CM1" = 9,
	# 		"CM2" = 25,
	# 		"CM3" = 0,
	# 		"CM4" = 0,
	# 		"CM5" = 2,
	# 		"CM1345" = 0,
	# 		"CM6" = 0),
	#   base_dir = base_dir,
	# 	output_dir = output_dir,
	# 	stagenum = stagenum,
	#   execute = (first_stage <= stagenum & last_stage >= stagenum))

	# archstabms_ddPCA_linear_model_datasets(
	# 	dataset_names = c(
	# 		"CM1",
	# 		"CM2",
	# 		"CM3",
	# 		"CM4",
	# 		"CM5",
	# 		"CM1345",
	# 		"CM6"),
	# 	seq_position_offset = list(
	# 		"CM1" = 9,
	# 		"CM2" = 25,
	# 		"CM3" = 0,
	# 		"CM4" = 0,
	# 		"CM5" = 2,
	# 		"CM1345" = 0,
	# 		"CM6" = 0),
	#   base_dir = base_dir,
	# 	output_dir = output_dir,
	# 	stagenum = stagenum,
	#   execute = (first_stage <= stagenum & last_stage >= stagenum))

	archstabms_ddPCA_biophysical_model_datasets(
		dataset_names = c(
			"CM1",
			"CM2",
			"CM3",
			"CM4",
			"CM5",
			"CM1345",
			"CM6"),
		seq_position_offset = list(
			"CM1" = 9,
			"CM2" = 25,
			"CM3" = 0,
			"CM4" = 0,
			"CM5" = 2,
			"CM1345" = 0,
			"CM6" = 0),
	  base_dir = base_dir,
		output_dir = output_dir,
		stagenum = stagenum,
	  execute = (first_stage <= stagenum & last_stage >= stagenum))

	### Evaluate thermo model results
	###########################

	stagenum <- 1
	archstabms_thermo_model_results_datasets(
		dataset_names = c(
			"CM1",
			"CM1binding",
			"CM2",
			"CM3",
			"CM4",
			"CM5",
			"CM1345",
			"CM6"),
		literature_free_energies = file.path(base_dir, "Data", "in_vitro", "GRB2_literature_free_energies.txt"),
		ddPCA_free_energies = file.path(base_dir, "Data", "mochi", "ddPCA", "mochi__fit_tmodel_3state_sparse_dimsum128", "model_weights_0.txt"),
		base_dir = base_dir,
		output_dir = output_dir,
		stagenum = stagenum,
		colour_scheme = colour_scheme,
		execute = (first_stage <= stagenum & last_stage >= stagenum))

	### Add 3D structure metrics
	###########################

	stagenum <- 2
	archstabms_structure_metrics_datasets(
		dataset_names = c(
			"CM1",
			"CM1binding",
			"CM2",
			"CM3",
			"CM4",
			"CM5",
			"CM1345",
			"CM6"),
		pdb_file = file.path(base_dir, "Data", "pdb", "2vwf.pdb"),
		base_dir = base_dir,
		output_dir = output_dir,
		stagenum = stagenum,
		execute = (first_stage <= stagenum & last_stage >= stagenum))

	### Coupling scatterplots
	###########################

	stagenum <- 3
	archstabms_coupling_scatter_datasets(
		dataset_names = c(
			"CM1",
			"CM1binding",
			"CM2",
			"CM3",
			"CM4",
			"CM5",
			"CM1345",
			"CM6"),
		base_dir = base_dir,
		output_dir = output_dir,
		stagenum = stagenum,
		colour_scheme = colour_scheme,
		execute = (first_stage <= stagenum & last_stage >= stagenum))

	### Fitness plots
	###########################

	stagenum <- 4
	archstabms_fitness_plots_datasets(
		dataset_names = c(
			"CM1",
			"CM1binding",
			# "CM2",
			# "CM3",
			# "CM4",
			# "CM5",
			# "CM1345",
			"CM6"),
		base_dir = base_dir,
		output_dir = output_dir,
		stagenum = stagenum,
		colour_scheme = colour_scheme,
		execute = (first_stage <= stagenum & last_stage >= stagenum))

	### Coupling heatmaps
	###########################

	stagenum <- 5
	archstabms_coupling_heatmap_datasets(
		dataset_names = c(
			"CM1",
			"CM1binding",
			"CM2",
			"CM3",
			"CM4",
			"CM5",
			"CM1345",
			"CM6"),
		base_dir = base_dir,
		output_dir = output_dir,
		stagenum = stagenum,
		execute = (first_stage <= stagenum & last_stage >= stagenum))

	### Coupling predictive decay
	###########################

	stagenum <- 6
	archstabms_coupling_decay_datasets(
		dataset_names = c(
			# "CM1",
			# "CM2",
			# "CM3",
			# "CM4",
			# "CM5",
			# "CM1345",
			"CM6"),
		base_dir = base_dir,
		output_dir = output_dir,
		stagenum = stagenum,
		colour_scheme = colour_scheme,
		execute = (first_stage <= stagenum & last_stage >= stagenum))

	### Plot fitness heatmaps
	###########################

	# stagenum <- 4
	# #GB1
	# archstabms_fitness_heatmaps(
	#   input_file = file.path(base_dir, paste0("002", "_archstabms_structure_metrics_GB1"), "dg_singles.txt"),
	#   input_file_fitness = file.path(base_dir, "Data", "fitness", "GB1"),
	#   domain_name = "GB1",
	#   outpath = archstabms__format_dir(dir_suffix="_archstabms_fitness_heatmaps_GB1", stagenum=stagenum, base_dir=output_dir),
	#   colour_scheme = colour_scheme,
	#   plot_width = 11,
	#   plot_traits = "Binding",
	#   execute = (first_stage <= stagenum & last_stage >= stagenum))
	# #PSD95-PDZ3
	# archstabms_fitness_heatmaps(
	#   input_file = file.path(base_dir, paste0("002", "_archstabms_structure_metrics_PSD95-PDZ3"), "dg_singles.txt"),
	#   input_file_fitness = file.path(base_dir, "Data", "fitness", "PSD95-PDZ3"),
	#   input_file_MSA = file.path(base_dir, "Data", "MSA", "PSD95-PDZ3", "frequencies.csv"),
	#   domain_name = "PSD95 PDZ3",
	#   outpath = archstabms__format_dir(dir_suffix="_archstabms_fitness_heatmaps_PSD95-PDZ3", stagenum=stagenum, base_dir=output_dir),
	#   colour_scheme = colour_scheme,
	#   plot_width = 15,
	#   execute = (first_stage <= stagenum & last_stage >= stagenum))
	# #GRB2-SH3
	# archstabms_fitness_heatmaps(
	#   input_file = file.path(base_dir, paste0("002", "_archstabms_structure_metrics_GRB2-SH3"), "dg_singles.txt"),
	#   input_file_fitness = file.path(base_dir, "Data", "fitness", "GRB2-SH3"),
	#   input_file_MSA = file.path(base_dir, "Data", "MSA", "GRB2-SH3", "frequencies.csv"),
	#   domain_name = "GRB2 SH3",
	#   outpath = archstabms__format_dir(dir_suffix="_archstabms_fitness_heatmaps_GRB2-SH3", stagenum=stagenum, base_dir=output_dir),
	#   colour_scheme = colour_scheme,
	#   plot_width = 11,
	#   execute = (first_stage <= stagenum & last_stage >= stagenum))

	### Plot free energy scatterplots
	###########################

	# stagenum <- 5
	# #All domains
	# archstabms_free_energy_scatterplots(
	#   input_list = list(
	#     "GB1" = file.path(base_dir, paste0("002", "_archstabms_structure_metrics_GB1"), "dg_singles.txt"),
	#     "PSD95-PDZ3" = file.path(base_dir, paste0("002", "_archstabms_structure_metrics_PSD95-PDZ3"), "dg_singles.txt"),
	#     "GRB2-SH3" = file.path(base_dir, paste0("002", "_archstabms_structure_metrics_GRB2-SH3"), "dg_singles.txt")),
	#   input_MSA_list = list(
	#     "GB1" = file.path(base_dir, "Data", "MSA", "GB1", "frequencies.csv"),
	#     "PSD95-PDZ3" = file.path(base_dir, "Data", "MSA", "PSD95-PDZ3", "frequencies.csv"),
	#     "GRB2-SH3" = file.path(base_dir, "Data", "MSA", "GRB2-SH3", "frequencies.csv")
	#   ),
	#   outpath = archstabms__format_dir(dir_suffix="_archstabms_free_energy_scatterplots", stagenum=stagenum, base_dir=output_dir),
	#   colour_scheme = colour_scheme,
	#   execute = (first_stage <= stagenum & last_stage >= stagenum))

	### Plot free energy heatmaps
	###########################

	# stagenum <- 6
	# #GB1
	# archstabms_free_energy_heatmaps(
	#   input_file = file.path(base_dir, paste0("002", "_archstabms_structure_metrics_GB1"), "dg_singles.txt"),
	#   domain_name = "GB1",
	#   outpath = archstabms__format_dir(dir_suffix="_archstabms_free_energy_heatmaps_GB1", stagenum=stagenum, base_dir=output_dir),
	#   colour_scheme = colour_scheme,
	#   plot_width = 11,
	#   execute = (first_stage <= stagenum & last_stage >= stagenum))
	# #PSD95-PDZ3
	# archstabms_free_energy_heatmaps(
	#   input_file = file.path(base_dir, paste0("002", "_archstabms_structure_metrics_PSD95-PDZ3"), "dg_singles.txt"),
	#   domain_name = "PSD95 PDZ3",
	#   outpath = archstabms__format_dir(dir_suffix="_archstabms_free_energy_heatmaps_PSD95-PDZ3", stagenum=stagenum, base_dir=output_dir),
	#   colour_scheme = colour_scheme,
	#   plot_width = 15,
	#   execute = (first_stage <= stagenum & last_stage >= stagenum))
	# #GRB2-SH3
	# archstabms_free_energy_heatmaps(
	#   input_file = file.path(base_dir, paste0("002", "_archstabms_structure_metrics_GRB2-SH3"), "dg_singles.txt"),
	#   domain_name = "GRB2 SH3",
	#   outpath = archstabms__format_dir(dir_suffix="_archstabms_free_energy_heatmaps_GRB2-SH3", stagenum=stagenum, base_dir=output_dir),
	#   colour_scheme = colour_scheme,
	#   plot_width = 11,
	#   execute = (first_stage <= stagenum & last_stage >= stagenum))

}