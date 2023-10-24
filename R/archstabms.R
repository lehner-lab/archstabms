
#' archstabms
#'
#' Main analysis script.
#'
#' @param startStage Start at a specified analysis stage (default:2)
#' @param stopStage Stop at a specified analysis stage (default:7)
#' @param base_dir Base directory for all input files (default:NB private CRG server path; change accordingly)
#' @param output_dir Output directory for all output files (default:same as base_dir)
#'
#' @return Nothing
#' @export
archstabms <- function(
	startStage = 2,
	stopStage = 7,
	base_dir = "/users/project/prj004631/afaure/DMS/Results/archstabms_proj",
	output_dir = ""
	){

	# startStage=5
	# stopStage=5
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

	### Library design
	###########################

	stagenum <- 0
	archstabms_library_design(
		base_dir = base_dir,
		output_dir = output_dir,
		stagenum = stagenum,
		colour_scheme = colour_scheme,
		execute = (first_stage <= stagenum & last_stage >= stagenum))

	### Fit linear models Fit linear models
	###########################

	stagenum <- 1
	archstabms_fit_linear_model_datasets(
		dataset_names = c(
			"CM1",
			"CM2",
			"CM6"),
		seq_position_offset = list(
			"CM1" = 9,
			"CM2" = 25,
			"CM6" = 0),
	  base_dir = base_dir,
		output_dir = output_dir,
		stagenum = stagenum,
	  execute = (first_stage <= stagenum & last_stage >= stagenum))

	archstabms_ddPCA_singlebackground_linear_model_datasets(
		dataset_names = c(
			"CM1",
			"CM2",
			"CM6"),
		seq_position_offset = list(
			"CM1" = 9,
			"CM2" = 25,
			"CM6" = 0),
	  base_dir = base_dir,
		output_dir = output_dir,
		stagenum = stagenum,
	  execute = (first_stage <= stagenum & last_stage >= stagenum))

	archstabms_ddPCA_linear_model_datasets(
		dataset_names = c(
			"CM1",
			"CM2",
			"CM6"),
		seq_position_offset = list(
			"CM1" = 9,
			"CM2" = 25,
			"CM6" = 0),
	  base_dir = base_dir,
		output_dir = output_dir,
		stagenum = stagenum,
	  execute = (first_stage <= stagenum & last_stage >= stagenum))

	archstabms_ddPCA_biophysical_model_datasets(
		dataset_names = c(
			"CM1",
			"CM2",
			"CM6"),
		seq_position_offset = list(
			"CM1" = 9,
			"CM2" = 25,
			"CM6" = 0),
	  base_dir = base_dir,
		output_dir = output_dir,
		stagenum = stagenum,
	  execute = (first_stage <= stagenum & last_stage >= stagenum))

	### Evaluate thermo model results
	###########################

	stagenum <- 2
	archstabms_thermo_model_results_datasets(
		dataset_names = c(
			"CM1",
			"CM1binding",
			"CM2",
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

	stagenum <- 3
	archstabms_structure_metrics_datasets(
		dataset_names = c(
			"CM1",
			"CM1binding",
			"CM2",
			"CM6"),
		pdb_file = file.path(base_dir, "Data", "pdb", "2vwf.pdb"),
		base_dir = base_dir,
		output_dir = output_dir,
		stagenum = stagenum,
		execute = (first_stage <= stagenum & last_stage >= stagenum))

	### Coupling scatterplots
	###########################

	stagenum <- 4
	archstabms_coupling_scatter_datasets(
		dataset_names = c(
			"CM1",
			"CM1binding",
			"CM2",
			"CM6"),
		base_dir = base_dir,
		output_dir = output_dir,
		stagenum = stagenum,
		colour_scheme = colour_scheme,
		execute = (first_stage <= stagenum & last_stage >= stagenum))

	### Fitness plots
	###########################

	stagenum <- 5
	archstabms_fitness_plots_datasets(
		dataset_names = c(
			"CM1",
			"CM1binding",
			"CM2",
			"CM6"),
		base_dir = base_dir,
		output_dir = output_dir,
		stagenum = stagenum,
		colour_scheme = colour_scheme,
		execute = (first_stage <= stagenum & last_stage >= stagenum))

	### Coupling heatmaps
	###########################

	stagenum <- 6
	archstabms_coupling_heatmap_datasets(
		dataset_names = c(
			"CM1",
			"CM1binding",
			"CM2",
			"CM6"),
		base_dir = base_dir,
		output_dir = output_dir,
		stagenum = stagenum,
		execute = (first_stage <= stagenum & last_stage >= stagenum))

	### Binding fitness plots
	###########################

	stagenum <- 7
	archstabms_binding_plots_datasets(
		dataset_names = c(
			"CM1binding"),
		base_dir = base_dir,
		output_dir = output_dir,
		stagenum = stagenum,
		colour_scheme = colour_scheme,
		execute = (first_stage <= stagenum & last_stage >= stagenum))

}

