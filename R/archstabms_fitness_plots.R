
#' archstabms_fitness_plots
#'
#' Plot fitness distributions and scatterplots.
#'
#' @param dataset_name character dataset name (required)
#' @param base_dir Base directory (required)
#' @param colour_scheme colour scheme file (required)
#' @param outpath output path for plots and saved objects (required)
#'
#' @return Nothing
#' @export
#' @import data.table
archstabms_fitness_plots <- function(
  dataset_name,
  base_dir,
  colour_scheme,
  outpath
  ){

  #Display status
  message(paste("\n\n*******", "running stage: archstabms_fitness_plots", "*******\n\n"))

  #Create output directory
  archstabms__create_dir(archstabms_dir = outpath)

  #Report path
  report_outpath <- outpath

  ### Plot abundance fitness distributions, replicate correlations, violins
  ###########################

  #Load abundance fitness data
  mochi_outpath <- file.path(base_dir, "Data", "mochi", paste0(dataset_name, "o1"))
  fitness_dt <- fread(file.path(mochi_outpath, 'task_1', 'predictions', "predicted_phenotypes_all.txt"))
  # #Save for supplement
  # write.table(fitness_dt[Nham_aa==Nmut_codons & STOP==F & STOP_readthrough==F,.SD,,.SDcols = c("protein", "pca_type", "aa_seq", "Nham_aa", "WT", "fitness", "sigma", "growthrate", "growthrate_sigma")], 
  #   file = file.path(report_outpath, "fitness_supp.txt"), 
  #   quote = F, sep = "\t", row.names = F)

  #Plot mutation distributions
  dataset_names <- names(fitness_dt)[grepl("Abundance_", names(fitness_dt))]
  for(i in dataset_names){
    archstabms__plot_mutation_distributions(
      input_dt = fitness_dt[get(i)==1,], 
      report_outpath = outpath,
      dataset_name = gsub("Abundance_", "", i),
      colour_scheme = colour_scheme)
  }

  #Plot fitness replicate correlation
  archstabms__plot_fitness_replicates_cor(
    input_dt = fitness_dt[phenotype==1],
    output_file = file.path(outpath, "replicate_scatter_binhex_abundance.pdf"),
    colour_scheme = colour_scheme)

  #Plot abundance growthrate violins with hamming distance
  dataset_names <- names(fitness_dt)[grepl("Abundance_|Binding_", names(fitness_dt))]
  for(i in dataset_names){
    archstabms__plot_growthrate_violins(
      input_dt = fitness_dt[get(i)==1,], 
      report_outpath = outpath,
      dataset_name = gsub("Abundance_|Binding_", "", i),
      trait = sapply(strsplit(i, "_"), '[', 1),
      colour_scheme = colour_scheme)
  }

  ### Plot abundance fitness distributions, replicate correlations, violins - lenient
  ###########################

  if(dir.exists(file.path(base_dir, "Data", "fitness", dataset_name))){
    #Load abundance fitness data
    fitness_filename <- list.files(file.path(base_dir, "Data", "fitness", dataset_name, ""), pattern = 'RData$')[1]
    load(file.path(base_dir, "Data", "fitness", dataset_name, fitness_filename))
    fitness_dt <- all_variants

    #Plot mutation distributions
    archstabms__plot_mutation_distributions(
      input_dt = fitness_dt, 
      report_outpath = outpath,
      dataset_name = paste0(gsub("Abundance_", "", dataset_name), "_lenient"),
      colour_scheme = colour_scheme)

    #Plot fitness replicate correlation
    archstabms__plot_fitness_replicates_cor(
      input_dt = fitness_dt,
      output_file = file.path(outpath, "replicate_scatter_binhex_abundance_lenient.pdf"),
      colour_scheme = colour_scheme)

    #Plot abundance growthrate violins with hamming distance
    dataset_names <- names(fitness_dt)[grepl("Abundance_|Binding_", names(fitness_dt))]
    for(i in dataset_names){
      archstabms__plot_growthrate_violins(
        input_dt = fitness_dt[get(i)==1,], 
        report_outpath = outpath,
        dataset_name = gsub("Abundance_", "", i),
        trait = sapply(strsplit(i, "_"), '[', 1),
        suffix = "_lenient",
        colour_scheme = colour_scheme)
    }
  }

  ### Plot binding fitness replicate correlations
  ###########################

  #Return if no binding data
  if(dataset_name!="CM1binding"){
    return()
  }

  #Load abundance and binding fitness data
  fitness_dt <- fread(file.path(mochi_outpath, 'task_1', 'predictions', "predicted_phenotypes_all.txt"))

  archstabms__plot_fitness_replicates_cor(
    input_dt = fitness_dt[phenotype==2],
    output_file = file.path(outpath, "replicate_scatter_binhex_binding.pdf"),
    colour_scheme = colour_scheme)

  ### Plot binding vs folding fitness
  ###########################

  #File paths
  coef_file <- file.path(base_dir, paste0("002_archstabms_structure_metrics_", dataset_name, "o1"), 'model_coefficients.txt')

  #Return if file doesn't exist
  if(!file.exists(coef_file)){
    return()
  }

  #Load fitness
  fitness_dt <- fread(file.path(mochi_outpath, 'task_1', 'predictions', "predicted_phenotypes_all.txt"))
  a_dt <- fitness_dt[phenotype==1]
  b_dt <- fitness_dt[phenotype==2]

  ab_dt <- merge(
    a_dt[,.(aa_seq, mut_order = Nham_aa, abundance_fitness = fitness, abundance_growthrate = growthrate)],
    b_dt[,.(aa_seq, mut_order = Nham_aa, binding_fitness = fitness, binding_growthrate = growthrate)],
    by = c('aa_seq', 'mut_order'))

  #Growth rate scatter
  plot_dt <- ab_dt[mut_order<15,.(mut_order, abundance_growthrate, binding_growthrate)]
  d <- ggplot2::ggplot(plot_dt[mut_order>0],ggplot2::aes(abundance_growthrate, binding_growthrate)) +
    # ggplot2::stat_binhex(bins = 100, size = 0, color = "lightgrey") +
    # ggplot2::scale_fill_gradientn(colours = c("white", "black"), trans = "log10") +
    ggplot2::geom_point(alpha = 1/10, ggplot2::aes(color = as.factor(mut_order))) +
    ggplot2::xlab("Abundance growth rate") +
    ggplot2::ylab("Binding growth rate") +
    ggplot2::geom_text(data = plot_dt[mut_order>0,.(label = paste("R-squared = ", round(cor(abundance_growthrate, binding_growthrate, use = "pairwise.complete")^2, 2), sep="")),mut_order], ggplot2::aes(label=label, x=-Inf, y=Inf, hjust = 0, vjust = 1)) +
    ggplot2::facet_wrap(~mut_order, nrow = 3) + 
    ggplot2::geom_hline(yintercept = plot_dt[mut_order==0,binding_growthrate], linetype = 2) +
    ggplot2::geom_vline(xintercept = plot_dt[mut_order==0,abundance_growthrate], linetype = 2) +
    ggplot2::theme_classic()
  d <- d + ggplot2::geom_abline(linetype = 2)
  ggplot2::ggsave(file.path(report_outpath, "growthrate_scatter_binding_vs_abundance.pdf"), d, width = 15, height = 9, useDingbats=FALSE)

  #Fitness scatter
  plot_dt <- ab_dt[mut_order<15,.(mut_order, abundance_fitness, binding_fitness)]
  d <- ggplot2::ggplot(plot_dt[mut_order>0],ggplot2::aes(abundance_fitness, binding_fitness)) +
    # ggplot2::stat_binhex(bins = 100, size = 0, color = "lightgrey") +
    # ggplot2::scale_fill_gradientn(colours = c("white", "black"), trans = "log10") +
    ggplot2::geom_point(alpha = 1/10, ggplot2::aes(color = as.factor(mut_order))) +
    ggplot2::xlab("Abundance fitness") +
    ggplot2::ylab("Binding fitness") +
    ggplot2::geom_text(data = plot_dt[mut_order>0,.(label = paste("R-squared = ", round(cor(abundance_fitness, binding_fitness, use = "pairwise.complete")^2, 2), sep="")),mut_order], ggplot2::aes(label=label, x=-Inf, y=Inf, hjust = 0, vjust = 1)) +
    ggplot2::facet_wrap(~mut_order, nrow = 3) + 
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::theme_classic()
  d <- d + ggplot2::geom_abline(linetype = 2)
  ggplot2::ggsave(file.path(report_outpath, "fitness_scatter_binding_vs_abundance.pdf"), d, width = 15, height = 9, useDingbats=FALSE)

  #Folding/binding energy scatter
  coef_dt <- fread(coef_file)
  plot_dt <- merge(
    coef_dt[trait_name=='Folding',.(id_ref, ddg = .SD[[2]], ddg_ci = .SD[[3]]),,.SDcols = c("id_ref", "mean_kcal/mol", "ci95_kcal/mol")],
    coef_dt[trait_name=='Binding',.(id_ref, ddg = .SD[[2]], ddg_ci = .SD[[3]]),,.SDcols = c("id_ref", "mean_kcal/mol", "ci95_kcal/mol")], 
    by = "id_ref", suffixes = c("_folding", "_binding"))
  d <- ggplot2::ggplot(plot_dt[id_ref!="WT"],ggplot2::aes(ddg_folding, ddg_binding)) +
    ggplot2::geom_point() +
    ggplot2::geom_linerange(ggplot2::aes(ymin = ddg_binding-ddg_ci_binding/2, ymax = ddg_binding+ddg_ci_binding/2), alpha = 1/4) +
    ggplot2::geom_linerange(ggplot2::aes(xmin = ddg_folding-ddg_ci_folding/2, xmax = ddg_folding+ddg_ci_folding/2), alpha = 1/4) +
    ggplot2::xlab("Folding free energy change") +
    ggplot2::ylab("Binding free energy change") +
    ggplot2::geom_text(data = plot_dt[id_ref!="WT",.(label = paste("r = ", round(cor(ddg_folding, ddg_binding, use = "pairwise.complete"), 2), sep=""))], ggplot2::aes(label=label, x=-Inf, y=Inf, hjust = 0, vjust = 1)) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_vline(xintercept = 0) +
    # ggplot2::coord_fixed() +
    ggplot2::theme_classic()
  d <- d + ggplot2::geom_abline(linetype = 2)
  ggplot2::ggsave(file.path(report_outpath, "free_energy_scatter_binding_vs_folding.pdf"), d, width = 4, height = 4, useDingbats=FALSE)

  #Binding free energy change vs. ligand distance
  coef_dt <- fread(coef_file)
  plot_dt <- coef_dt[trait_name=='Binding']
  plot_dt[, ddg_abs := abs(.SD[[1]]),,.SDcols = "mean_kcal/mol"]
  plot_dt[, ddg_abs_ci := .SD[[1]],,.SDcols = "ci95_kcal/mol"]
  d <- ggplot2::ggplot(plot_dt[id_ref!="WT"],ggplot2::aes(scHAmin_ligand, ddg_abs)) +
    ggplot2::geom_point() +
    ggplot2::geom_linerange(ggplot2::aes(ymin = ddg_abs-ddg_abs_ci/2, ymax = ddg_abs+ddg_abs_ci/2), alpha = 1/4) +
    ggplot2::xlab("Distance to ligand (Angstrom)") +
    ggplot2::ylab("|Binding free energy change|") +
    ggplot2::geom_text(data = plot_dt[id_ref!="WT",.(label = paste("r = ", round(cor(scHAmin_ligand, ddg_abs, use = "pairwise.complete"), 2), sep=""))], ggplot2::aes(label=label, x=-Inf, y=Inf, hjust = 0, vjust = 1)) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::geom_vline(xintercept = 5, linetype = 2) +
    # ggplot2::coord_fixed() +
    ggplot2::theme_classic()
  ggplot2::ggsave(file.path(report_outpath, "free_energy_scatter_binding_vs_ligand_distance.pdf"), d, width = 4, height = 4, useDingbats=FALSE)

}