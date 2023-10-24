
#' archstabms_binding_plots
#'
#' Plot fitness distributions and scatterplots.
#'
#' @param dataset_name character dataset name (required)
#' @param coef_order coefficient order (required)
#' @param base_dir Base directory (required)
#' @param colour_scheme colour scheme file (required)
#' @param outpath output path for plots and saved objects (required)
#'
#' @return Nothing
#' @export
#' @import data.table
archstabms_binding_plots <- function(
  dataset_name,
  coef_order,
  base_dir,
  colour_scheme,
  outpath
  ){

  #Display status
  message(paste("\n\n*******", "running stage: archstabms_binding_plots", "*******\n\n"))

  #Create output directory
  archstabms__create_dir(archstabms_dir = outpath)

  #Report path
  report_outpath <- outpath

  ### Plot binding fitness replicate correlations
  ###########################

  #Load abundance and binding fitness data
  mochi_outpath <- file.path(base_dir, "Data", "mochi", paste0(dataset_name, "o", coef_order))
  fitness_dt <- fread(file.path(mochi_outpath, 'task_1', 'predictions', "predicted_phenotypes_all.txt"))

  ### Plot binding vs folding fitness
  ###########################

  #File paths
  coef_file <- file.path(base_dir, paste0("003_archstabms_structure_metrics_", dataset_name, "o", coef_order), 'model_coefficients.txt')

  #Return if file doesn't exist
  if(!file.exists(coef_file)){
    return()
  }

  #Load fitness
  fitness_dt <- fread(file.path(mochi_outpath, 'task_1', 'predictions', "predicted_phenotypes_all.txt"))
  a_dt <- fitness_dt[phenotype==1]
  b_dt <- fitness_dt[phenotype==2]
  dataset_names <- names(fitness_dt)[grepl("Abundance_|Binding_", names(fitness_dt))]

  ab_dt <- merge(
    a_dt[,.(aa_seq, mut_order = Nham_aa, abundance_fitness = fitness, abundance_growthrate = growthrate)],
    b_dt[,.(aa_seq, mut_order = Nham_aa, binding_fitness = fitness, binding_growthrate = growthrate)],
    by = c('aa_seq', 'mut_order'))

  #Growth rate thresholds
  gr_thresholds <- list()
  for(i in dataset_names){
    linears_dt <- fread(file.path(mochi_outpath, 'task_1', 'weights', paste0("linears_weights_", i, ".txt")))
    gr_lm <- lm(growthrate ~ fitness, data = fitness_dt[get(i)==1])
    gr_thresholds[[i]] <- as.numeric(predict(gr_lm, newdata = data.frame(fitness = linears_dt[, mean(0.5*kernel+bias)])))
  }

  #9th order variants indistinguishable from WT abundance (nominal P=0.05)
  wtl9 <- a_dt[Nham_aa==9 & fitness+1.96*sigma>a_dt[WT==T,fitness],aa_seq]
  #...and bound (predicted fraction bound >50%)
  wtl9b <- b_dt[aa_seq %in% wtl9,aa_seq]
  wtl9bb <- b_dt[aa_seq %in% wtl9 & growthrate > gr_thresholds[["Binding_CM1"]],aa_seq]
  print(paste0("Percentage variants with >9 substitutions and WT-like fitness (", dataset_name, "): ", length(wtl9), "/", a_dt[Nham_aa==9,.N], "*100=", round(length(wtl9)/a_dt[Nham_aa==9,.N]*100)))
  print(paste0("Percentage variants with >9 substitutions and WT-like fitness also bound (", dataset_name, "): ", length(wtl9bb), "/", length(wtl9b), "*100=", round(length(wtl9bb)/length(wtl9b)*100)))

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
    ggplot2::stat_binhex(bins = 50, size = 0, color = "lightgrey") +
    ggplot2::scale_fill_viridis_c(trans = "log10") +
    ggplot2::xlab("Abundance fitness") +
    ggplot2::ylab("Binding fitness") +
    ggplot2::facet_wrap(~mut_order, nrow = 3) + 
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::theme_classic()
  d <- d + ggplot2::geom_abline(linetype = 2)
  ggplot2::ggsave(file.path(report_outpath, "fitness_scatter_binding_vs_abundance.pdf"), d, width = 15, height = 9, useDingbats=FALSE)

  #Fitness scatter - no facet
  plot_dt <- ab_dt[,.(mut_order, abundance_fitness, binding_fitness)]
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(abundance_fitness, binding_fitness)) +
    ggplot2::stat_binhex(bins = 100, size = 0, color = "lightgrey") +
    ggplot2::scale_fill_viridis_c(trans = "log10") +
    ggplot2::xlab("Abundance fitness") +
    ggplot2::ylab("Binding fitness") +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::theme_classic()
  d <- d + ggplot2::geom_abline(linetype = 2)
  ggplot2::ggsave(file.path(report_outpath, "fitness_scatter_binding_vs_abundance_nofacet.pdf"), d, width = 5, height = 4, useDingbats=FALSE)

  #Fitness scatter - 3,6,9 mutation orders
  plot_dt <- ab_dt[mut_order %in% c(3,6,9),.(mut_order, abundance_fitness, binding_fitness)]
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(abundance_fitness, binding_fitness)) +
    # ggplot2::stat_binhex(bins = 100, size = 0, color = "lightgrey") +
    # ggplot2::scale_fill_gradientn(colours = c("white", "black"), trans = "log10") +
    ggplot2::stat_binhex(bins = 50, size = 0, color = "lightgrey") +
    # ggplot2::scale_fill_gradientn(colours = c("white", "black"), trans = "log10") +
    ggplot2::scale_fill_viridis_c(trans = "log10") +
    ggplot2::xlab("Abundance fitness") +
    ggplot2::ylab("Binding fitness") +
    # ggplot2::geom_text(data = plot_dt[,.(label = paste("R-squared = ", round(cor(abundance_fitness, binding_fitness, use = "pairwise.complete")^2, 2), sep="")),mut_order], ggplot2::aes(label=label, x=-Inf, y=Inf, hjust = 0, vjust = 1)) +
    ggplot2::facet_wrap(~mut_order, ncol = 3) + 
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::theme_classic()
  d <- d + ggplot2::geom_abline(linetype = 2)
  ggplot2::ggsave(file.path(report_outpath, "fitness_scatter_binding_vs_abundance_order369.pdf"), d, width = 8, height = 3, useDingbats=FALSE)

  #Fitness scatter - 4,8,12 mutation orders
  plot_dt <- ab_dt[mut_order %in% c(4,8,12),.(mut_order, abundance_fitness, binding_fitness)]
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(abundance_fitness, binding_fitness)) +
    # ggplot2::stat_binhex(bins = 100, size = 0, color = "lightgrey") +
    # ggplot2::scale_fill_gradientn(colours = c("white", "black"), trans = "log10") +
    ggplot2::stat_binhex(bins = 50, size = 0, color = "lightgrey") +
    # ggplot2::scale_fill_gradientn(colours = c("white", "black"), trans = "log10") +
    ggplot2::scale_fill_viridis_c(trans = "log10") +
    ggplot2::xlab("Abundance fitness") +
    ggplot2::ylab("Binding fitness") +
    # ggplot2::geom_text(data = plot_dt[,.(label = paste("R-squared = ", round(cor(abundance_fitness, binding_fitness, use = "pairwise.complete")^2, 2), sep="")),mut_order], ggplot2::aes(label=label, x=-Inf, y=Inf, hjust = 0, vjust = 1)) +
    ggplot2::facet_wrap(~mut_order, ncol = 3) + 
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::theme_classic()
  d <- d + ggplot2::geom_abline(linetype = 2)
  ggplot2::ggsave(file.path(report_outpath, "fitness_scatter_binding_vs_abundance_order4812.pdf"), d, width = 8, height = 3, useDingbats=FALSE)

  #Folding/binding energy scatter
  coef_dt <- fread(coef_file)[coef_order<=1]
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
    ggrepel::geom_text_repel(ggplot2::aes(label = id_ref), show.legend = T, 
      max.overlaps = Inf) +
    # ggplot2::coord_fixed() +
    ggplot2::theme_classic()
  d <- d + ggplot2::geom_abline(linetype = 2)
  ggplot2::ggsave(file.path(report_outpath, "free_energy_scatter_binding_vs_folding.pdf"), d, width = 4, height = 4, useDingbats=FALSE)

  #Folding/binding energy scatter - order 2
  if(coef_order==2){
    coef_dt <- fread(coef_file)[coef_order==2]

    #Test whether folding ddGs larger than binding ddGs
    test_results <- archstabms__mann_whitney_U_wrapper(
      coef_dt[trait_name=="Folding", abs(.SD[[1]]),,.SDcols = "mean_kcal/mol"],
      coef_dt[trait_name=="Binding", abs(.SD[[1]]),,.SDcols = "mean_kcal/mol"])
    print(paste0("Magnitude of folding ddGs larger than binding ddGs (", dataset_name, "): Mann-Whitney U test P=", test_results[['p_value']], ", AUC=", test_results[['effect_size']], ", n=", coef_dt[,.N]))

    plot_dt <- merge(
      coef_dt[trait_name=='Folding',.(id_ref, ddg = .SD[[2]], ddg_ci = .SD[[3]]),,.SDcols = c("id_ref", "mean_kcal/mol", "ci95_kcal/mol")],
      coef_dt[trait_name=='Binding',.(id_ref, ddg = .SD[[2]], ddg_ci = .SD[[3]]),,.SDcols = c("id_ref", "mean_kcal/mol", "ci95_kcal/mol")], 
      by = "id_ref", suffixes = c("_folding", "_binding"))
    d <- ggplot2::ggplot(plot_dt,ggplot2::aes(ddg_folding, ddg_binding)) +
      ggplot2::geom_point() +
      ggplot2::geom_linerange(ggplot2::aes(ymin = ddg_binding-ddg_ci_binding/2, ymax = ddg_binding+ddg_ci_binding/2), alpha = 1/4) +
      ggplot2::geom_linerange(ggplot2::aes(xmin = ddg_folding-ddg_ci_folding/2, xmax = ddg_folding+ddg_ci_folding/2), alpha = 1/4) +
      ggplot2::xlab("Folding free energy change") +
      ggplot2::ylab("Binding free energy change") +
      ggplot2::geom_text(data = plot_dt[,.(label = paste("r = ", round(cor(ddg_folding, ddg_binding, use = "pairwise.complete"), 2), sep=""))], ggplot2::aes(label=label, x=-Inf, y=Inf, hjust = 0, vjust = 1)) +
      ggplot2::geom_hline(yintercept = 0) +
      ggplot2::geom_vline(xintercept = 0) +
      # ggplot2::coord_fixed() +
      ggplot2::theme_classic()
    d <- d + ggplot2::geom_abline(linetype = 2)
    ggplot2::ggsave(file.path(report_outpath, "free_energy_scatter_binding_vs_folding_order2.pdf"), d, width = 4, height = 4, useDingbats=FALSE)
    #Folding/binding energy scatter - order 2 coloured by secondary structure
    plot_dt <- merge(
      coef_dt[trait_name=='Folding',.(id_ref, ddg = .SD[[2]], ddg_ci = .SD[[3]], SS1, SS2),,.SDcols = c("id_ref", "mean_kcal/mol", "ci95_kcal/mol")],
      coef_dt[trait_name=='Binding',.(id_ref, ddg = .SD[[2]], ddg_ci = .SD[[3]], SS1, SS2),,.SDcols = c("id_ref", "mean_kcal/mol", "ci95_kcal/mol")], 
      by = c("id_ref", "SS1", "SS2"), suffixes = c("_folding", "_binding"))
    plot_dt[, SS := apply(.SD=="sheet", 1, sum, na.rm = T),,.SDcols = c("SS1", "SS2")]
    plot_dt[, SS_col := 'loop-loop']
    plot_dt[SS==1, SS_col := 'strand-loop']
    plot_dt[SS==2, SS_col := 'strand-strand']
    plot_dt[, SS_size := "0"]
    plot_dt[SS!=0, SS_size := "1"]
    plot_cols <- c('black', colour_scheme[["shade 0"]][2], colour_scheme[["shade 0"]][1])
    names(plot_cols) <- c("loop-loop", "strand-loop", "strand-strand")
    plot_sizes <- c(0.5, 1.5)
    names(plot_sizes) <- c(0, 1)
    d <- ggplot2::ggplot(plot_dt,ggplot2::aes(ddg_folding, ddg_binding, color = SS_col)) +
      ggplot2::geom_point(data = plot_dt[SS_size!=0], ggplot2::aes(size = SS_size)) +
      ggplot2::geom_point(data = plot_dt[SS_size==0], ggplot2::aes(size = SS_size)) +
      ggplot2::geom_linerange(ggplot2::aes(ymin = ddg_binding-ddg_ci_binding/2, ymax = ddg_binding+ddg_ci_binding/2), alpha = 1/4) +
      ggplot2::geom_linerange(ggplot2::aes(xmin = ddg_folding-ddg_ci_folding/2, xmax = ddg_folding+ddg_ci_folding/2), alpha = 1/4) +
      ggplot2::xlab("Folding free energy change") +
      ggplot2::ylab("Binding free energy change") +
      ggplot2::geom_text(data = plot_dt[,.(label = paste("r = ", round(cor(ddg_folding, ddg_binding, use = "pairwise.complete"), 2), sep=""), SS_col)], ggplot2::aes(label=label, x=-Inf, y=Inf, hjust = 0, vjust = 1)) +
      ggplot2::geom_hline(yintercept = 0) +
      ggplot2::geom_vline(xintercept = 0) +
      # ggplot2::coord_fixed() +
      ggplot2::scale_colour_manual(values = plot_cols) +
      ggplot2::scale_size_manual(values = plot_sizes) +
      ggplot2::theme_classic()
    d <- d + ggplot2::geom_abline(linetype = 2)
    ggplot2::ggsave(file.path(report_outpath, "free_energy_scatter_binding_vs_folding_order2_SScol.pdf"), d, width = 5, height = 3, useDingbats=FALSE)

    #Coupling scatter - scHAmin ligand distance - position class colour
    for(trait in c("Binding", "Folding")){
      coef_dt <- fread(coef_file)[coef_order==2 & trait_name==trait]
      #Weighted mean absolute value
      coef_dt[coef_order==2, dddg_abs_mut := abs(.SD[[1]]),.SDcols = c("mean_kcal/mol")]
      coef_dt[coef_order==2, dddg_abs := sum(abs(.SD[[1]])/.SD[[2]]^2, na.rm = T)/sum(1/.SD[[2]]^2, na.rm = T),Pos_ref,.SDcols = c("mean_kcal/mol", "std_kcal/mol")]
      coef_dt[coef_order==2, dddg_ci := sqrt(1/sum(1/.SD[[2]]^2, na.rm = T))*1.96*2,Pos_ref,.SDcols = c("mean_kcal/mol", "std_kcal/mol")]
      #Mean
      plot_dt <- coef_dt[!duplicated(Pos_ref),.(dddg_abs, dddg_ci, scHAmin, scHAmin_ligand1, scHAmin_ligand2, SS1, SS2)]
      plot_dt[is.na(dddg_abs), dddg_abs := 0]
      plot_dt[is.na(dddg_ci), dddg_ci := 0]
      plot_dt[, SS := apply(.SD=="sheet", 1, sum, na.rm = T),,.SDcols = c("SS1", "SS2")]
      plot_dt[, SS_col := 'loop-loop']
      plot_dt[SS==1, SS_col := 'strand-loop']
      plot_dt[SS==2, SS_col := 'strand-strand']
      plot_dt[, SS_size := "0"]
      plot_dt[SS!=0, SS_size := "1"]
      plot_cols <- c('black', colour_scheme[["shade 0"]][2], colour_scheme[["shade 0"]][1])
      names(plot_cols) <- c("loop-loop", "strand-loop", "strand-strand")
      plot_sizes <- c(0.5, 1.5)
      names(plot_sizes) <- c(0, 1)
      plot_dt[, scHAmin_ligmean := apply(.SD, 1, mean),,.SDcols = c("scHAmin_ligand1", "scHAmin_ligand2")]
      d <- ggplot2::ggplot(plot_dt,ggplot2::aes(scHAmin_ligmean, dddg_abs)) +
        ggplot2::geom_point(ggplot2::aes(color = SS_col, size = SS_size)) +
        ggplot2::geom_linerange(data = plot_dt, ggplot2::aes(ymin = dddg_abs-dddg_ci/2, ymax = dddg_abs+dddg_ci/2, color = SS_col), alpha = 1/4) +
        ggplot2::geom_vline(xintercept = 5, linetype = 2) + 
        ggplot2::geom_text(data = plot_dt[,.(label = paste("rho = ", round(cor(dddg_abs, scHAmin_ligmean, use = "pairwise.complete", method = 'spearman'), 2), sep=""))], ggplot2::aes(label=label, x=Inf, y=Inf, hjust = 1, vjust = 1)) +
        ggplot2::xlab("Mean ligand distance (Angstrom)") +
        ggplot2::ylab(bquote("|"*.(trait) ~ Delta*Delta*Delta*"G| (kcal/mol)")) +
        ggplot2::theme_classic() +
        ggplot2::scale_colour_manual(values = plot_cols) +
        ggplot2::scale_size_manual(values = plot_sizes)
      ggplot2::ggsave(file.path(report_outpath, paste0("coupling_scatter_scHAminligandmean_scHAmincol_coreshape_", trait, ".pdf")), d, width = 4, height = 3, useDingbats=FALSE)
      #Min
      plot_dt <- coef_dt[!duplicated(Pos_ref),.(dddg_abs, dddg_ci, scHAmin, scHAmin_ligand1, scHAmin_ligand2, SS1, SS2)]
      plot_dt[is.na(dddg_abs), dddg_abs := 0]
      plot_dt[is.na(dddg_ci), dddg_ci := 0]
      plot_dt[, SS := apply(.SD=="sheet", 1, sum, na.rm = T),,.SDcols = c("SS1", "SS2")]
      plot_dt[, SS_col := 'loop-loop']
      plot_dt[SS==1, SS_col := 'strand-loop']
      plot_dt[SS==2, SS_col := 'strand-strand']
      plot_dt[, SS_size := "0"]
      plot_dt[SS!=0, SS_size := "1"]
      plot_cols <- c('black', colour_scheme[["shade 0"]][2], colour_scheme[["shade 0"]][1])
      names(plot_cols) <- c("loop-loop", "strand-loop", "strand-strand")
      plot_sizes <- c(0.5, 1.5)
      names(plot_sizes) <- c(0, 1)      
      plot_dt[, scHAmin_ligmin := apply(.SD, 1, min),,.SDcols = c("scHAmin_ligand1", "scHAmin_ligand2")]
      d <- ggplot2::ggplot(plot_dt,ggplot2::aes(scHAmin_ligmin, dddg_abs)) +
        ggplot2::geom_point(ggplot2::aes(color = SS_col, size = SS_size)) +
        ggplot2::geom_linerange(data = plot_dt, ggplot2::aes(ymin = dddg_abs-dddg_ci/2, ymax = dddg_abs+dddg_ci/2, color = SS_col), alpha = 1/4) +
        ggplot2::geom_vline(xintercept = 5, linetype = 2) + 
        ggplot2::geom_text(data = plot_dt[,.(label = paste("rho = ", round(cor(dddg_abs, scHAmin_ligmin, use = "pairwise.complete", method = 'spearman'), 2), sep=""))], ggplot2::aes(label=label, x=Inf, y=Inf, hjust = 1, vjust = 1)) +
        ggplot2::xlab("Minimum ligand distance (Angstrom)") +
        ggplot2::ylab(bquote("|"*.(trait) ~ Delta*Delta*Delta*"G| (kcal/mol)")) +
        ggplot2::theme_classic() +
        ggplot2::scale_colour_manual(values = plot_cols) +
        ggplot2::scale_size_manual(values = plot_sizes)
      ggplot2::ggsave(file.path(report_outpath, paste0("coupling_scatter_scHAminligandmin_scHAmincol_coreshape_", trait, ".pdf")), d, width = 4, height = 3, useDingbats=FALSE)
    }

    #Binding free energy change vs. ligand distance
    coef_dt <- fread(coef_file)[coef_order<=1 & trait_name=='Binding' & id!="WT"]
    int_dt <- fread(coef_file)[coef_order==2 & trait_name=='Binding' & id!="WT",.(dddg_abs = abs(.SD[[1]]), Pos_ref),,.SDcols = "mean_kcal/mol"]
    int_dt <- int_dt[dddg_abs>0.1]
    int_dt <- rbind(
      int_dt[,.(dddg_abs, int_id = Pos_ref, Pos_ref = sapply(strsplit(Pos_ref, '_'), '[', 1))], 
      int_dt[,.(dddg_abs, int_id = Pos_ref, Pos_ref = sapply(strsplit(Pos_ref, '_'), '[', 2))])
    plot_dt <- merge(
      coef_dt,
      int_dt,
      by = "Pos_ref", all.x = T)
    plot_dt[, ddg := .SD[[1]],,.SDcols = "mean_kcal/mol"]
    plot_dt[, ddg_ci := .SD[[1]],,.SDcols = "ci95_kcal/mol"]
    plot_dt[, SS_plot := 'loop']
    plot_dt[SS=='sheet', SS_plot := 'strand']
    plot_cols <- c(colour_scheme[["shade 0"]][4], colour_scheme[["shade 0"]][2], 'black')
    names(plot_cols) <- c("second_shell", "adjacent", "distal")
    d <- ggplot2::ggplot(plot_dt,ggplot2::aes(scHAmin_ligand, ddg, color = anno_ligand, shape = SS_plot, group = int_id)) +
      ggplot2::geom_line(ggplot2::aes(size = dddg_abs), color = 'lightgrey') +
      ggplot2::geom_point(size = 3) +
      # ggplot2::geom_linerange(ggplot2::aes(ymin = ddg-ddg_ci/2, ymax = ddg+ddg_ci/2), alpha = 1/4) +
      ggplot2::xlab("Distance to ligand (Angstrom)") +
      ggplot2::ylab(bquote("Binding "*Delta*Delta*"G (kcal/mol)")) +
      ggplot2::geom_text(data = plot_dt[,.(anno_ligand, SS_plot, int_id, label = paste("rho = ", round(cor(scHAmin_ligand, ddg, use = "pairwise.complete", method = "spearman"), 2), sep=""))], ggplot2::aes(label=label, x=-Inf, y=Inf, hjust = 0, vjust = 1)) +
      # ggplot2::geom_hline(yintercept = 0) +
      ggrepel::geom_text_repel(data = plot_dt[!duplicated(id_ref)], ggplot2::aes(label = id_ref), show.legend = T, 
        max.overlaps = Inf, size = 3) +
      ggplot2::scale_colour_manual(values = plot_cols) +
      ggplot2::scale_size(range = c(0.05,2)) +
      ggplot2::theme_classic()
    ggplot2::ggsave(file.path(report_outpath, "free_energy_scatter_binding_vs_ligand_distance_anno.pdf"), d, width = 4, height = 3, useDingbats=FALSE)

    #Binding absolute free energy change vs. ligand distance - abs
    coef_dt <- fread(coef_file)[coef_order<=1 & trait_name=='Binding' & id!="WT"]
    int_dt <- fread(coef_file)[coef_order==2 & trait_name=='Binding' & id!="WT",.(dddg_abs = abs(.SD[[1]]), Pos_ref),,.SDcols = "mean_kcal/mol"]
    int_dt <- int_dt[dddg_abs>0.1]
    int_dt <- rbind(
      int_dt[,.(dddg_abs, int_id = Pos_ref, Pos_ref = sapply(strsplit(Pos_ref, '_'), '[', 1))], 
      int_dt[,.(dddg_abs, int_id = Pos_ref, Pos_ref = sapply(strsplit(Pos_ref, '_'), '[', 2))])
    plot_dt <- merge(
      coef_dt,
      int_dt,
      by = "Pos_ref", all = T)
    plot_dt[, ddg_abs := abs(.SD[[1]]),,.SDcols = "mean_kcal/mol"]
    plot_dt[, ddg_abs_ci := .SD[[1]],,.SDcols = "ci95_kcal/mol"]
    plot_dt[, SS_plot := 'loop']
    plot_dt[SS=='sheet', SS_plot := 'strand']
    plot_cols <- c(colour_scheme[["shade 0"]][4], colour_scheme[["shade 0"]][2], 'black')
    names(plot_cols) <- c("second_shell", "adjacent", "distal")
    d <- ggplot2::ggplot(plot_dt,ggplot2::aes(scHAmin_ligand, ddg_abs, color = anno_ligand, shape = SS_plot, group = int_id)) +
      ggplot2::geom_line(ggplot2::aes(size = dddg_abs), color = 'lightgrey') +
      ggplot2::geom_point(size = 3) +
      # ggplot2::geom_linerange(ggplot2::aes(ymin = ddg_abs-ddg_abs_ci/2, ymax = ddg_abs+ddg_abs_ci/2), alpha = 1/4) +
      ggplot2::xlab("Distance to ligand (Angstrom)") +
      ggplot2::ylab(bquote("|Binding "*Delta*Delta*"G| (kcal/mol)")) +
      ggplot2::geom_text(data = plot_dt[,.(anno_ligand, SS_plot, int_id, label = paste("rho = ", round(cor(scHAmin_ligand, ddg_abs, use = "pairwise.complete", method = "spearman"), 2), sep=""))], ggplot2::aes(label=label, x=-Inf, y=Inf, hjust = 0, vjust = 1)) +
      # ggplot2::geom_hline(yintercept = 0) +
      ggrepel::geom_text_repel(data = plot_dt[!duplicated(id_ref)], ggplot2::aes(label = id_ref), show.legend = T, 
        max.overlaps = Inf, size = 3) +
      ggplot2::scale_colour_manual(values = plot_cols) +
      ggplot2::scale_size(range = c(0.05,2)) +
      ggplot2::theme_classic()
    ggplot2::ggsave(file.path(report_outpath, "free_energy_scatter_binding_vs_ligand_distance_anno_abs.pdf"), d, width = 4, height = 3, useDingbats=FALSE)

    #Folding free energy change vs. ligand distance
    coef_dt <- fread(coef_file)[coef_order<=1 & trait_name=='Folding' & id!="WT"]
    int_dt <- fread(coef_file)[coef_order==2 & trait_name=='Folding' & id!="WT",.(dddg_abs = abs(.SD[[1]]), Pos_ref),,.SDcols = "mean_kcal/mol"]
    int_dt <- int_dt[dddg_abs>0.1]
    int_dt <- rbind(
      int_dt[,.(dddg_abs, int_id = Pos_ref, Pos_ref = sapply(strsplit(Pos_ref, '_'), '[', 1))], 
      int_dt[,.(dddg_abs, int_id = Pos_ref, Pos_ref = sapply(strsplit(Pos_ref, '_'), '[', 2))])
    plot_dt <- merge(
      coef_dt,
      int_dt,
      by = "Pos_ref", all = T)
    plot_dt[, ddg := .SD[[1]],,.SDcols = "mean_kcal/mol"]
    plot_dt[, ddg_ci := .SD[[1]],,.SDcols = "ci95_kcal/mol"]
    plot_dt[, SS_plot := 'loop']
    plot_dt[SS=='sheet', SS_plot := 'strand']
    plot_cols <- c(colour_scheme[["shade 0"]][4], colour_scheme[["shade 0"]][2], 'black')
    names(plot_cols) <- c("second_shell", "adjacent", "distal")
    d <- ggplot2::ggplot(plot_dt,ggplot2::aes(scHAmin_ligand, ddg, color = anno_ligand, shape = SS_plot, group = int_id)) +
      ggplot2::geom_line(ggplot2::aes(size = dddg_abs), color = 'lightgrey') +
      ggplot2::geom_point(size = 3) +
      # ggplot2::geom_linerange(ggplot2::aes(ymin = ddg-ddg_ci/2, ymax = ddg+ddg_ci/2), alpha = 1/4) +
      ggplot2::xlab("Distance to ligand (Angstrom)") +
      ggplot2::ylab(bquote("Folding "*Delta*Delta*"G (kcal/mol)")) +
      ggplot2::geom_text(data = plot_dt[,.(anno_ligand, SS_plot, int_id, label = paste("rho = ", round(cor(scHAmin_ligand, ddg, use = "pairwise.complete", method = "spearman"), 2), sep=""))], ggplot2::aes(label=label, x=-Inf, y=Inf, hjust = 0, vjust = 1)) +
      # ggplot2::geom_hline(yintercept = 0) +
      ggrepel::geom_text_repel(data = plot_dt[!duplicated(id_ref)], ggplot2::aes(label = id_ref), show.legend = T, 
        max.overlaps = Inf, size = 3) +
      ggplot2::scale_colour_manual(values = plot_cols) +
      ggplot2::scale_size(range = c(0.05,2)) +
      ggplot2::theme_classic()
    ggplot2::ggsave(file.path(report_outpath, "free_energy_scatter_folding_vs_ligand_distance_anno.pdf"), d, width = 4, height = 3, useDingbats=FALSE)

    #Folding absolute free energy change vs. ligand distance - abs
    coef_dt <- fread(coef_file)[coef_order<=1 & trait_name=='Folding' & id!="WT"]
    int_dt <- fread(coef_file)[coef_order==2 & trait_name=='Folding' & id!="WT",.(dddg_abs = abs(.SD[[1]]), Pos_ref),,.SDcols = "mean_kcal/mol"]
    int_dt <- int_dt[dddg_abs>0.1]
    int_dt <- rbind(
      int_dt[,.(dddg_abs, int_id = Pos_ref, Pos_ref = sapply(strsplit(Pos_ref, '_'), '[', 1))], 
      int_dt[,.(dddg_abs, int_id = Pos_ref, Pos_ref = sapply(strsplit(Pos_ref, '_'), '[', 2))])
    plot_dt <- merge(
      coef_dt,
      int_dt,
      by = "Pos_ref")
    plot_dt[, ddg_abs := abs(.SD[[1]]),,.SDcols = "mean_kcal/mol"]
    plot_dt[, ddg_abs_ci := .SD[[1]],,.SDcols = "ci95_kcal/mol"]
    plot_dt[, SS_plot := 'loop']
    plot_dt[SS=='sheet', SS_plot := 'strand']
    plot_cols <- c(colour_scheme[["shade 0"]][4], colour_scheme[["shade 0"]][2], 'black')
    names(plot_cols) <- c("second_shell", "adjacent", "distal")
    d <- ggplot2::ggplot(plot_dt,ggplot2::aes(scHAmin_ligand, ddg_abs, color = anno_ligand, shape = SS_plot, group = int_id)) +
      ggplot2::geom_line(ggplot2::aes(size = dddg_abs), color = 'lightgrey') +
      ggplot2::geom_point(size = 3) +
      # ggplot2::geom_linerange(ggplot2::aes(ymin = ddg_abs-ddg_abs_ci/2, ymax = ddg_abs+ddg_abs_ci/2), alpha = 1/4) +
      ggplot2::xlab("Distance to ligand (Angstrom)") +
      ggplot2::ylab(bquote("|Folding "*Delta*Delta*"G| (kcal/mol)")) +
      ggplot2::geom_text(data = plot_dt[,.(anno_ligand, SS_plot, int_id, label = paste("rho = ", round(cor(scHAmin_ligand, ddg_abs, use = "pairwise.complete", method = "spearman"), 2), sep=""))], ggplot2::aes(label=label, x=-Inf, y=Inf, hjust = 0, vjust = 1)) +
      # ggplot2::geom_hline(yintercept = 0) +
      ggrepel::geom_text_repel(data = plot_dt[!duplicated(id_ref)], ggplot2::aes(label = id_ref), show.legend = T, 
        max.overlaps = Inf, size = 3) +
      ggplot2::scale_colour_manual(values = plot_cols) +
      ggplot2::scale_size(range = c(0.05,2)) +
      ggplot2::theme_classic()
    ggplot2::ggsave(file.path(report_outpath, "free_energy_scatter_folding_vs_ligand_distance_anno_abs.pdf"), d, width = 4, height = 3, useDingbats=FALSE)

    #Interaction barplots
    coef_dt1 <- fread(coef_file)[trait_name=="Binding" & coef_order==1]
    coef_dt2 <- fread(coef_file)[trait_name=="Binding" & coef_order==2 & scHAmin<8]
    #Annotations
    anno_list <- as.list(unlist(coef_dt1[,anno_ligand]))
    names(anno_list) <- unlist(coef_dt1[,Pos_ref])
    #Interactions
    int_tab <- table(unlist(strsplit(coef_dt2[,Pos_ref], "_")))
    all_int_dt <- data.table(
      N = as.numeric(int_tab),
      Pos_ref = names(int_tab))
    #Interactions - second shell
    coef_dt2[,Pos_ref1 := sapply(strsplit(Pos_ref, "_"), '[', 1)]
    coef_dt2[,Pos_ref2 := sapply(strsplit(Pos_ref, "_"), '[', 2)]
    coef_dt2[,second_shell1 := anno_list[[Pos_ref1]]=='second_shell',Pos_ref]
    coef_dt2[,second_shell2 := anno_list[[Pos_ref2]]=='second_shell',Pos_ref]
    coef_dt2_redundant <- rbind(
      coef_dt2[,.(Pos_ref1, Pos_ref2, second_shell1, second_shell2)],
      coef_dt2[,.(Pos_ref1 = Pos_ref2, Pos_ref2 = Pos_ref1, second_shell1 = second_shell2, second_shell2 = second_shell1)])
    all_int_dt_ss <- coef_dt2_redundant[second_shell2==T,.N,.(Pos_ref = Pos_ref1)]
    #Combine
    all_int_dt <- merge(
      all_int_dt[,.(Pos_ref, num_all = N)], 
      all_int_dt_ss[,.(Pos_ref, num_ss = N)], all = T, fill = 0)
    all_int_dt[is.na(num_ss), num_ss := 0]
    all_int_dt[, num_nss := num_all-num_ss]
    all_int_dt <- rbind(
      all_int_dt,
      data.table(Pos_ref = 10:31, num_all = 0, num_ss = 0, num_nss = 0)[!Pos_ref %in% all_int_dt])
    #Plot
    plot_dt <- data.table(reshape2::melt(all_int_dt[,.(Pos_ref, num_ss, num_nss)], id = 'Pos_ref'))
    plot_dt[, variable := relevel(variable, ref = "num_nss")]
    plot_cols <- c(colour_scheme[["shade 0"]][4], 'lightgrey')
    names(plot_cols) <- c("num_ss", "num_nss")
    d <- ggplot2::ggplot(plot_dt,ggplot2::aes(Pos_ref, value, fill = variable)) +
      ggplot2::geom_bar(stat = 'identity') +
      # ggplot2::xlab("Distance to ligand (Angstrom)") +
      # ggplot2::ylab(bquote("|Folding "*Delta*Delta*"G| (kcal/mol)")) +
      ggplot2::scale_fill_manual(values = plot_cols) +
      # ggplot2::coord_fixed() +
      ggplot2::theme_classic()
    ggplot2::ggsave(file.path(report_outpath, "interaction_freq_barplot.pdf"), d, width = 7, height = 2, useDingbats=FALSE)
  }

}

