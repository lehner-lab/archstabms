
#' archstabms_coupling_scatter_datasets
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
archstabms_coupling_scatter_datasets <- function(
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
    #order 2
    archstabms_coupling_scatter(
      input_file = file.path(base_dir, paste0("002", "_archstabms_structure_metrics_", i, "o2"), "model_coefficients.txt"),
      outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_coupling_scatter_", i, "o2"), stagenum=stagenum, base_dir=output_dir),
      colour_scheme = colour_scheme)
    #order 2 -linear
    archstabms_coupling_scatter(
      input_file = file.path(base_dir, paste0("002", "_archstabms_structure_metrics_", i, "o2-linear"), "model_coefficients.txt"),
      outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_coupling_scatter_", i, "o2-linear"), stagenum=stagenum, base_dir=output_dir),
      colour_scheme = colour_scheme)
    # #order 3
    # archstabms_coupling_scatter(
    #   input_file = file.path(base_dir, paste0("002", "_archstabms_structure_metrics_", i, "o3"), "model_coefficients.txt"),
    #   outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_coupling_scatter_", i, "o3"), stagenum=stagenum, base_dir=output_dir))
    # #order 3 -linear
    # archstabms_coupling_scatter(
    #   input_file = file.path(base_dir, paste0("002", "_archstabms_structure_metrics_", i, "o3-linear"), "model_coefficients.txt"),
    #   outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_coupling_scatter_", i, "o3-linear"), stagenum=stagenum, base_dir=output_dir))

    # #order 2 - ensemble
    # archstabms_coupling_scatter(
    #   input_file = file.path(base_dir, paste0("002", "_archstabms_structure_metrics_", i, "o2-ensemble"), "model_coefficients.txt"),
    #   outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_coupling_scatter_", i, "o2-ensemble"), stagenum=stagenum, base_dir=output_dir))
    # #order 2 -linear - ensemble
    # archstabms_coupling_scatter(
    #   input_file = file.path(base_dir, paste0("002", "_archstabms_structure_metrics_", i, "o2-linear-ensemble"), "model_coefficients.txt"),
    #   outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_coupling_scatter_", i, "o2-linear-ensemble"), stagenum=stagenum, base_dir=output_dir))
    # #order 3 - ensemble
    # archstabms_coupling_scatter(
    #   input_file = file.path(base_dir, paste0("002", "_archstabms_structure_metrics_", i, "o3-ensemble"), "model_coefficients.txt"),
    #   outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_coupling_scatter_", i, "o3-ensemble"), stagenum=stagenum, base_dir=output_dir))
    # #order 3 -linear - ensemble
    # archstabms_coupling_scatter(
    #   input_file = file.path(base_dir, paste0("002", "_archstabms_structure_metrics_", i, "o3-linear-ensemble"), "model_coefficients.txt"),
    #   outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_coupling_scatter_", i, "o3-linear-ensemble"), stagenum=stagenum, base_dir=output_dir))

    # #order 3 - sparse
    # archstabms_coupling_scatter(
    #   input_file = file.path(base_dir, paste0("002", "_archstabms_structure_metrics_", i, "o3-sparse"), "model_coefficients.txt"),
    #   outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_coupling_scatter_", i, "o3-sparse"), stagenum=stagenum, base_dir=output_dir))
    # #order 3 -linear - sparse
    # archstabms_coupling_scatter(
    #   input_file = file.path(base_dir, paste0("002", "_archstabms_structure_metrics_", i, "o3-sparse-linear"), "model_coefficients.txt"),
    #   outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_coupling_scatter_", i, "o3-sparse-linear"), stagenum=stagenum, base_dir=output_dir))
    # #order 3 - ensemble - sparse
    # archstabms_coupling_scatter(
    #   input_file = file.path(base_dir, paste0("002", "_archstabms_structure_metrics_", i, "o3-sparse-ensemble"), "model_coefficients.txt"),
    #   outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_coupling_scatter_", i, "o3-sparse-ensemble"), stagenum=stagenum, base_dir=output_dir))
    # #order 3 -linear - ensemble - sparse
    # archstabms_coupling_scatter(
    #   input_file = file.path(base_dir, paste0("002", "_archstabms_structure_metrics_", i, "o3-sparse-linear-ensemble"), "model_coefficients.txt"),
    #   outpath = archstabms__format_dir(dir_suffix=paste0("_archstabms_coupling_scatter_", i, "o3-sparse-linear-ensemble"), stagenum=stagenum, base_dir=output_dir))
  }

  ### Coupling partial correlations
  ###########################

  outpath <- archstabms__format_dir(dir_suffix=paste0("_archstabms_coupling_scatter_datasets"), stagenum=stagenum, base_dir=output_dir)
  
  #Create output directory
  archstabms__create_dir(archstabms_dir = outpath)

  cor_list <- list()
  for(suffix in c("o2", "o2-linear")){
    for(i in dataset_names){
      if(!file.exists(file.path(base_dir, paste0("002", "_archstabms_structure_metrics_", i, suffix), "model_coefficients.txt"))){next}
      dg_dt <- fread(file.path(base_dir, paste0("002", "_archstabms_structure_metrics_", i, suffix), "model_coefficients.txt"))[coef_order==2]
      dg_dt[, contact := as.numeric(scHAmin<5)]
      #Weighted mean absolute value
      dg_dt[, dddg_abs := sum(abs(.SD[[1]])/.SD[[2]]^2, na.rm = T)/sum(1/.SD[[2]]^2, na.rm = T),Pos_ref,.SDcols = c("mean_kcal/mol", "std_kcal/mol")]
      dg_dt[, dddg_ci := sqrt(1/sum(1/.SD[[2]]^2, na.rm = T))*1.96*2,Pos_ref,.SDcols = c("mean_kcal/mol", "std_kcal/mol")]
      dg_dt <- dg_dt[!duplicated(Pos_ref)]
      #Correlation
      cor_obj_scHAmin <- cor.test(dg_dt[,dddg_abs], dg_dt[, scHAmin], method = 'spearman')
      cor_obj_backbone <- cor.test(dg_dt[,dddg_abs], dg_dt[, backbone], method = 'spearman')
      cor_dt <- data.table(
        dataset = rep(i, 2),
        model_type = rep(c("biophysical", "linear")[as.numeric(grepl("linear", suffix))+1], 2),
        metric = c('scHAmin', 'backbone'),
        type = rep('full', 2),
        cor = c(
          cor_obj_scHAmin$estimate,
          cor_obj_backbone$estimate),
        pvalue = c(
          cor_obj_scHAmin$p.value,
          cor_obj_backbone$p.value))
      #Partial correlation
      pcor_obj <- ppcor::pcor(as.data.frame(dg_dt[,.SD,,.SDcols = c("dddg_abs", "scHAmin", "backbone")]), method = 'spearman')
      pcor_dt <- data.table(
        dataset = rep(i, 2),
        model_type = rep(c("biophysical", "linear")[as.numeric(grepl("linear", suffix))+1], 2),
        metric = c('scHAmin', 'backbone'),
        type = rep('partial', 2),
        cor = c(
          pcor_obj$estimate['dddg_abs','scHAmin'],
          pcor_obj$estimate['dddg_abs','backbone']),
        pvalue = c(
          pcor_obj$p.value['dddg_abs','scHAmin'],
          pcor_obj$p.value['dddg_abs','backbone']))
      cor_list[[paste0(i, "_", suffix)]] <- rbind(cor_dt, pcor_dt)
    }
  }

  plot_dt <- rbindlist(cor_list)
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(dataset, cor)) +
    ggplot2::geom_bar(stat = "identity", ggplot2::aes(fill = type), position = 'dodge') +
    ggplot2::xlab("Dataset") +
    ggplot2::ylab("Correlation") +
    ggplot2::theme_bw() +
    ggplot2::facet_grid(model_type~metric, scales = "free")
  ggplot2::ggsave(file.path(outpath, "coupling_correlation.pdf"), d, width = 10, height = 7, useDingbats=FALSE)

  ### Coupling linear models - training data 
  ###########################

  dg_dt <- fread(file.path(base_dir, paste0("002", "_archstabms_structure_metrics_", "CM6o2"), "model_coefficients.txt"))
  dg_dt1 <- dg_dt[coef_order==1]
  dg_dt2 <- dg_dt[coef_order==2]

  #Get single mutant energies
  dg_dt1[, id_ref1 := id_ref]
  dg_dt1[, id_ref2 := id_ref]
  dg_dt2[, id_ref1 := sapply(strsplit(id_ref, "_"), '[', 1)]
  dg_dt2[, id_ref2 := sapply(strsplit(id_ref, "_"), '[', 2)]
  dg_dt <- merge(dg_dt2, dg_dt1[,.(id_ref1, energy1 = `mean_kcal/mol`)], by = "id_ref1") 
  dg_dt <- merge(dg_dt, dg_dt1[,.(id_ref2, energy2 = `mean_kcal/mol`)], by = "id_ref2") 

  #Other features
  dg_dt[, energy := `mean_kcal/mol`]
  dg_dt[, abs_energy := abs(`mean_kcal/mol`)]
  dg_dt[, core := as.numeric(Pos1_class=="core")+as.numeric(Pos2_class=="core")]
  dg_dt[, bi := as.numeric(Pos1_class=="binding_interface")+as.numeric(Pos2_class=="binding_interface")]
  dg_dt[, helix := sum(c(as.numeric(SS1=="helix"), as.numeric(SS2=="helix")), na.rm = T), id_ref]
  dg_dt[, sheet := sum(c(as.numeric(SS1=="sheet"), as.numeric(SS2=="sheet")), na.rm = T), id_ref]
  dg_dt[, energy_sum := energy1+energy2]
  dg_dt[, energy_absdiff := abs(energy1-energy2)]
  all_ids <- dg_dt[,id_ref]

  #Model - |coupling energy|
  mm <- dg_dt[,.(abs_energy, core, bi, sheet, scHAmin, backbone, WT_hbbb, WT_hbsb, WT_hbss, WT_pc, WT_ps, WT_sb, WT_vdw)]
  mm_lm <- lm(abs_energy~.,data = mm)

  #Model performance on training data
  plot_dt <- data.table(
    observed_energy = mm[,abs_energy],
    predicted_energy = mm_lm$fitted.values)
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(observed_energy, predicted_energy)) +
    ggplot2::stat_binhex(bins = 50, size = 0.2, color = "black") +
    ggplot2::geom_smooth(method = "lm", formula = 'y~x', color = colour_scheme[["shade 0"]][1], se = T) +
    ggplot2::scale_fill_gradientn(colours = c("white", "black")) +
    # ggplot2::geom_point(alpha = 1/10) +
    ggplot2::xlab(expression("Observed |Folding "*Delta*Delta*Delta*"G|")) +
    ggplot2::ylab(expression("Predicted |Folding "*Delta*Delta*Delta*"G|")) +
    ggplot2::geom_text(data = plot_dt[,.(label = paste("r = ", round(cor(observed_energy, predicted_energy), 2), sep=""))], ggplot2::aes(label=label, x=Inf, y=-Inf, hjust = 1, vjust = 0), color = "black") +
    ggplot2::theme_classic()
  # d <- d + ggplot2::geom_abline(linetype = 2)
  ggplot2::ggsave(file.path(outpath, paste0("CM6_coupling_lm_train.pdf")), d, width = 4, height = 3, useDingbats=FALSE)

  #Feature significance
  mm_lm_summary <- summary(mm_lm)$coefficients
  plot_dt <- data.table(
    feature = rownames(mm_lm_summary),
    nlog10pvalue = -log10(as.data.table(mm_lm_summary)[,`Pr(>|t|)`]))
  plot_dt <- plot_dt[-1]
  plot_dt <- plot_dt[order(nlog10pvalue, decreasing = T)]
  plot_dt[, feature := factor(feature, levels = feature)]
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(feature, nlog10pvalue)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::xlab("Feature") +
    ggplot2::ylab("Significance, -log10(P-value)") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1),
      axis.text.y = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1))
  ggplot2::ggsave(file.path(outpath, "CM6_coupling_lm_coefficients.pdf"), d, width = 4, height = 3, useDingbats=FALSE)

  ### Coupling linear models - test data 
  ###########################

  dg_dt <- fread(file.path(base_dir, paste0("002", "_archstabms_structure_metrics_", "CM1o2"), "model_coefficients.txt"))
  dg_dt1 <- dg_dt[coef_order==1]
  dg_dt2 <- dg_dt[coef_order==2]

  #Get single mutant energies
  dg_dt1[, id_ref1 := id_ref]
  dg_dt1[, id_ref2 := id_ref]
  dg_dt2[, id_ref1 := sapply(strsplit(id_ref, "_"), '[', 1)]
  dg_dt2[, id_ref2 := sapply(strsplit(id_ref, "_"), '[', 2)]
  dg_dt <- merge(dg_dt2, dg_dt1[,.(id_ref1, energy1 = `mean_kcal/mol`)], by = "id_ref1") 
  dg_dt <- merge(dg_dt, dg_dt1[,.(id_ref2, energy2 = `mean_kcal/mol`)], by = "id_ref2") 

  #Other features
  dg_dt[, energy := `mean_kcal/mol`]
  dg_dt[, abs_energy := abs(`mean_kcal/mol`)]
  dg_dt[, core := as.numeric(Pos1_class=="core")+as.numeric(Pos2_class=="core")]
  dg_dt[, bi := as.numeric(Pos1_class=="binding_interface")+as.numeric(Pos2_class=="binding_interface")]
  dg_dt[, helix := sum(c(as.numeric(SS1=="helix"), as.numeric(SS2=="helix")), na.rm = T), id_ref]
  dg_dt[, sheet := sum(c(as.numeric(SS1=="sheet"), as.numeric(SS2=="sheet")), na.rm = T), id_ref]
  dg_dt[, energy_sum := energy1+energy2]
  dg_dt[, energy_absdiff := abs(energy1-energy2)]

  #Model matrix - |coupling energy|
  mm <- dg_dt[!id_ref %in% all_ids,.(abs_energy, core, bi, sheet, scHAmin, backbone, WT_hbbb, WT_hbsb, WT_hbss, WT_pc, WT_ps, WT_sb, WT_vdw, energy_sum, energy_absdiff)]
  
  #Model performance on training data
  plot_dt <- data.table(
    observed_energy = mm[,abs_energy],
    predicted_energy = predict(mm_lm, mm))
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(observed_energy, predicted_energy)) +
    ggplot2::geom_point() +
    ggplot2::geom_smooth(method = "lm", formula = 'y~x', color = colour_scheme[["shade 0"]][1], se = T) +
    # ggplot2::scale_fill_gradientn(colours = c("white", "black")) +
    ggplot2::xlab(expression("Observed |Folding "*Delta*Delta*Delta*"G|")) +
    ggplot2::ylab(expression("Predicted |Folding "*Delta*Delta*Delta*"G|")) +
    ggplot2::geom_text(data = plot_dt[,.(label = paste("r = ", round(cor(observed_energy, predicted_energy), 2), sep=""))], ggplot2::aes(label=label, x=Inf, y=-Inf, hjust = 1, vjust = 0), color = "black") +
    ggplot2::theme_classic()
  # d <- d + ggplot2::geom_abline(linetype = 2)
  ggplot2::ggsave(file.path(outpath, paste0("CM6_coupling_lm_test.pdf")), d, width = 3, height = 3, useDingbats=FALSE)

}
