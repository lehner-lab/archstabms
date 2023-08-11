
#' archstabms_coupling_decay
#'
#' Coupling predictive decay plots.
#'
#' @param input_dir_list list of MoCHI directories (required)
#' @param min_mut_order minimum mutation order (required)
#' @param test_variants whether to use unseen test variants only: T/F (required)
#' @param outpath output path for plots and saved objects (required)
#' @param colour_scheme colour scheme file (required)
#'
#' @return Nothing
#' @export
#' @import data.table
archstabms_coupling_decay <- function(
	input_dir_list,
  min_mut_order,
  test_variants,
	outpath,
  colour_scheme
	){

	#Display status
	message(paste("\n\n*******", paste("running stage: archstabms_coupling_decay"), "*******\n\n"))

	#Create output directory
	archstabms__create_dir(archstabms_dir = outpath)

  ### Explainable variance - order 2
  ###########################

  result_list <- list()
  for(i in 1:3){
    pred_dt <- data.table()
    pred_file_train <- file.path(input_dir_list[[paste0("order", i)]], "task_1", "predictions", "predicted_phenotypes_all.txt")
    pred_file <- file.path(input_dir_list[[paste0("order", i)]], "task_1", "predictions", "predicted_phenotypes_supp.txt")
    if(test_variants){
      pred_dt_model <- fread(pred_file_train)
      pred_dt <- fread(pred_file)[!aa_seq %in% pred_dt_model[,aa_seq]]
    }else{
      pred_dt <- fread(pred_file_train)
    }
    result_list[[i]] <- pred_dt[Nham_aa>=min_mut_order,.(
      r2_fold1 = cor(fitness, fold_1)^2, 
      r2_fold2 = cor(fitness, fold_2)^2, 
      r2_fold3 = cor(fitness, fold_3)^2, 
      r2_fold4 = cor(fitness, fold_4)^2, 
      r2_fold5 = cor(fitness, fold_5)^2, 
      r2_fold6 = cor(fitness, fold_6)^2, 
      r2_fold7 = cor(fitness, fold_7)^2, 
      r2_fold8 = cor(fitness, fold_8)^2, 
      r2_fold9 = cor(fitness, fold_9)^2, 
      r2_fold10 = cor(fitness, fold_10)^2, 
      # r2_fold1 = 1 - var(fitness - fold_1)/var(fitness), 
      # r2_fold2 = 1 - var(fitness - fold_2)/var(fitness), 
      # r2_fold3 = 1 - var(fitness - fold_3)/var(fitness), 
      # r2_fold4 = 1 - var(fitness - fold_4)/var(fitness), 
      # r2_fold5 = 1 - var(fitness - fold_5)/var(fitness), 
      # r2_fold6 = 1 - var(fitness - fold_6)/var(fitness), 
      # r2_fold7 = 1 - var(fitness - fold_7)/var(fitness), 
      # r2_fold8 = 1 - var(fitness - fold_8)/var(fitness), 
      # r2_fold9 = 1 - var(fitness - fold_9)/var(fitness), 
      # r2_fold10 = 1 - var(fitness - fold_10)/var(fitness), 
      n = .N, 
      r2_max = mean(c(
        cor(fitness1_uncorr, fitness2_uncorr, use = "pairwise.complete")^2,
        cor(fitness1_uncorr, fitness3_uncorr, use = "pairwise.complete")^2,
        cor(fitness2_uncorr, fitness3_uncorr, use = "pairwise.complete")^2))), Nham_aa]
    names(result_list[[i]])[grepl("r2_fold", names(result_list[[i]]))] <- paste0("r2o", i, "_fold", 1:10)
    # result_list[[i]][, paste0("r2o", i) := rowMeans(.SD),,.SDcols = paste0("r2_fold", 1:10)]
    # result_list[[i]][, paste0("r2o", i, "_sd") := apply(.SD, 1, sd),,.SDcols = paste0("r2_fold", 1:10)]
  }

  #Merge
  common_names <- names(result_list[[1]])[names(result_list[[1]]) %in% names(result_list[[2]])]
  r2_dt <- merge(result_list[[1]], result_list[[2]], by = common_names)
  r2_dt <- merge(r2_dt, result_list[[3]], by = common_names)
  for(i in 1:10){
    r2_dt[, paste0("additive", i) := .SD[[1]]/r2_max,,.SDcols = paste0("r2o1_fold", i)]
    r2_dt[, paste0("r2o2p", i) := .SD[[1]]/r2_max,,.SDcols = paste0("r2o2_fold", i)]
    r2_dt[, paste0("pairwise", i) := .SD[[1]]-.SD[[2]],,.SDcols = c(paste0("r2o2p", i), paste0("additive", i))]
    # r2_dt[get(paste0("pairwise", i))<0, paste0("pairwise", i) := 0]
    r2_dt[, paste0("r2o3p", i) := .SD[[1]]/r2_max,,.SDcols = paste0("r2o3_fold", i)]
    r2_dt[, paste0("threeway", i) := .SD[[1]]-.SD[[2]]-.SD[[3]],,.SDcols = c(paste0("r2o3p", i), paste0("additive", i), paste0("pairwise", i))]
    # r2_dt[get(paste0("threeway", i))<0, paste0("threeway", i) := 0]
    r2_dt[, paste0("higher", i) := 1-.SD[[1]]-.SD[[2]]-.SD[[3]],,.SDcols = c(paste0("additive", i), paste0("pairwise", i), paste0("threeway", i))]
  }
  #Mean across 10 folds
  r2_dt[, "additive" := rowMeans(.SD),,.SDcols = paste0("additive", 1:10)]
  r2_dt[, "pairwise" := rowMeans(.SD),,.SDcols = paste0("pairwise", 1:10)]
  r2_dt[, "threeway" := rowMeans(.SD),,.SDcols = paste0("threeway", 1:10)]
  r2_dt[, "higher" := rowMeans(.SD),,.SDcols = paste0("higher", 1:10)]
  #Standard deviation across 10 folds
  r2_dt[, "additive_sd" := apply(.SD, 1, sd),,.SDcols = paste0("additive", 1:10)]
  r2_dt[, "pairwise_sd" := apply(.SD, 1, sd),,.SDcols = paste0("pairwise", 1:10)]
  r2_dt[, "threeway_sd" := apply(.SD, 1, sd),,.SDcols = paste0("threeway", 1:10)]
  r2_dt[, "higher_sd" := apply(.SD, 1, sd),,.SDcols = paste0("higher", 1:10)]

  #Barplot
  plot_dt <- melt(copy(r2_dt)[,.(Nham_aa, additive, pairwise, threeway, higher)], id = c('Nham_aa'))
  plot_dt_sd <- melt(copy(r2_dt)[,.(Nham_aa, additive = additive_sd, pairwise = pairwise_sd, threeway = threeway_sd)], id = c('Nham_aa'))
  names(plot_dt_sd)[names(plot_dt_sd)=="value"] = "value_sd"
  plot_dt <- merge(plot_dt, plot_dt_sd, by = c('Nham_aa', 'variable'), all = T)
  plot_dt[, variable := factor(variable, levels = c('higher', 'threeway',  'pairwise', 'additive'))]
  plot_dt[, mut_order := factor(Nham_aa)]
  plot_col <- c(colour_scheme[["shade 0"]][3], colour_scheme[["shade 0"]][1], colour_scheme[["shade 0"]][2], "lightgrey")
  names(plot_col) <- c("additive", "pairwise", "threeway", "higher")
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(mut_order, value, fill = variable)) +
    ggplot2::geom_col(position = 'stack') +
    ggplot2::theme_classic() +
    ggplot2::xlab("AA substitution order") +
    ggplot2::ylab("Proportion explainable variance") +
    ggplot2::scale_fill_manual(values = plot_col)
  ggplot2::ggsave(file.path(outpath, paste0("explainable_variance_Nham_aa.pdf")), d, width = 5, height = 3, useDingbats=FALSE)
  d <- d + ggplot2::coord_cartesian(xlim = c(1, 8))
  ggplot2::ggsave(file.path(outpath, paste0("explainable_variance_Nham_aa_xlim.pdf")), d, width = 5, height = 3, useDingbats=FALSE)

  #Line plot
  d <- ggplot2::ggplot(plot_dt[variable!='higher'],ggplot2::aes(Nham_aa, value, color = variable)) +
    ggplot2::geom_point() +
    ggplot2::geom_smooth(data = plot_dt[Nham_aa<=21 & variable!='higher'], method="lm", formula= (y ~ exp(-x))) +
    ggplot2::theme_classic() +
    ggplot2::xlab("AA substitution order") +
    ggplot2::ylab("Proportion explainable variance") +
    ggplot2::scale_colour_manual(values = plot_col)
  ggplot2::ggsave(file.path(outpath, paste0("explainable_variance_Nham_aa_delta.pdf")), d, width = 5, height = 3, useDingbats=FALSE)
  d <- d + ggplot2::geom_linerange(ggplot2::aes(ymin = value-value_sd*1.96, ymax = value+value_sd*1.96), alpha = 1/4)
  ggplot2::ggsave(file.path(outpath, paste0("explainable_variance_Nham_aa_delta_ci95.pdf")), d, width = 5, height = 3, useDingbats=FALSE)
  d <- d + ggplot2::coord_cartesian(xlim = c(14, 21))
  ggplot2::ggsave(file.path(outpath, paste0("explainable_variance_Nham_aa_delta_xlim.pdf")), d, width = 5, height = 3, useDingbats=FALSE)

  #Line plot - no additive
  d <- ggplot2::ggplot(plot_dt[variable %in% c('pairwise', 'threeway')],ggplot2::aes(Nham_aa, value, color = variable)) +
    ggplot2::geom_point() +
    ggplot2::geom_smooth(data = plot_dt[Nham_aa<=21 & variable %in% c('pairwise', 'threeway')], method="lm", formula= (y ~ exp(-x))) +
    ggplot2::theme_classic() +
    ggplot2::xlab("AA substitution order") +
    ggplot2::ylab("Proportion explainable variance") +
    ggplot2::coord_cartesian(ylim = c(0, 0.15)) +
    ggplot2::scale_colour_manual(values = plot_col)
  ggplot2::ggsave(file.path(outpath, paste0("explainable_variance_Nham_aa_delta_noadditive.pdf")), d, width = 5, height = 3, useDingbats=FALSE)
  d <- d + ggplot2::geom_linerange(ggplot2::aes(ymin = value-value_sd*1.96, ymax = value+value_sd*1.96), alpha = 1/4)
  ggplot2::ggsave(file.path(outpath, paste0("explainable_variance_Nham_aa_delta_noadditive_ci95.pdf")), d, width = 5, height = 3, useDingbats=FALSE)
  d <- d + ggplot2::coord_cartesian(xlim = c(14, 21))
  ggplot2::ggsave(file.path(outpath, paste0("explainable_variance_Nham_aa_delta_noadditive_xlim.pdf")), d, width = 5, height = 3, useDingbats=FALSE)




}
