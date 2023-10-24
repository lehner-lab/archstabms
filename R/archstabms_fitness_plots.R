
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

  ### Plot fitness distributions, replicate correlations, violins
  ###########################

  #Load fitness data
  mochi_outpath <- file.path(base_dir, "Data", "mochi", paste0(dataset_name, "o1"))
  fitness_dt <- fread(file.path(mochi_outpath, 'task_1', 'predictions', "predicted_phenotypes_all.txt"))
  # #Save for supplement
  # write.table(fitness_dt[Nham_aa==Nmut_codons & STOP==F & STOP_readthrough==F,.SD,,.SDcols = c("protein", "pca_type", "aa_seq", "Nham_aa", "WT", "fitness", "sigma", "growthrate", "growthrate_sigma")], 
  #   file = file.path(report_outpath, "fitness_supp.txt"), 
  #   quote = F, sep = "\t", row.names = F)
  dataset_names <- names(fitness_dt)[grepl("Abundance_|Binding_", names(fitness_dt))]

  #Supplementary files
  for(i in dataset_names){
    output_dt <- fitness_dt[get(i)==1,.(library = sapply(strsplit(i, "_"), '[', 2), pca_type = sapply(strsplit(i, "_"), '[', 1), aa_seq, Nham_aa, WT, fitness, sigma, growthrate, growthrate_sigma)]
    write.table(output_dt, file = file.path(outpath, paste0("Supplementary_fitness_", i, ".txt")), row.names = F, quote = F, sep = "\t")
  }

  #Growth rate thresholds
  gr_thresholds <- list()
  for(i in dataset_names){
    linears_dt <- fread(file.path(mochi_outpath, 'task_1', 'weights', paste0("linears_weights_", i, ".txt")))
    gr_lm <- lm(growthrate ~ fitness, data = fitness_dt[get(i)==1])
    gr_thresholds[[i]] <- as.numeric(predict(gr_lm, newdata = data.frame(fitness = linears_dt[, mean(0.5*kernel+bias)])))
  }

  #Plot fitness replicate correlation
  for(i in dataset_names){
    archstabms__plot_fitness_replicates_cor(
      input_dt = fitness_dt[get(i)==1,],
      output_file = file.path(outpath, paste0("replicate_scatter_binhex_", sapply(strsplit(i, "_"), '[', 1), ".pdf")),
      colour_scheme = colour_scheme)
  }

  #Plot growthrate violins with hamming distance
  for(i in dataset_names){
    archstabms__plot_growthrate_violins(
      input_dt = fitness_dt[get(i)==1,], 
      report_outpath = outpath,
      dataset_name = gsub("Abundance_|Binding_", "", i),
      trait = sapply(strsplit(i, "_"), '[', 1),
      threshold = gr_thresholds[[i]],
      colour_scheme = colour_scheme)
  }

  ### Plot abundance fitness distributions, replicate correlations, violins - lenient
  ###########################

  if(dir.exists(file.path(base_dir, "Data", "fitness", dataset_name))){
    #Load abundance fitness data
    fitness_filename <- list.files(file.path(base_dir, "Data", "fitness", dataset_name, ""), pattern = 'RData$')[1]
    load(file.path(base_dir, "Data", "fitness", dataset_name, fitness_filename))
    fitness_dt <- all_variants

    #Supplementary files
    output_dt <- fitness_dt[,.(library = paste0(sapply(strsplit(i, "_"), '[', 2), "_lenient"), pca_type = sapply(strsplit(i, "_"), '[', 1), aa_seq, Nham_aa, WT, fitness, sigma, growthrate, growthrate_sigma)]
    write.table(output_dt, file = file.path(outpath, paste0("Supplementary_fitness_", i, "_lenient.txt")), row.names = F, quote = F, sep = "\t")

    #Plot fitness replicate correlation
    archstabms__plot_fitness_replicates_cor(
      input_dt = fitness_dt,
      output_file = file.path(outpath, "replicate_scatter_binhex_Abundance_lenient.pdf"),
      colour_scheme = colour_scheme)

    #Plot abundance growthrate violins with hamming distance
    archstabms__plot_growthrate_violins(
      input_dt = fitness_dt, 
      report_outpath = outpath,
      dataset_name = dataset_name,
      trait = 'Abundance',
      suffix = "_lenient",
      threshold = gr_thresholds[[1]],
      colour_scheme = colour_scheme)

    #Number of variants with 20 substitutions and WT-like fitness
    fitness_dt[, pval := pnorm(fitness, fitness_dt[WT==T,fitness], sigma, lower.tail = T)]
    print(paste0("Num. variants with >20 substitutions and WT-like fitness (", dataset_name, "): ", fitness_dt[Nham_aa>20 & pval>0.05,.N]))

    #Number of variants with greater than 5 substitutions and WT-like fitness
    fitness_dt[, pval := pnorm(fitness, fitness_dt[WT==T,fitness], sigma, lower.tail = T)]
    print(paste0("Num. variants with >5 substitutions and WT-like fitness (", dataset_name, "): ", fitness_dt[Nham_aa>5 & pval>0.05,.N]))

  }

}