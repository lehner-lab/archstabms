
#' archstabms__random_library_design
#'
#' Random library design.
#'
#' @param ddpca_outpath path to ddPCA data (required)
#' @param position_offset residue position offset (default:0)
#' @param fitness_abundance_file path to abundance fitness data (required)
#' @param fitness_binding_file path to binding fitness data (required)
#' @param max_num_mutants maximum number of mutations (required)
#' @param sample_size_order number of samples per mutation order (default:1000) 
#' @param load_mochi_data python (required)
#' @param colour_scheme colour scheme file (required)
#' @param outpath output path for plots and saved objects (required)
#'
#' @return Nothing
#' @export
#' @import data.table
archstabms__random_library_design <- function(
  ddpca_outpath,
  position_offset = 0,
  fitness_abundance_file,
  fitness_binding_file,
  max_num_mutants, 
  sample_size_order = 10000,
  load_mochi_data,
  colour_scheme,
  outpath
  ){

  #Domain name
  domain_name <- rev(unlist(strsplit(basename(outpath), "_")))[1]

  #Display status
  message(paste("\n\n*******", paste("running stage: archstabms__random_library_design for", domain_name), "*******\n\n"))

  #Create output directory
  archstabms__create_dir(archstabms_dir = outpath)

  #Fitness - abundance
  load(fitness_abundance_file)
  gra_lm <- all_variants[,lm(growthrate~fitness)]
  wta_gr <- all_variants[WT==T,growthrate] #WT growth rate

  #Fitness - binding
  load(fitness_binding_file)
  grb_lm <- all_variants[,lm(growthrate~fitness)]
  wtb_gr <- all_variants[WT==T,growthrate] #WT growth rate

  #Growth rate thresholds
  modpar <- fread(file.path(ddpca_outpath, "model_parameters_0.txt"), header = F)[,V1]
  modpar_list <- as.list(as.numeric(modpar[seq(2, length(modpar), 2)]))
  names(modpar_list) <- modpar[seq(1, length(modpar), 2)]
  gr_thresholds <- list(
    "Abundance" = as.numeric(predict(gra_lm, newdata = data.frame(fitness = 0.5*modpar_list[["folding_linear_kernel"]]+modpar_list[["folding_linear_bias"]]))),
    "Binding" = as.numeric(predict(grb_lm, newdata = data.frame(fitness = 0.5*modpar_list[["binding_linear_kernel"]]+modpar_list[["binding_linear_bias"]]))))

  #Determine candidate mutations (confident binding and folding ddGs, reachable by single nt substitutions, at least 20 doubles with single)
  #WT nucleotide sequence
  load(fitness_abundance_file)
  nt_wtseq <- all_variants[WT==T,nt_seq]
  #Model results
  pred_dt <- fread(file.path(ddpca_outpath, "../model_results.txt"))
  ddg_dt <- fread(file.path(ddpca_outpath, "../dg_singles.txt"))
  pred_double1_dt <- copy(pred_dt[mut_order==2 & dataset_binding==0])
  pred_double1_dt[, id_ref_double := sapply(strsplit(id_ref, ","), "[", 1)]
  pred_double2_dt <- copy(pred_dt[mut_order==2 & dataset_binding==0])
  pred_double2_dt[, id_ref_double := sapply(strsplit(id_ref, ","), "[", 2)]
  pred_double_dt <- rbind(pred_double1_dt, pred_double2_dt)
  pred_double_dt[, num_doubles := .N,id_ref_double,]
  pred_double_dt <- pred_double_dt[!duplicated(id_ref_double),.(id_ref = id_ref_double, num_doubles)]
  ddg_dt <- merge(ddg_dt, pred_double_dt, by = "id_ref", all.x = T) #Merge with ddGs
  #Candidate mutations
  # candmut_dt <- ddg_dt[f_ddg_pred_conf==T & b_ddg_pred_conf==T & id_ref %in% archstabms__reachable_single_AA_mutants(nt_wtseq, position_offset)[["reachable_muts"]] & num_doubles>=20,]
  candmut_dt <- ddg_dt[f_ddg_pred_conf==T & b_ddg_pred_conf==T,]

  #Load energies
  energies <- fread(file.path(ddpca_outpath, "model_weights_0.txt"))[,.(id, folding_coefficient, binding_coefficient)]
  #Subset to candidate mutations
  energies <- energies[id %in% candmut_dt[,id] | id == "WT"]
  energies[, WT_AA := substr(id, 1, 1)]
  energies[id!="WT", Mut := substr(id, nchar(id), nchar(id))]
  energies[id!="WT", Pos := as.integer(substr(id, 2, nchar(id)-1))]
  energies[id!="WT", id := paste0(WT_AA, Pos, Mut)]

  #Energies dictionary
  folding_energies_dict <- list()
  binding_energies_dict <- list()
  for(i in energies[id!="WT",unique(Pos)]){
    folding_energies_dict[[as.character(i)]] <- energies[Pos==i,folding_coefficient]
    #Duplicate if only one list element
    if(length(folding_energies_dict[[as.character(i)]])==1){
      folding_energies_dict[[as.character(i)]] <- c(folding_energies_dict[[as.character(i)]], folding_energies_dict[[as.character(i)]])
    }
    binding_energies_dict[[as.character(i)]] <- energies[Pos==i,binding_coefficient]
    #Duplicate if only one list element
    if(length(binding_energies_dict[[as.character(i)]])==1){
      binding_energies_dict[[as.character(i)]] <- c(binding_energies_dict[[as.character(i)]], binding_energies_dict[[as.character(i)]])
    }
  }

  #Random variant data.table
  set.seed(123456)
  #Set variant ids
  rand_dt <- data.table(varid = 1:sample_size_order)
  #Select random position to be mutated at each step
  rand_dt <- rand_dt[, mutpos := paste0(sample(names(folding_energies_dict), max_num_mutants), collapse = ","),varid]
  #Select random mutation at each position
  rand_dt[, mutcoeff_0 := energies[id=="WT",folding_coefficient]]
  rand_dt[, mutcoefb_0 := energies[id=="WT",binding_coefficient]]
  for(i in 1:max_num_mutants){
    rand_dt[, paste0("mutcoeff_", i) := sample(folding_energies_dict[[as.character(unlist(strsplit(mutpos, ","))[i])]], 1), varid]
    rand_dt[, paste0("mutcoefb_", i) := sample(binding_energies_dict[[as.character(unlist(strsplit(mutpos, ","))[i])]], 1), varid]
  }
  #Predict growth rate
  for(i in 1:max_num_mutants){
    rand_dt[, paste0("mutenergyf_", i) := apply(.SD, 1, sum),,.SDcols = paste0("mutcoeff_", 0:i)]
    rand_dt[, paste0("mutenergyb_", i) := apply(.SD, 1, sum),,.SDcols = paste0("mutcoefb_", 0:i)]
    rand_dt[, paste0("mutfitnessa_", i) := archstabms__predict_fitness_ddPCA(mochi_outpath = ddpca_outpath, folding_energy = .SD[[1]], binding_energy = NULL, RT = 1)[["fitness_folding"]],,.SDcols = paste0("mutenergyf_", i)]
    rand_dt[, paste0("mutfitnessb_", i) := archstabms__predict_fitness_ddPCA(mochi_outpath = ddpca_outpath, folding_energy = .SD[[1]], binding_energy = .SD[[2]], RT = 1)[["fitness_binding"]],,.SDcols = paste0(c("mutenergyf_", "mutenergyb_"), i)]
    rand_dt[, paste0("mutgra_", i) := gra_lm[["coefficients"]][1] + gra_lm[["coefficients"]][2] * .SD[[1]],,.SDcols = paste0("mutfitnessa_", i)]
    rand_dt[, paste0("mutgrb_", i) := grb_lm[["coefficients"]][1] + grb_lm[["coefficients"]][2] * .SD[[1]],,.SDcols = paste0("mutfitnessb_", i)]
  }

  #Plot - abundance
  plot_dt <- melt(rand_dt[,.SD,,.SDcols = paste0("mutgra_", 1:max_num_mutants)], measure.vars = paste0("mutgra_", 1:max_num_mutants))
  plot_dt[, mut_order := factor(sapply(strsplit(as.character(variable), "_"), '[', 2), levels = as.character(c(1:max_num_mutants)))]
  #Percentage of folded variants
  print(paste0("Percentage of folded variants with ", max_num_mutants, " random mutations: ", plot_dt[mut_order==max_num_mutants,sum(value>gr_thresholds[["Abundance"]])]/sample_size_order*100))
  text_dt <- plot_dt[,.(label = round(sum(value>gr_thresholds[["Abundance"]])/.N*100, 1), value = 0.17),mut_order]
  cc <- scales::seq_gradient_pal(colour_scheme[["shade 0"]][3], colour_scheme[["shade 0"]][1], "Lab")(seq(0,1,length.out=plot_dt[, max(as.integer(as.character(mut_order)))]))
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(mut_order, value, fill = mut_order, color = mut_order)) +
    ggplot2::geom_violin(scale = 'width') +
    # ggplot2::geom_boxplot() +
    ggplot2::scale_colour_manual(values = cc) +
    ggplot2::scale_fill_manual(values = cc) +
    ggplot2::geom_hline(yintercept = wta_gr) +
    ggplot2::geom_hline(yintercept = gr_thresholds[["Abundance"]], linetype = 2) +
    ggplot2::theme_classic() +
    ggplot2::geom_text(data = plot_dt[,.(label = paste("n = ", .N, sep=""), mut_order),][1], ggplot2::aes(label=label, x=-Inf, y=-Inf, hjust = 0, vjust = 0), color = "black") +
    ggplot2::geom_text(data = text_dt, ggplot2::aes(label=label)) +
    ggplot2::xlab("AA substitution order") +
    ggplot2::ylab("Predicted growth rate")
  ggplot2::ggsave(file.path(outpath, paste0("growthrate_violins_", domain_name, "_abundance_Nham_aa.pdf")), d, width = 5, height = 3, useDingbats=FALSE)
  #Plot - abundance up to 5th order
  plot_dt <- plot_dt[as.integer(as.character(mut_order))<=5]
  #Percentage of folded variants
  print(paste0("Percentage of folded variants with ", 5, " random mutations: ", plot_dt[mut_order==5,sum(value>gr_thresholds[["Abundance"]])]/sample_size_order*100))
  text_dt <- plot_dt[,.(label = round(sum(value>gr_thresholds[["Abundance"]])/.N*100, 1), value = 0.17),mut_order]
  cc <- scales::seq_gradient_pal(colour_scheme[["shade 0"]][3], colour_scheme[["shade 0"]][1], "Lab")(seq(0,1,length.out=5))
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(mut_order, value, fill = mut_order, color = mut_order)) +
    ggplot2::geom_violin(scale = 'width') +
    # ggplot2::geom_boxplot() +
    ggplot2::scale_colour_manual(values = cc) +
    ggplot2::scale_fill_manual(values = cc) +
    ggplot2::geom_hline(yintercept = wta_gr) +
    ggplot2::geom_hline(yintercept = gr_thresholds[["Abundance"]], linetype = 2) +
    ggplot2::theme_classic() +
    ggplot2::geom_text(data = plot_dt[,.(label = paste("n = ", .N, sep=""), mut_order),][1], ggplot2::aes(label=label, x=-Inf, y=-Inf, hjust = 0, vjust = 0), color = "black") +
    ggplot2::geom_text(data = text_dt, ggplot2::aes(label=label)) +
    ggplot2::xlab("AA substitution order") +
    ggplot2::ylab("Predicted growth rate")
  ggplot2::ggsave(file.path(outpath, paste0("growthrate_violins_", domain_name, "_abundance_Nham_aa_max5.pdf")), d, width = 3, height = 3, useDingbats=FALSE)

  #Plot - binding
  plot_dt <- melt(rand_dt[,.SD,,.SDcols = paste0("mutgrb_", 1:max_num_mutants)], measure.vars = paste0("mutgrb_", 1:max_num_mutants))
  plot_dt[, mut_order := factor(sapply(strsplit(as.character(variable), "_"), '[', 2), levels = as.character(c(1:max_num_mutants)))]
  #Percentage of bound variants
  print(paste0("Percentage of bound variants with ", max_num_mutants, " random mutations: ", plot_dt[mut_order==max_num_mutants,sum(value>gr_thresholds[["Binding"]])]/sample_size_order*100))
  text_dt <- plot_dt[,.(label = round(sum(value>gr_thresholds[["Binding"]])/.N*100, 1), value = 0.17),mut_order]
  cc <- scales::seq_gradient_pal(colour_scheme[["shade 0"]][3], colour_scheme[["shade 0"]][1], "Lab")(seq(0,1,length.out=plot_dt[, max(as.integer(as.character(mut_order)))]))
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(mut_order, value, fill = mut_order, color = mut_order)) +
    ggplot2::geom_violin(scale = 'width') +
    # ggplot2::geom_boxplot() +
    ggplot2::scale_colour_manual(values = cc) +
    ggplot2::scale_fill_manual(values = cc) +
    ggplot2::geom_hline(yintercept = wtb_gr) +
    ggplot2::geom_hline(yintercept = gr_thresholds[["Binding"]], linetype = 2) +
    ggplot2::theme_classic() +
    ggplot2::geom_text(data = plot_dt[,.(label = paste("n = ", .N, sep=""), mut_order),][1], ggplot2::aes(label=label, x=-Inf, y=-Inf, hjust = 0, vjust = 0), color = "black") +
    ggplot2::geom_text(data = text_dt, ggplot2::aes(label=label)) +
    ggplot2::xlab("AA substitution order") +
    ggplot2::ylab("Predicted growth rate")
  ggplot2::ggsave(file.path(outpath, paste0("growthrate_violins_", domain_name, "_binding_Nham_aa.pdf")), d, width = 5, height = 3, useDingbats=FALSE)
  #Plot - binding up to 5th order
  plot_dt <- plot_dt[as.integer(as.character(mut_order))<=5]
  #Percentage of bound variants
  print(paste0("Percentage of bound variants with ", 5, " random mutations: ", plot_dt[mut_order==5,sum(value>gr_thresholds[["Binding"]])]/sample_size_order*100))
  text_dt <- plot_dt[,.(label = round(sum(value>gr_thresholds[["Binding"]])/.N*100, 1), value = 0.17),mut_order]
  cc <- scales::seq_gradient_pal(colour_scheme[["shade 0"]][3], colour_scheme[["shade 0"]][1], "Lab")(seq(0,1,length.out=5))
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(mut_order, value, fill = mut_order, color = mut_order)) +
    ggplot2::geom_violin(scale = 'width') +
    # ggplot2::geom_boxplot() +
    ggplot2::scale_colour_manual(values = cc) +
    ggplot2::scale_fill_manual(values = cc) +
    ggplot2::geom_hline(yintercept = wtb_gr) +
    ggplot2::geom_hline(yintercept = gr_thresholds[["Binding"]], linetype = 2) +
    ggplot2::theme_classic() +
    ggplot2::geom_text(data = plot_dt[,.(label = paste("n = ", .N, sep=""), mut_order),][1], ggplot2::aes(label=label, x=-Inf, y=-Inf, hjust = 0, vjust = 0), color = "black") +
    ggplot2::geom_text(data = text_dt, ggplot2::aes(label=label)) +
    ggplot2::xlab("AA substitution order") +
    ggplot2::ylab("Predicted growth rate")
  ggplot2::ggsave(file.path(outpath, paste0("growthrate_violins_", domain_name, "_binding_Nham_aa_max5.pdf")), d, width = 3, height = 3, useDingbats=FALSE)


}
