
#' archstabms_coupling_scatter
#'
#' Coupling scatterplots.
#'
#' @param input_file path to MoCHI thermo model fit results (required)
#' @param trait trait name (default:"Folding")
#' @param outpath output path for plots and saved objects (required)
#' @param colour_scheme colour scheme file (required)
#'
#' @return Nothing
#' @export
#' @import data.table
archstabms_coupling_scatter <- function(
	input_file,
  trait="Folding",
	outpath,
  colour_scheme
	){

  #Return if input_file doesn't exist
  if(!file.exists(input_file)){
    return()
  }

	#Load free energies
	dg_dt <- fread(input_file)

  #Return if trait doesn't exist
  if(!trait %in% dg_dt[,trait_name]){
    return()
  }

  #Domain name
  domain_name <- rev(unlist(strsplit(basename(outpath), "_")))[1]

  #Display status
  message(paste("\n\n*******", paste0("running stage: archstabms_coupling_scatter for ", domain_name, "(trait: ", trait, ")"), "*******\n\n"))

  #Create output directory
  archstabms__create_dir(archstabms_dir = outpath)

  #Subset to trait free energies
  dg_dt <- dg_dt[trait_name==trait]

  #Free energy densities
  plot_dt <- dg_dt[coef_order %in% c(1,2)]
  plot_dt <- plot_dt[, energy := .SD[[1]],,.SDcols = "mean_kcal/mol"]
  plot_dt <- plot_dt[, coef_order_plot := factor(coef_order, levels = c(1, 2))]
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(energy, fill = coef_order_plot)) +
    ggplot2::geom_density() +
    ggplot2::geom_vline(xintercept = 0) + 
    ggplot2::xlab(bquote(""*.(trait) ~ Delta*Delta*Delta*"G (kcal/mol)")) +
    ggplot2::ylab("Density") +
    ggplot2::theme_classic() +
    ggplot2::facet_wrap(coef_order_plot~., nrow = 2) +
    ggplot2::scale_fill_manual(values = c('grey', colour_scheme[["shade 0"]][1]))
  ggplot2::ggsave(file.path(outpath, "coefficient_densities.pdf"), d, width = 4, height = 3, useDingbats=FALSE)
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(energy, fill = coef_order_plot)) +
    ggplot2::geom_density() +
    ggplot2::geom_vline(xintercept = 0) + 
    ggplot2::xlab(bquote(""*.(trait) ~ Delta*Delta*Delta*"G (kcal/mol)")) +
    ggplot2::ylab("Density") +
    ggplot2::theme_classic() +
    ggplot2::facet_wrap(coef_order_plot~., scales = "free", nrow = 2) +
    ggplot2::scale_fill_manual(values = c('grey', colour_scheme[["shade 0"]][1]))
  ggplot2::ggsave(file.path(outpath, "coefficient_densities_scalesfree.pdf"), d, width = 4, height = 3, useDingbats=FALSE)

  #WT contacts less VDW
  dg_dt[coef_order==2, WT_contacts_nVDW := WT_contacts - WT_vdw]

	#Weighted mean absolute value
  dg_dt[coef_order==2, dddg_abs_mut := abs(.SD[[1]]),.SDcols = c("mean_kcal/mol")]
  dg_dt[coef_order==2, dddg_abs := sum(abs(.SD[[1]])/.SD[[2]]^2, na.rm = T)/sum(1/.SD[[2]]^2, na.rm = T),Pos_ref,.SDcols = c("mean_kcal/mol", "std_kcal/mol")]
	dg_dt[coef_order==2, dddg_ci := sqrt(1/sum(1/.SD[[2]]^2, na.rm = T))*1.96*2,Pos_ref,.SDcols = c("mean_kcal/mol", "std_kcal/mol")]

  #Enrichment for top 10 couplings in pairs where at least one is a core residue
  top10_int <- dg_dt[coef_order==2][order(abs(mean), decreasing = T)][1:10, Pos1_class=='core' | Pos2_class=='core']
  remainder_int <- dg_dt[coef_order==2][order(abs(mean), decreasing = T)][11:.N, Pos1_class=='core' | Pos2_class=='core']
  test_results <- fisher.test(matrix(c(sum(top10_int==T), sum(top10_int==F), sum(remainder_int==T), sum(remainder_int==F)), nrow = 2))
  print(paste0("Enrichment for core residues in top 10 coupling pairs: Fisher's Exact Test P=", test_results$p.value, ", odds ratio=", test_results$estimate))

  #Enrichment for physical contacts in pairs where at least one is a core residue
  contact_pairs <- dg_dt[coef_order==2][scHAmin<8, Pos1_class=='core' | Pos2_class=='core']
  remainder_pairs <- dg_dt[coef_order==2][scHAmin>=8, Pos1_class=='core' | Pos2_class=='core']
  test_results <- fisher.test(matrix(c(sum(contact_pairs==T), sum(contact_pairs==F), sum(remainder_pairs==T), sum(remainder_pairs==F)), nrow = 2))
  print(paste0("Enrichment for core residues in physical contacting pairs: Fisher's Exact Test P=", test_results$p.value, ", odds ratio=", test_results$estimate))


  #Number of contacts per residue
  contact_list <- as.list(table(dg_dt[scHAmin<8,unlist(strsplit(Pos_ref, "_"))]))
  contact_dt <- dg_dt[coef_order==1]
  contact_dt[, num_contacts := 0]
  contact_dt[Pos_ref %in% names(contact_list), num_contacts := contact_list[[Pos_ref]],Pos_ref]
  contact_dt[RSASA<0.25,summary(num_contacts)]
  contact_dt[RSASA>=0.25,summary(num_contacts)]

  #Number
  dg_dt1 <- copy(dg_dt[coef_order==2])
  dg_dt2 <- copy(dg_dt[coef_order==2])
  dg_dt1[, Pos := sapply(strsplit(Pos_ref, "_"), '[', 1)]
  dg_dt2[, Pos := sapply(strsplit(Pos_ref, "_"), '[', 2)]
  merge_dt <- rbind(dg_dt1, dg_dt2)
  merge_dt[, dddg_abs := sum(abs(.SD[[1]])/.SD[[2]]^2, na.rm = T)/sum(1/.SD[[2]]^2, na.rm = T),Pos,.SDcols = c("mean_kcal/mol", "std_kcal/mol")]
  merge_dt[, dddg_ci := sqrt(1/sum(1/.SD[[2]]^2, na.rm = T))*1.96*2,Pos,.SDcols = c("mean_kcal/mol", "std_kcal/mol")]
  plot_dt <- data.table(
    )

	#Coupling scatter - calpha distance
	plot_dt <- dg_dt[coef_order==2][!duplicated(Pos_ref),.(dddg_abs, dddg_ci, calphamin, scHAmin)]
	plot_dt[is.na(dddg_abs), dddg_abs := 0]
	plot_dt[is.na(dddg_ci), dddg_ci := 0]
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(calphamin, dddg_abs)) +
    ggplot2::geom_point(ggplot2::aes(color = scHAmin)) +
    ggplot2::geom_linerange(data = plot_dt, ggplot2::aes(ymin = dddg_abs-dddg_ci/2, ymax = dddg_abs+dddg_ci/2, color = scHAmin), alpha = 1/4) +
    ggplot2::geom_vline(xintercept = 5, linetype = 2) + 
    ggplot2::geom_text(data = plot_dt[,.(label = paste("rho = ", round(cor(dddg_abs, calphamin, use = "pairwise.complete", method = 'spearman'), 2), sep=""))], ggplot2::aes(label=label, x=Inf, y=Inf, hjust = 1, vjust = 1)) +
    ggplot2::xlab("Inter-residue distance (Angstrom)") +
    ggplot2::ylab(bquote("|"*.(trait) ~ Delta*Delta*Delta*"G| (kcal/mol)")) +
    ggplot2::theme_classic() +
		ggplot2::scale_colour_viridis_c(option = "plasma")
  ggplot2::ggsave(file.path(outpath, "coupling_scatter_calphamin.pdf"), d, width = 4, height = 3, useDingbats=FALSE)

	#Coupling scatter - scHAmin distance
	plot_dt <- dg_dt[coef_order==2][!duplicated(Pos_ref),.(dddg_abs, dddg_ci, scHAmin)]
	plot_dt[is.na(dddg_abs), dddg_abs := 0]
	plot_dt[is.na(dddg_ci), dddg_ci := 0]
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(scHAmin, dddg_abs)) +
    ggplot2::geom_point(ggplot2::aes(color = scHAmin)) +
    ggplot2::geom_linerange(data = plot_dt, ggplot2::aes(ymin = dddg_abs-dddg_ci/2, ymax = dddg_abs+dddg_ci/2, color = scHAmin), alpha = 1/4) +
    ggplot2::geom_vline(xintercept = 5, linetype = 2) + 
    ggplot2::geom_text(data = plot_dt[,.(label = paste("rho = ", round(cor(dddg_abs, scHAmin, use = "pairwise.complete", method = 'spearman'), 2), sep=""))], ggplot2::aes(label=label, x=Inf, y=Inf, hjust = 1, vjust = 1)) +
    ggplot2::xlab("Inter-residue distance (Angstrom)") +
    ggplot2::ylab(bquote("|"*.(trait) ~ Delta*Delta*Delta*"G| (kcal/mol)")) +
    ggplot2::theme_classic() +
		ggplot2::scale_colour_viridis_c(option = "plasma")
  ggplot2::ggsave(file.path(outpath, "coupling_scatter_scHAmin.pdf"), d, width = 4, height = 3, useDingbats=FALSE)

  #Coupling scatter - scHAmin distance - backbone colour
  plot_dt <- dg_dt[coef_order==2][!duplicated(Pos_ref),.(dddg_abs, dddg_ci, scHAmin, backbone)]
  plot_dt[is.na(dddg_abs), dddg_abs := 0]
  plot_dt[is.na(dddg_ci), dddg_ci := 0]
  plot_dt[backbone==1, backbone_col := "1"]
  plot_dt[backbone>1 & backbone<=5, backbone_col := "2-5"]
  plot_dt[backbone>5, backbone_col := "6-10"]
  plot_dt[backbone>10, backbone_col := "11-20"]
  plot_dt[backbone>20, backbone_col := ">20"]
  plot_dt[, backbone_col := factor(backbone_col, levels = c("1", "2-5", "6-10", "11-20", ">20"))]
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(scHAmin, dddg_abs)) +
    ggplot2::geom_point(ggplot2::aes(color = backbone_col)) +
    ggplot2::geom_linerange(data = plot_dt, ggplot2::aes(ymin = dddg_abs-dddg_ci/2, ymax = dddg_abs+dddg_ci/2, color = backbone_col), alpha = 1/4) +
    ggplot2::geom_vline(xintercept = 5, linetype = 2) + 
    ggplot2::geom_text(data = plot_dt[,.(label = paste("rho = ", round(cor(dddg_abs, scHAmin, use = "pairwise.complete", method = 'spearman'), 2), sep=""))], ggplot2::aes(label=label, x=Inf, y=Inf, hjust = 1, vjust = 1)) +
    ggplot2::xlab("Inter-residue distance (Angstrom)") +
    ggplot2::ylab(bquote("|"*.(trait) ~ Delta*Delta*Delta*"G| (kcal/mol)")) +
    ggplot2::theme_classic() +
    ggplot2::scale_colour_viridis_d(option = "plasma")
  ggplot2::ggsave(file.path(outpath, "coupling_scatter_scHAmin_bbcol.pdf"), d, width = 4, height = 3, useDingbats=FALSE)

  #Coupling scatter - scHAmin distance - scHAmin colour
  plot_dt <- dg_dt[coef_order==2][!duplicated(Pos_ref),.(dddg_abs, dddg_ci, scHAmin, Pos1_class, Pos2_class, backbone)]
  plot_dt[is.na(dddg_abs), dddg_abs := 0]
  plot_dt[is.na(dddg_ci), dddg_ci := 0]
  plot_dt[scHAmin<5, scHAmin_col := "<5A"]
  plot_dt[scHAmin>=5, scHAmin_col := "5-10"]
  plot_dt[scHAmin>=10, scHAmin_col := "10-20"]
  plot_dt[scHAmin>=20, scHAmin_col := ">20"]
  plot_dt[, scHAmin_col := factor(scHAmin_col, levels = c("<5A", "5-10", "10-20", ">20"))]
  plot_dt[, Pos_class_plot := "remainder"]
  plot_dt[Pos1_class=='core' | Pos2_class=='core', Pos_class_plot := "core"]
  plot_dt[, Pos_class_plot := factor(Pos_class_plot, levels = c("remainder", "core"))]
  plot_cols <- c(colour_scheme[["shade 0"]][1], colour_scheme[["shade 1"]][1], colour_scheme[["shade 1"]][3], colour_scheme[["shade 0"]][3])
  names(plot_cols) <- c("<5A", "5-10", "10-20", ">20")
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(scHAmin, dddg_abs)) +
    ggplot2::geom_point(ggplot2::aes(color = scHAmin_col, shape = Pos_class_plot)) +
    ggplot2::geom_linerange(data = plot_dt, ggplot2::aes(ymin = dddg_abs-dddg_ci/2, ymax = dddg_abs+dddg_ci/2, color = scHAmin_col), alpha = 1/4) +
    ggplot2::geom_vline(xintercept = 5, linetype = 2) + 
    ggplot2::geom_text(data = plot_dt[backbone>5,.(label = paste("rho = ", round(cor(dddg_abs, scHAmin, use = "pairwise.complete", method = 'spearman'), 2), sep=""))], ggplot2::aes(label=label, x=Inf, y=Inf, hjust = 1, vjust = 1)) +
    ggplot2::geom_text(data = plot_dt[,.(label = paste("rho = ", round(cor(dddg_abs, scHAmin, use = "pairwise.complete", method = 'spearman'), 2), sep=""))], ggplot2::aes(label=label, x=Inf, y=-Inf, hjust = 1, vjust = 0)) +
    ggplot2::xlab("Inter-residue distance (Angstrom)") +
    ggplot2::ylab(bquote("|"*.(trait) ~ Delta*Delta*Delta*"G| (kcal/mol)")) +
    ggplot2::theme_classic() +
    ggplot2::scale_colour_manual(values = plot_cols)
  ggplot2::ggsave(file.path(outpath, "coupling_scatter_scHAmin_scHAmincol_coreshape.pdf"), d, width = 4, height = 3, useDingbats=FALSE)

  #Coupling scatter - scHAmin distance - ligand scHAmin distance colour
  plot_dt <- dg_dt[coef_order==2][!duplicated(Pos_ref),.(dddg_abs, dddg_ci, scHAmin, Pos1_class, Pos2_class, backbone, scHAmin_ligand1, scHAmin_ligand2)]
  plot_dt[is.na(dddg_abs), dddg_abs := 0]
  plot_dt[is.na(dddg_ci), dddg_ci := 0]
  plot_dt[, scHAmin_ligmean := apply(.SD, 1, mean),,.SDcols = c("scHAmin_ligand1", "scHAmin_ligand2")]
  plot_dt[, Pos_class_plot := "remainder"]
  plot_dt[Pos1_class=='core' | Pos2_class=='core', Pos_class_plot := "core"]
  plot_dt[, Pos_class_plot := factor(Pos_class_plot, levels = c("remainder", "core"))]
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(scHAmin, dddg_abs)) +
    ggplot2::geom_point(ggplot2::aes(color = scHAmin_ligmean, shape = Pos_class_plot)) +
    ggplot2::geom_linerange(data = plot_dt, ggplot2::aes(ymin = dddg_abs-dddg_ci/2, ymax = dddg_abs+dddg_ci/2, color = scHAmin_ligmean), alpha = 1/4) +
    ggplot2::geom_vline(xintercept = 5, linetype = 2) + 
    ggplot2::geom_text(data = plot_dt[backbone>5,.(label = paste("rho = ", round(cor(dddg_abs, scHAmin, use = "pairwise.complete", method = 'spearman'), 2), sep=""))], ggplot2::aes(label=label, x=Inf, y=Inf, hjust = 1, vjust = 1)) +
    ggplot2::geom_text(data = plot_dt[,.(label = paste("rho = ", round(cor(dddg_abs, scHAmin, use = "pairwise.complete", method = 'spearman'), 2), sep=""))], ggplot2::aes(label=label, x=Inf, y=-Inf, hjust = 1, vjust = 0)) +
    ggplot2::xlab("Inter-residue distance (Angstrom)") +
    ggplot2::ylab(bquote("|"*.(trait) ~ Delta*Delta*Delta*"G| (kcal/mol)")) +
    ggplot2::theme_classic() +
    ggplot2::scale_colour_viridis_c()
  ggplot2::ggsave(file.path(outpath, "coupling_scatter_scHAmin_scHAminligandcol_coreshape.pdf"), d, width = 4, height = 3, useDingbats=FALSE)

  #Coupling scatter - scHAmin distance - backbone colour
  plot_dt <- dg_dt[coef_order==2][!duplicated(Pos_ref),.(dddg_abs, dddg_ci, scHAmin, Pos1_class, Pos2_class, backbone)]
  plot_dt[is.na(dddg_abs), dddg_abs := 0]
  plot_dt[is.na(dddg_ci), dddg_ci := 0]
  plot_dt[scHAmin<5, scHAmin_col := "<5A"]
  plot_dt[scHAmin>=5, scHAmin_col := "5-10"]
  plot_dt[scHAmin>=10, scHAmin_col := "10-20"]
  plot_dt[scHAmin>=20, scHAmin_col := ">20"]
  plot_dt[, scHAmin_col := factor(scHAmin_col, levels = c("<5A", "5-10", "10-20", ">20"))]
  plot_dt[, Pos_class_plot := "remainder"]
  plot_dt[Pos1_class=='core' | Pos2_class=='core', Pos_class_plot := "core"]
  plot_dt[, Pos_class_plot := factor(Pos_class_plot, levels = c("remainder", "core"))]
  plot_cols <- c(colour_scheme[["shade 0"]][1], colour_scheme[["shade 1"]][1], colour_scheme[["shade 1"]][3], colour_scheme[["shade 0"]][3])
  names(plot_cols) <- c("<5A", "5-10", "10-20", ">20")
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(backbone, dddg_abs)) +
    ggplot2::geom_point(ggplot2::aes(color = scHAmin_col, shape = Pos_class_plot)) +
    ggplot2::geom_linerange(data = plot_dt, ggplot2::aes(ymin = dddg_abs-dddg_ci/2, ymax = dddg_abs+dddg_ci/2, color = scHAmin_col), alpha = 1/4) +
    ggplot2::geom_vline(xintercept = 5, linetype = 2) + 
    ggplot2::geom_text(data = plot_dt[scHAmin>5,.(label = paste("rho = ", round(cor(dddg_abs, backbone, use = "pairwise.complete", method = 'spearman'), 2), sep=""))], ggplot2::aes(label=label, x=Inf, y=Inf, hjust = 1, vjust = 1)) +
    ggplot2::geom_text(data = plot_dt[,.(label = paste("rho = ", round(cor(dddg_abs, backbone, use = "pairwise.complete", method = 'spearman'), 2), sep=""))], ggplot2::aes(label=label, x=Inf, y=-Inf, hjust = 1, vjust = 0)) +
    ggplot2::xlab("Backbone distance (residues)") +
    ggplot2::ylab(bquote("|"*.(trait) ~ Delta*Delta*Delta*"G| (kcal/mol)")) +
    ggplot2::theme_classic() +
    ggplot2::scale_colour_manual(values = plot_cols)
  ggplot2::ggsave(file.path(outpath, "coupling_scatter_backbone_scHAmincol_coreshape.pdf"), d, width = 4, height = 3, useDingbats=FALSE)

  #Coupling scatter - scHAmin distance - position class colour
  plot_dt <- dg_dt[coef_order==2][!duplicated(Pos_ref),.(dddg_abs, dddg_ci, scHAmin, Pos1_class, Pos2_class)]
  plot_dt[is.na(dddg_abs), dddg_abs := 0]
  plot_dt[is.na(dddg_ci), dddg_ci := 0]
  plot_dt[, Pos_class_plot := "surface"]
  plot_dt[Pos1_class=='binding_interface' | Pos2_class=='binding_interface', Pos_class_plot := "binding_interface"]
  plot_dt[Pos1_class=='core' | Pos2_class=='core', Pos_class_plot := "core"]
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(scHAmin, dddg_abs)) +
    ggplot2::geom_point(ggplot2::aes(color = Pos_class_plot)) +
    ggplot2::geom_linerange(data = plot_dt, ggplot2::aes(ymin = dddg_abs-dddg_ci/2, ymax = dddg_abs+dddg_ci/2, color = Pos_class_plot), alpha = 1/4) +
    ggplot2::geom_vline(xintercept = 5, linetype = 2) + 
    ggplot2::geom_text(data = plot_dt[,.(label = paste("rho = ", round(cor(dddg_abs, scHAmin, use = "pairwise.complete", method = 'spearman'), 2), sep=""))], ggplot2::aes(label=label, x=Inf, y=Inf, hjust = 1, vjust = 1)) +
    ggplot2::xlab("Inter-residue distance (Angstrom)") +
    ggplot2::ylab(bquote("|"*.(trait) ~ Delta*Delta*Delta*"G| (kcal/mol)")) +
    ggplot2::theme_classic() +
    ggplot2::scale_colour_viridis_d(option = "plasma")
  ggplot2::ggsave(file.path(outpath, "coupling_scatter_scHAmin_posclasscol.pdf"), d, width = 4, height = 3, useDingbats=FALSE)

  #Coupling scatter - scHAmin distance - WT contacts
  plot_dt <- dg_dt[coef_order==2][!duplicated(Pos_ref),.(dddg_abs, dddg_ci, scHAmin, WT_contacts)]
  plot_dt[is.na(dddg_abs), dddg_abs := 0]
  plot_dt[is.na(dddg_ci), dddg_ci := 0]
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(scHAmin, dddg_abs)) +
    ggplot2::geom_point(ggplot2::aes(color = WT_contacts)) +
    ggplot2::geom_linerange(data = plot_dt, ggplot2::aes(ymin = dddg_abs-dddg_ci/2, ymax = dddg_abs+dddg_ci/2, color = WT_contacts), alpha = 1/4) +
    ggplot2::geom_vline(xintercept = 5, linetype = 2) + 
    ggplot2::geom_text(data = plot_dt[,.(label = paste("rho = ", round(cor(dddg_abs, scHAmin, use = "pairwise.complete", method = 'spearman'), 2), sep=""))], ggplot2::aes(label=label, x=Inf, y=Inf, hjust = 1, vjust = 1)) +
    ggplot2::xlab("Inter-residue distance (Angstrom)") +
    ggplot2::ylab(bquote("|"*.(trait) ~ Delta*Delta*Delta*"G| (kcal/mol)")) +
    ggplot2::theme_classic() +
    ggplot2::scale_colour_viridis_c(option = "plasma")
  ggplot2::ggsave(file.path(outpath, "coupling_scatter_scHAmin_wtcontactscol.pdf"), d, width = 4, height = 3, useDingbats=FALSE)

  #Coupling scatter - scHAmin distance - WT contacts
  plot_dt <- dg_dt[coef_order==2][!duplicated(Pos_ref),.(dddg_abs, dddg_ci, scHAmin, WT_contacts_nVDW)]
  plot_dt[is.na(dddg_abs), dddg_abs := 0]
  plot_dt[is.na(dddg_ci), dddg_ci := 0]
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(scHAmin, dddg_abs)) +
    ggplot2::geom_point(ggplot2::aes(color = WT_contacts_nVDW)) +
    ggplot2::geom_linerange(data = plot_dt, ggplot2::aes(ymin = dddg_abs-dddg_ci/2, ymax = dddg_abs+dddg_ci/2, color = WT_contacts_nVDW), alpha = 1/4) +
    ggplot2::geom_vline(xintercept = 5, linetype = 2) + 
    ggplot2::geom_text(data = plot_dt[,.(label = paste("rho = ", round(cor(dddg_abs, scHAmin, use = "pairwise.complete", method = 'spearman'), 2), sep=""))], ggplot2::aes(label=label, x=Inf, y=Inf, hjust = 1, vjust = 1)) +
    ggplot2::xlab("Inter-residue distance (Angstrom)") +
    ggplot2::ylab(bquote("|"*.(trait) ~ Delta*Delta*Delta*"G| (kcal/mol)")) +
    ggplot2::theme_classic() +
    ggplot2::scale_colour_viridis_c(option = "plasma")
  ggplot2::ggsave(file.path(outpath, "coupling_scatter_scHAmin_wtcontactsnVDWcol.pdf"), d, width = 4, height = 3, useDingbats=FALSE)

	#Coupling scatter - scHAmin distance no adjacent
	plot_dt <- dg_dt[coef_order==2 & backbone>=5][!duplicated(Pos_ref),.(dddg_abs, dddg_ci, scHAmin)]
	plot_dt[is.na(dddg_abs), dddg_abs := 0]
	plot_dt[is.na(dddg_ci), dddg_ci := 0]
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(scHAmin, dddg_abs)) +
    ggplot2::geom_point(ggplot2::aes(color = scHAmin)) +
    ggplot2::geom_linerange(data = plot_dt, ggplot2::aes(ymin = dddg_abs-dddg_ci/2, ymax = dddg_abs+dddg_ci/2, color = scHAmin), alpha = 1/4) +
    ggplot2::geom_smooth(data = plot_dt[scHAmin>=5], method = "lm", formula = 'y ~ x') +
    ggplot2::geom_vline(xintercept = 5, linetype = 2) + 
    ggplot2::geom_text(data = plot_dt[,.(label = paste("rho = ", round(cor(dddg_abs, scHAmin, use = "pairwise.complete", method = 'spearman'), 2), sep=""))], ggplot2::aes(label=label, x=Inf, y=Inf, hjust = 1, vjust = 1)) +
    ggplot2::xlab("Inter-residue distance (Angstrom)") +
    ggplot2::ylab(bquote("|"*.(trait) ~ Delta*Delta*Delta*"G| (kcal/mol)")) +
    ggplot2::theme_classic() +
		ggplot2::scale_colour_viridis_c(option = "plasma")
  ggplot2::ggsave(file.path(outpath, "coupling_scatter_scHAmin_geq5res.pdf"), d, width = 4, height = 3, useDingbats=FALSE)

	plot_dt <- dg_dt[coef_order==2][!duplicated(Pos_ref),.(dddg_abs, dddg_ci, backbone, scHAmin)]
	plot_dt[is.na(dddg_abs), dddg_abs := 0]
	plot_dt[is.na(dddg_ci), dddg_ci := 0]
	plot_dt[, dddg_abs_resid := lm(dddg_abs~backbone)$residuals]
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(scHAmin, dddg_abs_resid)) +
    ggplot2::geom_point(ggplot2::aes(color = scHAmin)) +
    ggplot2::geom_linerange(data = plot_dt, ggplot2::aes(ymin = dddg_abs_resid-dddg_ci/2, ymax = dddg_abs_resid+dddg_ci/2, color = scHAmin), alpha = 1/4) +
    ggplot2::geom_vline(xintercept = 5, linetype = 2) + 
    ggplot2::geom_text(data = plot_dt[,.(label = paste("rho = ", round(cor(dddg_abs_resid, scHAmin, use = "pairwise.complete", method = 'spearman'), 2), sep=""))], ggplot2::aes(label=label, x=Inf, y=Inf, hjust = 1, vjust = 1)) +
    ggplot2::xlab("Inter-residue distance (Angstrom)") +
    ggplot2::ylab(bquote("Residual |"*.(trait) ~ Delta*Delta*Delta*"G| (kcal/mol)")) +
    ggplot2::theme_classic() +
		ggplot2::scale_colour_viridis_c(option = "plasma")
  ggplot2::ggsave(file.path(outpath, "coupling_scatter_scHAmin_backboneresid.pdf"), d, width = 4, height = 3, useDingbats=FALSE)

	#Coupling scatter - backbone distance
	plot_dt <- dg_dt[coef_order==2][!duplicated(Pos_ref),.(dddg_abs, dddg_ci, backbone, scHAmin)]
	plot_dt[is.na(dddg_abs), dddg_abs := 0]
	plot_dt[is.na(dddg_ci), dddg_ci := 0]
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(backbone, dddg_abs)) +
    ggplot2::geom_point(ggplot2::aes(color = scHAmin)) +
    ggplot2::geom_linerange(data = plot_dt, ggplot2::aes(ymin = dddg_abs-dddg_ci/2, ymax = dddg_abs+dddg_ci/2, color = scHAmin), alpha = 1/4) +
    ggplot2::geom_vline(xintercept = 5, linetype = 2) + 
    ggplot2::geom_text(data = plot_dt[,.(label = paste("rho = ", round(cor(dddg_abs, backbone, use = "pairwise.complete", method = 'spearman'), 2), sep=""))], ggplot2::aes(label=label, x=Inf, y=Inf, hjust = 1, vjust = 1)) +
    ggplot2::xlab("Backbone distance (residues)") +
    ggplot2::ylab(bquote("|"*.(trait) ~ Delta*Delta*Delta*"G| (kcal/mol)")) +
    ggplot2::theme_classic() +
		ggplot2::scale_colour_viridis_c(option = "plasma")
  ggplot2::ggsave(file.path(outpath, "coupling_scatter_backbone.pdf"), d, width = 4, height = 3, useDingbats=FALSE)

	#Coupling scatter - backbone distance no contact
	plot_dt <- dg_dt[coef_order==2 & scHAmin>=5][!duplicated(Pos_ref),.(dddg_abs, dddg_ci, backbone, scHAmin)]
	plot_dt[is.na(dddg_abs), dddg_abs := 0]
	plot_dt[is.na(dddg_ci), dddg_ci := 0]
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(backbone, dddg_abs)) +
    ggplot2::geom_point(ggplot2::aes(color = scHAmin)) +
    ggplot2::geom_linerange(data = plot_dt, ggplot2::aes(ymin = dddg_abs-dddg_ci/2, ymax = dddg_abs+dddg_ci/2, color = scHAmin), alpha = 1/4) +
    ggplot2::geom_smooth(data = plot_dt[backbone>=5], method = "lm", formula = 'y ~ x') +
    ggplot2::geom_vline(xintercept = 5, linetype = 2) + 
    ggplot2::geom_text(data = plot_dt[,.(label = paste("rho = ", round(cor(dddg_abs, backbone, use = "pairwise.complete", method = 'spearman'), 2), sep=""))], ggplot2::aes(label=label, x=Inf, y=Inf, hjust = 1, vjust = 1)) +
    ggplot2::xlab("Backbone distance (residues)") +
    ggplot2::ylab(bquote("|"*.(trait) ~ Delta*Delta*Delta*"G| (kcal/mol)")) +
    ggplot2::theme_classic() +
		ggplot2::scale_colour_viridis_c(option = "plasma")
  ggplot2::ggsave(file.path(outpath, "coupling_scatter_backbone_geq5A.pdf"), d, width = 4, height = 3, useDingbats=FALSE)

	#Coupling scatter - backbone distance no contact
	plot_dt <- dg_dt[coef_order==2][!duplicated(Pos_ref),.(dddg_abs, dddg_ci, backbone, scHAmin)]
	plot_dt[is.na(dddg_abs), dddg_abs := 0]
	plot_dt[is.na(dddg_ci), dddg_ci := 0]
	plot_dt[, dddg_abs_resid := lm(dddg_abs~scHAmin)$residuals]
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(backbone, dddg_abs_resid)) +
    ggplot2::geom_point(ggplot2::aes(color = scHAmin)) +
    ggplot2::geom_linerange(data = plot_dt, ggplot2::aes(ymin = dddg_abs_resid-dddg_ci/2, ymax = dddg_abs_resid+dddg_ci/2, color = scHAmin), alpha = 1/4) +
    ggplot2::geom_vline(xintercept = 5, linetype = 2) + 
    ggplot2::geom_text(data = plot_dt[,.(label = paste("rho = ", round(cor(dddg_abs_resid, backbone, use = "pairwise.complete", method = 'spearman'), 2), sep=""))], ggplot2::aes(label=label, x=Inf, y=Inf, hjust = 1, vjust = 1)) +
    ggplot2::xlab("Backbone distance (residues)") +
    ggplot2::ylab(bquote("|"*.(trait) ~ Delta*Delta*Delta*"G| (kcal/mol)")) +
    ggplot2::theme_classic() +
		ggplot2::scale_colour_viridis_c(option = "plasma")
  ggplot2::ggsave(file.path(outpath, "coupling_scatter_backbone_scHAminresid.pdf"), d, width = 4, height = 3, useDingbats=FALSE)

  if(dg_dt[,max(table(Pos_ref))>10]){
    #Coupling scatter - scHAmin distance
    plot_dt <- dg_dt[coef_order==2][,.(dddg_abs_mut, scHAmin)]
    plot_dt[is.na(dddg_abs_mut), dddg_abs_mut := 0]
    plot_dt[, scHAmin_col := as.character(round(scHAmin, 1))]
    cc <- c(colour_scheme[["shade 0"]][3], colour_scheme[["shade 2"]][3], colour_scheme[["shade 1"]][3], colour_scheme[["shade 1"]][1], colour_scheme[["shade 2"]][1], colour_scheme[["shade 0"]][1])
    names(cc) <- unique(plot_dt[order(scHAmin_col, decreasing = T),scHAmin_col])
    d <- ggplot2::ggplot(plot_dt,ggplot2::aes(scHAmin_col, dddg_abs_mut)) +
      ggplot2::geom_boxplot(ggplot2::aes(fill = scHAmin_col), notch = T, outlier.shape = NA) +
      ggplot2::scale_fill_manual(values = cc) +
      ggplot2::geom_text(data = plot_dt[,.(label = paste("rho = ", round(cor(dddg_abs_mut, scHAmin, use = "pairwise.complete", method = 'spearman'), 2), sep=""))], ggplot2::aes(label=label, x=Inf, y=Inf, hjust = 1, vjust = 1)) +
      ggplot2::xlab("Inter-residue distance (Angstrom)") +
    ggplot2::ylab(bquote("|"*.(trait) ~ Delta*Delta*Delta*"G| (kcal/mol)")) +
      ggplot2::theme_classic() +
      ggplot2::coord_cartesian(ylim = c(0, 1))
    ggplot2::ggsave(file.path(outpath, "coupling_scatter_scHAmin_boxplot.pdf"), d, width = 4, height = 3, useDingbats=FALSE)

    #Coupling scatter - backbone distance
    plot_dt <- dg_dt[coef_order==2][,.(dddg_abs_mut, backbone, scHAmin)]
    plot_dt[is.na(dddg_abs_mut), dddg_abs_mut := 0]
    plot_dt[, scHAmin_col := as.character(round(scHAmin, 1))]
    cc <- c(colour_scheme[["shade 0"]][3], colour_scheme[["shade 2"]][3], colour_scheme[["shade 1"]][3], colour_scheme[["shade 1"]][1], colour_scheme[["shade 2"]][1], colour_scheme[["shade 0"]][1])
    names(cc) <- unique(plot_dt[order(scHAmin, decreasing = T),scHAmin_col])
    d <- ggplot2::ggplot(plot_dt,ggplot2::aes(as.factor(backbone), dddg_abs_mut)) +
      ggplot2::geom_boxplot(ggplot2::aes(fill = scHAmin_col), notch = T, outlier.shape = NA) +
      ggplot2::scale_fill_manual(values = cc) +
      ggplot2::geom_text(data = plot_dt[,.(label = paste("rho = ", round(cor(dddg_abs_mut, backbone, use = "pairwise.complete", method = 'spearman'), 2), sep=""))], ggplot2::aes(label=label, x=Inf, y=Inf, hjust = 1, vjust = 1)) +
      ggplot2::xlab("Backbone distance (residues)") +
    ggplot2::ylab(bquote("|"*.(trait) ~ Delta*Delta*Delta*"G| (kcal/mol)")) +
      ggplot2::theme_classic() +
      ggplot2::coord_cartesian(ylim = c(0, 1))
    ggplot2::ggsave(file.path(outpath, "coupling_scatter_backbone_boxplot.pdf"), d, width = 4, height = 3, useDingbats=FALSE)
  }

  return(dg_dt)
}
