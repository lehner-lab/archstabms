
#' archstabms__reachable_single_AA_mutants
#'
#' Determine single AA mutants reachable from single nucleotide substitutions.
#'
#' @param wt_coding_seq WT coding sequence string (required)
#' @param position_offset residue position offset (default:0)
#'
#' @return List of reachable mutant IDs and a dictionary mapping IDs to codons 
#' @export
archstabms__reachable_single_AA_mutants <- function(
  wt_coding_seq,
  position_offset = 0){

  #Mutate WT seq to all possible singles
  mut_list <- list()
  mut_dict <- list()
  #For all codons
  for(i in 1:(nchar(wt_coding_seq)/3)){
    this_codon <- archstabms__get_codons(wt_coding_seq)[i]
    #For all positions in this codon
    for(j in 1:3){
      for(k in c("A", "C", "G", "T")){
        mut_codon <- unlist(strsplit(this_codon, ""))
        wt_AA <- paste0(seqinr::translate(mut_codon), collapse="")
        mut_codon[j] <- k
        mut_AA <- paste0(seqinr::translate(mut_codon), collapse="")
        if(mut_AA!="*"){
          this_mut_code <- paste0(wt_AA, i+position_offset, mut_AA)
          if(!this_mut_code %in% names(mut_dict)){
            mut_dict[[this_mut_code]] <- paste0(mut_codon, collapse = "")
          }else{
            mut_dict[[this_mut_code]] <- unique(c(mut_dict[[this_mut_code]], paste0(mut_codon, collapse = "")))
          }
          mut_list <- append(mut_list, this_mut_code)
        }
      }
    }
  }
  reachable_muts <- unique(unlist(mut_list))
  return(list(
    reachable_muts = reachable_muts,
    reachable_muts_codons = mut_dict))
}
