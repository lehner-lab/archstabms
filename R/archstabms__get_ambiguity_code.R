
#' archstabms__get_ambiguity_code
#'
#' Get IUPAC ambiguity codes.
#'
#' @param input_nucs vector of nucleotides (required)
#'
#' @return ambiguity code 
#' @export
archstabms__get_ambiguity_code <- function(input_nucs){
  nuc_codes <- list(
    "A" = "A",
    "C" = "C",
    "G" = "G",
    "T" = "T",
    "R" = c("A", "G"),
    "Y" = c("C", "T"),
    "S" = c("C", "G"),
    "W" = c("A", "T"),
    "K" = c("G", "T"),
    "M" = c("A", "C"),
    "B" = c("C", "G", "T"),
    "D" = c("A", "G", "T"),
    "H" = c("A", "C", "T"),
    "V" = c("A", "C", "G"),
    "N" = c("A", "C", "G", "T")
    )
  input_nucs <- unlist(strsplit(input_nucs, ""))
  for(acode in names(nuc_codes)){
    if(sum(input_nucs %in% nuc_codes[[acode]])==length(nuc_codes[[acode]]) & length(nuc_codes[[acode]])==length(input_nucs)){
      return(acode)
    }
  }
}
