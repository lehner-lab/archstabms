
#' archstabms__get_codons
#'
#' Split nucleotide sequence into codons.
#'
#' @param nt_seq nucleotide sequence string (required)
#'
#' @return Vector of codons
#' @export
archstabms__get_codons <- function(
  nt_seq){
  return(apply(matrix(unlist(strsplit(nt_seq, "")), ncol = 3, byrow = T), 1, paste0, collapse = ""))
}
