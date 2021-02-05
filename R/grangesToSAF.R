grangesToSAF <- function(granges, add_cols=NA) {
  saf <- data.frame(granges)
  saf <- saf[,c(6, 1:3, 5, add_cols)]
  colnames(saf) <- c("GeneID", "Chr", "Start", "End", "Strand", colnames(saf)[add_cols])
  return(saf)
}
