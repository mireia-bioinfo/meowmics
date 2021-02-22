grangesToSAF <- function(granges, add_cols=NA) {
  saf <- data.frame(granges)

  if(length(unique(saf[,6]))!=nrow(saf)) stop("First mcol should be a unique identifier.")

  if (is.na(add_cols)) {
    saf <- saf[,c(6, 1:3, 5)]
    colnames(saf) <- c("GeneID", "Chr", "Start", "End", "Strand")
  } else {
    saf <- saf[,c(6, 1:3, 5, add_cols)]
    colnames(saf) <- c("GeneID", "Chr", "Start", "End", "Strand", colnames(saf)[add_cols])
  }

  saf$Strand[saf$Strand=="*"] <- "+"

  return(saf)
}
