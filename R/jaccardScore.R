#' Calculate Jaccard Scores between a list of peaks
#'
#' @param peaks List of GRanges objects representing different sets of peaks.
#' @import GenomicRanges
#' @return A data.frame containing the different comparisons between peaks and
#' the Jaccard Scores, as calculated by bedtools. As a result, the final
#' statistic ranges from 0.0 to 1.0, where 0.0 represents no overlap and 1.0
#' represent complete overlap.
#' @export
jaccardScore <- function(peaks) {

  # Make all possible permutations of samples
  comparisons <- expand.grid(names(peaks), names(peaks))

  comparisons$jaccard <- mapply(.calculateJaccardScore, x=peaks[comparisons$Var1], y=peaks[comparisons$Var2])
  return(comparisons)
}


.calculateJaccardScore <- function(x, y) {
  score <- sum(width(intersect(x, y))) /
    (sum(width(x)) + sum(width(y)) - sum(width(intersect(x, y))))
  return(score)
}
