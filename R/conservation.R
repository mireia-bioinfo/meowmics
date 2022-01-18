#' Get conservation scores
#'
#' Given a peak set, this function obtains the mean conservation scores in the given window and bins.
#' @param peak_file Path for the peak file or GRanges object.
#' @param id_col Number or name of `mcol` that contains unique identifier for each peak. Default: 1.
#' @param cons_score Annotation to use for getting the conservation scores. See \code{\link[GenomicScores]{gscores}}
#' for details.
#' @param bin_width Number of base pairs for binning the peak and obtaining the conservation scores. Default: 20
#' @param window_width Number of base pairs of a window centered in the center of the peak where to focus the
#' conservation analysis. Default = 2000
#' @param merge_fun Function for summarizing scores in the same region.
#' @param random_control Logical indicating whether to perform the same analysis using a randomized set of peaks
#' as a control. Default = TRUE
#' @param genome Character indicating the name of the genome where to randomize the regions.
#' See \code{\link[regioneR]{randomizeRegions}} for details.
#' @param per.chromosome Logical indicating if the randomization should be performed along the same chromosome as
#' the original set of peaks. See \code{\link[regioneR]{randomizeRegions}} for details.
#' @param summarise Logaical indicating whether to summarise the results, computing the mean of all peaks at a specific
#' position. Default: TRUE.
#' @return It returns a \code{data.frame} containing the summarized phastCons conservation scores for the file
#' and randomized control (if \code{random_control = TRUE}) at each relative position from the peak center.
#' The summarization is performed by computing the mean for all peaks in that specific position. Standard deviation
#' is also returned.
#' @export
get_conservation_scores <- function(peak_file,
                                    id_col=1,
                                    cons_score=phastCons7way.UCSC.hg38::phastCons7way.UCSC.hg38,
                                    bin_width=20,
                                    window_width=2e3,
                                    merge_fun="mean",
                                    random_control=TRUE,
                                    genome="hg38",
                                    mask=NULL,
                                    per.chromosome=FALSE,
                                    summarise=TRUE) {

  peaks <- regioneR::toGRanges(peak_file)
  peaks <- GenomicRanges::resize(peaks, width=window_width, fix="center")
  peaks_bin <- GenomicRanges::tile(peaks, width=bin_width)
  pos <- seq(-window_width/2, window_width/2-1, by=bin_width)

  peaks_unl <- unlist(peaks_bin)
  peaks_unl$id <- rep(paste0("sample_peak_", 1:length(peaks)), each=length(pos))
  peaks_unl$pos <- rep(pos, length(peaks_bin))

  scores <- GenomicScores::gscores(cons_score,
                                   peaks_unl,
                                   summaryFun=merge_fun)

  if (random_control) {
    rndm <- regioneR::randomizeRegions(peaks,
                                       genome=genome,
                                       mask=mask,
                                       per.chromosome=per.chromosome)
    rndm <- rndm[order(rndm),]
    rndm_bin <- GenomicRanges::tile(rndm, width=bin_width)
    rndm_unl <- unlist(rndm_bin)
    rndm_unl$id <- rep(paste0("random_peak_", 1:length(rndm)), each=length(pos))
    rndm_unl$pos <- rep(pos, length(rndm_bin))

    scores_rndm <- GenomicScores::gscores(cons_score,
                                          rndm_unl)

    scores <- c(scores, scores_rndm)
  }

  df <- data.frame(mcols(scores))
  df$type <- gsub("_peak_[[:digit:]]*", "", df$id)
  df$peakID <- rep(mcols(peaks)[,1], each=length(pos))
  df <- df[,-1]

  ## Summarise by type of sample
  if(summarise) {
    cons <- df %>%
      dplyr::group_by(type, pos) %>%
      dplyr::summarise(phastCons=mean(default, na.rm=TRUE),
                       sd=sd(default, na.rm=TRUE))
  } else {
    cons <- df
  }

  return(cons)
}

#' Get mean removing NAs
#'
#' @export
mean_rm <- function(x, ...) {
  mean(x, na.rm=TRUE, ...)
}
