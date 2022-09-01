#' Bin regions using a fixed number of bins
#'
#' @param gr `GRanges` containing the regions to be divided.
#' @param n_bins Numeric incating the number of bins to divide the regions into
#' (including the extended coordinates).
#' @param n_bins_exp Numeric indicating the value by which to expand the regions
#' (`width(gr)*n_bins_exp`).
#' @param id_col Numeric or character indicating the mcol number or name of the
#' unique identifier for the regions.
#' @return `GRanges` object with the peaks extended to `width(gr)*n_bins_exp`
#' and divided into a fixed number of bins.
#' @export
bin_regions_nbins <- function(gr,
                              n_bins=100,
                              n_bins_exp=2,
                              id_col="peakID") {
  ## Resize peaks to n_bins_exp
  peaks_rsz <- GenomicRanges::resize(gr, width=width(gr)*n_bins_exp, fix="center")

  ## Divide in a fixed number of bins
  peaks_tile <- GenomicRanges::tile(peaks_rsz, n=n_bins)

  ## Vector of positions
  # Obtain bin values from first instance
  bin_start <- queryHits(findOverlaps(peaks_tile[[1]], resize(gr[1], 1, fix="start")))
  bin_end <- queryHits(findOverlaps(peaks_tile[[1]], resize(gr[1], 1, fix="end")))
  pos_center <- seq(0, 1, length.out=bin_end-bin_start+1)

  # Value of each bin
  unit <- pos_center[2]

  # Value for bins upstream of peak
  pos_prev <- rev(-unit * 1:(bin_start-1))

  # Value for bins downstream of peak
  pos_post <- 1 + unit * 1:(length((bin_end+1):n_bins))

  # Vector of positions
  pos_vector <- c(pos_prev, pos_center, pos_post)

  ## Create final GRanges
  peaks_unl <- unlist(peaks_tile) %>%
    plyranges::mutate(PeakID=rep(mcols(gr)[,id_col], each=n_bins),
                      pos=rep(pos_vector, length(gr)),
                      id_pos=paste0(PeakID, "_", pos)) %>%
    plyranges::filter(start > 0 & end > 0) %>% ## Remove - bins
    plyranges::mutate(start = ifelse(start < 1, 1, start)) ## Resize bins outside scope

  return(peaks_unl)
}
