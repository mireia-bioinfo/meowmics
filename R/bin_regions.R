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
  # Resize peaks to n_bins_exp
  peaks_rsz <- GenomicRanges::resize(gr, width=width(gr)*n_bins_exp, fix="center")

  # Divide in a fixed number of bins``
  peaks_tile <- GenomicRanges::tile(peaks_rsz, n=n_bins)
  peaks_unl <- unlist(peaks_tile)
  peaks_unl$PeakID <- rep(mcols(gr)[,id_col], each=n_bins)

  # Identify start
  starts <- queryHits(findOverlaps(peaks_unl, promoters(gr, 0, 1)))
  first <- !duplicated(peaks_unl$PeakID[starts])
  starts <- starts[first]

  # Identify ends
  ends <- queryHits(findOverlaps(peaks_unl, resize(gr, 1, fix="end")))
  last <- rev(!duplicated(peaks_unl$PeakID[rev(ends)]))
  ends <- ends[last]

  # Position of start and end in vector
  bins_ini <- starts[1]-1
  bins_end <- starts[2] - ends[1] - bins_ini

  # Scales values to define start and end (center) and upstream and downstream
  # positions
  center <- seq(0, 1, length.out=unique(ends-starts+1)[1])
  prev <- seq(center[2]*bins_ini, center[2], length.out=bins_ini)
  post <- seq(1+center[2], (1+center[2]*bins_end), length.out=bins_end-1)

  # Final vector of positions
  pos <- c(
    -prev,
    center,
    post
  )

  # Add position value
  peaks_unl$pos <- rep(pos, length(peaks_tile))
  peaks_unl$id_pos <- paste0(peaks_unl$PeakID, "_", peaks_unl$pos)

  return(peaks_unl)
}
