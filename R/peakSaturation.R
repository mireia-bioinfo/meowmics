#' Peak saturation with increasing number of samples
#'
#' @param peakfiles Character vector of peak file paths or list of GRanges object
#' containing the peaks for each one of the different samples.
#' @param fc Numeric indicating a fold-change value to filter out peaks (from column
#' named "signalValue").
#' @param qval Numeric indicating a q-value threshold to filter out peaks (from column
#' named "qValue").
#' @param suffix Suffix to remove from the file paths to create sample names.
#'
#' @import GenomicRanges
#' @export
peakSaturation <- function(peakfiles,
                           fc=NA,
                           qval=NA,
                           suffix="_peaks.narrowPeak") {
  if (is.character(peakfiles)) {
    names <- gsub(suffix, "", basename(peakfiles))

    gr <- GRangesList(lapply(peakfiles, rtracklayer::import))
    names(gr) <- names

  } else {
    gr <- GRangesList(peakfiles)
  }

  if (!is.na(fc)) gr <- GRangesList(lapply(gr, function(x) x[x$signalValue>fc,]))

  if (!is.na(qval)) gr <- GRangesList(lapply(gr, function(x) x[x$qValue>=qval,]))

  ## Obtain all possible combinations
  comb <- sapply(1:length(gr),
                 function(x) utils::combn(names(gr), m=x))

  ## Obtain overlapping peaks & bp covered
  comb.df.all <- data.frame()
  for (i in 1:length(comb)) {
    message("Calculating combinations of ", i, " samples")

    tictoc::tic()
    rows <- BiocParallel::bplapply(1:ncol(comb[[i]]),
                                   function(j) .obtainOverlapMeasures(gr = gr,
                                                                      samples = comb[[i]][,j],
                                                                      num_samples = i)
                                   )

    rows_df <- do.call(rbind, rows)
    comb.df.all <- rbind(comb.df.all, rows_df)
    tictoc::toc()
  }

  ## Summarize stats
  stats <- comb.df.all %>%
    dplyr::group_by(num_samples) %>%
    dplyr::summarise(mean.peaks=mean(num_peaks),
              sd.peaks=sd(num_peaks),
              mean.bp=mean(bp_covered),
              sd.bp=sd(bp_covered))

  res <- list(df=comb.df.all,
              stats=stats)

  return(res)
}

#' Obtain number of overlaps and width of the overlapping regions
#'
#' @param gr GRangesList with the peaks for each of the different samples.
#' @param samples Name of the samples (should be names in the GRangesList object)
#' to collapse together.
.obtainOverlapMeasures <- function(gr, samples, num_samples) {
  # Merge regions
  sel <- gr[names(gr) %in% samples]
  sel <- unlist(sel)
  sel <- regioneR::joinRegions(sel)

  # DF
  comb.df <- data.frame(num_samples=as.factor(num_samples),
                        id_samples=paste0(samples, collapse=" "),
                        num_peaks=length(sel),
                        bp_covered=sum(width(sel)))

  return(comb.df)
}

#' Plot peak saturation
#'
#' @param res Results of the \code{peakSaturation} function.
#' @param comp Name of the comparison to plot. It can be either "peaks" (default),
#' which will plot the number of overlapping peaks or "bp" which will overlap the
#' number of overlapping bp.
#' @param title Character with the title of the plot
#' @param xlab Label for the x axis.
#' @import ggplot2
#' @export
plotPeakSaturation <- function(res,
                               comp="peaks", #or bp
                               title,
                               xlab="# of samples") {
  if (comp=="peaks") {
    yaxis <- "mean.peaks"
    ymin <- res$stats$mean.peaks - res$stats$sd.peaks
    ymax <- res$stats$mean.peaks + res$stats$sd.peaks
    dots <- "num_peaks"
    title.yaxis <- "# of Peaks"
  } else if (comp=="bp") {
    yaxis <- "mean.bp"
    ymin <- res$stats$mean.bp - res$stats$sd.bp
    ymax <- res$stats$mean.bp + res$stats$sd.bp
    dots <- "bp_covered"
    title.yaxis <- "Bp covered"
  }

  num_samples <- "num_samples"

  ggplot() +
    geom_line(data=res$stats,
              aes_string(num_samples, yaxis, group=1), size=0.7, linetype=2) +
    geom_jitter(data=res$df,
                aes_string(num_samples, dots,
                           fill=num_samples), width=0.1, pch=21,
                color="black", size=2, alpha=0.5) +
    geom_pointrange(data=res$stats,
                    aes_string(x=num_samples, y=yaxis, ymin=ymin, ymax=ymax),
                    size=1) +
    xlab(xlab) +
    scale_y_continuous(label=function(x) scales::comma(x),
                       breaks=scales::pretty_breaks(),
                       name=title.yaxis) +
    theme(legend.position="none") +
    ggtitle(title)
}
