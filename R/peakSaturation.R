peakfiles <- list.files("data/H3K27ac/Peaks/", pattern="NET.*.broadPeak", full.names=TRUE)
title <- "Saturation"
fc = NA
qval =NA
comp = "bp"
xlab="# of samples"
suffix="_peaks.broadPeak"

#' Peak saturation with increasing number of samples
#'
#' @param peakfiles
#' @param title
#' @param fc
#' @param qval
#' @param comp
#' @param xlab
#' @param suffix
#'
#' @import GenomicRanges
#' @import ggplot2
#' @export
peakSaturation <- function(peakfiles,
                           title,
                           fc=NA,
                           qval=NA,
                           comp="peaks", # or bp
                           xlab="# of samples",
                           suffix="_peaks.narrowPeak") {
  if (is.character(peakfiles)) {
    names <- gsub(suffix, "", basename(peakfiles))

    gr <- GRangesList(lapply(peakfiles, rtracklayer::import))
    names(gr) <- names

    if (!is.na(fc)) {
      gr <- lapply(gr,
                   function(x) x[x$signalValue>fc,])
    }

    if (!is.na(qval)) {
      gr <- lapply(gr,
                   function(x) x[x$qValue>=qval,])
    }

    gr <- GRangesList(gr)
  } else {
    gr <- GRangesList(peakfiles)
  }

  comb <- sapply(1:length(gr),
                 function(x) utils::combn(names(gr), m=x))

  comb.df.all <- data.frame()
  for (l in 1:length(comb)) {
    for (c in 1:ncol(comb[[l]])) {
      ## TODO Try parallelizing this section!!
      samples <- comb[[l]][,c]

      # Merge regions
      sel <- gr[names(gr) %in% samples]
      sel <- unlist(sel)
      sel <- regioneR::joinRegions(sel)

      # DF
      comb.df <- data.frame(num_samples=as.factor(l),
                            id_samples=paste0(samples, collapse=" "),
                            num_peaks=length(sel),
                            bp_covered=sum(width(sel)))
      comb.df.all <- rbind(comb.df.all, comb.df)
    }
  }


  stats <- comb.df.all %>%
    dplyr::group_by(num_samples) %>%
    dplyr::summarise(mean.peaks=mean(num_peaks),
              sd.peaks=sd(num_peaks),
              mean.bp=mean(bp_covered),
              sd.bp=sd(bp_covered))

  res <- list(df=comb.df.all,
              stats=stats)

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
                color="black", size=2) +
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
