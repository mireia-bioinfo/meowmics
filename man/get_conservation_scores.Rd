% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/conservation.R
\name{get_conservation_scores}
\alias{get_conservation_scores}
\title{Get conservation scores}
\usage{
get_conservation_scores(
  peak_file,
  cons_score = phastCons7way.UCSC.hg38::phastCons7way.UCSC.hg38,
  bin_width = 20,
  window_width = 2000,
  merge_fun = "mean",
  random_control = TRUE,
  genome = "hg38",
  per.chromosome = FALSE
)
}
\arguments{
\item{peak_file}{Path for the peak file or GRanges object.}

\item{cons_score}{Annotation to use for getting the conservation scores. See \code{\link[GenomicScores]{gscores}}
for details.}

\item{bin_width}{Number of base pairs for binning the peak and obtaining the conservation scores. Default: 20}

\item{window_width}{Number of base pairs of a window centered in the center of the peak where to focus the
conservation analysis. Default = 2000}

\item{merge_fun}{Function for summarizing scores in the same region.}

\item{random_control}{Logical indicating whether to perform the same analysis using a randomized set of peaks
as a control. Default = TRUE}

\item{genome}{Character indicating the name of the genome where to randomize the regions.
See \code{\link[regioneR]{randomizeRegions}} for details.}

\item{per.chromosome}{Logical indicating if the randomization should be performed along the same chromosome as
the original set of peaks. See \code{\link[regioneR]{randomizeRegions}} for details.}
}
\value{
It returns a \code{data.frame} containing the summarized phastCons conservation scores for the file
and randomized control (if \code{random_control = TRUE}) at each relative position from the peak center.
The summarization is performed by computing the mean for all peaks in that specific position. Standard deviation
is also returned.
}
\description{
Given a peak set, this function obtains the mean conservation scores in the given window and bins.
}
