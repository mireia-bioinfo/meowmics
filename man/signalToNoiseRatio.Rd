% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/signalToNoiseRatios.R
\name{signalToNoiseRatio}
\alias{signalToNoiseRatio}
\title{Signal-to-noise ratio}
\usage{
signalToNoiseRatio(
  bam_file,
  peak_file,
  genes = NULL,
  suffix = ".bam",
  build = "hg38",
  promoter_dist = 2000,
  nthreads = 5,
  paired_end = FALSE
)
}
\arguments{
\item{bam_file}{Character with the path of the BAM file.}

\item{peak_file}{Character with the path of the peak file.}

\item{genes}{Data.frame containing the coding genes to use as annotation,
as outputed by the function \code{\link{obtainCodingGenes}}. If set to NULL
(default) it will run the function to recover coding genes.}

\item{suffix}{Suffix to remove from the bam_file to generate the sample name.}

\item{build}{Name of the build to use to extract promoter coordinates. Can be
either hg38 (default) or hg19.}

\item{promoter_dist}{Distance from TSS to classify promoter and distal regions.}

\item{nthreads}{Number of threads to use for counting BAM reads.}
}
\value{
Data.frame containing the number and percentage of reads located in
promoter peaks, distal peaks and unassigned.
}
\description{
Given a bam file and peak file, it calculates the signal-to-ratio, that is,
the percentage of reads present in called peaks.
}
