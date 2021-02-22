#' Phantom Peak Quality Control
#'
#' Needs gawk unix package installed.
#' @param bam_file Character with the path of the BAM file.
#' @param suffix Suffix to remove from the bam_file to generate the sample name.
#' @param path_phantom Path to the phantom peak Rscript. Default assumes it is
#' in your PATH.
#' @param nthreads Number of threads to use for the analysis.
#' @return Data.frame containing the results of the QC.
#' @export
phantomPeakQC <- function(bam_file,
                          out_dir="phantomPeakQC/",
                          suffix=".bam",
                          path_phantom="run_spp.R",
                          nthreads=5) {
  name <- getNameFromPath(bam_file, suffix=suffix)

  ## Run Phantom Peak QC
  tmp <- tempdir()
  cmd <- paste(path_phantom,
               paste0("-c='", bam_file, "'"),
               paste0("-p=", nthreads),
               "-savp -rf",
               paste0("-odir='", out_dir, "'"),
               paste0("-out='", file.path(out_dir, paste0(name, ".txt'"))))
  system(cmd)

  ## Load files
  txt <- loadPhantomPeakQC(file_path=file.path(out_dir, paste0(name, ".txt")),
                             suffix = suffix)

  return(txt)
}

#' Load Phantom Peak Quality Control
#'
#' @inheritParams phantomPeakQC
#' @export
loadPhantomPeakQC <- function(file_path,
                                suffix = ".bam") {
  txt <- read.delim(file_path, stringsAsFactors = FALSE, header = FALSE)
  colnames(txt) <- c("sampleID", "numReads", "estFragLen", "corr_estFragLen", "phantomPeak",
                     "corr_phantomPeak", "argmin_corr", "min_corr", "NSC", "RSC", "QualityTag")

  txt$sampleID <- gsub(suffix, "", txt$sampleID)
  return(txt)
}
