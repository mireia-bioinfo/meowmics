#' Phantom Peak Quality Control
#'
#' @param bam_file Character with the path of the BAM file.
#' @param suffix Suffix to remove from the bam_file to generate the sample name.
#' @param path_phantom Path to the phantom peak Rscript. Default assumes it is
#' in your PATH.
#' @param nthreads Number of threads to use for the analysis.
#' @return Data.frame containing the results of the QC.
#' @export
phantomPeakQual <- function(bam_file,
                            suffix,
                            path_phantom="run_spp.R",
                            nthreads=5) {
  name <- getNameFromPath(bam_file, suffix=suffix)

  ## Run Phantom Peak QC
  tmp <- tempdir()
  cmd <- paste(path_phantom,
               paste0("-c='", bam_file, "'"),
               paste0("-p=", nthreads),
               "-savp -rf",
               paste0("-odir='", tmp, "'"),
               paste0("-out='", name, ".txt'"))
  system(cmd)

  ## Load files
  txt <- read.delim(paste0(tmp, "/test.txt"), stringsAsFactors = FALSE, header = FALSE)
  colnames(txt) <- c("sampleID", "numReads", "estFragLen", "corr_estFragLen", "phantomPeak",
                     "corr_phantomPeak", "argmin_corr", "min_corr", "NSC", "RSC", "QualityTag")

  txt$sampleID <- gsub(suffix, "", txt$sampleID)

  return(txt)
}
