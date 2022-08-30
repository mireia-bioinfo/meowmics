#' Whole-genome correlation between samples
#'
#' @param bam_files Character vector of bam file paths.
#' @param chr_sizes Text files with the sizes for each chormosome, as obtain from
#' @param suffix Suffix to remove from bam files to use as sample names. Default: ".bam".
#' @param nthreads Number of cores to use for the analysis.
#' @param paired_end Logical indicating whether the libraries are paired end or not. (Default: FALSE).
#' @export
wholeGenomeCor <- function(bam_files,
                           chr_sizes,
                           suffix=".bam",
                           bin_size=10e3,
                           nthreads = 3,
                           paired_end = FALSE) {

  ## Obtain binned genome
  genome_bins <- binGenome(chr_sizes, bin_size, out_type="saf")

  if(is(paired_end, "character")) paired_end <- paired_end == "PE"

  ## Count reads for each window
  counts <- Rsubread::featureCounts(bam_files,
                                    annot.ext = genome_bins,
                                    allowMultiOverlap = TRUE,
                                    nthreads = nthreads,
                                    isPairedEnd = paired_end)


  colnames(counts$counts) <- gsub(suffix, "", colnames(counts$counts))

  ## Obtain correlation matrix
  cormat <- cor(counts$counts)

  ## Return correlation matrix
  return(cormat)
}

#' Bin genome into windows
#'
#' @inheritParams wholeGenomeCor
#' @import GenomicRanges
binGenome <- function(chr_sizes,
                      bin_size=10e3,
                      out_type="saf") {
  chr_sizes <- read.delim(chr_sizes, header=FALSE, stringsAsFactors=FALSE)
  genome <- GRanges(seqnames=chr_sizes[,1],
                    ranges=IRanges::IRanges(start=1, end=chr_sizes[,2]))

  genome <- keepStandardChromosomes(genome, pruning.mode = "coarse")
  genome <- genome[seqnames(genome)!="chrM"]

  genome_bins <- unlist(tile(genome, width=bin_size))
  genome_bins$id <- paste0("genomeBin_", 1:length(genome_bins))

  if (out_type=="saf") genome_bins <- grangesToSAF(genome_bins)

  return(genome_bins)
}
