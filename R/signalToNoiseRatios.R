#' Signal-to-noise ratio
#'
#' Given a bam file and peak file, it calculated the signal-to-ratio, that is
#' the percentage of reads present in called peaks.
#' @param bam_file Character with the path of the BAM file.
#' @param peak_file Character with the path of the peak file.
#' @param suffix Suffix to remove from the bam_file to generate the sample name.
#' @param TxDb TxDb object to extract gene annotation from.
#' @param promoter_dist Distance from TSS to classify promoter and distal regions.
#' @param nthreads Number of threads to use for counting BAM reads.
#' @return Data.frame containing the number and percentage of reads located in
#' promoter peaks, distal peaks and unassigned.
#' @export
signalToNoiseRatio <- function(bam_file,
                               peak_file,
                               suffix=".bam",
                               TxDb,
                               promoter_dist=2e3,
                               nthreads=5) {
  name <- getNameFromPath(bam_file, suffix=suffix)
  ## Load and annotate peaks as promoter/distal
  peaks <- rtracklayer::import(peak_file)

  peaks$distanceToTSS <- mcols(distanceToNearest(peaks, promoters(TxDb),
                                                 ignore.strand = TRUE))$distance
  peaks$annotation <- "Promoter"
  peaks$annotation[abs(anno$distanceToTSS)>promoter_dist] <- "Distal"

  ## Convert to SAF
  peaks <- data.frame(peaks)
  saf <- peaks[,c(1:3,6,grep("annotation", colnames(peaks)))]
  colnames(saf) <- c("Chr", "Start", "End", "GeneID", "Annotation")
  saf$Strand <- "+"

  ## Obtain counts
  counts <- Rsubread::featureCounts(bam_file,
                                    annot.ext=saf,
                                    nthreads=nthreads)
  anno <- cbind(saf, counts$counts)
  colnames(anno)[ncol(anno)] <- "reads"

  sum <- anno %>%
    dplyr::group_by(Annotation) %>%
    dplyr::summarise(reads=sum(reads))
  sum[3,1] <- "Unassigned"
  sum[3,2] <- sum(counts$stat[grepl("Unas", counts$stat$Status),2])

  sum$percentage <- sum$reads / sum(sum$reads) * 100

  sum$sample <- name
  return(sum)
}
