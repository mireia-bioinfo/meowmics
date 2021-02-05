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
    group_by(Annotation) %>%
    summarise(reads=sum(reads))
  sum[3,1] <- "Unassigned"
  sum[3,2] <- sum(counts$stat[grepl("Unas", counts$stat$Status),2])

  sum$sample <- name
  df <- rbind(df, sum)
  return(df)
}
