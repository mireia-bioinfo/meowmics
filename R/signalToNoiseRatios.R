#' Signal-to-noise ratio
#'
#' Given a bam file and peak file, it calculates the signal-to-ratio, that is,
#' the percentage of reads present in called peaks.
#' @param bam_file Character with the path of the BAM file.
#' @param peak_file Character with the path of the peak file.
#' @param genes Data.frame containing the coding genes to use as annotation,
#' as outputed by the function \code{\link{obtainCodingGenes}}. If set to NULL
#' (default) it will run the function to recover coding genes.
#' @param suffix Suffix to remove from the bam_file to generate the sample name.
#' @param build Name of the build to use to extract promoter coordinates. Can be
#' either hg38 (default) or hg19.
#' @param promoter_dist Distance from TSS to classify promoter and distal regions.
#' @param nthreads Number of threads to use for counting BAM reads.
#' @return Data.frame containing the number and percentage of reads located in
#' promoter peaks, distal peaks and unassigned.
#' @export
signalToNoiseRatio <- function(bam_file,
                               peak_file,
                               genes = NULL,
                               suffix=".bam",
                               build="hg38",
                               promoter_dist=2e3,
                               nthreads=5,
                               paired_end = FALSE) {
  name <- getNameFromPath(bam_file, suffix=suffix)
  ## Load and annotate peaks as promoter/distal
  if (is(peak_file, "character")) {
    peaks <- rtracklayer::import(peak_file)
  } else if (is(peak_file, "GRanges")) {
    peaks <- peak_file
  }

  if (length(unique(mcols(peaks)[,1]))!=length(peaks)) stop("First mcol should be a unique identifier.")

  if (is.null(genes)) genes <- obtainCodingGenes(build=build)

  ## TODO: Add distances in another variable
  dist <- distanceToNearest(peaks, promoters(genes, upstream=0, downstream=1),
                            ignore.strand = TRUE)
  peaks$distanceToTSS <- NA
  peaks$distanceToTSS[queryHits(dist)] <- mcols(dist)$distance
  peaks$annotation <- "Distal"
  peaks$annotation[abs(peaks$distanceToTSS)<=promoter_dist] <- "Promoter"

  ## Convert to SAF
  peaks <- data.frame(peaks)
  saf <- peaks[,c(1:3,6,grep("annotation", colnames(peaks)))]
  colnames(saf) <- c("Chr", "Start", "End", "GeneID", "Annotation")
  saf$Strand <- "+"

  ## Obtain counts
  counts <- Rsubread::featureCounts(bam_file,
                                    annot.ext=saf,
                                    nthreads=nthreads,
                                    isPairedEnd = paired_end)
  anno <- cbind(saf, counts$counts)
  colnames(anno)[ncol(anno)] <- "reads"

  sum <- anno %>%
    dplyr::group_by(Annotation) %>%
    dplyr::summarise(reads=sum(reads))
  sum[nrow(sum)+1,1] <- "Unassigned"
  sum[nrow(sum),2] <- sum(counts$stat[grepl("Unas", counts$stat$Status),2])

  sum$percentage <- sum$reads / sum(sum$reads) * 100

  sum$sample <- name
  return(sum)
}

#' Obtain protein coding genes
#'
#' Obtains a list of protein coding genes from the human builds (hg38 or hg19).
#' @param build Name of the build to use to extract promoter coordinates. Can be
#' either hg38 (default) or hg19.
#' @return Data.frame with the coordinates, strand, ensembl_gene_id and
#' external_gene_name of protein coding genes.#'
#' @export
obtainCodingGenes <- function(build="hg38") {
  ## Select build
  if(tolower(build) %in% c("hg19", "grch37")) {
    host <- "https://grch37.ensembl.org"
  } else if (tolower(build) %in% c("hg38", "grch38")) {
    host <- "https://www.ensembl.org"
  } else {
    stop("Build not recognized. Choose from hg19 or hg38.")
  }

  ## Obtain gene annotation
  mart <- biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL",
                           host=host,
                           path="/biomart/martservice",
                           dataset="hsapiens_gene_ensembl")
  genes <- biomaRt::getBM(attributes=c("chromosome_name", "start_position", "end_position",
                                       "strand", "ensembl_gene_id", "external_gene_name"),
                          filters="biotype", values="protein_coding", useCache = FALSE,
                          mart=mart)

  genes$strand[genes$strand==-1] <- "-"
  genes$strand[genes$strand==1] <- "+"

  genes <- regioneR::toGRanges(genes)
  strand(genes) <- genes$strand
  mcols(genes) <- mcols(genes)[,-1]

  genes <- GenomeInfoDb::keepStandardChromosomes(genes, pruning.mode = "coarse")
  genes <- GenomeInfoDb::dropSeqlevels(genes, "MT", pruning.mode = "coarse")

  if (!grepl("grch", build))   GenomeInfoDb::seqlevelsStyle(genes) <- "UCSC"

  return(genes)
}
