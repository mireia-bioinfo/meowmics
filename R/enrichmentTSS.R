#' Get TSS enrichment
#'
#' @param bam_files Character vector of BAM files to obtain enrichment.
#' @param promoter_saf SAF promoter file generated with \code{binPromoterAnnotation}.
#' @param suffix Suffix to remove from bam files to use as sample names. Default: ".bam".
#' @param nthreads Number of cores to use for the analysis.
#' @param paired_end Logical indicating whether the libraries are paired end or not. (Default: FALSE).
#' @return Data.frame with the different promoter bins and the mean reads per sample
#' @export
enrichmentTSS <- function(bam_files,
                          promoter_saf,
                          suffix=".bam",
                          nthreads = 5,
                          paired_end = FALSE) {

  ## Obtain counts in binned promoters
  counts <- Rsubread::featureCounts(bam_files,
                                    annot.ext = promoter_saf,
                                    allowMultiOverlap = TRUE,
                                    nthreads = nthreads,
                                    isPairedEnd = paired_end)

  names <- getNameFromPath(bam_files, suffix=suffix)
  colnames(counts$counts) <- names
  colnames(counts$stat) <- c("Status", names)
  total_reads <- colSums(counts$stat[,-1])
  assigned_reads <- as.numeric(counts$stat[counts$stat$Status=="Assigned",-1])
  names(assigned_reads) <- names(total_reads)

  ## Join counts and create long df
  anno <- cbind(promoter_saf, counts$counts)
  anno.m <- reshape2::melt(anno,
                           id.vars=1:7,
                           value.vars=8:ncol(anno),
                           variable.name="sample",
                           value.name="reads")

  ## Fix Position of negative strand promoters
  anno.m$Position[anno.m$Strand=="-"] <- - anno.m$Position[anno.m$Strand=="-"]

  ## Get mean promoter enrichment for each sample
  tss <- anno.m %>%
    dplyr::left_join(tibble::rownames_to_column(data.frame(total_reads)),
                     by=c(sample="rowname")) %>%
    dplyr::left_join(tibble::rownames_to_column(data.frame(assigned_reads)),
                     by=c(sample="rowname")) %>%
    dplyr::group_by(sample, Position, total_reads, assigned_reads) %>%
    dplyr::summarise(mean=mean(reads*1e7/assigned_reads),
              sd=sd(reads*1e7/assigned_reads),
              median=median(reads*1e7/assigned_reads),
              mad=mad(reads*1e7/assigned_reads))

  return(tss)
}

#' Build binned promoter annotation for TSS enrichment
#'
#' @param scope Integer indicating how many base pairs upstream and downstream
#' to extend the promoter region. Returned window width will be scope*2.
#' Default: 2Kb.
#' @param bin Integer indicating the size of the bins. Default: 10 bp.
#' @param build Name of the build to use to extract promoter coordinates. Can be
#' either hg38 (default) or hg19.
#' @param out_format Output format, either "saf" or "granges".
#' @param genes Genes, as outputed by `obtainCodingGenes`, that is: a `GRanges`
#' object with one `mcol` named `GeneID` which is a unique identifyer for the
#' region.
#' @return Either a SAF of GRanges containing binned promoter regions.
#' @import GenomicRanges
#' @export
binPromoterAnnotation <- function(scope=2e3,
                                  bin=10,
                                  build="hg38",
                                  out_format="saf",
                                  genes=NULL) {
  ## Select build
  if(tolower(build) %in% c("hg19", "grch37")) {
    host <- "grch37.ensembl.org"
  } else if (tolower(build) %in% c("hg38", "grch38")) {
    host <- "www.ensembl.org"
  } else {
    stop("Build not recognized. Choose from hg19 or hg38.")
  }

  ## Obtain gene annotation
  if (is.null(genes)) genes <- obtainCodingGenes(build=build)
  if (!grepl("GeneID", colnames(mcols(genes)))) genes$GeneID <- genes$ensembl_gene_id

  # Extend regions to scope*2
  regions <- unique(promoters(genes, upstream=0, downstream=1))
  regions.bin <- binRegions(regions,
                            scope=scope,
                            bin=bin)

  # Output saf or granges
  if(out_format=="saf") {
    saf <- data.frame(regions.bin)[,c(1:3,5,6,8)]
    colnames(saf) <- c("Chr", "Start", "End", "Strand", "txID", "Position")
    saf$GeneID <- paste0(as.character(saf$txID), "_", saf$Position)
    return(saf)
  } else if (tolower(out_format) == "granges") {
    return(regions.bin)
  }
}

#' Bin a set of GRanges regions
#'
#' @param regions GRanges object containing the regions to bin. It assumes that
#' the first mcol contained the unique identifier of the region.
#' @inheritParams binPromoterAnnotation
#' @return Binned GRanges object.
#' @export
binRegions <- function(regions,
                       scope=2e3,
                       bin=10) {
  regions.ext <- resize(regions, width=scope*2, fix="center")
  regions.ext$center <- start(regions.ext) + scope

  # Bin regions
  regions.bin <- tile(regions.ext, width=bin)
  n <- unique(sapply(regions.bin, length))
  regions.bin <- unlist(regions.bin)
  regions.bin$GeneID <- rep(mcols(regions.ext)[,1], each=n) ## add identifier
  regions.bin$center <- rep(regions.ext$center, each=n) ## add center
  regions.bin$pos <- start(regions.bin) - regions.bin$center ## add pos relative to center
  return(regions.bin)
}
