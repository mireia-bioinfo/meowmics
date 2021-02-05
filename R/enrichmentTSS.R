#' Get TSS enrichment
#'
#' @param bam_files Character vector of BAM files to obtain enrichment.
#' @param promoter_saf SAF promoter file generated with \code{binPromoterAnnotation}.
#' @param suffix Suffix to remove from bam files to use as sample names. Default: ".bam".
#' @param nthreads Number of cores to use for the analysis.
#' @return Data.frame with the different promoter bins and the mean reads per sample
#' @export
enrichmentTSS <- function(bam_files,
                          promoter_saf,
                          suffix=".bam",
                          nthreads = 5) {
  ## Obtain counts in binned promoters
  counts <- Rsubread::featureCounts(files,
                                    annot.ext = promoter_saf,
                                    allowMultiOverlap = TRUE,
                                    nthreads = nthreads)

  names <- getNameFromPath(files, suffix=suffix)
  colnames(counts$counts) <- names
  colnames(counts$stat) <- c("Status", names)

  ## Join counts and create long df
  anno <- cbind(saf, counts$counts)
  anno.m <- reshape2::melt(anno,
                           id.vars=1:7,
                           value.vars=8:ncol(anno),
                           variable.name="sample",
                           value.name="reads")

  ## Fix Position of negative strand promoters
  anno.m$Position[anno.m$Strand=="-"] <- - anno.m$Position[anno.m$Strand=="-"]

  ## Get mean promoter enrichment for each sample
  tss <- anno.m %>%
    group_by(sample, Position) %>%
    summarise(mean=mean(reads),
              sd=sd(reads),
              median=median(reads),
              mad=mad(reads))

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
#' @param out_format
#' @return Either a SAF of GRanges containing binned promoter regions.
#' @import GenomicRanges
#' @export
binPromoterAnnotation <- function(scope=2e3,
                                  bin=10,
                                  build="hg38",
                                  out_format="saf") {
  ## Select build
  if(tolower(build) %in% c("hg19", "grch37")) {
    host <- "grch37.ensembl.org"
  } else if (tolower(build) %in% c("hg38", "grch38")) {
    host <- "www.ensembl.org"
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
                          filters="biotype", values="protein_coding",
                          mart=mart)

  genes$strand[genes$strand==-1] <- "-"
  genes$strand[genes$strand==1] <- "+"

  genes <- regioneR::toGRanges(genes)
  strand(genes) <- genes$strand
  mcols(genes) <- mcols(genes)[,-1]

  genes <- GenomeInfoDb::keepStandardChromosomes(genes, pruning.mode = "coarse")
  genes <- GenomeInfoDb::dropSeqlevels(genes, "MT", pruning.mode = "coarse")

  if (!grepl("grch", build))   GenomeInfoDb::seqlevelsStyle(genes) <- "UCSC"

  genes$GeneID <- genes$ensembl_gene_id

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
