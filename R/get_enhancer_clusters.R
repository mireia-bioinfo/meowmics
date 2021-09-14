#' Obtain enhancer clusters
#'
#' Given a list of peaks, obtain clusters.
#' @param gr Character with the name of the file containing peaks or GRanges object.
#' @param n_sites Minimum number of sites to consider a region enhancer cluster.
#' @param iterations Number of iterations for obtaining distance cutoff.
#' @param percentile Percentile of random distances to use for stitching peaks into the same cluster. Default: 0.25.
#' @param genome Character indicating the name of the genome to use. Default: hg38.
#' @param rm Names of chromosomes to remove. Default: chrX and chrY.
#' @export
#' @details To define enhancer clusters, input sites are first randomized \code{iterations}
#' times over the \code{genome} over idividual chromosomes). Then the \code{percentile}
#' percentile of inter-site distances of randomized is calculated for each chromosome. Next,
#' clusters are defined as any group of ≥ \code{n_sites} sites in which all adjacent sites
#' are separated by less than the abovementioned percentile distance.
#' @return GRanges object containing the coordinates for the enhancer clusters
#' and the following mcols:
#' \itemize{
#'     \item{clustID}{Unique identifier for the enhancer cluster.}
#'     \item{n_sites}{Number of sites in \code{gr} that are included in that specific cluster.}
#'     \item{distance_cutoff}{Distance used as cutoff for including regions in the same cluster.}
#' }
get_enhancer_clusters <- function(gr, n_sites=3, iterations=500, percentile=0.25, genome="hg38",
                                  rm=c("chrX", "chrY")) {

  if (is(gr, "character")) gr <- rtracklayer::import(gr)

  gr <- gr[order(gr)]
  gr <- dropSeqlevels(gr, rm, pruning.mode = "coarse")

  # 1) Randomize candidate enhancers + Calcuate 'percentile' of inter-site distances
  rndm_dist_chr <- BiocParallel::bplapply(as.character(unique(seqnames(gr))),
                                          function(chr) {
    lapply(1:iterations, function(i) interSiteDistances(regioneR::randomizeRegions(gr[seqnames(gr)==chr],
                                                                                   genome=genome,
                                                                                   per.chromosome=TRUE,
                                                                                   allow.overlaps=FALSE)))
  })

  # 2) Obtain percentile as cutoff - for each chromosome
  cutoff <- sapply(rndm_dist_chr, function(x) quantile(unlist(x), probs=percentile, na.rm=TRUE))

  # 3) Obtain distances betweeen enhancers
  distances <- split(interSiteDistances(gr), seqnames(gr))

  # 4) Create clusters while distances are shorter than cutoff
  clusters_list <- BiocParallel::bplapply(1:(length(distances)-1), function(i) {
    cluster <- cumsum(c(1, as.numeric(!(distances[[i]] < cutoff[[i]]))))
    gr_clust <- split(gr[seqnames(gr)==unique(seqnames(gr))[i]], cluster)
    num_sites <- sapply(gr_clust, length)
    gr_clust <- gr_clust[num_sites>=n_sites]

    ## Create output GRanges
    clusters <- lapply(gr_clust, function(x) GRanges(paste0(unique(seqnames(x)),
                                                            ":", min(start(x)),
                                                            "-", max(end(x))))
    )
    clusters <- unlist(GRangesList(clusters))
    clusters$n_sites <- num_sites[num_sites>=n_sites]
    clusters$distance_cutoff <- cutoff[[i]]
    clusters
  })

  clusters <- unlist(GRangesList(clusters_list))
  clusters$clustID <- paste0("cluster_", 1:length(clusters))

  return(clusters)
}

#' Obtain distances between consecutive regions
#'
#' @param gr Character with the name of the file containing peaks or GRanges object.
#' @export
#' @return A vector of the same length as gr containing the distances to the following
#' region. When there is no following region (i.e. end of a chromosome), returns NA.
interSiteDistances <- function(gr) {
  if (is(gr, "character")) gr <- rtracklayer::import(gr)

  gr <- gr[order(gr)]

  distances <- c(distance(gr[-length(gr)], gr[-1]), NA)
  return(distances)
}
