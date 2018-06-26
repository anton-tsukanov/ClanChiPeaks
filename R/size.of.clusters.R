#' Calculation sizes of clusters (faster than clustering function, becouse return only vector of sizes)
#'
#'
#'
#' @title Calculation sizes of clusters
#' @param peaks is a GenomeRanges object
#' @param cut.off is the maximal distance between centres of peaks that belongs to the same cluster
#' @return vector of cluster sizes
#' @import GenomicRanges GenomicFeatures fastcluster BiocParallel
#' @export
size.of.clusters <- function(peaks, cut.off = 25){

        peaks$BS <- round((start(peaks) + end(peaks))/2)
        chrs <- seqlevels(peaks)
        r.out <- vector()
        peaks <- as.data.frame(peaks)

        t <- function(i, peaks.df, cut.off){
                peaks.on.chr <- peaks.df[peaks.df$seqnames == i,]
                peaks.hc <- fastcluster::hclust.vector(peaks.on.chr$BS)
                cut.peaks <- cutree(peaks.hc, h = cut.off)
                sapply(unique(cut.peaks), function(x) base::sum(cut.peaks == x))
        }

        r.out <- unlist(bplapply(chrs, t, peaks, cut.off))
        return(r.out)
}
