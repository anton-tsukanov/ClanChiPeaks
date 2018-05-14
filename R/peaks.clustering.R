#' Find clusters based on distance between peaks
#'
#'
#'
#' @title peaks.clustering
#' @param peaks is a GenomeRanges object with matadata colum that is named 'TF' (TF - Transcroption factor) and contain name for each peaks
#' @param cut.off is the maximal distance between two centres of peaks that belongs to the same cluster
#' @param min.size is a minimal size of cluster
#' @return GRangesList object where each element contains cluster of peaks
#' @import IRanges GenomicRanges GenomicFeatures fastcluster BiocParallel
#' @export
peaks.clustering <- function(peaks, cut.off = 20, min.size = 1) {
        peaks$BS <- round((start(peaks) + end(peaks))/2)
        chrs <- seqlevels(peaks)
        chrs <- sapply(chrs, cheking, peaks)
        chrs <- chrs[!sapply(chrs, is.null)]

        peaks <- as.data.frame(peaks)
        out <- bplapply(chrs, function(i){
                peaks.on.chr <- peaks[peaks$seqnames == i,]
                peaks.hc <- fastcluster::hclust.vector(peaks.on.chr$BS)
                cut.peaks <- cutree(peaks.hc, h = cut.off)
                bigger.than <- sapply(unique(cut.peaks), function(x) sum(cut.peaks == x) >= min.size)
                clusters <- unique(cut.peaks)[bigger.than]
                lapply(clusters, function(j){
                        cl <- peaks.on.chr[cut.peaks == j,]
                        dplyr::arrange(cl, BS)
                })
        })

        new.out <- list()

        for(i in out){
                for(j in i){
                        new.out[['']] <- j
                }
        }
        return(new.out)
}

#Cheking on contening of peaks on each chromosome
cheking <- function(chr, peaks){
        if(length(peaks[seqnames(peaks) == chr]) >= 2){
                return(chr)
        }
}
