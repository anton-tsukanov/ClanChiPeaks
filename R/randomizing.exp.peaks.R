#' Random redistribution of ChIP-seq peaks
#'
#'
#'
#' @title randomizing.exp.peaks
#' @param peaks is GRangesObject with experimental ChIP-seq peaks
#' @param ref is a TxDb object (annotation)
#' @return GrangesObject
#' @import GenomicRanges GenomicFeatures
#' @export
randomizing.exp.peaks <- function(peaks, ref){
        df <- data.frame()
        chrs <- seqlevels(ref)
        chrs.length <- seqlengths(ref)
        for(i in chrs){
                sub.peaks <- peaks[seqnames(peaks) == i]
                num.of.peaks <- length(sub.peaks)
                l <- as.numeric(chrs.length[i])
                r.start <- sample(1:l, size = num.of.peaks, replace = T)
                r.width <- sample(width(sub.peaks), size = num.of.peaks, replace = F)
                sub.df <- data.frame(chr = i, start = r.start, width = r.width,
                                     strand = rep('*', num.of.peaks))
                df <- rbind(df, sub.df)
        }
        df$TF <- 'TF'
        df$end <- df$start + df$width
        r.peaks <- makeGRangesFromDataFrame(df, keep.extra.columns = T)
        return(r.peaks)
}
