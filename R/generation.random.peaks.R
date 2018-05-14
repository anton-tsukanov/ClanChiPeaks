#' Generation random peaks with certain width
#'
#'
#'
#' @title generation.random.peaks
#' @param number is amount of peaks te be generated
#' @param ref is a TxDb object (annotation)
#' @param width is a width of generated peaks
#' @return GrangesObject
#' @import GenomicRanges GenomicFeatures
#' @export
generation.random.peaks <- function(number, ref, width = 300){
        df <- data.frame()
        chrs <- seqlevels(ref)
        chrs.length <- seqlengths(ref)
        for(i in chrs){
                num.of.peaks <- round(number*(as.numeric(chrs.length[i])/sum(as.numeric(chrs.length))))
                l <- as.numeric(chrs.length[i])
                r.start <- sample(1:l, size = num.of.peaks, replace = T)
                r.width <- rep(300, num.of.peaks)
                r.chr <- rep(i, num.of.peaks)
                sub.df <- data.frame(chr = r.chr, start = r.start, width = r.width,
                                     strand = rep('*', num.of.peaks))
                df <- rbind(df, sub.df)
        }
        df$TF <- 'TF'
        df$end <- df$start + df$width
        r.peaks <- makeGRangesFromDataFrame(df, keep.extra.columns = T)
        return(r.peaks)
}
