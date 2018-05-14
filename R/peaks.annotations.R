#' Counting of peaks in genome
#'
#'
#'
#' @title peaks.annotaions
#' @param peaks is a GenomeRanges object with matadata colum that is named 'TF' (TF - Transcroption factor) and contain name for each peaks
#' @param ref is TxDb object
#' @param is.Ratio if False returns a number of peaks in diff regions, if False returns proportions of peaks in diif regions
#' @return data.frame
#' @import IRanges GenomicRanges GenomicFeatures ggplot2
#' @export
peaks.annotaions <- function(peaks, ref, is.Ratio){

        ref.promoters <- promoters(genes(ref))
        ref.exons <- unlist(exonsBy(ref, by = 'tx'))
        ref.introns <- unlist(intronsByTranscript(ref))
        ref.5UTR <- unlist(fiveUTRsByTranscript(ref))
        ref.3UTR <- unlist(threeUTRsByTranscript(ref))

        tfs <- unique(peaks$TF)

        #Object contained summary about TFs
        tfs.summary <- data.frame(TF = tfs,
                                  inPromoters = rep(integer(1), length(tfs)),
                                  in5UTR = rep(integer(1), length(tfs)),
                                  inExons = rep(integer(1), length(tfs)),
                                  inIntrons = rep(integer(1), length(tfs)),
                                  in3UTR = rep(integer(1), length(tfs)),
                                  inIntergenic = rep(integer(1), length(tfs)))

        #Calculating summory about TFs
        for (i in tfs){
                tf <- peaks[peaks$TF == i]
                in.promoters <- length(subsetByOverlaps(tf, ref.promoters))
                in.introns <- length(subsetByOverlaps(tf, ref.introns))
                in.exons <- length(subsetByOverlaps(tf, ref.exons))
                in.5UTR <- length(subsetByOverlaps(tf, ref.5UTR))
                in.3UTR <- length(subsetByOverlaps(tf, ref.3UTR))

                in.intergenic <- length(tf) - (in.promoters + in.introns + in.exons + in.5UTR + in.3UTR)
                tfs.summary[tfs.summary$TF == i, 2:7] <- c(in.promoters, in.5UTR, in.exons, in.introns,
                                                           in.3UTR, in.intergenic)

        }

        if(is.Ratio){
                #Summary in % (Where TFs are located)
                total.tfs <- apply(tfs.summary[,2:7], 1, sum)
                out <- apply(tfs.summary[,2:7], 2, function(x) x/total.tfs)
                out <- round(out, digits = 2)
                out <- data.frame(out)
                out$TF <- tfs
                out <- cbind(out[,7], out[,1:6])
                names(out)[1] <- 'TF'
        } else {
                out <- tfs.summary
        }
        return(out)
}
