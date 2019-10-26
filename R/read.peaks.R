#' Read .bed files in GRangesList
#'
#'
#'
#' @title Peaks reading
#' @param peaks is vector of filenames (shortnames, for reading files you should be in the directory where files are located) te be readed
#' @return GRangesList object where each element contains ChIP-seq peaks for certain transcription factors
#' @import GenomicRanges
#' @export
read.peaks <- function(files) {
        peaks <- lapply(files, read.csv, sep = '\t')
        names(peaks) <- gsub('.bed','', files)
        out <- list()
        for (i in names(peaks)){
                names(peaks[[i]]) <- c("chr", "start", "end")
                data <-  makeGRangesFromDataFrame(peaks[[i]], keep.extra.columns = T)
                data <- data[,-1]
                data$TF <- i
                data$BS <- round((start(data) + end(data)) / 2)
                out[[i]] <- data
        }
        out <- as(out, 'GRangesList')
        return(out)
}
