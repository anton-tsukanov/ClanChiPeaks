#' Calculation of matrix for co-localized ChIP-seq peaks in clusters
#'
#'
#'
#' @title cal.matrix
#' @param clusters is a list object with the clusters (maded by peaks.clustering)
#' @param tfs is vector of transcription factors that are contained in clusters
#' @return matrix
#' @export
cal.matrix <- function(clusters, tfs){
        m.all <- matrix(data = 0, nrow = length(tfs), ncol = length(tfs))
        tfs.list <- lapply(clusters, function(x) x$TF)
        for(x in tfs.list){
                y <- tfs %in% x
                m1 <- matrix(y, nrow = length(y))
                m2 <- matrix(y, ncol = length(y))
                m3 <- m1 %*% m2
                m.all <- m.all + m3

        }
        colnames(m.all) <- tfs
        rownames(m.all) <- tfs
        return(m.all)

}
