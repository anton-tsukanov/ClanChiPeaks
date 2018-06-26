#' Plot the density profile of peaks in certain sites
#'
#' Plot the density profile of peaks in certain sites
#'
#' @title peaks.near.TSS
#' @param peaks is a Granges object that contain TF (names of transcription factors) and Condition (colum with names of cell,tissue and so on) columns in metadata
#' @param sites is a number of sites (GRanges object) with the equal width, eg. promoter regions
#' @param range is a range with min and max values on x axis (abs sum is equale to width of sites)
#' @param ignor.strand is a boolean and used when strand of sites is '*'
#' @param vertical.facet TRUE for making vertical facet (useful when you have a lot of different TFs or Conditions)
#' @param horizontal.facet TRUE for making vertical facet (useful when you have TF and Condition colums)
#' @param wrap.facet TRUE for making wrap facet (useful if there are a lot of TFs)
#' @param x.label is a name of x axis
#' @param y.label is a name of y axis
#' @param legend.title is a title of legend
#' @param wide.of.line is width of curve on plot
#' @param axis.text.size is size of text on axis
#' @param axis.title.size is size of title text
#' @param legend.title.size is size of title text on legend
#' @param legend.text.size is size of text on legend
#' @param vertical.line is a coordinate of vertical line on plot
#' @return ggplot object
#' @import IRanges GenomicRanges ggplot2
#' @export
peaks.near.TSS <- function(peaks, sites, range = c(-2000, 2000), ignor.strand = F, vertical.facet = F,
                               horizonatl.facet = F, wrap.facet = T, x.label = 'Coordinates', y.label = 'Density',
                               legend.title = 'TF', wide.of.line = 0.7, axis.text.size = 8, axis.title.size = 12,
                               legend.title.size = 12, legend.text.size = 10, vertical.line = 0){

        if(ignor.strand){
                if(!is.null(peaks$Condition)){
                        TFs <- unique(peaks$TF)
                        Conditions <- unique(peaks$Condition)
                        df.with.cov <- data.frame(Coordinates = vector(),
                                                  Density = vector(),
                                                  TF = character(),
                                                  Condition = character())


                        for(i in Conditions){
                                sub.con.peaks <- peaks[peaks$Condition == i]
                                for(j in TFs){
                                        sub.tfs.peaks <- sub.con.peaks[sub.con.peaks$TF == j]
                                        cov <- coverage(sub.tfs.peaks)
                                        len <- abs(range[1]) + abs(range[2])

                                        norm.cov.all <- calc.cov.on.sites(cov = cov, sites = sites, len = len)
                                        norm.cov.all <- norm.cov.all/length(sub.tfs.peaks)/sum(width(sites))*1e6

                                        out <- data.frame(Coordinates = c(range[1]:-1, 1:range[2]),
                                                          Density = norm.cov.all,
                                                          TF = rep(j, len),
                                                          Condition = rep(i, len))
                                        df.with.cov <- rbind(df.with.cov, out)

                                }

                        }

                } else {
                        TFs <- unique(peaks$TF)
                        df.with.cov <- data.frame(Coordinates = vector(),
                                                  Density = vector(),
                                                  TF = character())

                        for(j in TFs){
                                sub.tfs.peaks <- peaks[peaks$TF == j]
                                cov <- coverage(sub.tfs.peaks)
                                len <- abs(range[1]) + abs(range[2])

                                norm.cov.all <- calc.cov.on.sites(cov = cov, sites = sites, len = len)
                                norm.cov.all <- norm.cov.all/length(sub.tfs.peaks)/sum(width(sites))*1e6

                                out <- data.frame(Coordinates = c(range[1]:-1, 1:range[2]),
                                                  Density = norm.cov.all,
                                                  TF = rep(j, len))
                                df.with.cov <- rbind(df.with.cov, out)
                        }
                }
        } else {
                if(!is.null(peaks$Condition)){
                        TFs <- unique(peaks$TF)
                        Conditions <- unique(peaks$Condition)
                        df.with.cov <- data.frame(Coordinates = vector(),
                                                  Density = vector(),
                                                  TF = character(),
                                                  Condition = character())


                        for(i in Conditions){
                                sub.con.peaks <- peaks[peaks$Condition == i]
                                for(j in TFs){
                                        sub.tfs.peaks <- sub.con.peaks[sub.con.peaks$TF == j]
                                        cov <- coverage(sub.tfs.peaks)
                                        len <- abs(range[1]) + abs(range[2])

                                        sites.1 <- sites[strand(sites) == '+']
                                        norm.cov.1 <- calc.cov.on.sites(cov = cov, sites = sites.1, len = len)

                                        sites.0 <- sites[strand(sites) == '-']
                                        norm.cov.0 <- calc.cov.on.sites(cov = cov, sites = sites.0, len = len)

                                        norm.cov.all <- norm.cov.1 + norm.cov.0[length(norm.cov.0):1]
                                        norm.cov.all <- norm.cov.all/length(sub.tfs.peaks)/sum(width(sites))*1e6

                                        out <- data.frame(Coordinates = c(range[1]:-1, 1:range[2]),
                                                          Density = norm.cov.all,
                                                          TF = rep(j, len),
                                                          Condition = rep(i, len))
                                        df.with.cov <- rbind(df.with.cov, out)

                                }

                        }

                } else {
                        TFs <- unique(peaks$TF)
                        df.with.cov <- data.frame(Coordinates = vector(),
                                                  Density = vector(),
                                                  TF = character())

                        for(j in TFs){
                                sub.tfs.peaks <- peaks[peaks$TF == j]
                                cov <- coverage(sub.tfs.peaks)
                                len <- abs(range[1]) + abs(range[2])

                                sites.1 <- sites[strand(sites) == '+']
                                norm.cov.1 <- calc.cov.on.sites(cov = cov, sites = sites.1, len = len)

                                sites.0 <- sites[strand(sites) == '-']
                                norm.cov.0 <- calc.cov.on.sites(cov = cov, sites = sites.0, len = len)

                                norm.cov.all <- norm.cov.1 + norm.cov.0[length(norm.cov.0):1]
                                norm.cov.all <- norm.cov.all/length(sub.tfs.peaks)/sum(width(sites))*1e6

                                out <- data.frame(Coordinates = c(range[1]:-1, 1:range[2]),
                                                  Density = norm.cov.all,
                                                  TF = rep(j, len))
                                df.with.cov <- rbind(df.with.cov, out)
                        }
                }
        }

        pic <- ggplot(df.with.cov, aes(x = Coordinates, y = Density))+
                geom_line(size = wide.of.line)+
                geom_vline(xintercept = vertical.line)+
                scale_fill_grey()+
                xlab(label = x.label)+
                ylab(label = y.label)+
                theme_linedraw()+
                theme(axis.text=element_text(size=axis.text.size),
                      axis.title=element_text(size=axis.title.size,face="bold"),
                      legend.title = element_text(size = legend.title.size),
                      legend.text = element_text(size = legend.text.size),
                      strip.text = element_text(size = 14))

        if(vertical.facet & length(unique(df.with.cov$TF)) == 1 & !is.null(peaks$Condition)){pic <- pic + facet_grid(Condition ~ .) + labs(colour = legend.title)}
        if(vertical.facet & length(unique(df.with.cov$TF)) == 1){pic <-  pic + facet_grid(Condition ~ .) + labs(colour = legend.title)}
        if(horizonatl.facet & !is.null(peaks$Condition)) {pic <- pic + facet_grid(. ~ Condition) + labs(colour = legend.title)}
        if(vertical.facet & !wrap.facet){pic <- pic + facet_grid(TF ~ .) + labs(colour = legend.title)}
        if(wrap.facet){pic <- pic + facet_wrap(~TF)}


        return(pic)
}

#' @import IRanges GenomicRanges
calc.cov.on.sites <- function(cov, sites, len){

        vi <- Views(cov, as(sites, 'IntegerRangesList'))
        cov.of.sites <- lapply(vi, function(x) viewApply(x,
                                                         function(el){
                                                                 el <- as.vector(el)
                                                                 names(el) <-  NULL
                                                                 if(length(el) != len){el <- c(el, rep(0, len - length(el)))}
                                                                 return(el)
                                                         }))

        norm.cov <- rep(0, len)
        for(i in 1:length(cov.of.sites)){
                sub.cov <- apply(cov.of.sites[[i]], 1, sum)
                norm.cov <- norm.cov + sub.cov
        }
        return(norm.cov)
}
