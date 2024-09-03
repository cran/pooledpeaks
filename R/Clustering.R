#' K-means Clustering
#'
#' @param RawData A data frame containing the raw data as read in by
#' [pooledpeaks::LoadData]
#' @param K An integer specifying the number of clusters.
#'
#' @importFrom stats kmeans
#' @importFrom utils capture.output
#'
#' @return A list containing the results of the K-means cluster analysis,
#' including cluster assignments and original data.
#' @export
#'
#' @examples
#' genetic_data <- data.frame(
#' Locus = c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2),
#' Locus_allele = c("Marker1", "n", 1, 2, 3, "Marker2", "n", 1, 2, 3),
#' Sample1 = c(NA, 10, 0.5, 0.5, 0, NA, 10, 0.2, 0.3, 0.5),
#' Sample2 = c(NA, 20, 0.1, 0.2, 0.7, NA, 20, 0.3, 0.4, 0.3),
#' Sample3 = c(NA, 30, 0.3, 0.4, 0.3, NA, 30, 0.4, 0.2, 0.4)
#' )
#' cluster(RawData=genetic_data, K=2)

cluster <- function(RawData=data.frame,K=2) {
  d <- RawData
  N <- TypedLoci(d)
  Ns <- 0
  for (i in 1:ncol(N)) Ns[i] <- N[i,i]
  w <- which(Ns==max(N))
  w <- w+2
  d <- d[,c(1:2,w)]
  dat <- d
  names <- colnames(d)
  names <- names[-(1:2)]
  PopLabel <- substr(names,1,1)
  PopFactor <- as.factor(PopLabel)
  d <- subset(d,(d[,2] != 'n')&(d[,3]!='NA'))
  d <- d[,-c(1:2)]

  fit <- stats::kmeans(t(d), centers = K, iter.max = 20, nstart = 50)
  fit_output <- utils::capture.output(print(fit))

  fit_output <- paste(fit_output, collapse = "\n")

  fit_cluster <- fit$cluster
  result <- list(PopFactor = PopFactor, clust = fit_cluster, dat = dat,
                 fit_output = fit_output)

  return(result)
}


#' Sample Of Loci
#'
#' An internal function that supports ClusterFromSamples. Sample loci from a
#' dataset based on the number of loci specified.
#'
#' @param aaax A data frame containing the input data must be in LoadData
#' style [pooledpeaks::LoadData].
#' @param NLoci An integer specifying the number of loci to sample.
#'
#' @return A data frame containing the sampled loci.
#'
SampleOfLoci <- function(aaax=data.frame,NLoci=max(aaax[,1])) {
  x <- sample(1:max(aaax[,1]),NLoci,replace = TRUE)
  z <- aaax[1,]
  for (i in x) {
    b <- subset(aaax,aaax[,1]==i)
    z <- rbind(z,b)
  }
  aaax <- z[-1,]
  return(aaax)
}



#' Cluster From Samples
#'
#' Perform clustering on samples of loci from a data frame and calculate
#' statistics.
#'
#' @param datafile A data frame containing the input data must be in LoadData
#' style [pooledpeaks::LoadData].
#' @param numloci An integer specifying the number of loci to sample.
#' @param reps An integer specifying the number of repetitions.
#'
#' @return A matrix containing statistics calculated from the clustering
#' results.
#' @export
#'
#' @examples
#' genetic_data <- data.frame(
#' Locus = c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2),
#' Locus_allele = c("Marker1", "n", 1, 2, 3, "Marker2", "n", 1, 2, 3),
#' Sample1 = c(NA, 10, 0.5, 0.5, 0, NA, 10, 0.2, 0.3, 0.5),
#' Sample2 = c(NA, 20, 0.1, 0.2, 0.7, NA, 20, 0.3, 0.4, 0.3),
#' Sample3 = c(NA, 30, 0.3, 0.4, 0.3, NA, 30, 0.4, 0.2, 0.4)
#' )
#'
#' ClusterFromSamples(datafile=genetic_data, numloci=5, reps=10)

ClusterFromSamples <- function(datafile=data.frame,numloci=5,reps=100) {
  for (j in 1:reps) {
    A <- SampleOfLoci(datafile,numloci)
    names <- colnames(datafile)
    names <- names[-(1:2)]
    PopLabel <- substr(names,1,1)
    PopFactor <- as.factor(PopLabel)
    H <- cluster(A,K=length(levels(PopFactor)))
    R <- table(H$clust,H$PopFactor)
    S <- colSums(R)
    m <-0
    for (i in 1:ncol(R)) {
      m[i] <- max(R[,i]/S[i])
    }

    if (j==1)
    {M <- m} else
    {M <- rbind(M,m)}
  }

  M <- round(M,4)
  return(M)
}



#' Multi Dimensional Scaling (MDS) Plot
#'
#' Generate a multidimensional scaling (MDS) plot from genetic distance data.
#'
#' @param distance A matrix containing the genetic distance data.
#' @param pcs A numeric vector specifying the principal coordinates to plot.
#' @param PF A factor vector specifying population labels.
#' @param y A character vector specifying colors for population labels.
#'
#' @importFrom stats cmdscale
#' @importFrom graphics par
#' @importFrom graphics title
#'
#' @return The output is the MDS plot for the samples for the specified principal
#' coordinates.
#'
#' @export
#'
#' @examples
#' genetic_distance_matrix <- matrix(c(
#' 0, 0.2836333, 0.2760485, 0.2685221, 0.2797302,0.3202661,
#' 0.2836333, 0, 0.2867215, 0.2687472, 0.2596309, 0.2957862,
#' 0.2760485,0.2867215, 0, 0.297918, 0.3057039, 0.3153261,
#' 0.2685221, 0.2687472, 0.297918,0, 0.2753477, 0.3042383,
#' 0.2797302, 0.2596309, 0.3057039, 0.2753477, 0,0.3398558,
#' 0.3202661, 0.2957862, 0.3153261, 0.3042383, 0.3398558, 0),
#'  nrow = 6, byrow = TRUE,dimnames = list(c("Sample1", "Sample2",
#'  "Sample3", "Ind1", "Ind2", "Ind3"),
#'  c("Sample1", "Sample2", "Sample3", "Ind1", "Ind2", "Ind3")))
#'
#'  MDSplot(distance=genetic_distance_matrix, pcs=c(1,3))
#'

MDSplot<- function(distance=matrix,pcs=c(1,2),PF=NULL,
                   y= c('dodgerblue','red','turquoise3','purple','olivedrab3')
                   ) {
  K <- nrow(distance)
  if (K < 11){
    K <- K - 1
  } else {
    K <- 10
  }
  names <- colnames(distance)
  PopLabel <- substr(names,1,1)

  if (is.null(PF)) PF <- as.factor(PopLabel)

  if (is.null(y)){ y <- c('dodgerblue','red','turquoise3','purple','olivedrab3')
  } else {y<- y}


  z <- rep(1,length(PF))

  for (i in 1:length(PF)) z[i] <- y[PF[i]]

  E <- stats::cmdscale(distance,eig=TRUE, k=K)# multidimensional scaling from R

  W <- 0
  for (i in 1:length(E$eig)) if (E$eig[i]>0) W <- W + E$eig[i]

  percent <- round(E$eig/W,3)*100
  LtextX <-  paste('Principal Coordinate (',pcs[1],') :  ',percent[pcs[1]],
                   '% Dispersion')
  LtextY <-paste('Principal Coordinate (',pcs[2],') :  ',percent[pcs[2]],
                 '% Dispersion')

  plot(E$points[,pcs[1]],E$points[,pcs[2]],pch=21,col='gray25',bg=z,cex=1,
       xlab=NA,ylab=NA,las=1)
  legend('topright', legend = c(substr(levels(PF), 1,6)), xjust=1, yjust=1,
         inset=c(-0.2, 0), fill = y, bty = "n", cex = 1)
  graphics::title(xlab=LtextX)
  graphics::title(ylab=LtextY)
  graphics::title(main='Principal Coordinates of Genetic Distances',
                  font.main=1, cex.main=1.25)
}
