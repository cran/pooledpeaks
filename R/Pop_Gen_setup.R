#' Load Genetic Data
#'
#' This function imports data for genetic analysis.
#'
#' @param datafile The path to your datafile. The format of your data should be
#'  .txt or .csv.
#' @importFrom utils read.table
#'
#' @return A data frame containing the imported data formatted in the way
#' necessary for downstream population genetic functions.
#' @export
#'
#' @examples
#' file<-system.file("extdata", "Multiplex_frequencies.txt",
#' package = "pooledpeaks")
#' LoadData(file)
LoadData <- function(datafile=NULL)
{
  if (is.null(datafile))
    datafile <- readline(" Please enter the path to your data file :    ")
  N1453 <- utils::read.table(datafile,header=TRUE)

  Loci <- 0
  Locus <- 1

  for (i in 1:(nrow(N1453)-1)) {
    if (N1453[i+1,1]=='n') {
      Loci <- Loci+1
    }
    Locus[i] <- Loci
  }
  Locus[nrow(N1453)]<- Loci

  N1453 <- cbind(Locus,N1453)

  return(N1453)
}


#' Typed Loci
#'
#' This function calculates the number of loci successfully genotyped by each
#' individual included in our data set
#'
#' @param datafile A data frame containing the input data must be in LoadData
#'  style [pooledpeaks::LoadData].
#'
#' @return A matrix representing processed data.
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
#' TypedLoci(datafile=genetic_data)

TypedLoci <- function(datafile=data.frame) {

  b <- subset(datafile,datafile[,3]!='NA')

  for (i in 1:nrow(b))
    for (j in 3:ncol(b))
      if (b[i,j]>1) b[i,j]<-1

  b <- subset(b,b[,2]=='n')
  b <- b[,-(1:2)]
  X <- as.matrix(b)
  N <- t(X)%*%X

}

#' Gene Identity Matrix
#'
#' Using the number of typed loci, this function calculates the gene identity
#'  between all possible pairwise combinations between individuals for all
#'  markers creating a matrix.
#'
#' @param RawData A data frame containing the input data must be in LoadData
#' style [pooledpeaks::LoadData].
#' @param LociGenotyped The Output from the TypedLoci function
#'
#' @return The Gene Identity Matrix
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
#' n_alleles <- matrix(c(
#' 3, 3, 3,
#' 3, 3, 3,
#' 3, 3, 3
#' ), nrow = 3, byrow = TRUE,
#' dimnames = list(paste0("Sample", 1:3), paste0("Sample", 1:3)))
#'
#' GeneIdentityMatrix(RawData=genetic_data,LociGenotyped=n_alleles)

GeneIdentityMatrix <- function(RawData=data.frame,LociGenotyped=matrix) {
  d <- subset(RawData,(RawData[,2] != 'n')&(RawData[,3]!='NA'))
  d <- d[,-(1:2)]
  X <- as.matrix(d)
  J <- t(X)%*%X
  J <- J/LociGenotyped
  return(J)
}


#' Genetic Distance Matrix
#'
#' This function calculates the genetic distance matrix from a given gene
#' identity matrix.
#'
#' @param J The Gene Identity Matrix created using
#' [pooledpeaks::GeneIdentityMatrix]
#'
#' @return The Genetic Distance Matrix
#' @export
#'
#' @examples
#' gene_identity_matrix <- matrix(c(
#' 0.3164550, 0.2836333, 0.2760485,
#' 0.2836333, 0.3106084, 0.2867215,
#' 0.2760485, 0.2867215, 0.3338663
#' ), nrow = 3, byrow = TRUE,
#' dimnames = list(paste0("Sample", 1:3), paste0("Sample", 1:3)))
#'
#' GeneticDistanceMatrix(gene_identity_matrix)

GeneticDistanceMatrix <- function(J=matrix) {
  D <- rep(0,nrow(J)*ncol(J))
  D <- array(D,dim=c(nrow(J),ncol(J)))
  rownames(D) <- rownames(J)
  colnames(D) <- colnames(J)

  for (i in 1:nrow(J))
    for (j in 1:ncol(J))
      D[i,j] <- (J[i,i]+J[j,j]-2*J[i,j])/2
  return(D)
}

#' Random Walk Covariance Distance Matrix
#'
#' This function calculates the RWC (Random Walk Covariance) distance matrix
#' from a given matrix of genetic distances.
#'
#' @param J The Genetic Distance Matrix calculated using
#' [pooledpeaks::GeneticDistanceMatrix]
#'
#' @return A matrix representing the distance matrix calculated using the
#' Random Walk Covariance method.
#' @export
#'
#' @examples
#' genetic_distance_matrix <- matrix(c(0.316455, 0.2836333, 0.2760485,
#' 0.2685221, 0.2797302,0.3202661,0.2836333, 0.3106084, 0.2867215, 0.2687472,
#'  0.2596309, 0.2957862,0.2760485,0.2867215, 0.3338663, 0.297918, 0.3057039,
#'   0.3153261,0.2685221, 0.2687472, 0.297918,0.3107094, 0.2753477, 0.3042383,
#'   0.2797302, 0.2596309, 0.3057039, 0.2753477, 0.3761386,0.3398558,0.3202661,
#'    0.2957862, 0.3153261, 0.3042383, 0.3398558, 0.4402125),
#'  nrow = 6, byrow = TRUE, dimnames = list(c("Sample1", "Sample2", "Sample3",
#'  "Ind1", "Ind2", "Ind3"),
#'  c("Sample1", "Sample2", "Sample3", "Ind1", "Ind2", "Ind3")))
#'  RWCDistanceMatrix(genetic_distance_matrix)

RWCDistanceMatrix <- function(J=matrix){
  D <- rep(0,nrow(J)*ncol(J))
  D <- array(D,dim=c(nrow(J),ncol(J)))
  rownames(D) <- rownames(J)
  colnames(D) <- colnames(J)

  for (i in 1:nrow(J))
    for (j in 1:ncol(J)) {
      D[i,j] <- (J[i,i]+J[j,j]-2*J[i,j])/2
      D[i,j] <- D[i,j]/(1-J[i,j])
    }
  return(D)
}
