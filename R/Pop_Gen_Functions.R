#' Distance Correlation
#'
#' Calculate the correlation between expected and realized genetic distances
#' and plot them.
#'
#' @param GD A matrix containing the genetic distance data.
#'
#' @importFrom stats cophenetic
#' @importFrom stats cor
#' @importFrom ape nj
#' @importFrom graphics par
#' @importFrom graphics abline
#' @importFrom graphics title
#'
#' @return A plot showing the Expected Genetic Distance versus Realized Genetic
#'  Distance
#' @export
#'
#' @examples
#' genetic_distance_matrix <- matrix(c(0.316455, 0.2836333, 0.2760485,
#' 0.2685221, 0.2797302,0.3202661,0.2836333, 0.3106084, 0.2867215, 0.2687472,
#' 0.2596309, 0.2957862,0.2760485,0.2867215, 0.3338663, 0.297918, 0.3057039,
#' 0.3153261,0.2685221, 0.2687472, 0.297918,0.3107094, 0.2753477, 0.3042383,
#' 0.2797302, 0.2596309, 0.3057039, 0.2753477, 0.3761386, 0.3398558,0.3202661,
#' 0.2957862, 0.3153261, 0.3042383, 0.3398558, 0.4402125),
#'  nrow = 6, byrow = TRUE,dimnames = list(c("Sample1", "Sample2", "Sample3",
#'   "Ind1", "Ind2", "Ind3"),
#'  c("Sample1", "Sample2", "Sample3", "Ind1", "Ind2", "Ind3")))
#'
#' DC<- DistCor(genetic_distance_matrix)
#'
DistCor <- function(GD=matrix) {
  R <- ape::nj(GD)
  ED <- stats::cophenetic(R)
  plot(GD,ED)
  graphics::abline(0,1,col=2)

  k<-0
  x<-0
  y<-0
  for (i in 1:(nrow(GD)-1))
    for (j in 2:ncol(GD))
    {k<-k+1
    x[k]<- GD[i,j]
    y[k]<- ED[i,j] }
  message(stats::cor(x,y),las=1)
  graphics::title(main='Expected Genetic Distance \n versus Realized Genetic
                  Distance',
        font.main=1,cex.main=1.25)
}




#' Calculate Empirical Standard Error
#'
#' This function calculates the empirical standard error based on repeated
#' sampling.
#'
#' @param datafile A data frame containing genetic data from the
#' [pooledpeaks::LoadData]
#' @param NLoci Number of loci to sample in each iteration.
#'
#' @return A numeric vector containing the empirical standard error estimates.
#' @export
#'
#' @examples
#'
#' genetic_data <- data.frame(
#' Locus = c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2),
#' Locus_allele = c("Marker1", "n", 1, 2, 3, "Marker2", "n", 1, 2, 3),
#' Sample1 = c(NA, 10, 0.5, 0.5, 0, NA, 10, 0.2, 0.3, 0.5),
#' Sample2 = c(NA, 20, 0.1, 0.2, 0.7, NA, 20, 0.3, 0.4, 0.3),
#' Sample3 = c(NA, 30, 0.3, 0.4, 0.3, NA, 30, 0.4, 0.2, 0.4)
#' )
#'
#'
#' EmpiricalSE(datafile=genetic_data, NLoci=3)

EmpiricalSE <- function(datafile=data.frame,NLoci=10) {
  G<-0
  for ( i in 1:100)
  { if(i %% 10==0) message(i)
    A <- SampleOfLoci(datafile,NLoci)
    names <- colnames(A)
    names <- names[-(1:2)]

    PopLabel <- substr(names,1,1)
    PopFactor <- as.factor(PopLabel)

    Nx <- TypedLoci(A)
    Jx <- GeneIdentityMatrix(A,Nx)
    G[i] <- GST(Jx)
  }
  return(G)
}




#' Calculate Pre-Jost's D
#'
#' This function calculates the pre-Jost's D measure from a gene identity
#' matrix.
#'
#' @param G A square matrix representing a gene identity matrix.
#'
#' @return The Jost's D value.

preJostD <- function(G=matrix) {

  JS <- 0
  for (i in 1:nrow(G)) {
    JS <- JS+G[i,i]
  }
  JS <- JS/nrow(G)

  JT <-0
  for (i in 1:nrow(G)){
    for (j in 1:nrow(G)){
      JT <- JT+G[i,j]
    }
  }
  JT <- JT/nrow(G)/nrow(G)
  D <- (JT/JS -1)/(1/nrow(G)-1)
  return(D)
}




#' Calculate Jost's D
#'
#' This function calculates Jost's D measure from a gene identity matrix.
#'
#' @param J A gene identity matrix.
#' @param pairwise Logical indicating whether to calculate pairwise Jost's D.
#' If pairwise=FALSE, must not have any missing data.
#'
#' @return If pairwise = TRUE, returns a matrix of pairwise Jost's D values.
#' If pairwise = FALSE, returns the overall Jost's D value.
#'
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
#' JostD(J=gene_identity_matrix, pairwise=TRUE)
#' JostD(J=gene_identity_matrix, pairwise=FALSE)

JostD <- function(J=matrix, pairwise=TRUE) {
  if (pairwise==TRUE){
    N <- nrow(J)
    PJostD <- array(dim=c(N,N))
    for (i in 1:N)
      for (j in 1:N) {
        G <- J[c(i,j),c(i,j)]
        PJostD[i,j] <- preJostD(G)
      }
    return(PJostD)}

  if (pairwise==FALSE)
  {preJostD(J)}

}


#' Pre GST Calculation
#'
#' This function calculates the GST from a gene identity matrix.
#'
#' @param G The gene identity matrix
#'
#' @return The GST value.

preGST <- function(G=matrix) {
  JS <- 0
  for (i in 1:nrow(G)) {
    JS <- JS+G[i,i]
  }
    JS <- JS/nrow(G)

  JT <-0
  for (i in 1:nrow(G)) {
    for (j in 1:nrow(G)) {
      JT <- JT+G[i,j]
    }
  }
  JT <- JT/nrow(G)/nrow(G)

  GST <- (JS-JT)/(1-JT)
  return(GST)
}

#' Nei's GST
#'
#' This function calculates GST (Nei's standard genetic distance) measure from
#' a gene identity matrix.
#'
#' @param J A square matrix representing a gene identity matrix.
#' @param pairwise Logical indicating whether to calculate pairwise GST. If set
#'  to FALSE, must not contain any missing data.
#'
#' @return If pairwise = TRUE, returns a matrix of pairwise GST values.
#' If pairwise = FALSE, returns the overall GST value.
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
#' GST(J=gene_identity_matrix, pairwise=TRUE)
#' GST(J=gene_identity_matrix, pairwise=FALSE)

GST <- function(J=matrix, pairwise=TRUE) {
  if (pairwise==TRUE){
    N <- nrow(J)
    PGST <- array(dim=c(N,N))
    for (i in 1:N)
      for (j in 1:N) {
        G <- J[c(i,j),c(i,j)]
        PGST[i,j] <- preGST(G)
      }
    return(PGST)}

  if (pairwise==FALSE){
    OGST<-preGST(J)

  return(OGST)
  }

}

#' Calculate Two-Level GST
#'
#' This function calculates two-level GST (Nei's standard gene identity)
#' measure from a gene identity matrix.
#'
#' @param G A square matrix representing a gene identity matrix.
#'
#' @return A list containing the components of two-level GST including
#' within-group gene identity, between-group gene identity, and GST values.
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
#' TwoLevelGST(G=gene_identity_matrix)
#'

TwoLevelGST <- function(G=matrix) {
  names <- colnames(G)
  PopLabel <- substr(names,1,1)
  PF <- as.factor(PopLabel)
  z <- levels(PF)
  Gsc <- rep(0,length(z))
  Jc <- rep(0,length(z))
  n <- 0
  Jx <- 0  # Gene Identity Matrix for Component
  Js <- rep(0,length(z)) # Gene Identity Within
  k<-1

  for (k in 1:length(z)) {
    group <- which(PF==z[k])
    Jx <- G[group,group]
    n[k] <- Jx[k]
    Js[k]<-0
    for (i in 1:length(group)){
      Js[k] <- Js[k]+Jx[i,i]
    }
    Js[k] <- Js[k]/length(group)

    for (i in 1:length(group)){
      for (j in 1:length(group)){
        Jc[k] <- Jc[k]+Jx[i,j]
      }
    }

    Jc[k] <- Jc[k]/length(group)/length(group)

    Gsc[k] <- (Js[k]-Jc[k])/(1-Jc[k])
  }

  JS <- 0
  for (i in 1:nrow(G)){
    JS <- JS+G[i,i]
  }
  JS <- JS/nrow(G)

  JC <-0
  for (k in 1:length(z)){
    JC <- JC + n[k]*Jc[k]
  }
  JC <- JC/sum(n)

  JT <-0
  for (i in 1:nrow(G)){
    for (j in 1:nrow(G)){
      JT <- JT + G[i,j]
    }
  }
  JT <- JT/nrow(G)/nrow(G)

  GG <- list(Js=Js,Jc=Jc,Gsc=Gsc,JS=JS,JC=JC,JT=JT,
             GSC=(JS-JC)/(1-JC),GCT=(JC-JT)/(1-JT),GST=(JS-JT)/(1-JT))
  return(GG)
}


#' Calculate Allelic Richness
#'
#' This function calculates allelic richness based on provided genetic data.
#'
#' @param datafile A data frame containing the data as read in by
#' [pooledpeaks::LoadData]
#' @param n A matrix representing the number of markers successfully genotyped
#'  like the output of the [pooledpeaks::TypedLoci] function.
#'
#' @return A vector containing the allelic richness for each locus.
#' @export
#'
#' @examples
#'
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
#' AlRich(datafile=genetic_data,n=n_alleles)

AlRich <- function(datafile=data.frame,n=matrix) {

  loci <- rep(0,length(n))
  loci <- diag(n)

  g <- subset(datafile,(datafile[,2] != 'n')&(datafile[,3]!='NA'))
  g <- g[,-(1:2)]
  g <- as.matrix(g)

  allnum <- rep(0,ncol(g))

  i<-1
  h<-c()

  for (i in 1:ncol(g)) {
    h <- g[,i]
    h <- subset(h,h>0)
    allnum[i] <- length(h)
  }

  allnum <- allnum/loci

  return(allnum)
}

#' Perform Bootstrap Analysis
#'
#' This function performs bootstrap analysis on genetic data.
#'
#' @param A Data frame containing data as read in by [pooledpeaks::LoadData]
#' @param Rep Number of bootstrap replicates.
#' @param Stat Type of statistic to compute (1 for AlRich, 2 for TwoLevelGST)
#'
#' @return Either a matrix of AlRich statistics or a list containing various
#' statistics computed using TwoLevelGST.
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
#' BootStrap3(A=genetic_data, Rep=10, Stat=1)
#' BootStrap3(A=genetic_data, Rep=10, Stat=2)

BootStrap3 <- function(A=data.frame,Rep=20,Stat=1) {
  A <- A
  Rep <- Rep
  Stat <- Stat

  Loci <- max(A[,1])

  names <- colnames(A)
  names <- names[-(1:2)]

  PopLabel <- substr(names,1,1)
  PopFactor <- as.factor(PopLabel)
  K <- length(levels(PopFactor))

  gst <- 0

  if (Stat==1) {
    w <- rep(0,length(A[1,])-2)
  }

  if (Stat==2) {
    js <- rep(0,Rep*K)
    js <- array(js,dim=c(Rep,K))
    jc <- rep(0,Rep*K)
    jc <- array(jc,dim=c(Rep,K))
    gsc <- rep(0,Rep*K)
    gsc <- array(gsc,dim=c(Rep,K))
    JS <- rep(0,Rep)
    JC <- 0
    JT <- 0
    GSC <-0
    GCT <- 0
    GST <- 0
  }

  for (j in 1:Rep) {
    if (j%%10==0){
      message(j)
    }
    z <-A[1,]
    x <- sample(1:Loci,Loci,replace=TRUE)
    for (i in x) {
      b <- subset(A,A[,1]==i)
      z <- rbind(z,b)
    }
    z <- z[-1,]


    if (Stat==1) {
      N <- TypedLoci(z)
      v <- AlRich(z,N)
      if (j==1) {
        w <- v
      }
      else if (j>1)  {
        w <- rbind(w,v)
      }
    }  # End of Stat = 1 loop

    if (Stat==2) {
      nn <- TypedLoci(z)
      Jx<-GeneIdentityMatrix(z,nn)
      G <- TwoLevelGST(Jx)
      for (k in 1:length(G$Js)) {
        js[j,k] <- G$Js[k]
      }
      for (k in 1:length(G$Jc)) {
        jc[j,k] <- G$Jc[k]
      }
      for (k in 1:length(G$Gsc)){
        gsc[j,k] <- G$Gsc[k]
      }
      JS[j] <- G$JS
      JC[j] <- G$JC
      JT[j] <- G$JT
      GST[j] <-G$GST
      GSC[j] <- (JS[j]-JC[j])/(1-JC[j])
      GCT[j] <- (JC[j] -JT[j] )/(1-JT[j])

    }  # end of stat=2 loop
  }   # end of Replicates loop

  if (Stat==1) send_out <- w
  if (Stat == 2) send_out <- list(Js=js,Jc=jc,Gsc=gsc,JS=JS,JC=JC,JT=JT,
                                  GSC=GSC,GCT=GCT,GST=GST)

  return(send_out)

}  # End of Function Bootstrap3
