#' Pairwise Jost D between replicates
#'
#' This function calculates Jost's D between two columns, specifically designed
#'  for comparing duplicate samples based on allele frequencies.
#'
#' @param Ni1 Vector containing the allele frequencies of the first duplicate
#' sample.
#' @param Ni2 Vector containing the allele frequencies of the second duplicate
#' sample.
#'
#' @return The calculated Jost's D value.
JostD_KK <- function(Ni1, Ni2) {
  # Define Ni1 and Ni2 as the peak heights in the columns to be compared
  i <- 1
  Ni1_1 <- c()
  Ni2_1 <- c()
  pi1 <- c()
  pi2 <- c()
  pi1sq2 <- c()
  pi2sq2 <- c()
  pisq2 <- c()
  ai <- c()
  a <- c()
  bi <- c()
  b <- c()
  H1 <- c()
  H2 <- c()
  Ha <- c()
  Da <- c()
  Dg <- c()
  Db <- c()
  Dnominal <- c()
  Dest <- c()
  # Calculate Ni1-1 and Ni2-1
  for (i in 1:length(Ni1)) {
    if (Ni1[i] != 0) {
      Ni1_1 <- c(Ni1_1, Ni1[i] - 1)
    } else {
      Ni1_1 <- c(Ni1_1, 0)
    }
  }

  for (i in 1:length(Ni2)) {
    if (Ni2[i] != 0) {
      Ni2_1 <- c(Ni2_1, Ni2[i] - 1)
    } else {
      Ni2_1 <- c(Ni2_1, 0)
    }
  }
  # Calculate N1 and N2
  N1 <- sum(Ni1)
  N2 <- sum(Ni2)
  # Calculate pi1 and pi2

  for (i in 1:length(Ni1)) {
    if (Ni1[i] != 0) {
      pi1 <- c(pi1, Ni1[i] / N1)
    } else {
      pi1 <- c(pi1, 0)
    }
  }

  for (i in 1:length(Ni2)) {
    if (Ni2[i] != 0) {
      pi2 <- c(pi2, Ni2[i] / N2)
    } else {
      pi2 <- c(pi2, 0)
    }
  }

  # Calculate pi1sq2 and pi2sq2

  for (i in 1:length(pi1)) {
    if (pi1[i] != 0) {
      pi1sq2 <- c(pi1sq2, pi1[i]^2)
    } else {
      pi1sq2 <- c(pi1sq2, 0)
    }
  }

  for (i in 1:length(pi2)) {
    if (pi2[i] != 0) {
      pi2sq2 <- c(pi2sq2, pi2[i]^2)
    } else {
      pi2sq2 <- c(pi2sq2, 0)
    }
  }

  # Calculate pisq2
  for (i in 1:length(pi1)) {
    pisq2 <- c(pisq2, ((pi1[i] + pi2[i]) / 2)^2)
  }

  # Calculate ai
  for (i in 1:length(pi1)) {
    ai <- c(ai, ((pi1[i] + pi2[i])^2 - (pi1[i])^2 - (pi2[i])^2))
  }
  # Calculate a
  a <- sum(ai)

  # Calculate bi
  for (i in 1:length(pi1)) {
    bi <- bi + c((Ni1[i] * (Ni1_1[i]) /
                    (N1 * (N1 - 1))) + (Ni2[i] * (Ni2_1[i] / (N2 * (N2 - 1)))))
  }
  # Calculate b
  b <- sum(bi)
  # Calculate H1
  H1 <- 1 - sum(pi1sq2)
  # Calculate H2
  H2 <- 1 - sum(pi2sq2)
  # Calculate Ha
  Ha <- 0.5 * (H1 + H2)
  # Calculate Da
  Da <- 1 / (1 - Ha)
  # Calculate Dg
  Dg <- 1 / sum(pisq2)
  # Calculate Db
  Db <- Dg / Da
  # Calculate Dnominal
  Dnominal <- -2 * ((1 / Db) - 1)
  # Calculate D est
  Dest <- 1 - a / b

  return(Dnominal)
}


#' Replicate Check for Duplicate Samples
#'
#' This function checks for duplicate samples in the input data frame and
#' calculates the average peak heights for each sample. If the Jost's D
#' between duplicate samples exceeds 0.05, it flags those samples.
#'
#' @param df The input data frame containing peak heights for each sample.
#'
#' @importFrom dplyr select
#' @importFrom dplyr select_if
#' @importFrom magrittr %>%
#'
#' @return A data frame containing the average peak heights for each sample,
#' with flagged samples where duplicates have a Jost's D exceeding 0.05.
#' @export
#'
#' @examples
#' marker_data <- data.frame(
#' Sample.1a = c(400, 600, 700),
#' Sample.1b = c(420, 606, 710),
#' Sample.2a = c(450, 550, 480),
#' Sample.2b = c(500, 540, 480),
#' Sample.3a = c(300, 200, 500),
#' Sample.3b = c(290, 100, 400),
#' row.names=c(185,188,191)
#' )
#'
#' Rep_check(marker_data)
#'
Rep_check <- function(df) {
  if (ncol(df) < 2) {
    stop("Input data frame must have at least two columns")
  }
  # Define placeholder variables
  newdf <- df
  repD <- c()
  redflag <- c()
  temp_df <- data.frame(matrix(nrow = nrow(newdf), ncol = 2))
  newmarker <- data.frame(matrix(nrow = nrow(df), ncol = ncol(df)))
  sampleID <- character()

  # Sort df so duplicate columns are next to eachother
  newdf <- newdf %>%
    dplyr::select(sort(names(newdf)))

  # Remove a and b demarcation from column names

  colnames(newdf) <- gsub("(?<=\\d)(a|b|c|d|e|f|ab|ac|bb|bc)", "",
                          colnames(newdf), perl = TRUE)

  # Define Columnames vector
  Columnames <- colnames(newdf)

  # Define counters
  j <- 1 # for vectors
  k <- 1 # newmarker df
  sample_counter <- 1

  # For Loop to calculate the average and JostD of duplicate columns
  i <- 1
  while (i <= ncol(newdf)) {
    count_dup <- sum(Columnames == Columnames[i])

    if (count_dup == 1) {
      newmarker[, k] <- newdf[, i]
      sampleID[sample_counter] <- Columnames[i]
      repD[j] <- "NA" # Indicate no duplicates for this column
      redflag[j] <- "NA" # Indicate no duplicates for this column

      sample_counter <- sample_counter + 1
      k <- k + 1
    } else if (count_dup > 1) {
      if (i < ncol(newdf) && Columnames[i] == Columnames[i + 1]) {
        repD[j] <- JostD_KK(newdf[, i], newdf[, i + 1])
        temp_df <- data.frame(newdf[, i], newdf[, i + 1])

        if (repD[j] <= 0.05) {
          newmarker[, k] <- rowMeans(temp_df, na.rm = TRUE)
          sampleID[sample_counter] <- Columnames[i]
        } else {
          newmarker[, k] <- rowMeans(temp_df, na.rm = TRUE)
          redflag[j] <- Columnames[i]
          sampleID[sample_counter] <- Columnames[i]
        }
        j <- j + 1
        i <- i + 1
      } else {
        newmarker[, k] <- newdf[, i]
        sampleID[sample_counter] <- Columnames[i]
        repD[j] <- "NA" # Indicate no duplicates for this column
      }

      sample_counter <- sample_counter + 1
      k <- k + 1
    }
    i <- i + 1
  }

  sampleID <- sampleID[!is.na(sampleID)]
  avg_height <- newmarker %>% dplyr::select_if(~ !all(is.na(.)))
  colnames(avg_height) <- sampleID
  rownames(avg_height) <- rownames(newdf)
  repD <- repD[!is.na(repD)]
  redflag <- redflag[!is.na(redflag)]


  message("Jost D between duplicate samples")
  message(paste(sampleID, collapse = ", "))
  message(paste(repD, collapse = ", "))
  message("Samples where duplicates have a Jost D exceeding 0.05")
  message(paste(redflag, collapse = ", "))
  return(avg_height)
}
