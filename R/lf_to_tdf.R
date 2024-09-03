#' Clean Scores Data
#'
#' This function cleans the score_markers_rev3 data by applying specified
#' patterns and replacements to the ID and filename columns.
#'
#' @param scores_data The list containing the output scores data from the
#' score_markers_rev3.
#' @param pattern1 The first pattern to replace in the ID.This is intended to
#' clean up the ID names for when the machine adds substrings to the names. For
#'  example 104.1a_FA060920_2020-06-09_C05.fsa.1 becomes 104.1a using
#'  pattern1="_FA.*" and replacement1= ""
#' @param replacement1 Replacement for the first pattern.
#' @param pattern2 The second pattern to replace in the ID. See pattern1 for
#' more details.
#' @param replacement2 Replacement for the second pattern.
#' @param pattern3 The pattern to replace in the file name.This is intended to
#' clean up the file names for when the machine adds sub strings to the names.
#' For example 104.1a_FA060920_2020-06-09_C05.fsa.1 becomes
#' 104.1a_FA060920_2020-06-09_C05.fsa using pattern3= "\\.1*$" and
#' replacement3= ""
#' @param replacement3 Replacement for the file name pattern.
#'
#' @importFrom dplyr distinct
#'
#' @return A cleaned long format data frame
#' @export
#'
#' @examples
#' scores_data <- list(
#' data.frame(Score = c(90, 85, 70), stringsAsFactors = FALSE),
#' data.frame(Score = c(80, 75, 60), stringsAsFactors = FALSE)
#' )
#' rownames(scores_data[[1]]) <- c("104.1a_FA060920_2020-06-09_C05.fsa_Sa.1",
#'                                 "105.2b_FA060920_2020-06-09_C05.fsa_Sa.1",
#'                                 "106.3c_FA060920_2020-06-09_C05.fsa_Fa.1")
#' rownames(scores_data[[2]]) <- c("107.4d_FA060920_2020-06-09_C05.fsa_Sa.1",
#'                                 "108.5e_FA060920_2020-06-09_C05.fsa_Sa.1",
#'                                 "109.6f_SA060920_2020-06-09_C05.fsa_Fa.1")
#' clean_scores(scores_data,pattern1= "_SA.*", replacement1="",
#' pattern2= "_FA.*",replacement2="")
#'

clean_scores <- function(scores_data, pattern1 = NULL, replacement1 = NULL,
                         pattern2 = NULL, replacement2 = NULL, pattern3 = NULL,
                         replacement3 = NULL) {
  scores_df <- do.call(rbind.data.frame, scores_data)
  scores_df$ID <- rownames(scores_df)
  if (!is.null(pattern1) && !is.null(replacement1)) {
    scores_df$ID <- gsub(pattern1, replacement1, scores_df$ID)
  }
  if (!is.null(pattern2) && !is.null(replacement2)) {
    scores_df$ID <- gsub(pattern2, replacement2, scores_df$ID)
  }
  scores_df$filename <- rownames(scores_df)
  if (!is.null(pattern3) && !is.null(replacement3)) {
    scores_df$filename <- gsub(pattern3, replacement3, scores_df$filename)
  }
  scores_df <- scores_df %>% dplyr::distinct()
  return(scores_df)
}




#' Transform LF to TDF
#'
#' This function transforms a data frame from LF (long format) to TDF (table
#' format),performing various data manipulation steps such as spreading data
#'  across columns,removing NA and/or 0 columns, merging ID allele heights
#'  within each replicate,transposing the table, converting from character to
#'  numeric class, and replacing empty data with "0".
#'
#' @param x A data frame in LF format ideally coming out of the clean_scores
#' function.
#'
#' @importFrom dplyr coalesce
#' @importFrom dplyr select
#' @importFrom dplyr relocate
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise_all
#' @importFrom tidyr pivot_wider
#' @importFrom tibble rowid_to_column
#'
#' @return A transformed data frame in TDF format.
#' @export
#'
#' @examples
#' scores<- data.frame(ID=c("104.1a","105.2b","106.3c","107.4d","108.5e",
#' "109.6f"),
#' filename=c("104.1a_FA060920_2020-06-09_C05.fsa_Sa.1",
#' "105.2b_FA060920_2020-06-09_C05.fsa_Sa.1",
#' "106.3c_FA060920_2020-06-09_C05.fsa_Fa.1",
#' "107.4d_FA060920_2020-06-09_C05.fsa_Sa.1" ,
#' "108.5e_FA060920_2020-06-09_C05.fsa_Sa.1" ,
#' "109.6f_SA060920_2020-06-09_C05.fsa_Fa.1"),
#' hei=c(2000,3000,4000,5000,2500, 1000),
#' pos=c(2000,3000,4000,5000,2500, 1000),
#' wei=c(290,285,280,275,270,260),
#' row.names= c("104.1a_FA060920_2020-06-09_C05.fsa_Sa.1",
#' "105.2b_FA060920_2020-06-09_C05.fsa_Sa.1",
#' "106.3c_FA060920_2020-06-09_C05.fsa_Fa.1",
#' "107.4d_FA060920_2020-06-09_C05.fsa_Sa.1" ,
#' "108.5e_FA060920_2020-06-09_C05.fsa_Sa.1" ,
#' "109.6f_SA060920_2020-06-09_C05.fsa_Fa.1"))
#'
#' lf_to_tdf(scores)

lf_to_tdf <- function(x) {
  if (!all(c("ID", "filename", "hei", "pos", "wei") %in% colnames(x))) {
    stop("Input data frame does not contain required columns: ID, filename,
         hei, pos, wei")
  }
  ID <- x$ID
  filename <- x$filename
  hei <- x$hei
  pos <- x$pos
  wei <- x$wei
  rowid<-x$rowid
  score <- x
  coalesce_by_column <- function(df) {
    return(dplyr::coalesce(!!!as.list(df)))
  }
  # Remove trace position column now that there is bp weight info
  out <- dplyr::select(score, -pos, -filename) %>%
    dplyr::relocate(ID)
  rownames(out) <- NULL
  ## ID trim and table manipulation

  # Spread data across columns - weight by height
  out <- out %>%
    tibble::rowid_to_column() %>%
    tidyr::pivot_wider(names_from = wei, values_from = hei) %>%
    dplyr::select(-rowid)

  # Remove NA and/or 0 columns and replace na with 0
  out <- out[, colnames(out) != "<NA>" & colnames(out) > 0]

  out <- out[, c("ID", sort(setdiff(names(out), "ID")))]

  # Merge ID allele heights within each replicate
  out <- out %>%
    dplyr::group_by(ID) %>%
    dplyr::summarise_all(coalesce_by_column)

  # Transpose table
  tout <- t(out) %>% data.frame()

  # Replace column names after transposing
  colnames(tout) <- tout[1, ]
  tout <- tout[-1, ]

  # Convert from character to numeric class - calculating frequencies
  tout[-1, ] <- sapply(tout[-1, ], as.numeric)

  # Replace empty data with "0"
  tout[is.na(tout)] <- 0

  return(tout)
}



#' Data Manipulation for Marker Data
#'
#' This function ensures that at least one peak for each sample is greater than
#'  a specified threshold (default: 500) and then formats the data frame for
#'  the next steps in the analysis.
#'
#' @param marker A data frame containing marker data, where each row represents
#'  a marker and each column represents a sample.
#' @param threshold The threshold value for peak height. Peaks below this
#' threshold will be replaced with 0.
#'
#' @importFrom dplyr mutate_all
#' @importFrom dplyr mutate
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom magrittr %>%
#' @importFrom tibble rownames_to_column
#' @importFrom tibble column_to_rownames
#' @importFrom rlang .data
#'
#' @return A formatted data frame where at least one peak for each sample is
#' greater than the specified threshold.
#' @export
#'
#' @examples
#'
#' marker_data <- data.frame(
#' Sample1 = c(400, 600, 700,0),
#' Sample2 = c(450, 550, 480,0),
#' Sample3 = c(300, 200, 400,200),
#' Sample4 = c(0,0,0,0),
#' row.names=c(185,188,191,194)
#' )
#' data_manipulation(marker_data,threshold=500)

data_manipulation<-function(marker, threshold=500){

  marker<- marker
  threshold<- threshold
  if(is.na(threshold)){threshold <- 500}

  #Ensure at least one peak is >500
  marker <- as.data.frame(marker[, apply(marker, 2, function(marker)
    any(abs(marker)>threshold)), drop=FALSE])

  marker<- as.data.frame(marker%>%
                           dplyr::mutate_all(~replace(., is.na(.), 0))%>%
                           sapply( as.integer), row.names = rownames(marker))

  #drop the rows that only contain zeros,
  # because these are alleles that weren't truly observed in this set

  marker<- marker%>%
    tibble::rownames_to_column()

  marker$total<-c(rowSums(marker[-1]))

  marker<- marker%>%
    tibble::column_to_rownames()%>%
    dplyr::filter(marker$total != 0)

  marker$total<-NULL
  return(marker)
}


#' Post-consolidation Data Manipulation
#'
#' This function manipulates consolidated marker data and egg count data to
#' prepare them for further analysis.
#'
#' @param consolidated_marker A data frame containing consolidated marker data.
#' @param eggcount A data frame containing egg count data.
#' @param marker_name A string specifying the marker name.
#'
#' @importFrom dplyr left_join
#' @importFrom dplyr relocate
#' @importFrom dplyr rename
#' @importFrom tibble add_row
#' @importFrom tibble rownames_to_column
#'
#' @return A dataframe containing the allele frequencies and eggcounts for each
#'  sample.
#' @export
#'
#' @examples
#' marker_data <- data.frame(
#' Sample1 = c(400, 600, 700),
#' Sample2 = c(450, 550, 480),
#' Sample3 = c(300, 200, 500),
#' row.names=c(185,188,191)
#' )
#'
#' eggs<-data.frame(
#' ID=c("Sample1","Sample2","Sample3"),n=c(3000,400,50))
#'
#' PCDM(consolidated_marker=marker_data, eggcount= eggs,"SMMS2")
#'

PCDM<- function(consolidated_marker=data.frame, eggcount=data.frame,
                marker_name){
  #Ensure input data is in the correct format
  if (!is.data.frame(consolidated_marker) || !is.data.frame(eggcount)) {
    stop("The 'consolidated_marker' and 'eggcount' arguments must be data
         frames.")
  }

  if (!is.character(marker_name) || marker_name == "") {
    stop("The 'marker_name' argument must be a non-empty string.")
  }


  ID <- eggcount$ID
  n <- eggcount$n

  marker<- consolidated_marker
  eggcount$ID<- as.character(eggcount$ID)
  marker_name<-as.vector(marker_name)
  alleles<-row.names(marker)
  alleles_2<-c("Marker_Name", alleles)

  #Convert Heights to Frequencies

  marker <- as.data.frame(sapply(marker, function(x) {
    if (is.numeric(x)) {
      x[is.na(x) | !is.numeric(x)] <- 0
      return(x / sum(x))
    } else {
      return(x)
    }
  }))

  #Adding the alleles back
  rownames(marker)<-alleles

  #Transpose the data frame to match eggcount
  marker<- as.data.frame(t(marker))
  #transposes the df so that participant allele loci are now column names
  marker$ID <- rownames(marker) #adds a column of the ID (aka rownames)
  marker<- marker %>% dplyr::relocate(ID)
  #move the ID column to the first column

  #Match Eggcount and fill the egg count row

  if (!"ID" %in% names(eggcount)) {
    stop("The 'eggcount' data frame must contain an 'ID' column.")
  }

  if (colnames(eggcount)[1] != "ID") {
    eggcount <- eggcount %>%
      dplyr::rename(ID = colnames(eggcount)[1])
  }

  marker<- dplyr::left_join(marker,eggcount, by = "ID")
  #match ID name to ID of eggs
  marker<- marker %>%
    dplyr::relocate(n, .after = 1) # move the n column to the 2nd column

  #Reformat Data
  marker<- as.data.frame(t(marker))
  #transposes the data frame to the original format
  colnames(marker)<- marker[1,]
  marker<-marker[-1,]
  marker<-tibble::rownames_to_column(marker, var = "Locus_allele")

  marker<- marker %>%
    tibble::add_row(.before = 1)
  marker[1,1]<-c(marker_name)

  return(marker)
}
