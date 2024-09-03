#' Associate Dye Names in Batch Import Output
#'
#' This function associates dye info with fragman channel names. It was
#' designed to be performed on any fsa formats after final columns are
#' correctly imported.
#' @param x The Output list of data frames from fsa_batch_imp.
#' @param y The path to the folder from the current directory where the
#' .fsa files that will be analyzed are stored.
#'
#' @importFrom Fragman read.abif
#' @return The input dataframe with an added column assigning fluorescent
#' dye colors.
#' @export
#'
#' @examples
#' y <- system.file("extdata", package = "pooledpeaks")
#' x <- fsa_batch_imp(y, channels = 5, fourier = FALSE, saturated = FALSE ,
#' lets.pullup = FALSE, plotting = FALSE, rawPlot = FALSE,
#' llength = 3000, ulength = 80000 )
#' associate_dyes(x,y)

associate_dyes <- function(x,y){

  frag_obj <- x
  folder <- y

  fsalist <- paste0(folder, "/", dir(folder, "*.fsa$"))

  # Import fsa files with read.abif
  fsaFile<-lapply(fsalist, function(x) suppressWarnings(Fragman::read.abif(x)))

  # Add names to the fsa list, from the absolute paths
  names(fsaFile) <- gsub(".*/","",as.character(fsalist))

  # Check presence of names in specified directory
  chknm <- unique((names(frag_obj)%in% names(fsaFile)))

  if (FALSE %in% chknm ){
    mimatchnm <-  table(names(frag_obj)%in% names(fsaFile))
    message(paste0(mimatchnm["FALSE"],
    " Imported filenames do not match files in the directory.
            Check the specified directory path."))
  }

  if (!FALSE %in% chknm ) {

    for (i in names(fsaFile)) {

      fsa <- fsaFile[[i]]

      # Extract dye names for each channel and create column name list
      dyedf <- unlist(fsa$Data[grepl("DyeN",names(fsa$Data))])

      dyecolnm <- paste0("Dye_ch", 1:length(dyedf),"_", dyedf)

      # Replace column names in imported object

      if (length(colnames(frag_obj[[i]])) > length(dyecolnm)){
        message("This import object has more columns than Dye channels.
                If the .fsa files are v3 (SeqStudio), make sure you import
                with the 'storing_inds_rev2' command.")
      }

      colnames(frag_obj[[i]]) <- dyecolnm

    }

    message("Columns have been renamed with associated fluoresecent dye.")

    class(frag_obj) <- "fsa_stored"
    return(frag_obj)
  }
}
