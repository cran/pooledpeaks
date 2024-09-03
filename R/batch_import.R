#' Batch Import of .fsa files
#'
#' This function imports and extracts all of the information out of the .fsa
#' files and combines them into one list type object.`fsa_batch_imp` is a
#' modification of the original Fragman import script function, `storing.inds`,
#' This revised script accommodates ABI's .fsa file format up to version 3.
#' It retains Fragman functions for Fourier transformation, saturated peaks,
#' and pull-up correction. Notable adjustments include updating channel
#' parameters, utilizing Dyechannel count from the file directory, and
#' streamlining the script by extracting data only from "DATA" tags.
#' Major changes involve column selection for v3 formats and modifications to
#' the "channel" parameter. Minor changes include allowing relative paths for
#' the data directory, importing only .fsa files, and renaming channels with
#' dye names. This revision ensures successful execution for any format version
#' up to 3.
#'
#'
#' @param folder The path to the folder from the current directory where the
#' .fsa files that will be analyzed are stored.
#' @param channels The number of dye channels expected, including the ladder.
#' @param fourier True/False Should fourier transformation be applied.
#' @param saturated True/False whether to Check and correct for saturated peaks.
#' @param lets.pullup True/False Applying pull up correction to the samples to
#' decrease noise from channel to channel. The default is FALSE, please do not
#' change this.
#' @param plotting True/False Should plots be drawn of all channels after data
#' cleaning.
#' @param rawPlot True/False indicating whether a plot should be drawn of all
#' vectors.
#' @param llength A numeric value for the minimum number of indexes in each
#' channel.
#' @param ulength A numeric value for the maximum number fo indexes in each
#' channel.
#'
#' @importFrom Fragman read.abif
#' @importFrom Fragman transfft
#' @importFrom Fragman saturate
#' @importFrom Fragman pullup
#' @importFrom graphics layout
#'
#' @return Output is a LIST where each element of the list is a DATAFRAME with
#' the channels in columns for each FSA file
#' @export
#'
#' @examples
#' file_path <- system.file("extdata", package = "pooledpeaks")
#' fsa_batch_imp(file_path, channels = 5, fourier = FALSE, saturated = FALSE ,
#' lets.pullup = FALSE,plotting = FALSE, rawPlot = FALSE)

fsa_batch_imp <- function(folder, channels = NULL, fourier = TRUE,
                          saturated = TRUE,
                          lets.pullup = FALSE, plotting = FALSE,
                          rawPlot = FALSE, llength = 3000, ulength = 80000) {

  listp2 <- dir(folder, "*.fsa$")

  if (length(listp2) == 0) {
    stop(paste(
      "We have not found files with extension .fsa. Please \nmake sure this is
      the right folder:",
      folder, "\n"
    ), call. = FALSE)
  }

  all.inds.mats <- list(NA)

  message("Reading FSA files")

  tot <- length(listp2)

  for (i in 1:length(listp2)) {
    fsaFile <- suppressWarnings(
      Fragman::read.abif(paste0(folder, "/", listp2[i]))
    )

    fsaFile.data <- fsaFile$Data[(grepl("DATA|Dye#|DyeN", names(fsaFile$Data)))]
    fsaFile.dir <- fsaFile$Directory[(grepl("DATA|Dye#|DyeN",
                                            fsaFile$Directory$name)), ]

    lens <- lapply(fsaFile.data, length)
    aaa <- table(unlist(lens))

    # Get channel number from fsa directory
    channels_d <- fsaFile.data$`Dye#.1`

    if (is.numeric(channels) && (channels != channels_d)) {
      message(paste0("Your specified channel number of (", channels, ") is not
                     the same as\n the number of dye channels recorded in the
                     fsa file (", channels_d, ").\n Results will use ",
                     channels_d, " channels."))
    }

#Set channel param for each loop and accommodate double entries in v3 format
    if (fsaFile$Header$version == 300 && !is.null(channels_d)) {
      channels_l <- 2 * channels_d
    } else {
      channels_l <- channels_d
    }

    # Search for channels candidates, even if channel parameter is set
    cfound <- as.vector(aaa[which(
      as.numeric(names(aaa)) > llength &
        as.numeric(names(aaa)) < ulength
    )])

    # For files where location of indexes is not clear,
    #as seen in some files with shorter runtimes
    if ((length(cfound) > 1) & !(channels_l %in% cfound)) {

      warning(paste("\nYour data for file", listp2[1], "has multiple possible
                places where\nrun indexes could be stored and we don't
                know which is the correct one.\n"))
      prov <- aaa[which(as.numeric(names(aaa)) > llength &
        as.numeric(names(aaa)) < ulength)]

      if (fsaFile$Header$version == 300) {
        prov2 <- matrix(prov, nrow = 1) / 2
      } else {
        prov2 <- matrix(prov, nrow = 1)
      }

      rownames(prov2) <- "number.of.channels.found"
      colnames(prov2) <- paste("Run_Length", names(prov), "indexes", sep = "_")
      prov2 <- rbind(prov2, prov2)
      rownames(prov2)[2] <- "number.to.type.if.selected"
      prov2[2, 1:ncol(prov2)] <- 1:ncol(prov2)

      warning("Please tell us which option has AT LEAST the number of expected
          channels\n\n")

      message(prov2)
      inut <- as.numeric(readline(prompt = "Enter one of the number.to.type: "))
      channels_l <- cfound[inut]
    }

    # Get length of runs (select length with entry # equal to # of channels)
    real.len <- as.numeric(
      names(aaa)[which(
        aaa == channels_l &
          as.numeric(names(aaa)) > llength &
          as.numeric(names(aaa)) < ulength
      )]
    )

    # Identify which DATA tags are of correct length to be runs
    v0 <- as.vector(which(unlist(lens) == real.len))

    # Filter to only those which have identical data size,
    #and are multiple of dye channel number
    chsize <- as.numeric(
      names(which((table(fsaFile.dir$datasize[v0]) %% channels_d) == 0))
    )

    v <- intersect(v0, which(fsaFile.dir$datasize == chsize))

    # Subset list to tags of correct length
    reads <- fsaFile.data[v]

    # Create dataframe of all selected tags
    prov0 <- as.data.frame(do.call(cbind, reads))

    if (fsaFile$Header$version == 300) {
      # Get vector of indexes for columns to keep
      keepcol <- c(
        # Sample channels_l
        (((channels_l - 2) / 2) + 1):(channels_l - 2),
        # Ladder channel
        channels_l
      )

      # Select columns to keep
      prov <- prov0[, keepcol]
    } else {
      # If not v3
      prov <- prov0
    }

    # Add column and row names

    dyenames <- fsaFile.data[grepl("DyeN", names(fsaFile.data))] %>% unlist()

    if (length(dyenames) == ncol(prov)) {
      colnames(prov) <- paste0("ch_", 1:ncol(prov), "__", dyenames)
      rownames(prov) <- paste("index_", 1:nrow(prov), sep = "")
    } else {
      # Add basic column and row names
      colnames(prov) <- paste("channel_", 1:ncol(prov), sep = "")
      rownames(prov) <- paste("index_", 1:nrow(prov), sep = "")
    }

    # Add data frame to list of batch imports
    all.inds.mats[[i]] <- prov

    prov <- NULL # Reset variable for next loop of import

    # Name list entry with original file name
    names(all.inds.mats)[i] <- as.character(listp2[i])
  }


  if (fourier == TRUE) {
    message("Applying Fourier tranformation for smoothing...")
    all.inds.mats <- lapply(all.inds.mats, function(x) {
      apply(x, 2, Fragman::transfft)
    })
  }

  if (saturated == TRUE) {
    message("Checking and correcting for saturated peaks...")
    all.inds.mats <- lapply(all.inds.mats, function(x) {
      apply(x, 2, Fragman::saturate)
    })
  }

  if (lets.pullup == TRUE) {
    message("Applying pull up correction to the samples to decrease noise
            from channel to channel")
    if (plotting == TRUE) {
      all.inds.mats <- lapply(all.inds.mats, Fragman::pullup,
        channel = channels_l, plotting = TRUE
      )
    }
    all.inds.mats <- lapply(all.inds.mats, Fragman::pullup,
      channel = channels_l
    )
  }

  if (rawPlot == TRUE) {
    oldpar<-par(no.readonly = TRUE)
    on.exit(par(oldpar))
    graphics::layout(matrix(1:2, 2, 1))
    coli <- c(
      "blue", "darkgreen", "yellow3",
      "red", "orange", "purple"
    )
    naname <- c("FAM?", "HEX?", "NED?", "ROX?", "LIZ?")
    message("Plotting raw data")

    tot <- length(listp2)
    for (i in 1:dim(all.inds.mats[[1]])[2]) {
      plot(all.inds.mats[[1]][, i],
        col = Fragman::transp(coli[i], 0.6),
        type = "l", ylab = "RFU", main = paste0(
          "Channel ", i," of ", dim(all.inds.mats[[1]])[2]," (", naname[i], ")"
        ),
        cex.main = 1, las = 2
      )
      graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3],
                     graphics::par("usr")[2],
                     graphics::par("usr")[4],
        col = "white"
      )
      if (length(all.inds.mats) > 1) {
        for (j in 1:length(all.inds.mats)) {
          graphics::lines(all.inds.mats[[j]][, i],
            col = Fragman::transp(
              coli[i],
              0.2
            ), lwd = 0.6
          )
        }
      }
    }
  }

  graphics::layout(matrix(1, 1, 1))
  class(all.inds.mats) <- c("fsa_stored")
  return(all.inds.mats)
}
