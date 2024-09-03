#' Check .fsa Version and Batch Information
#'
#' This function analyzes .fsa files in a specified folder, providing a summary
#' of their version and batch information.
#' @param x The path to the folder from the current directory where the .fsa
#' files that will be analyzed are stored.
#' @importFrom magrittr %>%
#' @importFrom Fragman read.abif
#' @return A written summary of how many .fsa files are in the folder and which
#' version they are.
#' @export
#'
#' @examples
#' file_path <- system.file("extdata", package = "pooledpeaks")
#' check_fsa_v_batch(x = file_path)
#'
check_fsa_v_batch <- function(x) {
  folder <- x
  fsalist <- paste0(folder, "/", dir(folder, "*\\.fsa$"))

  if (length(fsalist[grepl("\\.fsa", fsalist)]) > 0) {
    fsaFile <- lapply(fsalist, function(x)
      suppressWarnings(Fragman::read.abif(x)) %>% try())

    # Get format versions
    vers <- unlist(
      lapply(fsaFile, function(df) {
        vs <- df$Header$version
        v <- unique(vs)
        v / 100
      })
    )

    # Get batch (container) names
    ctnms <- unlist(
      lapply(fsaFile, function(df) {
        ctnm <- df$Data$CTNM.1
      })
    )

    message("-- Number of .fsa files found in batch:", length(fsaFile), "\n")

    message("\n-- Number of .fsa file formats present in batch:",
        paste(paste0("v", unique(vers)), collapse = ", "), "\n")
    if (length(unique(vers)) > 1) {
      message("-- Multiple version types found in directory,
      indicating multiple machine runs.
        Be aware of possible batch-related peak artifacts.", "\n")
    }

    message("\n-- Batch names found in directory:", paste(unique(ctnms),
                                                      collapse = ", "), "\n")

    if (length(unique(ctnms)) > 1) {
      message("-- Multiple batch names found in directory.
        Be aware of possible batch-related peak artifacts.", "\n")
    }
  } else {
    message("-- No fsa files are present in the specified directory")
  }
}





#' Retrieve Metadata
#'
#' Retrieves basic info from .fsa files about the sample and run,and aggregates
#'  multiple samples in a single object.
#' @param x The path to the folder from the current directory where the .fsa
#' files that will be analyzed are stored.
#' @importFrom tibble tibble
#' @import rlang
#' @return A data frame that contains the metadata of the machine and run
#'  extracted from the .fsa file.
#' One row for each .fsa file in directory x and the following columns:
#' retrieved_sample_name, batch_container_name,  fsa_version, user,
#' run_start_date, run_start_time, machine_type,machineN_serial.
#' @export
#'
#' @examples
#' file_path <- system.file("extdata", package = "pooledpeaks")
#' fsa_metadata(x = file_path)
#'
fsa_metadata <- function(x) {
  fsalist <- paste0(x, "/", dir(x, "*\\.fsa$"))
  . <- c()
  fsaFile <- lapply(fsalist, function(x)
    suppressWarnings(Fragman::read.abif(x)) %>% try())

  fsa_metadata <- tibble(
    file_name = gsub(".*/", "", fsalist),
    retrieved_sample_name = NA,
    batch_container_name = NA,
    fsa_version = NA,
    user = NA,
    run_start_date = NA,
    run_start_time = NA,
    machine_type = NA,
    machineN_serial = NA
  )

  # Get sample_name
  fsa_metadata$retrieved_sample_name <- unlist(
    lapply(fsaFile, function(df) {
      df$Data$SpNm.1
    })
  )

  # Get batch_container_name
  fsa_metadata$batch_container_name <- unlist(
    lapply(fsaFile, function(df) {
      df$Data$CTNM.1
    })
  )

  # Get format versions
  fsa_metadata$fsa_version <- unlist(
    lapply(fsaFile, function(df) {
      vs <- df$Header$version
      v <- unique(vs)
      v / 100
    })
  )

  # Get user
  fsa_metadata$user <- unlist(
    lapply(fsaFile, function(df) {
      df$Data$User.1
    })
  )

  # Get run_start_date
  fsa_metadata$run_start_date <- unlist(
    lapply(fsaFile, function(df) {
      df$Data$RUND.1 %>%
        unlist() %>%
        paste0(., collapse = ".")
    }),
    use.names = FALSE
  )

  # Get run_start_time
  fsa_metadata$run_start_time <- unlist(
    lapply(fsaFile, function(df) {
      df$Data$RUNT.1 %>%
        unlist() %>%
        paste0(., collapse = ":")
    }),
    use.names = FALSE
  )

  # Get machine_type
  fsa_metadata$machine_type <- unlist(
    lapply(fsaFile, function(df) {
      df$Data$HCFG.2
    })
  )

  # Get machine
  fsa_metadata$machineN_serial <- unlist(
    lapply(fsaFile, function(df) {
      df$Data$MCHN.1
    })
  )

  return(fsa_metadata)
}
