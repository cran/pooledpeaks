test_that("Function processes .fsa files correctly", {
  # Get the path to the .fsa file within the package
  fsa_file_1 <- system.file("extdata",
                "Multiplex_set_I_Shaem.1a_1_Sample_20221028_215632.fsa",
                package = "pooledpeaks")
  fsa_file_2 <- system.file("extdata",
                "Multiplex_set_I_Shaem.1b_1_Sample_20221028_232301.fsa",
                package = "pooledpeaks")
  fsa_file_3 <- system.file("extdata",
                "Multiplex_set_I_Shaem.3a_2_Sample_20221028_215633.fsa",
                package = "pooledpeaks")
  fsa_file_4 <- system.file("extdata",
                "Multiplex_set_I_Shaem.3b_2_Sample_20221028_232302.fsa",
                package = "pooledpeaks")
  fsa_file_5 <- system.file("extdata",
                "Multiplex_set_I_Shaem.4a_3_Sample_20221028_215634.fsa",
                package = "pooledpeaks")
  fsa_file_6 <- system.file("extdata",
                "Multiplex_set_I_Shaem.4b_3_Sample_20221028_232303.fsa",
                package = "pooledpeaks")
  fsa_file_7 <- system.file("extdata", "23.2a_I_A01_2012-07-18.fsa",
                            package = "pooledpeaks")
  fsa_file_8 <- system.file("extdata", "23.2b_I_A07_2012-07-18.fsa",
                            package = "pooledpeaks")
  fsa_file_9 <- system.file("extdata", "30.3a_I_B01_2012-07-18.fsa",
                            package = "pooledpeaks")
  fsa_file_10 <- system.file("extdata", "30.3b_I_B07_2012-07-18.fsa",
                             package = "pooledpeaks")
  fsa_file_11 <- system.file("extdata", "33.1a_I_C01_2012-07-18.fsa",
                             package = "pooledpeaks")
  fsa_file_12 <- system.file("extdata", "33.1b_I_C07_2012-07-18.fsa",
                             package = "pooledpeaks")


  # Ensure the files exist
  expect_true(file.exists(fsa_file_1))
  expect_true(file.exists(fsa_file_2))
  expect_true(file.exists(fsa_file_3))
  expect_true(file.exists(fsa_file_4))
  expect_true(file.exists(fsa_file_5))
  expect_true(file.exists(fsa_file_6))
  expect_true(file.exists(fsa_file_7))
  expect_true(file.exists(fsa_file_8))
  expect_true(file.exists(fsa_file_9))
  expect_true(file.exists(fsa_file_10))
  expect_true(file.exists(fsa_file_11))
  expect_true(file.exists(fsa_file_12))
})


test_that("fsa_batch_imp processes .fsa files correctly", {
  # Get the path to the folder containing .fsa files within the package
  fsa_folder <- system.file("extdata", package = "pooledpeaks")

  # Ensure the folder and files exist
  expect_true(file.exists(fsa_folder))
  expect_true(length(list.files(fsa_folder, pattern = "\\.fsa$")) > 0)

  # Run the function with default parameters
  result <- fsa_batch_imp(fsa_folder, channels = 5, fourier = TRUE,
                          saturated = TRUE,
                          lets.pullup = FALSE, plotting = FALSE,
                          rawPlot = FALSE,
                          llength = 3000, ulength = 80000)

  # Check the result is a list
  expect_type(result, "list")

  # Check that each element of the list is of type double
  expect_true(all(sapply(result, is.double)))

  # Check that each data frame has the expected number of columns
  # We assume that each .fsa file should result in a df with 'channels' columns
  expected_columns <- 5
  expect_true(all(sapply(result, function(df) ncol(df) == expected_columns)))

  # Check that the Fourier transformation, saturation, and pullup correction
  #(if applied) do not return errors
  # Specific checks for data integrity can be added here
})

test_that("fsa_batch_imp handles missing or incorrect files gracefully", {
  # Test with an empty directory
  temp_dir <- tempdir()
  expect_error(fsa_batch_imp(temp_dir, channels = 5, fourier = TRUE,
                             saturated = TRUE,
                             lets.pullup = FALSE, plotting = FALSE,
                             rawPlot = FALSE,
                             llength = 3000, ulength = 80000),
               "We have not found files with extension .fsa")

  # Test with a non-existent directory
  expect_error(fsa_batch_imp("non_existent_directory", channels = 5,
                             fourier = TRUE, saturated = TRUE,
                             lets.pullup = FALSE, plotting = FALSE,
                             rawPlot = FALSE,
                             llength = 3000, ulength = 80000),
               "We have not found files with extension .fsa")
})

test_that("fsa_batch_imp correctly applies data transformations", {
  fsa_folder <- system.file("extdata", package = "pooledpeaks")

  # Run the function with different transformation options
  result_fourier <- fsa_batch_imp(fsa_folder, channels = 5, fourier = TRUE,
                                  saturated = FALSE,
                                  lets.pullup = FALSE, plotting = FALSE,
                                  rawPlot = FALSE,
                                  llength = 3000, ulength = 80000)

  result_saturated <- fsa_batch_imp(fsa_folder, channels = 5, fourier = FALSE,
                                    saturated = TRUE,
                                    lets.pullup = FALSE, plotting = FALSE,
                                    rawPlot = FALSE,
                                    llength = 3000, ulength = 80000)


  # Check the results are lists
  expect_type(result_fourier, "list")
  expect_type(result_saturated, "list")

})
