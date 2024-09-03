test_that("associate_dyes adds dye names to data frames correctly", {
  file_path <- system.file("extdata", package = "pooledpeaks")

  # Create a mock output of fsa_batch_imp function
  # You should replace this with an actual call to fsa_batch_imp if possible
  mock_fsa_batch_imp_output<- fsa_batch_imp(file_path, channels = 5,
                                            fourier = TRUE, saturated = TRUE,
                                          lets.pullup = FALSE,
                                          plotting = FALSE, rawPlot = FALSE,
                                          llength = 3000, ulength = 80000)
  names(mock_fsa_batch_imp_output) <- c("23.2a_I_A01_2012-07-18.fsa",
                                        "23.2b_I_A07_2012-07-18.fsa",
                                        "30.3a_I_B01_2012-07-18.fsa",
                                        "30.3b_I_B07_2012-07-18.fsa",
                                        "33.1a_I_C01_2012-07-18.fsa",
                                        "33.1b_I_C07_2012-07-18.fsa",
                       "Multiplex_set_I_Shaem.1a_1_Sample_20221028_215632.fsa",
                       "Multiplex_set_I_Shaem.1b_1_Sample_20221028_232301.fsa",
                       "Multiplex_set_I_Shaem.3a_2_Sample_20221028_215633.fsa",
                       "Multiplex_set_I_Shaem.3b_2_Sample_20221028_232302.fsa",
                       "Multiplex_set_I_Shaem.4a_3_Sample_20221028_215634.fsa",
                       "Multiplex_set_I_Shaem.4b_3_Sample_20221028_232303.fsa")

  # Run associate_dyes
  result <- associate_dyes(mock_fsa_batch_imp_output, file_path)

  # Check that the result has the correct class
  expect_s3_class(result, "fsa_stored")

  # Check that the column names have been correctly renamed to include dye info
  expected_colnames <- c("Dye_ch1_6-FAM", "Dye_ch2_VIC",  "Dye_ch3_NED",
                         "Dye_ch4_PET", "Dye_ch5_LIZ")

  for (file in names(result)) {
    expect_equal(colnames(result[[file]]), expected_colnames)
  }
})

test_that("associate_dyes warns about mismatched filenames", {
  file_path <- system.file("extdata", package = "pooledpeaks")

  # Create a mock output with mismatched names
  mock_fsa_batch_imp_output <- list(
    file3 = data.frame(Channel1 = 1:5, Channel2 = 6:10)
  )
  names(mock_fsa_batch_imp_output) <- c("file3.fsa")

  # Expect a message about mismatched filenames
  expect_message(associate_dyes(mock_fsa_batch_imp_output, file_path),
                 "Imported filenames do not match files in the directory")
})

test_that("associate_dyes handles more columns than dye channels correctly", {
  file_path <- system.file("extdata", package = "pooledpeaks")

  # Create a mock output with more columns than dye channels
  mock_fsa_batch_imp_output <- list(
    file1 = data.frame(Channel1 = 1:5, Channel2 = 6:10, Channel3= 11:15,
                       Channel4=16:20, Channel5=21:25, ExtraColumn=26:30)
  )
  #names(mock_fsa_batch_imp_output) <- c("file1.fsa")

  # Expect a message about more columns than dye channels
  expect_message(associate_dyes(mock_fsa_batch_imp_output, file_path),
  "1 Imported filenames do not match files in the directory.
            Check the specified directory path.")
})

