test_that("clean_scores function cleans data as expected", {
  # Create some sample data
  scores_data <- list(
    data.frame(Score = c(90, 85, 70), stringsAsFactors = FALSE),
    data.frame(Score = c(80, 75, 60), stringsAsFactors = FALSE)
  )
  rownames(scores_data[[1]]) <- c("104.1a_FA060920_2020-06-09_C05.fsa.1",
                                  "105.2b_FA060920_2020-06-09_C05.fsa.1",
                                  "106.3c_FA060920_2020-06-09_C05.fsa.1")
  rownames(scores_data[[2]]) <- c("107.4d_FA060920_2020-06-09_C05.fsa.1",
                                  "108.5e_FA060920_2020-06-09_C05.fsa.1",
                                  "109.6f_FA060920_2020-06-09_C05.fsa.1")

  # Call the function
  cleaned_data <- pooledpeaks::clean_scores(scores_data,
                                            pattern1 = "_FA.*",
                                            replacement1 = "",
                                            pattern2 = "_.*",
                                            replacement2 = "",
                                            pattern3 = "\\.1*$",
                                            replacement3 = "")

  # Test that the function returns a data frame
  expect_true(isTRUE(is.data.frame(cleaned_data)))

  # Test that the cleaned data frame has expected number of rows
  expect_equal(nrow(cleaned_data), 6)

  # Test that the ID column is cleaned as expected
  expect_equal(unique(cleaned_data$ID), c("104.1a", "105.2b", "106.3c",
                                          "107.4d", "108.5e", "109.6f"))

  # Test that the filename column is cleaned as expected
  expect_equal(unique(cleaned_data$filename),
               c("104.1a_FA060920_2020-06-09_C05.fsa",
                 "105.2b_FA060920_2020-06-09_C05.fsa",
                 "106.3c_FA060920_2020-06-09_C05.fsa",
                 "107.4d_FA060920_2020-06-09_C05.fsa",
                 "108.5e_FA060920_2020-06-09_C05.fsa",
                 "109.6f_FA060920_2020-06-09_C05.fsa"))
})


test_that("lf_to_tdf function transforms data as expected", {
  # Create a sample LF data frame resembling cleaned scores data
  lf_data <- data.frame(ID = c("104.1a", "104.1a", "105.2b", "105.2b"),
                        filename = c("104.1a_FA060920_2020-06-09_C05.fsa",
                                     "104.1a_FA060920_2020-06-09_C05.fsa",
                                     "105.2b_FA060920_2020-06-09_C05.fsa",
                                     "105.2b_FA060920_2020-06-09_C05.fsa"),
                        hei = c(100, 120, 90, 110),
                        pos = c(1, 2, 1, 2),
                        wei = c(117, 120, 123, 126),
                        # Adjust weights to increments of 3
                        stringsAsFactors = FALSE)

  # Call the function
  tdf_data <- pooledpeaks::lf_to_tdf(lf_data)

  # Test that the function returns a data frame
  expect_true(isTRUE(is.data.frame(tdf_data)))

  # Test that the transformed data frame has expected number of rows
  expect_equal(nrow(tdf_data), 4)

  # Test that the transformed data frame has expected number of columns
  expect_equal(ncol(tdf_data), 2) # 3 weights + ID column

  # Test that the row names are correct
  expect_equal(rownames(tdf_data), c("117", "120", "123", "126"))

  # Test that the data values are correct
  expect_equal(tdf_data[ ,1], c("100", "120", "0", "0"))
  # Assuming "117" column is empty for the first row
  expect_equal(tdf_data[ ,2], c("0", "0", "90", "110"))
  # Assuming "120", "123", "126" columns are empty for the second row
})


test_that("data_manipulation function manipulates data as expected", {
  # Create a sample marker data frame
  marker_data <- data.frame(
    Sample1 = c(400, 600, 700),
    Sample2 = c(450, 550, 480),
    Sample3 = c(300, 200, 400)
  )

  # Call the function with default threshold
  manipulated_data <- pooledpeaks::data_manipulation(marker_data)

  # Test that the function returns a data frame
  expect_true(isTRUE(is.data.frame(manipulated_data)))

  # Test that the manipulated data frame has the correct dimensions
  expect_equal(nrow(manipulated_data), 3) # 3 markers
  expect_equal(ncol(manipulated_data), 2) # 2 samples

  # Test that at least one peak for each sample is greater than 500
  expect_true(all(apply(manipulated_data, 2, function(x) any(x > 500))))

  # Call the function with custom threshold
  manipulated_data_custom <- pooledpeaks::data_manipulation(marker_data,
                                                            threshold = 600)

  # Test that >= 1 peak for each sample is greater than the custom threshold
  expect_true(all(apply(manipulated_data_custom, 2, function(x) any(x > 600))))
})


test_that("PCDM function manipulates data as expected", {
  # Create sample consolidated marker data frame
  consolidated_marker <- data.frame(
    Sample1 = c(0.2, 0.3, 0.5),
    Sample2 = c(0.1, 0.2, 0.7),
    Sample3 = c(0.3, 0.4, 0.3)
  )

  # Create sample egg count data frame
  eggcount <- data.frame(
    ID = c("Sample1", "Sample2", "Sample3"),
    n = c(10, 20, 30)
  )

  # Call the function
  manipulated_data <- pooledpeaks::PCDM(consolidated_marker, eggcount,
                                        "Marker1")

  # Test that the function returns a data frame
  expect_true(isTRUE(is.data.frame(manipulated_data)))

  # Test that the manipulated data frame has the correct dimensions
  expect_equal(nrow(manipulated_data), 5) # 3 alleles + 2 header row
  expect_equal(ncol(manipulated_data), 4) # 3 samples + 1 ID column

  # Test that the column names are correct
  expect_equal(colnames(manipulated_data), c("Locus_allele", "Sample1",
                                             "Sample2", "Sample3"))

  # Test that the first column contains the marker name
  expect_equal(manipulated_data[1, "Locus_allele"], "Marker1")

  # Test that the data values are correct
  expect_equal(as.numeric(manipulated_data[3, "Sample1"]), 0.2)
  expect_equal(as.numeric(manipulated_data[4, "Sample2"]), 0.2)
  expect_equal(as.numeric(manipulated_data[5, "Sample3"]), 0.3)
})
