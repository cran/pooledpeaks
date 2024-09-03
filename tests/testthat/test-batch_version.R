test_that("check_fsa_v_batch correctly identifies number of .fsa files", {
  file_path <- system.file("extdata", package = "pooledpeaks")

  output <- capture_messages(check_fsa_v_batch(file_path))

  expect_true(any(grepl("-- Number of .fsa files found in batch:", output)))
  expect_true(any(grepl("12", output)))
  # Assuming there are 12 .fsa files in the test directory
})

test_that("check_fsa_v_batch correctly identifies version of .fsa files", {
  file_path <- system.file("extdata", package = "pooledpeaks")

  output <- capture_messages(check_fsa_v_batch(file_path))

  expect_true(any(grepl("-- Number of .fsa file formats present in batch:",
                        output)))
  expect_true(any(grepl("v1.01, v3", output)))
  # Adjust the versions based on your test files
})

test_that("check_fsa_v_batch correctly identifies batch names", {
  file_path <- system.file("extdata", package = "pooledpeaks")

  output <- capture_messages(check_fsa_v_batch(file_path))

  expect_true(any(grepl("-- Batch names found in directory:", output)))
  expect_true(any(grepl("PXA2012-07-17, Kathleen_76 Samples_Frag_10-28-22",
                        output)))
  # Adjust the batch names based on your test files
})



test_that("fsa_metadata returns a data frame with correct columns", {
  file_path <- system.file("extdata", package = "pooledpeaks")

  metadata <- fsa_metadata(file_path)

  expected_columns <- c("file_name", "retrieved_sample_name",
                        "batch_container_name", "fsa_version",
                        "user", "run_start_date", "run_start_time",
                        "machine_type", "machineN_serial")

  expect_true(all(expected_columns %in% colnames(metadata)))
})

test_that("fsa_metadata correctly retrieves metadata", {
  file_path <- system.file("extdata", package = "pooledpeaks")

  metadata <- fsa_metadata(file_path)

  # Ensure that there are as many rows as .fsa files
  expect_equal(nrow(metadata), length(list.files(file_path,
                                                 pattern = "\\.fsa$")))

  # Example checks for specific metadata (adjust according to your test files)
  expect_true(all(metadata$fsa_version[1:6] == 1.01))
  expect_true(all(metadata$fsa_version[7:12] == 3.00))
  expect_true(all(metadata$batch_container_name[1:6] == "PXA2012-07-17"))
  expect_true(all(metadata$batch_container_name[7:12] ==
                    "Kathleen_76 Samples_Frag_10-28-22"))
})

test_that("fsa_metadata handles empty directory correctly", {
  empty_dir <- tempdir() # Creating a temporary empty directory

  expect_error(metadata <- fsa_metadata(empty_dir))

})
