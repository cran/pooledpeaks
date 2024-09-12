df <- data.frame(
  Multiplex_set_I_Shaem.1a = c(0, 0, 0, 0, 1523, 250, 0, 0),
  Multiplex_set_I_Shaem.1b = c(0, 0, 0, 0, 1410, 216, 0, 0),
  Multiplex_set_I_Shaem.2a = c(0, 0, 0, 0, 0, 0, 24713, 29684),
  Multiplex_set_I_Shaem.3a = c(0, 515, 258, 339, 0, 0, 0, 0),
  Multiplex_set_I_Shaem.3b = c(2897, 4326, 4232, 6122, 6192, 5008, 0, 0),
  Multiplex_set_I_Shaem.4a = c(4876, 6587, 11638, 17054, 23975, 20009, 0, 0),
  Multiplex_set_I_Shaem.4b = c(2839, 4088, 12267, 18243, 31101, 25679, 0, 0),
  Multiplex_set_I_Smans.1a = c(6574, 9802, 8511, 12657, 17247, 12899, 0, 0),
  Multiplex_set_I_Smans.1b = c(4465, 7223, 4801, 7868, 3608, 3485, 0, 0)
)

# Modify row names to match the expected format
rownames(df) <- c(164, 173, 176, 179, 182, 185, 188, 191)

test_that("JostD_KK calculates Jost's D correctly for duplicate samples", {
  # Define the expected values for the test cases
  expected_value_1 <- 8.724602e-05
  expected_value_2 <- 0.372616
  expected_value_3 <- 0.01202364
  expected_value_4 <- 0.10885475


   # Pair 1: Multiplex_set_I_Shaem.1a and Multiplex_set_I_Shaem.1b
  expect_equal(JostD_KK(df$Multiplex_set_I_Shaem.1a,
                        df$Multiplex_set_I_Shaem.1b),
               expected_value_1, tolerance = 1e-5)

  # Pair 2: Multiplex_set_I_Shaem.3a and Multiplex_set_I_Shaem.3b
  expect_equal(JostD_KK(df$Multiplex_set_I_Shaem.3a,
                        df$Multiplex_set_I_Shaem.3b),
               expected_value_2, tolerance = 1e-5)

  # Pair 3: Multiplex_set_I_Shaem.4a and Multiplex_set_I_Shaem.4b
  expect_equal(JostD_KK(df$Multiplex_set_I_Shaem.4a,
                        df$Multiplex_set_I_Shaem.4b),
               expected_value_3, tolerance = 1e-5)

  # Pair 4: Multiplex_set_I_Smans.1a and Multiplex_set_I_Smans.1b
  expect_equal(JostD_KK(df$Multiplex_set_I_Smans.1a,
                        df$Multiplex_set_I_Smans.1b),
               expected_value_4, tolerance = 1e-5)
})







# Test case for Rep_check
test_that("Rep_check correctly processes and flags duplicate samples, and
          handles missing duplicates", {
  result <- Rep_check(df)

  # Check if the result is a data frame
  expect_true(is.data.frame(result))

  # Check the column names of the result
  expected_colnames <- c(
    "Multiplex_set_I_Shaem.1", "Multiplex_set_I_Shaem.2",
    "Multiplex_set_I_Shaem.3", "Multiplex_set_I_Shaem.4",
    "Multiplex_set_I_Smans.1"
  )
  expect_equal(colnames(result), expected_colnames)

  # Check the row names of the result
  expected_rownames <- c(164, 173, 176, 179, 182, 185, 188, 191)
  expect_equal(rownames(result), as.character(expected_rownames))

  # Check the values in the result (average peak heights)
  expected_values <- data.frame(
    Multiplex_set_I_Shaem.1 = c(0, 0, 0, 0, 1466.5, 233.0, 0, 0),
    Multiplex_set_I_Shaem.2 = c(0, 0, 0, 0, 0, 0, 24713, 29684),
    Multiplex_set_I_Shaem.3 = c(1448.5, 2420.5, 2245.0, 3230.5, 3096.0,
                                2504.0, 0, 0),
    Multiplex_set_I_Shaem.4 = c(3857.5, 5337.5, 11952.5, 17648.5, 27538,
                                22844, 0, 0),
    Multiplex_set_I_Smans.1 = c(5519.5, 8512.5, 6656, 10262.5, 10427.5,
                                8192, 0, 0)
  )
  rownames(expected_values) <- expected_rownames
  expect_equal(result, expected_values, tolerance = 1e-5)

  # Check if the function prints the expected messages
  expect_message(Rep_check(df), "Jost D between duplicate samples")
  expect_message(Rep_check(df), "Multiplex_set_I_Shaem.1, Multiplex_set_I_Shaem.2, Multiplex_set_I_Shaem.3, Multiplex_set_I_Shaem.4, Multiplex_set_I_Smans.1")
  #expect_message(Rep_check(df), "8.72460211127635e-05, 0.372616004803687, 0.0120236418304402, 0.108854747323611")
  expect_message(Rep_check(df), "Samples where duplicates have a Jost D exceeding 0.05")
  expect_message(Rep_check(df), "Multiplex_set_I_Shaem.3, Multiplex_set_I_Smans.1")
})



