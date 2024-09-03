test_that("DistCor function calculates and plots the correlation correctly", {
  # Create a sample genetic distance matrix
  GD <- matrix(c(
    0, 0.1, 0.2, 0.3,
    0.1, 0, 0.1, 0.2,
    0.2, 0.1, 0, 0.1,
    0.3, 0.2, 0.1, 0
  ), nrow = 4, byrow = TRUE)
  rownames(GD) <- colnames(GD) <- c("Sample1", "Sample2", "Sample3", "Sample4")

  # Capture the output of the print statement
  output <- capture_messages({
    DistCor(GD)
  })

  # Print captured output for debugging
  print(output)

  # Check if correlation is printed
  expect_true(any(grepl("[0-9.]+", output)), info = paste("Output was:",
                                                          output))

  # Check for plot elements by opening a mock graphical device
  temp_file <- tempfile(fileext = ".pdf")
  pdf(temp_file)
  DistCor(GD)
  dev.off()

  # Read the contents of the PDF file
  pdf_contents <- pdftools::pdf_text(temp_file)
  unlink(temp_file)

  # Print PDF contents for debugging
  print(pdf_contents)

  # Check if the title is present in the PDF content
  expect_true(any(grepl("Expected Genetic Distance", pdf_contents)),
              info = paste("PDF contents were:", pdf_contents))
})

DistCor <- function(GD=matrix) {
  T <- ape::nj(GD)
  ED <- stats::cophenetic(T)
  graphics::par(mar=c(5,5,5,5))
  plot(GD, ED, xlab="GD", ylab="ED")
  graphics::abline(0, 1, col = 2)

  k <- 0
  x <- numeric()
  y <- numeric()
  for (i in 1:(nrow(GD) - 1)) {
    for (j in (i + 1):ncol(GD)) {
      k <- k + 1
      x[k] <- GD[i, j]
      y[k] <- ED[i, j]
    }
  }

  correlation <- stats::cor(x, y)
  print(correlation, las = 1)
  graphics::title(main = 'Expected Genetic Distance \n versus Realized
                  Genetic Distance', font.main = 1, cex.main = 1.25)
}



# Sample data for TypedLoci, GeneIdentityMatrix, and SampleOfLoci
sample_data <- data.frame(
  Locus = c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2),
  Locus_allele = c("Marker1", "n", 1, 2, 3, "Marker2", "n", 1, 2, 3),
  Sample1 = c(NA, 0, 0, 0, 0, NA, 10, 0.2, 0.3, 0.5),
  Sample2 = c(NA, 20, 0.1, 0.2, 0.7, NA, 20, 0.3, 0.4, 0.3),
  Sample3 = c(NA, 30, 0.3, 0.4, 0.3, NA, 0, 0, 0, 0)
)

# Sample data for preGST and GST
gene_identity_matrix <- matrix(c(
  0.3164550, 0.2836333, 0.2760485,
  0.2836333, 0.3106084, 0.2867215,
  0.2760485, 0.2867215, 0.3338663
), nrow = 3, byrow = TRUE,
dimnames = list(paste0("Sample", 1:3), paste0("Sample", 1:3)))

#### 1. Test for TypedLoci

test_that("TypedLoci function calculates the number of loci correctly", {
  result <- pooledpeaks::TypedLoci(sample_data)

  # Check if the result is a matrix
  expect_true(is.matrix(result))

  # Check if the dimensions of the result are correct
  expect_equal(dim(result), c(3, 3))

  # Check if the values in the result are correct
  expected_result <- matrix(c(1, 1, 0, 1, 2, 1, 0, 1, 1), nrow = 3,
                            byrow = TRUE,dimnames = list(paste0("Sample", 1:3),
                                                         paste0("Sample", 1:3))
                            )
  expect_equal(result, expected_result)
})

#### 2. Test for GeneIdentityMatrix

test_that("GeneIdentityMatrix function calculates the gene identity matrix
          correctly", {
  typed_loci_result <- pooledpeaks::TypedLoci(sample_data)
  result <- pooledpeaks::GeneIdentityMatrix(sample_data, typed_loci_result)

  # Check if the result is a matrix
  expect_true(is.matrix(result))

  # Check if the dimensions of the result are correct
  expect_equal(dim(result), c(3, 3))

  # Example expected result
  expected_result <- matrix(c(0.38, 0.33, NaN, 0.33, 0.44, 0.32,
                              NaN, 0.32, 0.34), nrow = 3,
                            dimnames = list(paste0("Sample", 1:3),
                                            paste0("Sample", 1:3)))
  expect_equal(result, expected_result)
})

#### 3. Test for SampleOfLoci

test_that("SampleOfLoci function samples loci correctly", {
  result <- SampleOfLoci(sample_data, NLoci = 2)

  # Check if the result is a data frame
  expect_true(is.data.frame(result))

  # Check if the number of loci sampled is correct
  expect_equal(nrow(result), 10)
})

#### 4. Test for GST

test_that("GST function calculates GST correctly", {
  result_pairwise <- pooledpeaks::GST(gene_identity_matrix, pairwise = TRUE)

  # Check if the result is a matrix
  expect_true(is.matrix(result_pairwise))

  # Check if the dimensions of the result are correct
  expect_equal(dim(result_pairwise), c(3, 3))

  # Example expected result
  expected_pairwise_result <- matrix(c(0.00000000, 0.02131284, 0.03511043,
                                       0.02131284, 0, 0.02553185, 0.03511043,
                                       0.02553185, 0), nrow = 3)
  expect_equal(result_pairwise, expected_pairwise_result, tolerance = 1e-6)

  result_overall <- pooledpeaks::GST(gene_identity_matrix, pairwise = FALSE)
  expect_true(is.numeric(result_overall))
  expected_overall_result <- 0.03609254
  expect_equal(result_overall, expected_overall_result, tolerance = 1e-5)
})


gene_identity_matrix <- matrix(c(
  0.3164550, 0.2836333, 0.2760485,
  0.2836333, 0.3106084, 0.2867215,
  0.2760485, 0.2867215, 0.3338663
), nrow = 3, byrow = TRUE,
dimnames = list(paste0("Sample", 1:3), paste0("Sample", 1:3)))

#### 1. Test for preJostD

test_that("preJostD function calculates Jost's D correctly", {
  result <- preJostD(gene_identity_matrix)

  # Check if the result is a numeric value
  expect_true(is.numeric(result))

  # Example expected result
  expected_result <- 0.11918291
  expect_equal(result, expected_result, tolerance = 1e-5)
})

#### 2. Test for JostD

test_that("JostD function calculates Jost's D correctly", {
  result_pairwise <- pooledpeaks::JostD(gene_identity_matrix, pairwise = TRUE)

  # Check if the result is a matrix
  expect_true(is.matrix(result_pairwise))

  # Check if the dimensions of the result are correct
  expect_equal(dim(result_pairwise), c(3, 3))

  # Example expected result
  expected_pairwise_result <- matrix(c(0, 0.09536005, 0.1510396, 0.09536005, 0,
                                       0.1102164, 0.15103965, 0.11021643, 0),
                                     nrow = 3)
  expect_equal(result_pairwise, expected_pairwise_result, tolerance = 1e-5)

  result_overall <- pooledpeaks::JostD(gene_identity_matrix, pairwise = FALSE)
  expect_true(is.numeric(result_overall))
  expected_overall_result <- 0.1191829
  expect_equal(result_overall, expected_overall_result, tolerance = 1e-5)
})


gene_identity_matrix <- matrix(c(
  0.3164550, 0.2836333, 0.2760485,
  0.2836333, 0.3106084, 0.2867215,
  0.2760485, 0.2867215, 0.3338663
), nrow = 3, byrow = TRUE,
dimnames = list(paste0("Sample", 1:3), paste0("Sample", 1:3)))

#### 1. Test for preGST

test_that("preGST function calculates GST correctly", {
  result <- preGST(gene_identity_matrix)

  # Check if the result is a numeric value
  expect_true(is.numeric(result))

  # Example expected result
  expected_result <- 0.03609254
  expect_equal(result, expected_result, tolerance = 1e-5)
})

#### 2. Test for GST

test_that("GST function calculates GST correctly", {
  result_pairwise <- GST(gene_identity_matrix, pairwise = TRUE)

  # Check if the result is a matrix
  expect_true(is.matrix(result_pairwise))

  # Check if the dimensions of the result are correct
  expect_equal(dim(result_pairwise), c(3, 3))

  # Example expected result
  expected_pairwise_result <- matrix(c(
    0, 0.02131284, 0.03511043,
    0.02131284, 0, 0.02553185,
    0.03511043, 0.02553185, 0
  ), nrow = 3)
  expect_equal(result_pairwise, expected_pairwise_result, tolerance = 1e-5)

  result_overall <- GST(gene_identity_matrix, pairwise = FALSE)
  expect_true(is.numeric(result_overall))
  expected_overall_result <- 0.03609254
  expect_equal(result_overall, expected_overall_result, tolerance = 1e-5)
})



genetic_distance_matrix <- matrix(c(0.316455, 0.2836333, 0.2760485, 0.2685221,
                                    0.2797302, 0.3202661,0.2836333, 0.3106084,
                                    0.2867215, 0.2687472, 0.2596309, 0.2957862,
                                    0.2760485, 0.2867215, 0.3338663, 0.297918,
                                    0.3057039, 0.3153261,0.2685221, 0.2687472,
                                    0.297918, 0.3107094, 0.2753477, 0.3042383,
                                    0.2797302, 0.2596309, 0.3057039, 0.2753477,
                                    0.3761386, 0.3398558,0.3202661, 0.2957862,
                                    0.3153261, 0.3042383, 0.3398558, 0.4402125
                                    ),nrow = 6, byrow = TRUE,
                                  dimnames = list(c("Sample1", "Sample2",
                                                    "Sample3", "Ind1", "Ind2",
                                                    "Ind3"),
                                                  c("Sample1", "Sample2",
                                                    "Sample3", "Ind1", "Ind2",
                                                    "Ind3")))


genetic_data <- data.frame(
  Locus = c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2),
  Locus_allele = c("Marker1", "n", 1, 2, 3, "Marker2", "n", 1, 2, 3),
  Sample1 = c(NA, 10, 0.5, 0.5, 0, NA, 10, 0.2, 0.3, 0.5),
  Sample2 = c(NA, 20, 0.1, 0.2, 0.7, NA, 20, 0.3, 0.4, 0.3),
  Sample3 = c(NA, 30, 0.3, 0.4, 0.3, NA, 30, 0.4, 0.2, 0.4)
)

n_alleles <- matrix(c(
       2, 2, 2,
      2, 2, 2,
       2, 2, 2
   ), nrow = 3, byrow = TRUE,
   dimnames = list(paste0("Sample", 1:3), paste0("Sample", 1:3)))

#### 1. Test for TwoLevelGST

test_that("TwoLevelGST function calculates GST correctly", {
  result <- TwoLevelGST(genetic_distance_matrix)

  # Check if the result is a list
  expect_true(is.list(result))

  # Check the names of the list
  expected_names <- c("Js", "Jc", "Gsc", "JS", "JC", "JT",
                      "GSC", "GCT", "GST")
  expect_equal(names(result), expected_names)

  # Example expected result
  expected_result <- list(
    Js = c(0.375686833333333, 0.3203099),
    Jc = c(0.329549344444444, 0.294859588888889),
    Gsc = c(0.0688156369250739, 0.0360925437119797),
    JS = 0.3479984,
    JC = 0.3129946,
    JT = 0.3011928,
    GSC = 0.05095117,
    GCT = 0.01688851,
    GST = 0.06697919
  )
  expect_equal(result$Js, expected_result$Js, tolerance = 1e-5)
  expect_equal(result$Jc, expected_result$Jc, tolerance = 1e-5)
  expect_equal(result$Gsc, expected_result$Gsc, tolerance = 1e-5)
  expect_equal(result$JS, expected_result$JS, tolerance = 1e-5)
  expect_equal(result$JC, expected_result$JC, tolerance = 1e-5)
  expect_equal(result$JT, expected_result$JT, tolerance = 1e-5)
  expect_equal(result$GSC, expected_result$GSC, tolerance = 1e-5)
  expect_equal(result$GCT, expected_result$GCT, tolerance = 1e-5)
  expect_equal(result$GST, expected_result$GST, tolerance = 1e-5)
})

#### 2. Test for AlRich

test_that("AlRich function calculates allelic richness correctly", {
  result <- AlRich(genetic_data, n_alleles)

  # Check if the result is a numeric vector
  expect_true(is.numeric(result))

  # Example expected result
  expected_result <- c(2.5, 3, 3)
  names(expected_result)<-c('Sample1', 'Sample2', 'Sample3')
  expect_equal(result, expected_result, tolerance = 1e-5)
})


genetic_data <- data.frame(
  Locus = c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2),
  Locus_allele = c("Marker1", "n", 1, 2, 3, "Marker2", "n", 1, 2, 3),
  Sample1 = c(NA, 10, 0.5, 0.5, 0, NA, 10, 0.2, 0.3, 0.5),
  Sample2 = c(NA, 20, 0.1, 0.2, 0.7, NA, 20, 0.3, 0.4, 0.3),
  Sample3 = c(NA, 30, 0.3, 0.4, 0.3, NA, 30, 0.4, 0.2, 0.4)
)

# Test 1: Check if the output is a matrix when Stat=1
test_that("Output is a matrix when Stat=1", {
  output <- BootStrap3(A = genetic_data, Rep = 20, Stat = 1)
  expect_type(output, "double")
})

# Test 2: Check if the output is a list when Stat=2
test_that("Output is a list when Stat=2", {
  output <- BootStrap3(A = genetic_data, Rep = 20, Stat = 2)
  expect_type(output, "list")
})

# Test 3: Check if the number of replicates in the output matches Rep
test_that("Number of replicates in the output matches Rep", {
  Rep <- 20
  output <- BootStrap3(A = genetic_data, Rep = Rep, Stat = 1)
  expect_equal(nrow(output), Rep)
})

#Test 4: Check if the dimensions of the output list match expected dimensions
test_that("Dimensions of output list match the expected dimensions", {
  Rep <- 20
  K <- 3 # Assuming 3 levels in the factor variable
  output <- BootStrap3(A = genetic_data, Rep = Rep, Stat = 2)
  expect_equal(length(output), 9) # Expected number of elements in the list
})
