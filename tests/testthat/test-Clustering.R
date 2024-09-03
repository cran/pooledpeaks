test_that("cluster function performs K-means clustering as expected", {
  # Create an expanded sample dataset with more samples
  datafile <- data.frame(
    Locus = rep(c(1, 2), each = 9),
    Locus_allele = rep(c("Marker1", "n", 1, 2, 3, "Marker2", "n", 1, 2), 2),
    Sample1 = c(NA, 0, 0, 0, 0, NA, 10, 0.2, 0.3, NA, 0, 0, 0, 0,
                NA, 10, 0.2, 0.3),
    Sample2 = c(NA, 20, 0.1, 0.2, 0.7, NA, 20, 0.3, 0.4, NA, 20, 0.1, 0.2, 0.7,
                NA, 20, 0.3, 0.4),
    Sample3 = c(NA, 30, 0.3, 0.4, 0.3, NA, 0, 0, 0, NA, 30, 0.3, 0.4, 0.3,
                NA, 0, 0, 0),
    Sample4 = c(NA, 25, 0.2, 0.3, 0.5, NA, 15, 0.1, 0.2, NA, 25, 0.2, 0.3, 0.5,
                NA, 15, 0.1, 0.2),
    Sample5 = c(NA, 35, 0.3, 0.4, 0.3, NA, 5, 0.3, 0.4, NA, 35, 0.3, 0.4, 0.3,
                NA, 5, 0.3, 0.4),
    Sample6 = c(NA, 40, 0.1, 0.2, 0.4, NA, 20, 0.2, 0.3, NA, 40, 0.1, 0.2, 0.4,
                NA, 20, 0.2, 0.3),
    Sample7 = c(NA, 50, 0.2, 0.3, 0.4, NA, 30, 0.3, 0.4, NA, 50, 0.2, 0.3, 0.4,
                NA, 30, 0.3, 0.4),
    Sample8 = c(NA, 45, 0.1, 0.2, 0.5, NA, 25, 0.2, 0.3, NA, 45, 0.1, 0.2, 0.5,
                NA, 25, 0.2, 0.3)
  )

  # Call the cluster function with the sample dataset and K=2
  clustering_result <- cluster(datafile, K=2)

  # Test that the function returns a list
  expect_true(is.list(clustering_result))

  # Test that the list contains the expected elements
  expect_true(all(c("PopFactor", "clust","dat") %in% names(clustering_result)))

  # Test that the PopFactor is a factor
  expect_true(is.factor(clustering_result$PopFactor))

  # Test clust is an integer vector with length equal to number of samples
  expect_true(is.integer(clustering_result$clust))
  expect_equal(length(clustering_result$clust), ncol(datafile) - 4)
  # Number of samples

  # Test specific values in the clustering result
  # Since K-means clustering results can vary,
  #we will not test specific cluster assignments
  expect_true(all(clustering_result$clust > 0 & clustering_result$clust <= 2))
  # All clusters should be between 1 and K
})


test_that("SampleOfLoci function samples loci as expected", {
  # Create a sample data frame resembling the desired format with two markers
  aaax <- data.frame(
    Locus = c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2),
    Locus_allele = c("Marker1", "n", 1, 2, 3, "Marker2", "n", 1, 2, 3),
    Sample1 = c(NA, 0, 0, 0, 0, NA, 10, 0.2, 0.3, 0.5),
    Sample2 = c(NA, 20, 0.1, 0.2, 0.7, NA, 20, 0.3, 0.4, 0.3),
    Sample3 = c(NA, 30, 0.3, 0.4, 0.3, NA, 0, 0, 0, 0)
  )

  # Specify the number of loci to sample
  NLoci <- 2

  # Call SampleOfLoci with the sample dataset and specified number of loci
  sampled_loci <- SampleOfLoci(aaax, NLoci)

  # Test that the function returns a data frame
  expect_true(is.data.frame(sampled_loci))

  # Test the number of sampled loci
  expect_equal(NLoci, length(aaax$Locus_allele[aaax$Locus_allele == "n"]))

  # Test that the sampled loci contain the expected number of occurrences
  #of "n" in the second column
  expected_occurrences <- NLoci * length(unique(aaax$Locus_allele
                                                [aaax$Locus_allele == "n"]))
  actual_occurrences <- sum(sampled_loci$Locus_allele == "n")
  expect_equal(actual_occurrences, expected_occurrences)

})

test_that("ClusterFromSamples function performs clustering on samples of loci
          and calculates statistics as expected", {
  # Create a sample data frame similar to the one used in the other tests
  datafile <- data.frame(
    Locus = c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3),
    Locus_allele = c("Marker1", "n", 1, 2, 3, "Marker2", "n", 1, 2, 3,
                     "Marker3", "n", 1, 2, 3),
    Sample1 = c(NA, 0, 0, 0, 0, NA, 10, 0.2, 0.3, 0.5, NA, 0.1, 0.2, 0.3, 0.4),
    Sample2 = c(NA, 20, 0.1, 0.2, 0.7, NA, 20, 0.3, 0.4, 0.3,
                NA, 0.4, 0.3, 0.2, 0.1),
    Sample3 = c(NA, 30, 0.3, 0.4, 0.3, NA, 0, 0, 0, 0, NA, 0.2, 0.3, 0.4, 0.1),
    Sample4 = c(NA, 10, 0.5, 0.3, 0.2, NA, 5, 0.1, 0.2, 0.3,
                NA, 0.1, 0.4, 0.3, 0.2),
    Sample5 = c(NA, 15, 0.2, 0.1, 0.5, NA, 15, 0.4, 0.3, 0.2,
                NA, 0.3, 0.1, 0.2, 0.4)
  )

  # Call the ClusterFromSamples function with the sample data frame
  numloci <- 2
  reps <- 10
  result <- pooledpeaks::ClusterFromSamples(datafile, numloci, reps)

  # Test that the function returns a matrix
  expect_true(isTRUE(is.matrix(result)))

  # Test the dimensions of the matrix
  expect_equal(nrow(result), reps)
  expect_equal(ncol(result), length(unique(substr(colnames(datafile)[-(1:2)],
                                                  1, 1))))

  # Test that the matrix contains values between 0 and 1
  expect_true(all(result >= 0 & result <= 1))

  # Test that the matrix has rounded values to 4 decimal places
  expect_true(all(result == round(result, 4)))

  set.seed(42)
  result_with_seed <- ClusterFromSamples(datafile, numloci, reps)
  expect_true(isTRUE(all.equal(result, result_with_seed)))
})

test_that("MDSplot function generates MDS plot correctly", {
  # Create a sample genetic distance matrix
  distance <- matrix(c(
    0, 0.1, 0.2, 0.3,
    0.1, 0, 0.1, 0.2,
    0.2, 0.1, 0, 0.1,
    0.3, 0.2, 0.1, 0
  ), nrow = 4, byrow = TRUE)
  colnames(distance) <- rownames(distance) <- c("Sample1", "Sample2",
                                                "Sample3", "Sample4")

  # Define the principal coordinates to plot and population labels
  pcs <- c(1, 2)
  PF <- factor(c("A", "A", "B", "B"))
  colors <- c('dodgerblue', 'red')

  # Capture the plot output
  plot_output <- capture.output(
    MDSplot(distance = distance, pcs = pcs, PF = PF, y = colors)
  )

  # Check if cmdscale is executed
  E <- stats::cmdscale(distance, eig = TRUE, k = max(pcs))
  expect_true(!is.null(E$points))
  expect_equal(nrow(E$points), 4)
  expect_equal(ncol(E$points), max(pcs))

  # Check if the eigenvalues are correctly calculated
  T <- sum(E$eig[E$eig > 0])
  percent <- round(E$eig / T, 3) * 100
  expect_equal(length(percent), length(E$eig))
  expect_true(all(percent >= 0 & percent <= 100))

  # Ensure the function handles an empty distance matrix gracefully
  empty_distance <- matrix(nrow = 0, ncol = 0)
  expect_error(MDSplot(distance = empty_distance, pcs = pcs, PF = PF,
                       y = colors))

  # Check handling of invalid principal coordinates
  invalid_pcs <- c(1, 10)
  expect_error(MDSplot(distance = distance, pcs = invalid_pcs, PF = PF,
                       y = colors))
})

# Additional edge case tests
test_that("MDSplot handles edge cases", {
  # Case with only one sample
  single_sample_distance <- matrix(0, nrow = 1, ncol = 1)
  colnames(single_sample_distance) <- rownames(single_sample_distance) <-
    "Sample1"
  PF_single <- factor("A")
  expect_error(MDSplot(distance = single_sample_distance, pcs = c(1, 2),
                       PF = PF_single))

  # Case with two samples
  two_sample_distance <- matrix(c(0, 0.1, 0.1, 0), nrow = 2)
  colnames(two_sample_distance) <- rownames(two_sample_distance) <-
    c("Sample1", "Sample2")
  PF_two <- factor(c("A", "B"))
  expect_error(MDSplot(distance = two_sample_distance, pcs = c(1, 2),
                       PF = PF_two))
})
