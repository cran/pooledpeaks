test_that("LoadData function imports data as expected", {
  # Create a temporary file with sample data
  temp_file <- tempfile()
  write.table(data.frame(Locus_allele=c("Marker1","n","1","2","3"),
                         Sample1= c(NA, "10","0.2","0.3","0.5"),
                         Sample2= c(NA, "20","0.1","0.2","0.7"),
                         Sample3= c(NA, "30","0.3","0.4","0.3")),
              file = temp_file, row.names = FALSE)

  # Call the function with the temporary file path
  imported_data <- pooledpeaks::LoadData(temp_file)

  # Test that the function returns a data frame
  expect_true(isTRUE(is.data.frame(imported_data)))

  # Test that the imported data frame has the correct dimensions
  expect_equal(nrow(imported_data), 5) # 5 rows
  expect_equal(ncol(imported_data), 5) # 5 columns

  # Test that the column names are correct
  expect_equal(colnames(imported_data),
               c("Locus","Locus_allele", "Sample1","Sample2","Sample3"))

  # Test that the Locus column values are correct
  expect_equal(imported_data$Locus_allele, c('Marker1', 'n', '1', '2', '3'))

  # Cleanup the temporary file
  unlink(temp_file)
})


test_that("TypedLoci function calculates typed loci as expected", {
  # Create a sample data frame resembling the desired format with two markers
  datafile <- data.frame(
    Locus = c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2),
    Locus_allele = c("Marker1", "n", 1, 2, 3, "Marker2", "n", 1, 2, 3),
    Sample1 = c(NA, 0, 0, 0, 0, NA, 10, 0.2, 0.3, 0.5),
    Sample2 = c(NA, 20, 0.1, 0.2, 0.7, NA, 20, 0.3, 0.4, 0.3),
    Sample3 = c(NA, 30, 0.3, 0.4, 0.3, NA, 0, 0, 0, 0)
  )

  # Call the function with the sample data frame
  processed_data <- pooledpeaks::TypedLoci(datafile)

  # Test that the function returns a matrix
  expect_true(isTRUE(is.matrix(processed_data)))

  # Test that the processed data matrix has the correct dimensions
  expect_equal(nrow(processed_data), 3) # 3 rows (number of individuals)
  expect_equal(ncol(processed_data), 3) # 3 columns

  row<-c("Sample1", "Sample2","Sample3")
  col<-c("Sample1", "Sample2","Sample3")

  rowcol<- list(row,col)

  # Test that the processed data matrix has the correct values
  expect_equal(processed_data, matrix(c(1, 1, 0, 1, 2, 1, 0, 1, 1), nrow = 3,
                                      ncol = 3, dimnames = rowcol))
})


test_that("GeneIdentityMatrix function calculates gene identity matrix as
          expected", {
  # Create a sample dataset similar to the one used in the TypedLoci function
  datafile <- data.frame(
    Locus = c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2),
    Locus_allele = c("Marker1", "n", 1, 2, 3, "Marker2", "n", 1, 2, 3),
    Sample1 = c(NA, 0, 0, 0, 0, NA, 10, 0.2, 0.3, 0.5),
    Sample2 = c(NA, 20, 0.1, 0.2, 0.7, NA, 20, 0.3, 0.4, 0.3),
    Sample3 = c(NA, 30, 0.3, 0.4, 0.3, NA, 0, 0, 0, 0)
  )

  # Call the TypedLoci function with the sample dataset
  typed_loci <- pooledpeaks::TypedLoci(datafile)

  # Call GeneIdentityMatrix function with sample dataset and typed loci matrix
  gene_identity <- pooledpeaks::GeneIdentityMatrix(datafile, typed_loci)

  # Test that the function returns a matrix
  expect_true(isTRUE(is.matrix(gene_identity)))

  # Test the dimensions of the gene identity matrix
  expect_equal(nrow(gene_identity), nrow(typed_loci))
  expect_equal(ncol(gene_identity), nrow(typed_loci))

  # Test some values to ensure correctness
  # Here you can include some specific checks for the expected values in the
  # gene identity matrix
  # based on your knowledge of the input data and the calculations performed
  #by the function
})


test_that("GeneIdentityMatrix function calculates gene identity matrix as
          expected", {
  # Create a sample dataset similar to the one used in the TypedLoci function
  datafile <- data.frame(
    Locus = c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2),
    Locus_allele = c("Marker1", "n", 1, 2, 3, "Marker2", "n", 1, 2, 3),
    Sample1 = c(NA, 0, 0, 0, 0, NA, 10, 0.2, 0.3, 0.5),
    Sample2 = c(NA, 20, 0.1, 0.2, 0.7, NA, 20, 0.3, 0.4, 0.3),
    Sample3 = c(NA, 30, 0.3, 0.4, 0.3, NA, 0, 0, 0, 0)
  )

  # Call the TypedLoci function with the sample dataset
  typed_loci <- pooledpeaks::TypedLoci(datafile)

  # Call GeneIdentityMatrix function with sample dataset and typed loci matrix
  gene_identity <- pooledpeaks::GeneIdentityMatrix(datafile, typed_loci)

  # Test that the function returns a matrix
  expect_true(isTRUE(is.matrix(gene_identity)))

  # Test the dimensions of the gene identity matrix
  expect_equal(nrow(gene_identity), nrow(typed_loci))
  expect_equal(ncol(gene_identity), nrow(typed_loci))

  # Test specific values in the gene identity matrix
  expect_equal(gene_identity["Sample1", "Sample1"], 0.38)
  expect_equal(gene_identity["Sample1", "Sample2"], 0.33)
  expect_true(is.nan(gene_identity["Sample1", "Sample3"]))
  # Check for NaN value
  expect_equal(gene_identity["Sample2", "Sample1"], 0.33)
  expect_equal(gene_identity["Sample2", "Sample2"], 0.44)
  expect_equal(gene_identity["Sample2", "Sample3"], 0.32)
  expect_true(is.nan(gene_identity["Sample3", "Sample1"]))
  # Check for NaN value
  expect_equal(gene_identity["Sample3", "Sample2"], 0.32)
  expect_equal(gene_identity["Sample3", "Sample3"], 0.34)

})

test_that("GeneticDistanceMatrix function calculates genetic distance matrix as
          expected", {
  # Define a sample gene identity matrix
  J <- matrix(c(0.38, 0.33, NA,
                0.33, 0.44, 0.32,
                NA, 0.32, 0.34), nrow = 3, byrow = TRUE,
              dimnames = list(paste0("Sample", 1:3), paste0("Sample", 1:3)))

  # Call GeneticDistanceMatrix function with the sample gene identity matrix
  genetic_distance <- pooledpeaks::GeneticDistanceMatrix(J)

  # Test that the function returns a matrix
  expect_true(isTRUE(is.matrix(genetic_distance)))

  # Test the dimensions of the genetic distance matrix
  expect_equal(nrow(genetic_distance), nrow(J))
  expect_equal(ncol(genetic_distance), ncol(J))

  # Test specific values in the genetic distance matrix
  expect_equal(genetic_distance["Sample1", "Sample1"], 0)
  expect_equal(genetic_distance["Sample1", "Sample2"], 0.08)
  expect_true(is.na(genetic_distance["Sample1", "Sample3"]))
  # Check for NaN value
  expect_equal(genetic_distance["Sample2", "Sample1"], 0.08)
  expect_equal(genetic_distance["Sample2", "Sample2"], 0)
  expect_equal(genetic_distance["Sample2", "Sample3"], 0.07)
  expect_equal(genetic_distance["Sample3", "Sample2"], 0.07)
  expect_equal(genetic_distance["Sample3", "Sample3"], 0)
  expect_true(is.na(genetic_distance["Sample3", "Sample1"]))
  # Check for NaN value
})


test_that("RWCDistanceMatrix function calculates RWC distance matrix as
          expected", {
  # Define a sample genetic distance matrix
  J <- matrix(c(0, 0.08, NA,
                0.08, 0, 0.07,
                NA, 0.07, 0), nrow = 3, byrow = TRUE,
              dimnames = list(paste0("Sample", 1:3), paste0("Sample", 1:3)))

  # Call the RWCDistanceMatrix function with the sample genetic distance matrix
  rwc_distance <- pooledpeaks::RWCDistanceMatrix(J)

  # Test that the function returns a matrix
  expect_true(isTRUE(is.matrix(rwc_distance)))

  # Test the dimensions of the RWC distance matrix
  expect_equal(nrow(rwc_distance), nrow(J))
  expect_equal(ncol(rwc_distance), ncol(J))

  # Test specific values in the RWC distance matrix
  expect_equal(rwc_distance["Sample1", "Sample1"], 0)
  expect_equal(rwc_distance["Sample1", "Sample2"], -0.0869565217391304)
  expect_true(is.na(rwc_distance["Sample1", "Sample3"]))
  # Check for NaN value
  expect_equal(rwc_distance["Sample2", "Sample1"], -0.0869565217391304)
  expect_equal(rwc_distance["Sample2", "Sample2"], 0)
  expect_equal(rwc_distance["Sample2", "Sample3"], -0.0752688172043011)
  expect_equal(rwc_distance["Sample3", "Sample2"], -0.0752688172043011)
  expect_equal(rwc_distance["Sample3", "Sample3"], 0)
  expect_true(is.na(rwc_distance["Sample3", "Sample1"]))
  # Check for NaN value
})


