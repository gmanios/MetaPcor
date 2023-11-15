library(MetaPcor)



# Test Case 1: Check if the function throws an error when an invalid GEO_names is provided.
test_that("Invalid GEO_names should throw an error", {
  expect_error(meta_pcor(GEO_names = c("invalid_name"), option=2, target_namespace = c('ILLUMINA_HUMANREF_8_V3')))
})


# Test Case 2: Check if the function throws an error when an invalid target_namespace is provided.
test_that("Invalid target_namespace should throw an error", {
  expect_error(meta_pcor(GEO_names=c("GSE76427"), option=2, target_namespace = c('blablablablablalbala'),method="sparse", meta_method= "random",pvalue_thres = 0.01,l1  = 0.7 ,l2 = 0))
})



# Test Case 3: Check if mismatch in arguments should throw an error.

test_that("Type mismatch in arguments should throw an error", {
  expect_error(meta_pcor(GEO_names = 123, option=2, target_namespace = c('ILLUMINA_HUMANREF_8_V3')))
  expect_error( meta_pcor(GEO_names=c("GSE76427") ,target_namespace = c(12344), option=2, method="sparse", meta_method= "random", pvalue_thres = 0.01,l1  = 0.7 ,l2 = 0))

  })


# Test Case 4: Check if output is in the correct format.
test_that("Output is in the correct format", {
  result <- meta_pcor(GEO_names=c("GSE76427"), target_namespace = c('ILLUMINA_HUMANREF_8_V3'), option=2, method="sparse", meta_method="random", pvalue_thres = 0.01, l1 = 0.7, l2 = 0)
  expect_true(is.data.frame(result))
})


#Test Case 5: Check if the function gives every time the same result
pcor <- meta_pcor(GEO_names=c("GSE76427") ,target_namespace = c('ILLUMINA_HUMANREF_8_V3'), option=2,
                  method="sparse", meta_method= "random",
                  pvalue_thres = 0.01,l1  = 0.7 ,l2 = 0)

test_that("MetaPcor function", {
  expect_equal(meta_pcor(GEO_names=c("GSE76427") ,target_namespace = c('ILLUMINA_HUMANREF_8_V3'), option=2,
                         method="sparse", meta_method= "random",
                         pvalue_thres = 0.01,l1  = 0.7 ,l2 = 0),pcor)
})


