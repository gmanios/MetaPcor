pcor <- meta_pcor(GEO_names=c("GSE76427") ,target_namespace = c('ILLUMINA_HUMANREF_8_V3'), option=2, 
          method="sparse", meta_method= "random", 
          pvalue_thres = 0.01,l1  = 0.7 ,l2 = 0)

test_that("MetaPcor function", {
  expect_equal(meta_pcor(GEO_names=c("GSE76427") ,target_namespace = c('ILLUMINA_HUMANREF_8_V3'), option=2, 
            method="sparse", meta_method= "random", 
            pvalue_thres = 0.01,l1  = 0.7 ,l2 = 0),pcor)
})

