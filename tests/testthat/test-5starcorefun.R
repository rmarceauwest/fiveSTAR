context("test-5starcorefun")

test_that("node names clean correctly", {
  expect_equal(cleanNodeNames('A %in% c(\"0\") & B < 5'),"A = 0, B < 5")
  expect_equal(cleanNodeNames('A %in% c(\"0, 1, 2\") & B %in% c(\"0\")'),
               'A in {0, 1, 2}, B = 0')
  expect_equal(cleanNodeNames('A %in% c(\"0\") & B %in% c(\"0\")'),
               'A = 0, B = 0')
})
