context("test-5starcorefun")

# test_that("joint family/measure assignment is handled properly", {
#
#   #check for error if wrong measure is set for corresponding family
#   #check for proper handling of null default measure input
#
# })

test_that("node names clean correctly", {
  expect_equal(cleanNodeNames('A %in% c(\"0\") & B < 5'),"A = 0, B < 5")
  expect_equal(cleanNodeNames('A %in% c(\"0, 1, 2\") & B %in% c(\"0\")'),
               'A in {0, 1, 2}, B = 0')
  expect_equal(cleanNodeNames('A %in% c(\"0\") & B %in% c(\"0\")'),
               'A = 0, B = 0')
})


# test_that("cilevel gives correct coverage of tests and CIs",{
#
#   #(think about this behavior -> currently ALWAYS default to give 95% CIs)
#   #test that we get the 95% interval we're expecting
#   #test that nothing is halved, etc.
#   #test behavior is correct when not using default value
#
# })
