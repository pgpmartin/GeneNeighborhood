context("Test getTES function")
library(GenomicRanges)

test_that("the output of getTES function returns the TES", {

  testgn <- Genegr[1:10]
  res <- getTES(testgn)
  expect_is(res, "GRanges")
  expect_equal(start(res),
               end(res))
  expect_equal(as.character(strand(res)),
               as.character(strand(testgn)))
  expect_equal(start(res)[as.character(strand(testgn))=="+"],
               end(testgn)[as.character(strand(testgn))=="+"])
  expect_equal(start(res)[as.character(strand(testgn))=="-"],
               start(testgn)[as.character(strand(testgn))=="-"])

})

test_that("getTES throws a message for '*' strands", {
  gg <- Genegr[1:5]
  strand(gg)[1] <- "*"
  expect_message(getTES(gg))

})
