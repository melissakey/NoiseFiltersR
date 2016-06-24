context("Tests for GE filter")

data(iris)
out <- GE(Species~., data = iris)

test_that("the result has the correct class", {
      expect_is(out,"filter")
})

test_that("the new labels are acceptable levels for the class variable", {
      expect_true(length(out$repLab)==0 || all(as.character(out$repLab) %in% levels(iris$Species)))
})

test_that("the clean dataset can be correctly reconstructed from the filter object", {
      newData <- iris
      newData[out$repIdx,5] <- out$repLab
      cleanData <- newData[setdiff(1:nrow(iris),out$remIdx),]
      expect_equal(cleanData,out$cleanData)
})
