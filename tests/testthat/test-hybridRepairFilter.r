context("Tests for hybridRepairFilter")

data(iris)
trainData <- iris[c(1:20,51:70,101:120),]
out <- hybridRepairFilter(Species~., data = trainData, noiseAction = "hybrid")

test_that("the result has the correct class", {
      expect_is(out,"filter")
})

test_that("the new labels are acceptable levels for the class variable", {
      expect_true(length(out$repLab)==0 || all(as.character(out$repLab) %in% levels(iris$Species)))
})

test_that("the clean dataset can be correctly reconstructed from the filter object", {
      newData <- trainData
      newData[out$repIdx,5] <- out$repLab
      cleanData <- newData[setdiff(1:nrow(trainData),out$remIdx),]
      expect_equal(cleanData,out$cleanData)
})
