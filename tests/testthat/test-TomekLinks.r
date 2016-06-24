context("Tests for TomekLinks filter")

test_that("error is thrown when there are more than two classes", {
      expect_error(TomekLinks(Species~., data =iris))
})

test_that("the filter works when there are two classes", {
      skip_on_cran()
      twoClassIris <- iris[1:100,]
      twoClassIris$Species <- factor(twoClassIris$Species)
      out <- TomekLinks(twoClassIris)
      expect_is(out,"filter")
      expect_equal(twoClassIris[setdiff(1:100,out$remIdx),],out$cleanData)
})
