# formu DONE
# library(caret) #to use 'createFolds'
# library(RWeka) #for performing C4.5 and 1-NN
# library(MASS) #for performing LDA

#' Ensemble Filter
#'
#' Ensemble-based filter for removing label noise from a dataset as a
#' preprocessing step of classification. For more information, see 'Details' and
#' 'References' sections.
#'
#'Full description of the method can be looked up in the provided references.
#'Dataset is split in \code{nfolds} folds, an ensemble of three different base classifiers (C4.5, 1-KNN, LDA) is
#'built over every combination of \code{nfolds}-1 folds, and then tested on the other one. Finally, consensus
#'or majority voting scheme is applied to remove noisy instances.
#'
#' @param formula A formula describing the classification variable and the attributes to be used.
#' @param data,x data frame containing the tranining dataset to be filtered.
#' @param nfolds number of folds in which the dataset is split.
#' @param consensus logical. If TRUE, consensus voting scheme is used. If
#' FALSE, majority voting scheme is applied.
#' @param classColumn positive integer indicating the column which contains the
#' (factor of) classes. By default, the last column is considered.
#' @param ... Optional parameters to be passed to other methods.
#'
#' @return An object of class \code{filter}, which is a list with seven components:
#' \itemize{
#'    \item \code{cleanData} is a data frame containing the filtered dataset.
#'    \item \code{remIdx} is a vector of integers indicating the indexes for
#'    removed instances (i.e. their row number with respect to the original data frame).
#'    \item \code{repIdx} is a vector of integers indicating the indexes for
#'    repaired/relabelled instances (i.e. their row number with respect to the original data frame).
#'    \item \code{repLab} is a factor containing the new labels for repaired instances.
#'    \item \code{parameters} is a list containing the argument values.
#'    \item \code{call} contains the original call to the filter.
#'    \item \code{extraInf} is a character that includes additional interesting
#'    information not covered by previous items.
#' }
#'
#' @references Brodley C. E., Friedl M. A. (1996, May): Improving automated
#' land cover mapping by identifying and eliminating mislabeled observations
#' from training data. In \emph{Geoscience and Remote Sensing Symposium, 1996.
#' IGARSS'96.'Remote Sensing for a Sustainable Future.', International} (Vol.
#' 2, pp. 1379-1381). IEEE.
#'
#' Brodley C. E., Friedl M. A. (1996, August): Identifying and eliminating
#' mislabeled training instances. In \emph{AAAI/IAAI}, Vol. 1 (pp. 799-805).
#'
#' Brodley C. E., Friedl M. A. (1999): Identifying mislabeled training data.
#' \emph{Journal of Artificial Intelligence Research}, 131-167.
#' @examples
#' data(iris)
#' # We fix a seed since there exists a random partition for the ensemble
#' set.seed(1)
#' out <- EF(Species~., data = iris, consensus = FALSE)
#' summary(out, explicit = TRUE)
#' identical(out$cleanData, iris[setdiff(1:nrow(iris),out$remIdx),])
#' @name EF
#' @importFrom stats predict
#' @importFrom stats model.frame
#' @importFrom stats as.formula
#'
NULL

#' @export
EF <- function(x, ...)
{
      UseMethod("EF")
}

#' @rdname EF
#' @export
EF.formula <- function(formula,
                         data,
                         ...)
{
      if(!is.data.frame(data)){
            stop("data argument must be a data.frame")
      }
      modFrame <- model.frame(formula,data) # modFrame is a data.frame built from 'data' using the variables indicated in 'formula'. The first column of 'modFrame' is the response variable, thus we will indicate 'classColumn=1' when calling the HARF.default method in next line.
      attr(modFrame,"terms") <- NULL

      ret <- EF.default(x=modFrame,...,classColumn = 1)
      ret$call <- match.call(expand.dots = TRUE)
      ret$call[[1]] <- as.name("EF")
      # Next, we reconstruct the 'cleanData' from the removed and repaired indexes. Otherwise, the 'cleanData' would only contain those columns passed to the default method (for example imagine when running HARF(Species~Petal.Width+Sepal.Length,iris)).
      cleanData <- data
      if(!is.null(ret$repIdx)){
            cleanData[ret$repIdx,which(colnames(cleanData)==colnames(modFrame)[1])] <- ret$repLab  # This is not necessary in HARF because it only removes instances, it does not relabel. However, it must be used when the algorithm relabels instances (in our part there are some of them).
      }
      ret$cleanData <- cleanData[setdiff(1:nrow(cleanData),ret$remIdx),]
      return(ret)
}

#' @rdname EF
#' @export
EF.default <- function(x,
                       nfolds=4,
                       consensus=TRUE,
                       classColumn=ncol(x),
                       ...)
{
      if(!is.data.frame(x)){
            stop("data argument must be a data.frame")
      }
      if(!classColumn%in%(1:ncol(x))){
            stop("class column out of range")
      }
      if(!is.factor(x[,classColumn])){
            stop("class column of data must be a factor")
      }

      if(consensus){
            threshold <- 3
      }else{
            threshold <- 2
      }

      formu <- as.formula(paste(names(x)[classColumn],"~.",sep = ""))

      folds <- caret::createFolds(x[,classColumn],nfolds)

      votes <- vector("integer",nrow(x))
      for(i in 1:nfolds){
            #C4.5 algorithm predictions and votes
            votes[folds[[i]]] <- votes[folds[[i]]]+(predict(RWeka::J48(formu,x[-folds[[i]],]),newdata=x[folds[[i]],])!=x[,classColumn][folds[[i]]])
            #1-NN algorithm predictions and votes
            votes[folds[[i]]] <- votes[folds[[i]]]+(predict(RWeka::IBk(formu,x[-folds[[i]],]),newdata=x[folds[[i]],])!=x[,classColumn][folds[[i]]])
            #LDA algorithm predictions and votes
            votes[folds[[i]]] <- votes[folds[[i]]]+(predict(MASS::lda(x=x[-folds[[i]],-classColumn],grouping=x[-folds[[i]],classColumn]),newdata=x[folds[[i]],-classColumn])$class!=x[,classColumn][folds[[i]]])
      }
      ##### Building the 'filter' object ###########
      cleanData <- x[votes < threshold,]
      remIdx  <- which(votes>=threshold)
      repIdx <- NULL
      repLab <- NULL
      parameters <- list(nfolds=nfolds,
                         consensus=consensus)
      call <- match.call()
      call[[1]] <- as.name("EF")

      ret <- list(cleanData = cleanData,
                  remIdx = remIdx,
                  repIdx=repIdx,
                  repLab=repLab,
                  parameters=parameters,
                  call = call,
                  extraInf = NULL
      )
      class(ret) <- "filter"
      return(ret)
}
