#formu DONE
#library(randomForest)
#library(caret)

#'@title High Agreement Random Forest
#'
#'@description Ensemble-based filter for removing label noise from a dataset as a
#' preprocessing step of classification. For more information, see 'Details' and
#' 'References' sections.
#'
#'@param formula A formula describing the classification variable and the attributes to be used.
#'@param data,x Data frame containing the tranining dataset to be filtered.
#'@param nfolds Number of folds for the cross voting scheme.
#'@param agreementLevel Real number between 0.5 and 1. An instance is identified as
#'  noise when the classification confidences provided by the random forest to the
#'  classes that are not the actual class of the instance add up at least
#'  \code{agreementLevel}. Authors obtain the best performance in (Sluban et al., 2010)
#'  when setting it between 0.7 and 0.8.
#'@param ntrees Number of trees for the random forest.
#'@param classColumn Positive integer indicating the column which contains the (factor
#'  of) classes. By default, the last column is considered.
#'@param ... Optional parameters to be passed to other methods.
#'
#'@details Making use of a \code{nfolds}-folds cross validation scheme, instances are
#'identified as noise and removed when a random forest provides little confidence for
#'the actual instance's label (namely, less than 1-\code{agreementLevel}). The value of
#'\code{agreementLevel} allows to tune the precision and recall of the filter, getting
#'the best trade-off when moving between 0.7 and 0.8 (Sluban et al., 2010).
#'
#'@return An object of class \code{filter}, which is a list with seven components:
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
#'@references Sluban B., Gamberger D., Lavrac N. (2010, August): Advances in Class
#'Noise Detection. In \emph{ECAI} (pp. 1105-1106).
#'@examples
#'# Next example is not run in order to save time
#'\dontrun{
#'data(iris)
#'# We fix a seed since there exists a random partition for the ensemble
#'set.seed(1)
#'out <- HARF(Species~., data = iris, ntrees = 100)
#'print(out)
#'identical(out$cleanData, iris[setdiff(1:nrow(iris),out$remIdx),])
#'}
#'@name HARF
NULL


#' @export
HARF <- function (x, ...){
      UseMethod ("HARF")
}

#' @rdname HARF
#' @export
HARF.formula <- function(formula,
                         data,
                         ...)
{
      if (!is.data.frame(data)){
            stop("data argument must be a data.frame")
      }
      modFrame <- model.frame(formula,data) # modFrame is a data.frame built from 'data' using the variables indicated in 'formula'. The first column of 'modFrame' is the response variable, thus we will indicate 'classColumn=1' when calling the HARF.default method in next line.
      attr(modFrame,"terms") <- NULL

      ret <- HARF.default(x=modFrame,...,classColumn = 1)
      ret$call <- match.call(expand.dots = TRUE)
      ret$call[[1]] <- as.name("HARF")
      # Next, we reconstruct the 'cleanData' from the removed and repaired indexes. Otherwise, the 'cleanData' would only contain those columns passed to the default method (for example imagine when running HARF(Species~Petal.Width+Sepal.Length,iris)).
      cleanData <- data
      if(!is.null(ret$repIdx)){
            cleanData[ret$repIdx,which(colnames(cleanData)==colnames(modFrame)[1])] <- ret$repLab  # This is not necessary in HARF because it only removes instances, it does not relabel. However, it must be used when the algorithm relabels instances (in our part there are some of them).
      }
      ret$cleanData <- cleanData[setdiff(1:nrow(cleanData),ret$remIdx),]
      return(ret)
}

#' @rdname HARF
#' @export
HARF.default <- function(x,
                         nfolds = 10,
                         agreementLevel = 0.7,
                         ntrees = 500,
                         classColumn = ncol(x),
                         ...)
{
      if(agreementLevel<0.5 || agreementLevel>1){
            stop("the agreement level must range between 0.5 and 1")
      }
      if(!is.data.frame(x)){
            stop("data argument must be a data.frame")
      }
      if(!classColumn%in%(1:ncol(x))){
            stop("class column out of range")
      }
      if(!is.factor(x[,classColumn])){
            stop("class column of data must be a factor")
      }

      classes <- x[,classColumn]
      folds <- caret::createFolds(classes,nfolds)

      isNoise <- logical(nrow(x))
      for(i in 1:nfolds){
            probs <- predict(randomForest::randomForest(x=x[-folds[[i]],-classColumn],y=x[-folds[[i]],classColumn],ntrees=ntrees),
                             newdata = x[folds[[i]],],
                             type = "prob")
            isNoise[folds[[i]]] <- sapply(1:length(folds[[i]]),function(j){probs[j,classes[folds[[i]][j]]] <= 1-agreementLevel})
      }
      ##### Building the 'filter' object ###########
      cleanData <- x[!isNoise,]
      remIdx  <- which(isNoise)
      repIdx <- NULL
      repLab <- NULL
      parameters <- list(nfolds=nfolds,
                         ntrees=ntrees,
                         agreementLevel=agreementLevel)
      call <- match.call()
      call[[1]] <- as.name("HARF")

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
