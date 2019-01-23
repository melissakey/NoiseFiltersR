#' Decremental Reduction Optimization Procedures
#'
#' Similarity-based filters for removing label noise from a dataset as a
#' preprocessing step of classification. For more information, see 'Details' and
#' 'References' sections.
#'
#' \code{DROP1} goes over the dataset in the provided order, and removes those
#' instances whose removal does not decrease the accuracy of the 1-NN rule in
#' the remaining dataset.
#'
#' \code{DROP2} introduces two modifications against \code{DROP1}. Regarding the
#' order of processing instances, \code{DROP2} starts with those which are
#' furthest from their nearest "enemy" (two instances are said to be "enemies"
#' if they belong to different classes). Moreover, \code{DROP2} removes an
#' instance if its removal does not decrease the accuracy of the 1-NN rule in
#' the \emph{original} dataset (rather than the \emph{remaining} dataset as in
#' \code{DROP1}).
#'
#' \code{DROP3} is identical to \code{DROP2}, but it includes a preprocessing
#' step to clean the borders between classes. It consists of applying the
#' \code{\link{ENN}} method: any instance misclassified by its k nearest
#' neighbors is removed.
#'
#' @param formula A formula describing the classification variable and the attributes to be used.
#' @param data,x Data frame containing the tranining dataset to be filtered.
#' @param k Number of nearest neighbors to be used.
#' @param classColumn positive integer indicating the column which contains the
#'   (factor of) classes. By default, the last column is considered.
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
#' @references Wilson D. R., Martinez T. R. (2000): Reduction techniques for
#' instance-based learning algorithms. \emph{Machine learning}, 38(3), 257-286.
#' Wilson D. R., Martinez T. R. (1997, July): Instance pruning techniques. In
#' \emph{ICML} (Vol. 97, pp. 403-411).
#' @examples
#' # Next example is not run in order to save time
#' \dontrun{
#' data(iris)
#' trainData <- iris[c(1:20,51:70,101:120),]
#' out1 <- DROP1(Species~ Petal.Length + Petal.Width, data = trainData)
#' summary(out1, explicit = TRUE)
#' identical(out1$cleanData, trainData[setdiff(1:nrow(trainData),out1$remIdx),])
#' }
#' @name DROP
#' @aliases DROP1 DROP2 DROP3
NULL

#' @export
DROP1 <- function(x, ...)
{
      UseMethod("DROP1")
}

#' @rdname DROP
#' @export
DROP1.formula <- function(formula,
                          data,
                          ...){
      if(!is.data.frame(data)){
            stop("data argument must be a data.frame")
      }
      modFrame <- model.frame(formula,data) # modFrame is a data.frame built from 'data' using the variables indicated in 'formula'. The first column of 'modFrame' is the response variable, thus we will indicate 'classColumn=1' when calling the HARF.default method in next line.
      attr(modFrame,"terms") <- NULL

      ret <- DROP1.default(x=modFrame,...,classColumn = 1)
      ret$call <- match.call(expand.dots = TRUE)
      ret$call[[1]] <- as.name("DROP1")
      # Next, we reconstruct the 'cleanData' from the removed and repaired indexes. Otherwise, the 'cleanData' would only contain those columns passed to the default method (for example imagine when running HARF(Species~Petal.Width+Sepal.Length,iris)).
      cleanData <- data
      if(!is.null(ret$repIdx)){
            cleanData[ret$repIdx,which(colnames(cleanData)==colnames(modFrame)[1])] <- ret$repLab  # This is not necessary in HARF because it only removes instances, it does not relabel. However, it must be used when the algorithm relabels instances (in our part there are some of them).
      }
      ret$cleanData <- cleanData[setdiff(1:nrow(cleanData),ret$remIdx),]
      return(ret)
}

#' @rdname DROP
#' @export
DROP1.default <- function(x,
                          k = 1,
                          classColumn=ncol(x),
                          ...){
      if(!is.data.frame(x)){
            stop("data argument must be a data.frame")
      }
      if(!classColumn%in%(1:ncol(x))){
            stop("class column out of range")
      }
      if(!is.factor(x[,classColumn])){
            stop("class column of data must be a factor")
      }

      formu <- as.formula(paste(names(x)[classColumn],"~.",sep = ""))

      preds <- sapply(1:nrow(x),function(i){
            kknn::kknn(formula = formu, train = x[-i,], test = x[i,], k = k, kernel = "rectangular")$fitted.values
      })
      currentAcc <- sum(preds==x[,classColumn])

      toRemove <- logical(nrow(x))
      for(j in 1:nrow(x)){
            predsIn <- sapply(setdiff(1:nrow(x),which(toRemove)),function(i){
                  kknn::kknn(formula = formu, train = x[-c(i,j,which(toRemove)),], test = x[i,], k = k, kernel = "rectangular")$fitted.values
            })
            newAcc <- sum(predsIn==x[setdiff(1:nrow(x),which(toRemove)),classColumn])
            if(newAcc>=currentAcc){
                  currentAcc <- newAcc+ifelse(predsIn[j-length(toRemove)]==x[j,classColumn],-1,0)
                  toRemove[j] <- TRUE
            }
      }

      ##### Building the 'filter' object ###########
      remIdx  <- which(toRemove)
      cleanData <- x[setdiff(1:nrow(x),remIdx),]
      repIdx <- NULL
      repLab <- NULL
      parameters <- list(k=k)
      call <- match.call()
      call[[1]] <- as.name("DROP1")

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


#' @export
DROP2 <- function(x, ...)
{
      UseMethod("DROP2")
}

#' @rdname DROP
#' @export
DROP2.formula <- function(formula,
                          data,
                          ...)
{
      if(!is.data.frame(data)){
            stop("data argument must be a data.frame")
      }
      modFrame <- model.frame(formula,data) # modFrame is a data.frame built from 'data' using the variables indicated in 'formula'. The first column of 'modFrame' is the response variable, thus we will indicate 'classColumn=1' when calling the HARF.default method in next line.
      attr(modFrame,"terms") <- NULL

      ret <- DROP2.default(x=modFrame,...,classColumn = 1)
      ret$call <- match.call(expand.dots = TRUE)
      ret$call[[1]] <- as.name("DROP2")
      # Next, we reconstruct the 'cleanData' from the removed and repaired indexes. Otherwise, the 'cleanData' would only contain those columns passed to the default method (for example imagine when running HARF(Species~Petal.Width+Sepal.Length,iris)).
      cleanData <- data
      if(!is.null(ret$repIdx)){
            cleanData[ret$repIdx,which(colnames(cleanData)==colnames(modFrame)[1])] <- ret$repLab  # This is not necessary in HARF because it only removes instances, it does not relabel. However, it must be used when the algorithm relabels instances (in our part there are some of them).
      }
      ret$cleanData <- cleanData[setdiff(1:nrow(cleanData),ret$remIdx),]
      return(ret)
}

#' @rdname DROP
#' @export
DROP2.default <- function(x,
                          k = 1,
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

      formu <- as.formula(paste(names(x)[classColumn],"~.",sep = ""))

      distsEnemy <- sapply(1:nrow(x),function(i){distEnemy(i,x,classColumn)})
      removalOrder <- order(distsEnemy, decreasing = TRUE)

      preds <- sapply(1:nrow(x),function(i){
            kknn::kknn(formula = formu, train = x[-i,], test = x[i,], k = k, kernel = "rectangular")$fitted.values
      })
      currentAcc <- sum(preds==x[,classColumn])

      toRemove <- logical(length(removalOrder))
      for(j in removalOrder){
            predsIn <- sapply(1:nrow(x),function(i){
                  kknn::kknn(formula = formu, train = x[-c(i,j,which(toRemove)),], test = x[i,], k = k, kernel = "rectangular")$fitted.values
            })
            newAcc <- sum(predsIn==x[,classColumn])
            if(newAcc>=currentAcc){
                  currentAcc <- newAcc
                  toRemove[j] <- TRUE
            }
      }

      ##### Building the 'filter' object ###########
      remIdx  <- which(toRemove)
      cleanData <- x[setdiff(1:nrow(x),remIdx),]
      repIdx <- NULL
      repLab <- NULL
      parameters <- list(k=k)
      call <- match.call()
      call[[1]] <- as.name("DROP2")

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

#' @export
DROP3 <- function(x, ...)
{
      UseMethod("DROP3")
}

#' @rdname DROP
#' @export
DROP3.formula <- function(formula,
                          data,
                          ...)
{
      if(!is.data.frame(data)){
            stop("data argument must be a data.frame")
      }
      modFrame <- model.frame(formula,data) # modFrame is a data.frame built from 'data' using the variables indicated in 'formula'. The first column of 'modFrame' is the response variable, thus we will indicate 'classColumn=1' when calling the HARF.default method in next line.
      attr(modFrame,"terms") <- NULL

      ret <- DROP3.default(x=modFrame,...,classColumn = 1)
      ret$call <- match.call(expand.dots = TRUE)
      ret$call[[1]] <- as.name("DROP3")
      # Next, we reconstruct the 'cleanData' from the removed and repaired indexes. Otherwise, the 'cleanData' would only contain those columns passed to the default method (for example imagine when running HARF(Species~Petal.Width+Sepal.Length,iris)).
      cleanData <- data
      if(!is.null(ret$repIdx)){
            cleanData[ret$repIdx,which(colnames(cleanData)==colnames(modFrame)[1])] <- ret$repLab  # This is not necessary in HARF because it only removes instances, it does not relabel. However, it must be used when the algorithm relabels instances (in our part there are some of them).
      }
      ret$cleanData <- cleanData[setdiff(1:nrow(cleanData),ret$remIdx),]
      return(ret)
}

#' @rdname DROP
#' @export
DROP3.default <- function(x,
                          k = 1,
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

      formu <- as.formula(paste(names(x)[classColumn],"~.",sep = ""))

      initialPreds <- sapply(1:nrow(x),function(i){
            kknn::kknn(formula = formu, train = x[-i,], test = x[i,], k = k, kernel = "rectangular")$fitted.values
      })
      initiallyKept <- which(initialPreds==x[,classColumn])
      remIdx <- setdiff(1:nrow(x),initiallyKept)
      x <- x[initiallyKept,]



      distsEnemy <- sapply(1:nrow(x),function(i){distEnemy(i,x,classColumn)})
      removalOrder <- order(distsEnemy,decreasing = TRUE)

      preds <- sapply(1:nrow(x),function(i){
            kknn::kknn(formula = formu, train = x[-i,], test = x[i,], k = k, kernel = "rectangular")$fitted.values
      })
      currentAcc <- sum(preds==x[,classColumn])

      toRemove <- logical(removalOrder)
      for(j in removalOrder){
            predsIn <- sapply(1:nrow(x),function(i){
                  kknn::kknn(formula = formu, train = x[-c(i,j,which(toRemove)),], test = x[i,], k = k, kernel = "rectangular")$fitted.values
            })
            newAcc <- sum(predsIn==x[,classColumn])
            if(newAcc>=currentAcc){
                  currentAcc <- newAcc
                  toRemove[j] <- TRUE
            }
      }

      ##### Building the 'filter' object ###########
      remIdx  <- c(remIdx,initiallyKept[which(toRemove)])
      cleanData <- x[setdiff(1:nrow(x),which(toRemove)),]
      repIdx <- NULL
      repLab <- NULL
      parameters <- list(k=k)
      call <- match.call()
      call[[1]] <- as.name("DROP3")

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

