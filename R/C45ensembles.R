# formu DONE
# library(RWeka) J48
# library(caret) createFolds

#' Classical Filters based on C4.5
#'
#' Ensembled-based filters that use C4.5 classifier to remove label noise from a dataset as a
#' preprocessing step of classification. For more information, see 'Details' and
#' 'References' sections.
#'
#'@details Full description of the methods can be looked up in the provided reference. Notice that C4.5 is used as
#'base classifier instead of TILDE, since a standard attribute-value classification framework
#'is considered (instead of the ILP classification approach of the reference).
#'
#'\code{C45robustFilter} builds a C4.5 decision tree from the training data, and then
#'removes those instances misclassfied by this tree. The process is repeated until no instances are removed.
#'
#'\code{C45votingFilter} splits the dataset into \code{nfolds} folds, building and testing a C4.5 tree on every
#'combination of \code{nfolds}-1 folds. Thus \code{nfolds}-1 votes are gathered
#'for each instance. Removal is carried out by majority or consensus voting schemes.
#'
#'\code{C45iteratedVotingFilter} somehow combines the two previous filter, since
#'it iterates \code{C45votingFilter} until no more noisy instances are removed.
#'
#' @param formula A formula describing the classification variable and the attributes to be used.
#' @param data,x Data frame containing the tranining dataset to be filtered.
#' @param nfolds Number of folds in which the dataset is split.
#' @param consensus Logical. If TRUE, consensus voting scheme is used. If
#' FALSE, majority voting scheme is applied.
#' @param classColumn Positive integer indicating the column which contains the
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
#' @note By means of a message, the number of noisy
#' instances removed is displayed in the console.
#' @references
#' Verbaeten S. (2002, December): Identifying mislabeled training examples in ILP classification
#' problems, in \emph{Proc. 12th Belgian-Dutch Conf. Mach. Learn.}, Utrecht, The Netherlands, pp. 71-78.
#' @examples
#' # Next example is not run in order to save time
#' \dontrun{
#' data(iris)
#' out1 <- C45robustFilter(Species~.-Sepal.Length, iris)
#' # We fix a seed since next two functions create partitions of data for the ensemble
#' set.seed(1)
#' out2 <- C45votingFilter(iris, consensus = TRUE)
#' out3 <- C45iteratedVotingFilter(Species~., iris, nfolds = 5)
#' print(out1)
#' print(out2)
#' print(out3)
#' identical(out1$cleanData,iris[setdiff(1:nrow(iris),out1$remIdx),])
#' identical(out2$cleanData,iris[setdiff(1:nrow(iris),out2$remIdx),])
#' identical(out3$cleanData,iris[setdiff(1:nrow(iris),out3$remIdx),])
#' }
#' @name C45ensembles
#' @aliases C45robustFilter C45votingFilter C45iteratedVotingFilter
NULL

#' @export
C45robustFilter <- function(x, ...)
{
      UseMethod("C45robustFilter")
}

#' @rdname C45ensembles
#' @export
C45robustFilter.formula <- function(formula,
                               data,
                               ...)
{
      if(!is.data.frame(data)){
            stop("data argument must be a data.frame")
      }
      modFrame <- model.frame(formula,data) # modFrame is a data.frame built from 'data' using the variables indicated in 'formula'. The first column of 'modFrame' is the response variable, thus we will indicate 'classColumn=1' when calling the HARF.default method in next line.
      attr(modFrame,"terms") <- NULL

      ret <- C45robustFilter.default(x=modFrame,...,classColumn = 1)
      ret$call <- match.call(expand.dots = TRUE)
      ret$call[[1]] <- as.name("C45robustFilter")
      # Next, we reconstruct the 'cleanData' from the removed and repaired indexes. Otherwise, the 'cleanData' would only contain those columns passed to the default method (for example imagine when running HARF(Species~Petal.Width+Sepal.Length,iris)).
      cleanData <- data
      if(!is.null(ret$repIdx)){
            cleanData[ret$repIdx,which(colnames(cleanData)==colnames(modFrame)[1])] <- ret$repLab  # This is not necessary in HARF because it only removes instances, it does not relabel. However, it must be used when the algorithm relabels instances (in our part there are some of them).
      }
      ret$cleanData <- cleanData[setdiff(1:nrow(cleanData),ret$remIdx),]
      return(ret)
}

#' @rdname C45ensembles
#' @export
C45robustFilter.default <- function(x,
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

      origSize <- nrow(x)
      rownamesOrig <- attr(x,"row.names")
      formu <- as.formula(paste(names(x)[classColumn],"~.",sep = ""))
      row.names(x) <- 1:nrow(x)

      KeepOn <- TRUE
      counter <- 0
      while(KeepOn){
            counter <- counter+1
            pred <- predict(RWeka::J48(formu, x), x)
            IndexesToRemove <- which(pred!=x[,classColumn])
            if(length(IndexesToRemove)==0){
                  KeepOn <- FALSE
            }else{
                  x <- x[-IndexesToRemove,]
            }
            message("Iteration ",counter,": ",length(IndexesToRemove)," instances removed")
      }
      message("Summary: ", origSize - nrow(x), " instances removed in ", counter, " iterations")
      ##### Building the 'filter' object ###########
      remIdx  <- setdiff(1:origSize,as.integer(row.names(x)))
      row.names(x) <- rownamesOrig[as.integer(row.names(x))]
      cleanData <- x
      repIdx <- NULL
      repLab <- NULL
      parameters <- NULL
      call <- match.call()
      call[[1]] <- as.name("C45robustFilter")

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
C45votingFilter <- function(x, ...)
{
      UseMethod("C45votingFilter")
}

#' @rdname C45ensembles
#' @export
C45votingFilter.formula <- function(formula,
                                     data,
                                     ...)
{
      if(!is.data.frame(data)){
            stop("data argument must be a data.frame")
      }
      modFrame <- model.frame(formula,data) # modFrame is a data.frame built from 'data' using the variables indicated in 'formula'. The first column of 'modFrame' is the response variable, thus we will indicate 'classColumn=1' when calling the HARF.default method in next line.
      attr(modFrame,"terms") <- NULL

      ret <- C45votingFilter.default(x=modFrame,...,classColumn = 1)
      ret$call <- match.call(expand.dots = TRUE)
      ret$call[[1]] <- as.name("C45votingFilter")
      # Next, we reconstruct the 'cleanData' from the removed and repaired indexes. Otherwise, the 'cleanData' would only contain those columns passed to the default method (for example imagine when running HARF(Species~Petal.Width+Sepal.Length,iris)).
      cleanData <- data
      if(!is.null(ret$repIdx)){
            cleanData[ret$repIdx,which(colnames(cleanData)==colnames(modFrame)[1])] <- ret$repLab  # This is not necessary in HARF because it only removes instances, it does not relabel. However, it must be used when the algorithm relabels instances (in our part there are some of them).
      }
      ret$cleanData <- cleanData[setdiff(1:nrow(cleanData),ret$remIdx),]
      return(ret)
}

#' @rdname C45ensembles
#' @export
C45votingFilter.default <- function(x,
                                    nfolds=10,
                                    consensus=FALSE,
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

      if(consensus){
            threshold <- nfolds-1
      }else{
            threshold <- floor((nfolds-1)/2)+1
      }

      folds <- caret::createFolds(x[,classColumn],nfolds)
      votes <- vector("integer",nrow(x))
      for(i in 1:nfolds){
            classPreds <- predict(RWeka::J48(formu,x[-folds[[i]],]),newdata=x[-folds[[i]],])
            votes[-folds[[i]]] <- votes[-folds[[i]]]+(classPreds!=x[,classColumn][-folds[[i]]])
      }
      ##### Building the 'filter' object ###########
      cleanData <- x[votes < threshold,]
      remIdx  <- which(votes>=threshold)
      repIdx <- NULL
      repLab <- NULL
      parameters <- list(nfolds=nfolds,
                         consensus=consensus)
      call <- match.call()
      call[[1]] <- as.name("C45votingFilter")

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
C45iteratedVotingFilter <- function(x, ...)
{
      UseMethod("C45iteratedVotingFilter")
}

#' @rdname C45ensembles
#' @export
C45iteratedVotingFilter.formula <- function(formula,
                                             data,
                                             ...)
{
      if(!is.data.frame(data)){
            stop("data argument must be a data.frame")
      }
      modFrame <- model.frame(formula,data) # modFrame is a data.frame built from 'data' using the variables indicated in 'formula'. The first column of 'modFrame' is the response variable, thus we will indicate 'classColumn=1' when calling the HARF.default method in next line.
      attr(modFrame,"terms") <- NULL

      ret <- C45iteratedVotingFilter.default(x=modFrame,...,classColumn = 1)
      ret$call <- match.call(expand.dots = TRUE)
      ret$call[[1]] <- as.name("C45iteratedVotingFilter")
      # Next, we reconstruct the 'cleanData' from the removed and repaired indexes. Otherwise, the 'cleanData' would only contain those columns passed to the default method (for example imagine when running HARF(Species~Petal.Width+Sepal.Length,iris)).
      cleanData <- data
      if(!is.null(ret$repIdx)){
            cleanData[ret$repIdx,which(colnames(cleanData)==colnames(modFrame)[1])] <- ret$repLab  # This is not necessary in HARF because it only removes instances, it does not relabel. However, it must be used when the algorithm relabels instances (in our part there are some of them).
      }
      ret$cleanData <- cleanData[setdiff(1:nrow(cleanData),ret$remIdx),]
      return(ret)
}

#' @rdname C45ensembles
#' @export
C45iteratedVotingFilter.default <- function(x,
                                            nfolds=10,
                                            consensus=FALSE,
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

      origSize <- nrow(x)
      rownamesOrig <- attr(x,"row.names")
      formu <- as.formula(paste(names(x)[classColumn],"~.",sep = ""))
      row.names(x) <- 1:nrow(x)

      if(consensus){
            threshold <- nfolds-1
      }else{
            threshold <- floor((nfolds-1)/2)+1
      }

      KeepOn <- TRUE
      counter <- 0
      while(KeepOn){
            counter <- counter+1
            folds <- caret::createFolds(x[,classColumn],nfolds)
            votes <- vector("integer",nrow(x))
            for(i in 1:nfolds){
                  classPreds <- predict(RWeka::J48(formu,x[-folds[[i]],]),newdata=x[-folds[[i]],])
                  votes[-folds[[i]]] <- votes[-folds[[i]]]+(classPreds!=x[,classColumn][-folds[[i]]])
            }
            IndexesToRemove <- which(votes>=threshold)
            message("Iteration ",counter,": ",length(IndexesToRemove)," noisy instances removed")
            if(length(IndexesToRemove)==0){
                  KeepOn <- FALSE
            }else{
                  x <- x[-IndexesToRemove,]
            }
      }
      message("Summary: ", origSize - nrow(x), " instances removed in ", counter, " iterations")
      ##### Building the 'filter' object ###########
      remIdx  <- setdiff(1:origSize,as.integer(row.names(x)))
      row.names(x) <- rownamesOrig[as.integer(row.names(x))]
      cleanData <- x
      repIdx <- NULL
      repLab <- NULL
      parameters <- list(nfolds=nfolds,
                         consensus=consensus)
      call <- match.call()
      call[[1]] <- as.name("C45iteratedVotingFilter")

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
