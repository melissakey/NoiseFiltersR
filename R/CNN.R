#' Condensed Nearest Neighbors
#'
#' Similarity-based method designed to select the most relevant instances for
#' subsequent classification with a \emph{nearest neighbor} rule. For more
#' information, see 'Details' and 'References' sections.
#'
#' \code{CNN} searches for a 'consistent subset' of the provided dataset, i.e. a subset that is enough for
#' correctly classifying the rest of instances by means of 1-NN. To do so, \code{CNN} stores the first instance and
#' goes for a first sweep over the dataset, adding to the stored bag those instances which are not correctly classified by 1-NN taking the stored bag as training set.
#' Then, the process is iterated until all non-stored instances are correctly classified.
#'
#'Although \code{CNN} is not strictly a label noise filter, it is included here for completeness, since
#' the origins of noise filters are connected with instance selection algorithms.
#'
#' @param formula A formula describing the classification variable and the attributes to be used.
#' @param data,x Data frame containing the tranining dataset to be filtered.
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
#' @references
#' Hart P. (May, 1968): The condensed nearest neighbor rule, \emph{IEEE Trans. Inf.
#' Theory}, vol. 14, no. 5, pp. 515-516.
#' @seealso \code{\link{RNN}}
#' @examples
#' # Next example is not run in order to save time
#' \dontrun{
#' data(iris)
#' out <- CNN(iris)
#' print(out)
#' length(out$remIdx)
#' identical(out$cleanData, iris[setdiff(1:nrow(iris),out$remIdx),])
#' }
#' @name CNN
NULL

#' @export
CNN <- function(x, ...)
{
      UseMethod("CNN")
}

#' @rdname CNN
#' @export
CNN.formula <- function(formula,
                        data,
                        ...)
{
      if(!is.data.frame(data)){
            stop("data argument must be a data.frame")
      }
      modFrame <- model.frame(formula,data) # modFrame is a data.frame built from 'data' using the variables indicated in 'formula'. The first column of 'modFrame' is the response variable, thus we will indicate 'classColumn=1' when calling the HARF.default method in next line.
      attr(modFrame,"terms") <- NULL

      ret <- CNN.default(x=modFrame,...,classColumn = 1)
      ret$call <- match.call(expand.dots = TRUE)
      ret$call[[1]] <- as.name("CNN")
      # Next, we reconstruct the 'cleanData' from the removed and repaired indexes. Otherwise, the 'cleanData' would only contain those columns passed to the default method (for example imagine when running HARF(Species~Petal.Width+Sepal.Length,iris)).
      cleanData <- data
      if(!is.null(ret$repIdx)){
            cleanData[ret$repIdx,which(colnames(cleanData)==colnames(modFrame)[1])] <- ret$repLab  # This is not necessary in HARF because it only removes instances, it does not relabel. However, it must be used when the algorithm relabels instances (in our part there are some of them).
      }
      ret$cleanData <- cleanData[setdiff(1:nrow(cleanData),ret$remIdx),]
      return(ret)
}

#' @rdname CNN
#' @export
CNN.default <- function(x,
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

      firstDif <- which(x[,classColumn]!=x[1,classColumn])[1]
      store <- c(1,firstDif)
      grabBag <- setdiff(1:firstDif,store)
      for(i in (firstDif+1):nrow(x)){
            if(kknn::kknn(formula = formu,
                          train = x[store,],
                          test = x[i,],k = 1)$fitted.values==x[i,classColumn]){
                  grabBag <- c(grabBag,i)
            }else{
                  store <- c(store,i)
            }
      }

      KeepOn <- TRUE
      while(KeepOn){
            KeepOn <- FALSE
            for(i in grabBag){
                  if(kknn::kknn(formula = formu,
                                train = x[store,],
                                test = x[i,],k=1)$fitted.values!=x[i,classColumn]){
                        store <- c(store,i)
                        grabBag <- setdiff(grabBag,i)
                        KeepOn <- TRUE
                  }
            }
      }

      ##### Building the 'filter' object ###########
      remIdx  <- sort(grabBag)
      cleanData <- x[sort(store),]
      repIdx <- NULL
      repLab <- NULL
      parameters <- NULL
      call <- match.call()
      call[[1]] <- as.name("CNN")

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
