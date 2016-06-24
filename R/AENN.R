#' All-k Edited Nearest Neighbors
#'
#' Similarity-based filter for removing label noise from a dataset as a
#' preprocessing step of classification. For more information, see 'Details' and
#' 'References' sections.
#'
#' \code{AENN} applies the Edited Nearest Neighbor algorithm \code{\link{ENN}} for all integers between 1 and \code{k}
#' on the whole dataset. At the end, any instance considered noisy by some \emph{ENN} is removed.
#'
#' @param formula A formula describing the classification variable and the attributes to be used.
#' @param data,x Data frame containing the tranining dataset to be filtered.
#' @param k Total number of nearest neighbors to be used.
#' @param classColumn Positive integer indicating the column which contains the
#' (factor of) classes. By default, the last column is considered.
#' @param ... Optional parameters to be passed to other methods.
#' @return An object of class \code{filter}, which is a list with seven components:
#' \itemize{
#'    \item \code{cleanData} is a data frame containing the filtered dataset.
#'    \item \code{remIdx} is a vector of integers indicating the indexes for
#'    removed instances (i.e. their row number with respect to the original data frame).
#'    \item \code{repIdx} is a vector of integers indicating the indexes for
#'    repaired/relabelled instances (i.e. their row number with respect to the original data frame).
#'    \item \code{repLab} contains the new labels for repaired instances.
#'    \item \code{parameters} is a list containing the argument values.
#'    \item \code{call} contains the original call to the filter.
#'    \item \code{extraInf} is a character that includes additional interesting
#'    information not covered by previous items.
#' }
#' @references
#' Tomek I. (1976, June): An Experiment with the Edited Nearest-Neighbor Rule, in
#' \emph{Systems, Man and Cybernetics, IEEE Transactions on}, vol.SMC-6, no.6, pp. 448-452.
#' @examples
#' # Next example is not run in order to save time
#' \dontrun{
#' data(iris)
#' out <- AENN(Species~.-Petal.Length,iris)
#' print(out)
#' identical(out$cleanData, iris[setdiff(1:nrow(iris),out$remIdx),])
#' }
#' @name AENN
NULL

#' @export
AENN <- function(x, ...){
      UseMethod("AENN")
}

#' @rdname AENN
#' @export
AENN.formula <- function(formula,
                        data,
                        ...){
      if(!is.data.frame(data)){
            stop("data argument must be a data.frame")
      }
      modFrame <- model.frame(formula,data) # modFrame is a data.frame built from 'data' using the variables indicated in 'formula'. The first column of 'modFrame' is the response variable, thus we will indicate 'classColumn=1' when calling the HARF.default method in next line.
      attr(modFrame,"terms") <- NULL

      ret <- AENN.default(x = modFrame,...,classColumn = 1)
      ret$call <- match.call(expand.dots = TRUE)
      ret$call[[1]] <- as.name("AENN")
      # Next, we reconstruct the 'cleanData' from the removed and repaired indexes. Otherwise, the 'cleanData' would only contain those columns passed to the default method (for example imagine when running HARF(Species~Petal.Width+Sepal.Length,iris)).
      cleanData <- data
      if(!is.null(ret$repIdx)){
            cleanData[ret$repIdx,which(colnames(cleanData)==colnames(modFrame)[1])] <- ret$repLab  # This is not necessary in HARF because it only removes instances, it does not relabel. However, it must be used when the algorithm relabels instances (in our part there are some of them).
      }
      ret$cleanData <- cleanData[setdiff(1:nrow(cleanData),ret$remIdx),]
      return(ret)
}

#' @rdname AENN
#' @export
AENN.default <- function(x,
                         k=5,
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

      isNoise <- logical(nrow(x))
      for(j in 1:k){
            preds <- sapply(which(!isNoise),function(i){kknn::kknn(formu,x[-i,],x[i,],k=j)$fitted.values})
            isNoise[!isNoise] <- preds!=x[!isNoise,classColumn]
      }

      ##### Building the 'filter' object ###########
      remIdx  <- which(isNoise)
      cleanData <- x[!isNoise,]
      repIdx <- NULL
      repLab <- NULL
      parameters <- list(k=k)
      call <- match.call()
      call[[1]] <- as.name("AENN")

      ret <- list(cleanData = cleanData,
                  remIdx = remIdx,
                  repIdx=repIdx,
                  repLab=repLab,
                  parameters=parameters,
                  call = call,
                  extraInf = NULL)
      class(ret) <- "filter"
      return(ret)
}
