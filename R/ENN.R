#' Edited Nearest Neighbors
#'
#' Similarity-based filter for removing label noise from a dataset as a
#' preprocessing step of classification. For more information, see 'Details' and
#' 'References' sections.
#'
#' \code{ENN} finds the \code{k} nearest neighbors for each instance, which is removed
#' if the majority class in this neighborhood is different from its class.
#'
#' @param formula A formula describing the classification variable and the attributes to be used.
#' @param data,x Data frame containing the tranining dataset to be filtered.
#' @param k Number of nearest neighbors to be used.
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
#' Wilson D. L. (1972): Asymptotic properties of nearest neighbor rules using edited data.
#' \emph{Systems, Man and Cybernetics, IEEE Transactions on}, (3), 408-421.
#' @examples
#' data(iris)
#' out <- ENN(Species~., data = iris, k = 5)
#' summary(out)
#' identical(out$cleanData, iris[setdiff(1:nrow(iris),out$remIdx),])
#' @name ENN
NULL

#' @export
ENN <- function(x, ...)
{
      UseMethod("ENN")
}

#' @rdname ENN
#' @export
ENN.formula <- function(formula,
                         data,
                         ...)
{
      if(!is.data.frame(data)){
            stop("data argument must be a data.frame")
      }
      modFrame <- model.frame(formula,data) # modFrame is a data.frame built from 'data' using the variables indicated in 'formula'. The first column of 'modFrame' is the response variable, thus we will indicate 'classColumn=1' when calling the HARF.default method in next line.
      attr(modFrame,"terms") <- NULL

      ret <- ENN.default(x=modFrame,...,classColumn = 1)
      ret$call <- match.call(expand.dots = TRUE)
      ret$call[[1]] <- as.name("ENN")
      # Next, we reconstruct the 'cleanData' from the removed and repaired indexes. Otherwise, the 'cleanData' would only contain those columns passed to the default method (for example imagine when running HARF(Species~Petal.Width+Sepal.Length,iris)).
      cleanData <- data
      if(!is.null(ret$repIdx)){
            cleanData[ret$repIdx,which(colnames(cleanData)==colnames(modFrame)[1])] <- ret$repLab  # This is not necessary in HARF because it only removes instances, it does not relabel. However, it must be used when the algorithm relabels instances (in our part there are some of them).
      }
      ret$cleanData <- cleanData[setdiff(1:nrow(cleanData),ret$remIdx),]
      return(ret)
}

#' @rdname ENN
#' @export
ENN.default <- function(x,
                        k=3,
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

      preds <- sapply(1:nrow(x),function(i){
            kknn::kknn(formula = formu,
                       train = x[-i,],
                       test = x[i,],
                       k = k,
                       kernel = "rectangular")$fitted.values
      })
      ##### Building the 'filter' object ###########
      remIdx  <- which(preds!=x[,classColumn])
      cleanData <- x[setdiff(1:nrow(x),remIdx),]
      repIdx <- NULL
      repLab <- NULL
      parameters <- list(k=k)
      call <- match.call()
      call[[1]] <- as.name("ENN")

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
