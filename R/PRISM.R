#' PReprocessing Instances that Should be Misclassified
#'
#' Similarity-based filter for removing label noise from a dataset as a
#' preprocessing step of classification. For more information, see 'Details' and
#' 'References' sections.
#'
#' \code{PRISM} identifies \emph{ISMs} (Instances that Should be Misclassified) and removes them from the dataset.
#' In order to do so, it combines five heuristics based on varied approaches by means of a formula.
#' One heuristic relies on class distribution among nearest neighbors, two heuristics are based on the class
#' distribution in a leaf node of a C4.5 tree (either pruned or unpruned), and the other two are based on
#' the class likelihood for an instance, assuming gaussian distribution for continuous variables when necessary.
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
#'
#' @references
#' Smith M. R., Martinez T. (2011, July): Improving classification accuracy by identifying
#' and removing instances that should be misclassified.
#' In \emph{Neural Networks (IJCNN), The 2011 International Joint Conference on} (pp. 2690-2697). IEEE.
#' @examples
#' data(iris)
#' out <- PRISM(Species~., data = iris)
#' print(out)
#' identical(out$cleanData, iris[setdiff(1:nrow(iris),out$remIdx),])
#' @name PRISM
NULL

#' @export
PRISM <- function(x, ...)
{
      UseMethod("PRISM")
}

#' @rdname PRISM
#' @export
PRISM.formula <- function(formula,
                          data,
                          ...)
{
      if(!is.data.frame(data)){
            stop("data argument must be a data.frame")
      }
      modFrame <- model.frame(formula,data) # modFrame is a data.frame built from 'data' using the variables indicated in 'formula'. The first column of 'modFrame' is the response variable, thus we will indicate 'classColumn=1' when calling the HARF.default method in next line.
      attr(modFrame,"terms") <- NULL

      ret <- PRISM.default(x=modFrame,...,classColumn = 1)
      ret$call <- match.call(expand.dots = TRUE)
      ret$call[[1]] <- as.name("PRISM")
      # Next, we reconstruct the 'cleanData' from the removed and repaired indexes. Otherwise, the 'cleanData' would only contain those columns passed to the default method (for example imagine when running HARF(Species~Petal.Width+Sepal.Length,iris)).
      cleanData <- data
      if(!is.null(ret$repIdx)){
            cleanData[ret$repIdx,which(colnames(cleanData)==colnames(modFrame)[1])] <- ret$repLab  # This is not necessary in HARF because it only removes instances, it does not relabel. However, it must be used when the algorithm relabels instances (in our part there are some of them).
      }
      ret$cleanData <- cleanData[setdiff(1:nrow(cleanData),ret$remIdx),]
      return(ret)
}

#' @rdname PRISM
#' @export
PRISM.default <- function(x,
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

      toRemove <- sapply(1:nrow(x),function(i){
            if(cld(i,x,classColumn)>=0){
                  return(FALSE)
            }
            if(dn(i,x,classColumn)>0.8){
                  return(TRUE)
            }
            if(dcp(i,x,classColumn)>=0.5){
                  return(FALSE)
            }
            if(ds(i,x,classColumn)==0){
                  return(TRUE)
            }
            return(FALSE)
      })


      ##### Building the 'filter' object ###########
      remIdx  <- which(toRemove)
      cleanData <- x[setdiff(1:nrow(x),remIdx),]
      repIdx <- NULL
      repLab <- NULL
      call <- match.call()
      call[[1]] <- as.name("PRISM")

      ret <- list(cleanData = cleanData,
                  remIdx = remIdx,
                  repIdx=repIdx,
                  repLab=repLab,
                  parameters=NULL,
                  call = call,
                  extraInf = NULL
      )
      class(ret) <- "filter"
      return(ret)
}
