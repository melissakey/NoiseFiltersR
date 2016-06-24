#' Generalized Edition
#'
#' Similarity-based filter for removing or repairing label noise from a dataset as a
#' preprocessing step of classification. For more information, see 'Details' and
#' 'References' sections.
#'
#' \code{GE} is a generalization of \code{\link{ENN}} that integrates the possibility of 'repairing'
#' or 'relabeling' instances rather than only 'removing'. For each instance, \code{GE} considers
#' its \code{k-1} neighbors and the instance itself. If there are at least \code{kk} examples from the same class,
#' the instance is relabeled with that class (which could be its own). Otherwise, it is removed.
#'
#' @param formula A formula describing the classification variable and the attributes to be used.
#' @param data,x Data frame containing the tranining dataset to be filtered.
#' @param k Number of nearest neighbors to be considered.
#' @param kk Minimum size for local majority class in order to relabel an instance.
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
#' Koplowitz J., Brown T. A. (1981): On the relation of performance to editing
#' in nearest neighbor rules. \emph{Pattern Recognition}, 13(3), 251-255.
#' @examples
#' # Next example is not run in order to save time
#' \dontrun{
#' data(iris)
#' out <- GE(iris)
#' summary(out, explicit = TRUE)
#' # We check that the process was correct
#' irisCopy <- iris
#' irisCopy[out$repIdx,5] <- out$repLab
#' cleanData <- irisCopy[setdiff(1:nrow(iris),out$remIdx),]
#' identical(out$cleanData,cleanData)
#' }
#' @name GE
NULL

#' @export
GE <- function(x, ...)
{
      UseMethod("GE")
}

#' @rdname GE
#' @export
GE.formula <- function(formula,
                        data,
                        ...)
{
      if(!is.data.frame(data)){
            stop("data argument must be a data.frame")
      }
      modFrame <- model.frame(formula,data) # modFrame is a data.frame built from 'data' using the variables indicated in 'formula'. The first column of 'modFrame' is the response variable, thus we will indicate 'classColumn=1' when calling the HARF.default method in next line.
      attr(modFrame,"terms") <- NULL

      ret <- GE.default(x=modFrame,...,classColumn = 1)
      ret$call <- match.call(expand.dots = TRUE)
      ret$call[[1]] <- as.name("GE")
      # Next, we reconstruct the 'cleanData' from the removed and repaired indexes. Otherwise, the 'cleanData' would only contain those columns passed to the default method (for example imagine when running HARF(Species~Petal.Width+Sepal.Length,iris)).
      cleanData <- data
      if(!is.null(ret$repIdx)){
            cleanData[ret$repIdx,which(colnames(cleanData)==colnames(modFrame)[1])] <- ret$repLab  # This is not necessary in HARF because it only removes instances, it does not relabel. However, it must be used when the algorithm relabels instances (in our part there are some of them).
      }
      ret$cleanData <- cleanData[setdiff(1:nrow(cleanData),ret$remIdx),]
      return(ret)
}

#' @rdname GE
#' @export
GE.default <- function(x,
                       k=5,
                       kk=ceiling(k/2),
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

      localClasses <- kknn::kknn(formula = formu,
                                 train = x,
                                 test = x,
                                 k = k,
                                 kernel = "rectangular")$CL

      toRepair <- apply(localClasses,1,function(v){
            w <- table(v)
            if(max(w)>=kk){
                  TRUE
                  }
            else{
                  FALSE
                  }
            })

      majClass <- apply(localClasses,1,function(v){
            w <- table(v)
            names(w)[nnet::which.is.max(w)]
      })

      ##### Building the 'filter' object ###########
      remIdx  <- which(!toRepair)
      repIdx <- which(toRepair & majClass!=x[,classColumn])
      repLab <- factor(majClass[repIdx], levels = levels(x[,classColumn]))
      cleanData <- x
      cleanData[repIdx,classColumn] <- repLab
      cleanData <- cleanData[setdiff(1:nrow(cleanData),remIdx),]

      parameters <- list(k=k,
                         kk=kk)
      call <- match.call()
      call[[1]] <- as.name("GE")

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
