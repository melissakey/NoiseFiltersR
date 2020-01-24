#' Blame Based Noise Reduction
#'
#' Similarity-based filter for removing label noise from a dataset as a
#' preprocessing step of classification. For more information, see 'Details' and
#' 'References' sections.
#'
#' \code{BBNR} removes an instance 'X' if: (i) it participates in the misclassification of other instance
#' (i.e. 'X' is among the \code{k} nearest neighbors of a misclassified instance and has a different class);
#' and (ii) its removal does not produce a misclassification in instances that, initially, were correctly
#' classified by 'X' (i.e. 'X' was initially among the \code{k} nearest neighbors and had the same class).
#'
#' @param formula A formula describing the classification variable and the attributes to be used.
#' @param data,x Data frame containing the tranining dataset to be filtered.
#' @param k Number of nearest neighbors to be used.
#' @param classColumn positive integer indicating the column which contains the
#' (factor of) classes. By default, the last column is considered.
#' @param ... Optional parameters to be passed to other methods.
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
#' Delany S. J., Cunningham P. (2004): An analysis of case-base editing in a spam filtering system.
#' In \emph{Advances in Case-Based Reasoning} (pp. 128-141). Springer Berlin Heidelberg.
#' @examples
#' # Next example is not run in order to save time
#' \dontrun{
#' data(iris)
#' out <- BBNR(iris, k = 5)
#' print(out)
#' identical(out$cleanData, iris[setdiff(1:nrow(iris),out$remIdx),])
#' }
#' @name BBNR
NULL

#' @export
BBNR <- function(x, ...)
{
  UseMethod("BBNR")
}

#' @rdname BBNR
#' @export
BBNR.formula <- function(formula,
  data,
  ...)
{
  if(!is.data.frame(data)){
    stop("data argument must be a data.frame")
  }
  modFrame <- model.frame(formula,data) # modFrame is a data.frame built from 'data' using the variables indicated in 'formula'. The first column of 'modFrame' is the response variable, thus we will indicate 'classColumn=1' when calling the HARF.default method in next line.
  attr(modFrame,"terms") <- NULL
  
  ret <- BBNR.default(x=modFrame,...,classColumn = 1)
  ret$call <- match.call(expand.dots = TRUE)
  ret$call[[1]] <- as.name("BBNR")
  # Next, we reconstruct the 'cleanData' from the removed and repaired indexes. Otherwise, the 'cleanData' would only contain those columns passed to the default method (for example imagine when running HARF(Species~Petal.Width+Sepal.Length,iris)).
  cleanData <- data
  if(!is.null(ret$repIdx)){
    cleanData[ret$repIdx,which(colnames(cleanData)==colnames(modFrame)[1])] <- ret$repLab  # This is not necessary in HARF because it only removes instances, it does not relabel. However, it must be used when the algorithm relabels instances (in our part there are some of them).
  }
  ret$cleanData <- cleanData[setdiff(1:nrow(cleanData),ret$remIdx),]
  return(ret)
}

#' @rdname BBNR
#' @export
BBNR.default <- function(x,
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
  n <- nrow(x)

  
  #Instead of using LOOCV in which the model is fit nrow(x) times,
  # We are going to run it once on all the data, and ignore the closet point (which will be the data point itself).
  knn_result <- kknn::kknn(formula = formu,
    train = x,
    test = x,
    k = k + 1,
    kernel = 'rectangular')

  # ignore the first data point -- it's always the point itself  
  knnNeigh <- knn_result$C[,-1]
  
  # calculate the prediction by hand since the model prediction is no longer valid.
  pred_set <- sapply( 1:nrow(x), function(.x) {
    assigned_class <- x[.x,classColumn]
    predicted_classes <- x[knnNeigh[.x, ], classColumn]
    
    predicted_classes == assigned_class
  })
  
  isMisclassified <- sapply( 1:n, function(.x) {
    assigned_class <- x[.x,classColumn]
    predicted_classes <- x[knnNeigh[.x, ], classColumn]
    
    sum(predicted_classes == assigned_class) <= k / 2
  })
  CoverageSets <- lapply(1:n, function(.x) {
    nn_to <- which(.x == knnNeigh)
    
    CS <- nn_to[pred_set[nn_to]] %% n
    CS <- CS[!isMisclassified[CS]]
    if(!length(CS)) return(NULL)
    
    CS
  })
  LiabilitySets <- lapply(which(isMisclassified), function(.x) knnNeigh[.x, !pred_set[, .x]])
  # DissimilaritySets <- lapply(1:n, function(.x) {
  #   nn_to <- which(.x == knnNeigh)
  #   DS <- nn_to[!pred_set[nn_to]] %% n
  #   DS <- DS[isMisclassified[DS]]
  #   if(!length(DS)) return(NULL)
  #   
  #   DS
  # })
  LS_counts <- table(unlist(LiabilitySets))
  toExamine <- as.integer(names(sort(LS_counts, decreasing = TRUE)))
  
  
  
  # this is BBNRv1
  #Examine in order
  toRemove <- sapply(toExamine, function(i) {
    if(!is.null(CoverageSets[[i]])) {
      training <- x[-i, ]
      affectsRemoval <- kknn::kknn(formula = formu,
        train = training,
        test = x[CoverageSets[[i]],],
        k = k,
        kernel = "rectangular")$fitted.values != x[CoverageSets[[i]],classColumn]
      if(all(!affectsRemoval)){
        return(TRUE)
      }
      return(FALSE)
    }
    else{
      return(TRUE)
    }
  })
  
  # toRemove <- rep(NA,length(toExamine))
  # for(j in 1:length(toExamine)){
  #   i <- toExamine[j]
  #   if(!is.null(coverageSets[[i]])){
  #     training <- x[setdiff(1:nrow(x),c(toRemove,i)),]
  #     affectsRemoval <- kknn::kknn(formula = formu,
  #       train = training,
  #       test = x[coverageSets[[i]],],
  #       k = k,
  #       kernel = "rectangular")$fitted.values != x[coverageSets[[i]],classColumn]
  #     if(all(!affectsRemoval)){
  #       toRemove[j] <- i
  #     }
  #   }
  #   else{
  #     toRemove[j] <- i
  #   }
  # }
  
  ##### Building the 'filter' object ###########
  remIdx  <- toExamine[toRemove]
  cleanData <- x[setdiff(1:nrow(x),remIdx),]
  repIdx <- NULL
  repLab <- NULL
  parameters <- list(k=k)
  call <- match.call()
  call[[1]] <- as.name("BBNR")
  
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
