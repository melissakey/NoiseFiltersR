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

      knnInf <- sapply(1:nrow(x),function(i){
            model <- kknn::kknn(formula = formu,
                                train = x[-i,],
                                test = x[i,],
                                k = k,
                                kernel = "rectangular")
            isMisclassified <- model$fitted.values == x[i,classColumn]
            nearestNeigh <- setdiff(1:nrow(x),i)[model$C]
            c(isMisclassified,nearestNeigh)
      })

      isMisclassified <- as.logical(knnInf[1,])
      knnNeigh <- knnInf[-1,]

      classifiesWellOther <- lapply(which(!isMisclassified),function(i){
            knnNeigh[x[knnNeigh[,i],classColumn]==x[i,classColumn],i]
      })
      coverageSets <- list()
      length(coverageSets) <- nrow(x)
      for(i in length(classifiesWellOther)){
            items <- classifiesWellOther[[i]]
            for(j in items){
                  coverageSets[[j]] <- c(coverageSets[[j]],i)
            }
      }

      misclassifyOther <- lapply(which(isMisclassified),function(i){
            knnNeigh[x[knnNeigh[,i],classColumn]!=x[i,classColumn],i]
      })
      toExamine <- unlist(misclassifyOther)
      counts <- table(toExamine)
      toExamine <- as.integer(names(counts)[order(counts,decreasing = TRUE)])

      #Examine in order
      toRemove <- integer(0)
      for(i in toExamine){
            if(!is.null(coverageSets[[i]])){
                  training <- x[setdiff(1:nrow(x),c(toRemove,i)),]
                  affectsRemoval <- kknn::kknn(formula = formu,
                                               train = training,
                                               test = x[coverageSets[[i]],],
                                               k = k,
                                               kernel = "rectangular")$fitted.values != x[coverageSets[[i]],classColumn]
                  if(all(!affectsRemoval)){
                        toRemove <- c(toRemove,i)
                  }
            }
            else{
                  toRemove <- c(toRemove,i)
            }
      }

      ##### Building the 'filter' object ###########
      remIdx  <- toRemove
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
