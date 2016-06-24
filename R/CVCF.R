# formu DONE
# J48 from RWeka
# createFolds from caret

#' Cross-Validated Committees Filter
#'
#' Ensemble-based filter for removing label noise from a dataset as a
#' preprocessing step of classification. For more information, see 'Details' and
#' 'References' sections.
#'
#'Full description of the method can be looked up in the provided references.
#'Dataset is split in \code{nfolds} folds, a base classifiers (C4.5 in this implementation) is
#'built over every combination of \code{nfolds}-1 folds, and then tested on the whole dataset. Finally, consensus
#'or majority voting scheme is applied to remove noisy instances.
#'
#' @param formula A formula describing the classification variable and the attributes to be used.
#' @param data,x data frame containing the tranining dataset to be filtered.
#' @param nfolds number of folds in which the dataset is split.
#' @param consensus logical. If TRUE, consensus voting scheme is used. If
#' FALSE, majority voting scheme is applied.
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
#' Verbaeten S., Van Assche A. (2003, June): Ensemble methods for noise elimination in classification
#' problems. \emph{Proc. 4th Int. Conf. Multiple Classifier Syst.}, Guildford, U.K., pp. 317-325.
#' @examples
#' # Next example is not run in order to save time
#' \dontrun{
#' data(iris)
#' # We fix a seed since there exists a random partition for the ensemble
#' set.seed(1)
#' out <- CVCF(Species~.-Sepal.Width, data = iris)
#' print(out)
#' identical(out$cleanData, iris[setdiff(1:nrow(iris),out$remIdx),])
#' }
#' @name CVCF
NULL

#' @export
CVCF <- function(x, ...)
{
      UseMethod("CVCF")
}

#' @rdname CVCF
#' @export
CVCF.formula <- function(formula,
                       data,
                       ...)
{
      if(!is.data.frame(data)){
            stop("data argument must be a data.frame")
      }
      modFrame <- model.frame(formula,data) # modFrame is a data.frame built from 'data' using the variables indicated in 'formula'. The first column of 'modFrame' is the response variable, thus we will indicate 'classColumn=1' when calling the HARF.default method in next line.
      attr(modFrame,"terms") <- NULL

      ret <- CVCF.default(x=modFrame,...,classColumn = 1)
      ret$call <- match.call(expand.dots = TRUE)
      ret$call[[1]] <- as.name("CVCF")
      # Next, we reconstruct the 'cleanData' from the removed and repaired indexes. Otherwise, the 'cleanData' would only contain those columns passed to the default method (for example imagine when running HARF(Species~Petal.Width+Sepal.Length,iris)).
      cleanData <- data
      if(!is.null(ret$repIdx)){
            cleanData[ret$repIdx,which(colnames(cleanData)==colnames(modFrame)[1])] <- ret$repLab  # This is not necessary in HARF because it only removes instances, it does not relabel. However, it must be used when the algorithm relabels instances (in our part there are some of them).
      }
      ret$cleanData <- cleanData[setdiff(1:nrow(cleanData),ret$remIdx),]
      return(ret)
}

#' @rdname CVCF
#' @export
CVCF.default <- function(x,
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

      if(consensus){
            threshold <- nfolds
      }else{
            threshold <- floor(nfolds/2)+1
      }

      formu <- as.formula(paste(names(x)[classColumn],"~.",sep = ""))

      folds <- caret::createFolds(x[,classColumn],nfolds)

      votes <- vector("integer",nrow(x))
      for(i in 1:nfolds){
            votes <- votes+(predict(RWeka::J48(formu,x[-folds[[i]],]),newdata=x)!=x[,classColumn])
      }
      ##### Building the 'filter' object ###########
      cleanData <- x[votes < threshold,]
      remIdx  <- which(votes>=threshold)
      repIdx <- NULL
      repLab <- NULL
      parameters <- list(nfolds=nfolds,
                         consensus=consensus)
      call <- match.call()
      call[[1]] <- as.name("CVCF")

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
