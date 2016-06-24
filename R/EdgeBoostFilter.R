# formu DONE
# library(RWeka), C4.5

#'@title Edge Boosting Filter
#'
#'@description Ensemble-based filter for removing label noise from a dataset as a
#' preprocessing step of classification. For more information, see 'Details' and
#' 'References' sections.
#'
#'@param formula A formula describing the classification variable and the attributes to be used.
#'@param data,x Data frame containing the tranining dataset to be filtered.
#'@param m Number of boosting iterations
#'@param percent Real number between 0 and 1. It sets the percentage of instances to be removed (as long as
#'their edge value exceeds the parameter \code{threshold}).
#'@param threshold Real number between 0 and 1. It sets the minimum edge value required
#'by an instance in order to be removed.
#'@param classColumn Positive integer indicating the column which contains the (factor of) classes.
#'By default, the last column is considered.
#'@param ... Optional parameters to be passed to other methods.
#'
#'@details The full description of the method can be looked up in the provided reference.
#'
#'An AdaBoost scheme (Freund & Schapire) is applied with a default C4.5 tree as weak classifier.
#'After \code{m} iterations, those instances with larger (according to the constraints
#'\code{percent} and \code{threshold}) edge values (Wheway, Freund & Schapire) are considered noisy
#'and thus removed.
#'
#'Notice that making use of extreme values (i.e. \code{percent=1} or \code{threshold=0}) any
#''removing constraints' can be ignored.
#'
#'@return An object of class \code{filter}, which is a list with seven components:
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
#'@references
#'Freund Y., Schapire R. E. (1997): A decision-theoretic generalization of on-line learning and
#'an application to boosting. \emph{Journal of computer and system sciences}, 55(1), 119-139.
#'
#'Wheway V. (2001, January): Using boosting to detect noisy data. In \emph{Advances in Artificial Intelligence}.
#'PRICAI 2000 Workshop Reader (pp. 123-130). Springer Berlin Heidelberg.
#'
#'@examples
#'# Next example is not run in order to save time
#'\dontrun{
#' data(iris)
#' out <- edgeBoostFilter(Species~., data = iris, m = 10, percent = 0.05, threshold = 0)
#' print(out)
#' identical(out$cleanData, iris[setdiff(1:nrow(iris),out$remIdx),])
#' }
#'@name edgeBoostFilter
NULL

#' @export
edgeBoostFilter <- function(x, ...)
{
      UseMethod("edgeBoostFilter")
}

#' @rdname edgeBoostFilter
#' @export
edgeBoostFilter.formula <- function(formula,
                                    data,
                                    ...)
{
      if(!is.data.frame(data)){
            stop("data argument must be a data.frame")
      }
      modFrame <- model.frame(formula,data) # modFrame is a data.frame built from 'data' using the variables indicated in 'formula'. The first column of 'modFrame' is the response variable, thus we will indicate 'classColumn=1' when calling the HARF.default method in next line.
      attr(modFrame,"terms") <- NULL

      ret <- edgeBoostFilter.default(x=modFrame,...,classColumn = 1)
      ret$call <- match.call(expand.dots = TRUE)
      ret$call[[1]] <- as.name("edgeBoostFilter")
      # Next, we reconstruct the 'cleanData' from the removed and repaired indexes. Otherwise, the 'cleanData' would only contain those columns passed to the default method (for example imagine when running HARF(Species~Petal.Width+Sepal.Length,iris)).
      cleanData <- data
      if(!is.null(ret$repIdx)){
            cleanData[ret$repIdx,which(colnames(cleanData)==colnames(modFrame)[1])] <- ret$repLab  # This is not necessary in HARF because it only removes instances, it does not relabel. However, it must be used when the algorithm relabels instances (in our part there are some of them).
      }
      ret$cleanData <- cleanData[setdiff(1:nrow(cleanData),ret$remIdx),]
      return(ret)
}

#' @rdname edgeBoostFilter
#' @export
edgeBoostFilter.default <- function(x,
                                    m = 15,
                                    percent = 0.05,
                                    threshold = 0,
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

      initialInstances <- nrow(x)
      mOrig <- m

      weights <- rep(1/initialInstances,initialInstances)
      weightsLearners <- vector("numeric",m)
      predictions <- matrix(nrow = initialInstances, ncol = m)
      edges <- vector("numeric",initialInstances)

      for(i in 1:m){
            preds <- predict(weightedC45(x,weights,formu),x)
            error <- sum(weights[which(preds!=x[,classColumn])])
            error <- ifelse(error>1e-6,error,1e-6)
            if(error>1/2){
                  warning("AdaBoost process stopped at iteration ", i, " because error exceeded 0.5")
                  m <- i-1
                  break
            }
            predictions[,i] <- preds==x[,classColumn]
            weightsLearners[i] <- log((1-error)/error)
            weights <- weights*ifelse(preds==x[,classColumn],1/(2*(1-error)),1/(2*error))
      }
      weightsLearners <- weightsLearners/sum(weightsLearners)
      for(i in 1:m){
            edges <- edges + ifelse(predictions[,i],0,weightsLearners[i])
      }
      indexesToKeep <- union(which(edges<stats::quantile(edges,1-percent)),which(edges<threshold))
      ##### Building the 'filter' object ###########
      extraInf <- paste("Highest edge value kept:",max(edges[indexesToKeep]))
      cleanData <- x[indexesToKeep,]
      remIdx  <- setdiff(1:nrow(x),indexesToKeep)
      repIdx <- NULL
      repLab <- NULL
      parameters <- list(m = mOrig,
                         percent = percent,
                         threshold = threshold)
      call <- match.call()
      call[[1]] <- as.name("edgeBoostFilter")

      ret <- list(cleanData = cleanData,
                  remIdx = remIdx,
                  repIdx=repIdx,
                  repLab=repLab,
                  parameters=parameters,
                  call = call,
                  extraInf = extraInf
      )
      class(ret) <- "filter"
      return(ret)

}

weightedC45 <- function(data,weights,formu){
      N <- 1/min(weights[which(weights>1e-5)])
      dataTrain <- data[rep(seq_len(nrow(data)),times=round(weights*N)),]
      return(RWeka::J48(formu,dataTrain))
}
