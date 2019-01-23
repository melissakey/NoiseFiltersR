# library(kknn) # 3-NN
# library(RWeka) # C4.5
# library(nnet) # multinom: multinomial (more than two classes) logistic regression.

#'@title Iterative Noise Filter based on the Fusion of Classifiers
#'
#'@description Ensemble-based filter for removing label noise from a dataset as a
#' preprocessing step of classification. For more information, see 'Details' and
#' 'References' sections.
#'
#'@param formula A formula describing the classification variable and the attributes to be used.
#'@param data,x Data frame containing the tranining dataset to be filtered.
#'@param p Real number between 0 and 1. It sets the minimum proportion of original
#'instances which must be tagged as noisy in order to go for another iteration.
#'@param s Positive integer setting the stop criterion together with \code{p}. The filter stops
#'after \code{s} iterations with not enough noisy instances removed (according to the proportion \code{p}).
#'@param consensus Logical. If FALSE, majority voting scheme is used for 'preliminary filtering' and 'noise free filtering' (see 'Details' and References' section). If TRUE, consensus
#'voting scheme is applied.
#'@param k Parameter for the k-nearest neighbors algorithm used for the 'noise score' stage (see 'Details' and 'References').
#'@param threshold Real number between -1 and 1. It sets the noise score value above which an instance is removed.
#'@param classColumn Positive integer indicating the column which contains the (factor of) classes.
#'By default, the last column is considered.
#'@param ... Optional parameters to be passed to other methods.

#'@details The full description of the method can be looked up in the provided reference.
#'A 'preliminary filtering' is carried out with a fusion of classifiers (FC), including C4.5, 3NN, and logistic regression. Then,
#'potentially noisy instances are identified in a 'noise free filtering' process building the FC on the (preliminary) filtered
#'instances. Finally, a 'noise score' is computed on these potentially noisy instances, removing those exceeding the \code{threshold} value.
#'The process stops after \code{s} iterations with not enough (according to the proportion \code{p}) noisy
#'instances removed.
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
#'S\'{a}ez J. A., Galar M., Luengo J., Herrera F. (2016): INFFC: An iterative class noise filter
#'based on the fusion of classifiers with noise sensitivity control. \emph{Information Fusion}, 27, 19-32.
#'@note
#'By means of a message, the number of noisy instances removed
#'in each iteration is displayed in the console.
#'@examples
#'# Next example is not run because it might be time-consuming
#'\dontrun{
#' data(iris)
#' out <- INFFC(Species~., data = iris)
#' summary(out)
#' identical(out$cleanData, iris[setdiff(1:nrow(iris),out$remIdx),])
#'}
#'@name INFFC
NULL

#' @export
INFFC <- function(x, ...)
{
      UseMethod("INFFC")
}

#' @rdname INFFC
#' @export
INFFC.formula <- function(formula,
                          data,
                          ...)
{
      if(!is.data.frame(data)){
            stop("data argument must be a data.frame")
      }
      modFrame <- model.frame(formula,data) # modFrame is a data.frame built from 'data' using the variables indicated in 'formula'. The first column of 'modFrame' is the response variable, thus we will indicate 'classColumn=1' when calling the HARF.default method in next line.
      attr(modFrame,"terms") <- NULL

      ret <- INFFC.default(x=modFrame,...,classColumn = 1)
      ret$call <- match.call(expand.dots = TRUE)
      ret$call[[1]] <- as.name("INFFC")
      # Next, we reconstruct the 'cleanData' from the removed and repaired indexes. Otherwise, the 'cleanData' would only contain those columns passed to the default method (for example imagine when running HARF(Species~Petal.Width+Sepal.Length,iris)).
      cleanData <- data
      if(!is.null(ret$repIdx)){
            cleanData[ret$repIdx,which(colnames(cleanData)==colnames(modFrame)[1])] <- ret$repLab  # This is not necessary in HARF because it only removes instances, it does not relabel. However, it must be used when the algorithm relabels instances (in our part there are some of them).
      }
      ret$cleanData <- cleanData[setdiff(1:nrow(cleanData),ret$remIdx),]
      return(ret)
}

#' @rdname INFFC
#' @export
INFFC.default <- function(x,
                          consensus=FALSE,
                          p=0.01,
                          s=3,
                          k=5,
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

      origSize <- nrow(x)
      namesOrig <- names(x)
      rownamesOrig <- attr(x,"row.names")
      names(x)[classColumn] <- "class"
      if(any(names(x)[-classColumn]=="class")){
            v <- names(x)[-classColumn]
            v[v=="class"] <- paste("classss",1:sum(v=="class"),sep="")
            names(x)[-classColumn] <- v
      }
      row.names(x) <- 1:nrow(x)

      #setting some parameters and auxiliary variables
      if(consensus)
            majThreshold <- 3
      else
            majThreshold <- 2
      stopThreshold <- floor(origSize*p)
      KeepOn <- TRUE #will control the while loop
      counter <- 0 #will control how many consecutive times we filter few nosiy data
      countIter <- 0 # counts the total iterations

      while(KeepOn){
            countIter <- countIter+1

            ##Preliminary Filtering
            PreFiltIndexes <- FusionClassifiers(x,trainingIndexes=1:nrow(x),majThreshold)
            ##Noise free Filtering
            NoisyIndexes <- FusionClassifiers(x,trainingIndexes=PreFiltIndexes,majThreshold,returnNoisy = TRUE)
            ##Noise scoring
            scores <- sapply(NoisyIndexes,function(i){NoiseScore(x,NoisyIndexes,k,i)})

            IndexesToRemove <- NoisyIndexes[which(scores > threshold)]
            if(length(IndexesToRemove)>0){
                  x <- x[-IndexesToRemove,]
            }
            #refreshing stopping conditions
            if( length(IndexesToRemove) <= stopThreshold & counter+1==s) KeepOn <- FALSE
            if( length(IndexesToRemove) <= stopThreshold & counter+1<s) counter <- counter+1
            if(length(IndexesToRemove) > stopThreshold) counter <- 0

            message("Iteration ", countIter,": ", length(IndexesToRemove), " noisy instances removed \n")
      }
      ##### Building the 'filter' object ###########
      remIdx  <- setdiff(1:origSize,as.integer(row.names(x)))
      names(x) <- namesOrig
      row.names(x) <- rownamesOrig[as.integer(row.names(x))]
      cleanData <- x
      repIdx <- NULL
      repLab <- NULL
      parameters <- list(consensus=consensus,
                         p=p,
                         s=s,
                         k=k,
                         threshold=threshold)
      call <- match.call()
      call[[1]] <- as.name("INFFC")

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
