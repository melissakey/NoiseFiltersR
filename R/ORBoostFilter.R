# library(RWeka) # DecisionStump
# library(e1071) # naiveBayes

#'@title Outlier Removal Boosting Filter
#'
#'@description Ensemble-based filter for removing label noise from a dataset as a
#' preprocessing step of classification. For more information, see 'Details' and
#' 'References' sections.
#'
#'@param formula A formula describing the classification variable and the attributes to be used.
#'@param data,x Data frame containing the tranining dataset to be filtered.
#'@param N Number of boosting iterations.
#'@param d Threshold for removing noisy instances. Authors recommend to set it between 3 and 20. If it is set to \code{NULL},
#'the optimal threshold is chosen according to the procedure described in Karmaker & Kwek. However, this can be
#'very time-consuming, and in most cases is little relevant for the final result.
#'@param Naux Number of boosting iterations for AdaBoost when computing the optimal threshold 'd'.
#'@param useDecisionStump If \code{TRUE}, a decision stump is used as weak classifier.
#'Otherwise (default), naive-Bayes is applied. Recall decision stumps are not appropriate for multi-class problems.
#'@param classColumn Positive integer indicating the column which contains the (factor of) classes.
#'By default, the last column is considered.
#'@param ... Optional parameters to be passed to other methods.
#'
#'@details The full description of \code{ORBoostFilter} method can be looked up in Karmaker & Kwek.
#'In general terms, a weak classifier is built in each iteration, and misclassified instances have their weight
#'increased for the next round. Instances are removed when their weight exceeds the
#'threshold \code{d}, i.e. they have been misclassified in consecutive rounds.
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
#'Karmaker A., Kwek S. (2005, November): A boosting approach to remove class label noise.
#'In \emph{Hybrid Intelligent Systems}, 2005. HIS'05. Fifth International Conference on (pp. 6-pp). IEEE.
#'
#'Freund Y., Schapire R. E. (1997): A decision-theoretic generalization of on-line learning and
#'an application to boosting. \emph{Journal of computer and system sciences}, 55(1), 119-139.
#'@note
#'By means of a message, the number of noisy instances
#'removed in each iteration is displayed in the console.
#'@examples
#'# Next example is not run in order to save time
#'\dontrun{
#'data(iris)
#'out <- ORBoostFilter(Species~., data = iris, N = 10)
#'summary(out)
#'identical(out$cleanData, iris[setdiff(1:nrow(iris),out$remIdx),])
#'}
#'@name ORBoostFilter
NULL

#' @export
ORBoostFilter <- function(x, ...)
{
      UseMethod("ORBoostFilter")
}

#' @rdname ORBoostFilter
#' @export
ORBoostFilter.formula <- function(formula,
                                  data,
                                  ...)
{
      if(!is.data.frame(data)){
            stop("data argument must be a data.frame")
      }
      modFrame <- model.frame(formula,data) # modFrame is a data.frame built from 'data' using the variables indicated in 'formula'. The first column of 'modFrame' is the response variable, thus we will indicate 'classColumn=1' when calling the HARF.default method in next line.
      attr(modFrame,"terms") <- NULL

      ret <- ORBoostFilter.default(x=modFrame,...,classColumn = 1)
      ret$call <- match.call(expand.dots = TRUE)
      ret$call[[1]] <- as.name("ORBoostFilter")
      # Next, we reconstruct the 'cleanData' from the removed and repaired indexes. Otherwise, the 'cleanData' would only contain those columns passed to the default method (for example imagine when running HARF(Species~Petal.Width+Sepal.Length,iris)).
      cleanData <- data
      if(!is.null(ret$repIdx)){
            cleanData[ret$repIdx,which(colnames(cleanData)==colnames(modFrame)[1])] <- ret$repLab  # This is not necessary in HARF because it only removes instances, it does not relabel. However, it must be used when the algorithm relabels instances (in our part there are some of them).
      }
      ret$cleanData <- cleanData[setdiff(1:nrow(cleanData),ret$remIdx),]
      return(ret)
}

#' @rdname ORBoostFilter
#' @export
ORBoostFilter.default <- function(x,
                                  N = 20,
                                  d = 11,
                                  Naux = max(20,N),
                                  useDecisionStump = FALSE,
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

      namesOrig <- names(x)
      rownamesOrig <- attr(x,"row.names")
      names(x)[classColumn] <- "class"
      if(any(names(x)[-classColumn]=="class")){
            v <- names(x)[-classColumn]
            v[v=="class"] <- paste("classss",1:sum(v=="class"),sep="")
            names(x)[-classColumn] <- v
      }
      row.names(x) <- 1:nrow(x)
      origSize <- nrow(x)
      dNULL <- FALSE

      if(is.null(d)){
            d <- optimumThreshold(x,Naux,useDecisionStump)
            message("The optimum threshold is d = ",d,"\n")
            extraInf <- paste("The chosen threshold is d =",d)
            dNULL <- TRUE
      }
      else if(d<1){
            stop("the threshold must be a positive integer")
      }
      d <- as.integer(d)

      distributionWeights <- rep(1/origSize, origSize)
      weights <- rep(1, origSize)
      for(t in 1:N){
            if(nrow(x)==0){
                  stop("data set empty!")
            }
            predictions <- predict(weightedWeakLearner(x,distributionWeights,useDecisionStump),x)
            error <- sum(distributionWeights[which(predictions!=x$class)])
            error <- if(error!=0 & error!=1){#for values 0 or 1 we would have later problems when updating weights
                  error
            }else if(error==0){
              1e-10
            }else{
              1-1e-10
            }

            weights <- weights*ifelse(predictions==x$class,sqrt(error/(1-error)),sqrt((1-error)/error))
            indexesToKeep <- which(weights < d)
            distributionWeights <- distributionWeights[indexesToKeep]*(weights[indexesToKeep]^ifelse(predictions[indexesToKeep]==x$class[indexesToKeep],1,-1))
            distributionWeights <- distributionWeights/sum(distributionWeights)
            weights <- weights[indexesToKeep]

            message("Iteration ", t,": ", nrow(x)-length(indexesToKeep), " noisy instances removed.")
            x <- x[indexesToKeep,]
      }

      ##### Building the 'filter' object ###########
      remIdx  <- setdiff(1:origSize,as.integer(row.names(x)))
      names(x) <- namesOrig
      row.names(x) <- rownamesOrig[as.integer(row.names(x))]
      cleanData <- x
      repIdx <- NULL
      repLab <- NULL
      parameters <- list(N = N,
                         d = d,
                         Naux = Naux,
                         useDecisionStump = useDecisionStump)
      call <- match.call()
      call[[1]] <- as.name("ORBoostFilter")

      ret <- list(cleanData = cleanData,
                  remIdx = remIdx,
                  repIdx=repIdx,
                  repLab=repLab,
                  parameters=parameters,
                  call = call,
                  extraInf = NULL)
      if(dNULL){
            ret$extraInf <- extraInf
      }
      class(ret) <- "filter"
      return(ret)
}

# AdaBoostFilter <- function(data, N = 100, d = 11, useDecisionStump = FALSE, classColumn=ncol(data)){
#       data <- as.data.frame(data)
#       names(data)[classColumn] <- "class"
#       initialInstances <- nrow(data)
#
#       if(!is.null(d)){
#             return(AdaBoostFilterAux(data,N,d,useDecisionStump)$cleanData)
#       }
#       d <- optimumThreshold(data,N,useDecisionStump)
#       cat("The optimum threshold is d =",d,"\n\n")
#       AdaBoostFilterAux(data,N,d,useDecisionStump)$cleanData
# }

#next function is AdaBoostFilter for which 'd' must be explicitly specified.
#It returns also the 'AdaBoost-finalclassifier' (i.e. the 'N' weak classifiers together with their associated weights), so
#that 'optimumThreshold' can compute the accuracy for differents 'd' values.
AdaBoostFilterAux <- function(data, N, d, useDecisionStump){
      d <- as.integer(d)
      initialInstances <- nrow(data)
      weights <- rep(1/initialInstances,initialInstances)

      out <- list(classifiers=list(), classifWeights=vector("numeric",N), cleanData=NULL)
      length(out$classifiers) <- N #we initialize the size of the objects in the output, for greater efficiency in the 'for' loop
      for(t in 1:N){
            classif <- weightedWeakLearner(data,weights,useDecisionStump)
            predictions <- predict(classif,data)

            error <- sum(weights[which(predictions!=data$class)])
            beta <- if(error!=0 & error!=1){#for values 0 or 1 we would have later problems when updating weights
                  error/(1-error)
            }else if(error==0){
                1e-10
            }else{
                1e10
            }
            if(error>1/2){#This is a condition specified in AdaBoost paper
                  warning("error greater than 1/2 prevents from reaching the specified number of iterations")
                  N <- t-1
                  length(out$classifiers) <- t-1
                  length(out$classifWeights) <- t-1
                  break
            }
            #updating weights
            weights <- weights*ifelse(predictions==data$class,beta,1)
            weights <- weights/sum(weights)
            indexesToKeep <- which(weights <= d/nrow(data))
            weights <- weights[indexesToKeep]
            weights <- weights/sum(weights)
            #adding weak classifier and its weight to the output
            out$classifiers[[t]] <- classif
            out$classifWeights[t] <- log(1/beta)

            data <- data[indexesToKeep,]
      }
      out$cleanData <- data
      out
}

optimumThreshold <- function(data, N, useDecisionStump){
      trainIndexes <- caret::createDataPartition(data$class,p=0.8)[[1]]
      trainData <- data[trainIndexes,]
      validationData <- data[-trainIndexes,]

      bestAccuracy <- 0
      bestThreshold <- 0

      for(d in 20:3){
            results <- AdaBoostFilterAux(trainData,N,d,useDecisionStump)
            accuracy <- getAccuracy(results$classifiers,results$classifWeights,validationData)
            if(accuracy > bestAccuracy){
                  bestAccuracy <- accuracy
                  bestThreshold <- d
            }
            if(bestAccuracy==1){
                  break
            }
      }
      bestThreshold
}

#next function returns the accuracy when testing on 'testData' the classifier consisting of 'classifiers' weighted
#by means of 'weights'.
getAccuracy <- function(classifiers,weights,testData){
      classifications <- sapply(classifiers,function(cl){unclass(predict(cl,testData))})
      totalWeightForClass <- matrix(0,nrow = nrow(testData),ncol = nlevels(testData$class))
      for(i in 1:nrow(classifications)){
            for(j in 1:ncol(classifications)){
                  totalWeightForClass[i,classifications[i,j]] <- totalWeightForClass[i,classifications[i,j]]+weights[j]
            }
      }
      predictions <- apply(totalWeightForClass,1,which.max)
      sum(predictions==unclass(testData$class))/length(predictions)
}

#in order to obtain a 'weighted classifier', in next function we create an appropriately weighted training dataset from the original dataset and the associated weights
weightedWeakLearner <- function(data, weights, useDecisionStump){
      N <- 1/min(weights[which(weights>1e-5)])
      dataTrain <- data[rep(seq_len(nrow(data)),times=sapply(round(weights*N),function(x){max(x,1)})),] #each instance is repeated according to its weight (and appears at least once)
      if(useDecisionStump){
            return(RWeka::DecisionStump(class~.,dataTrain))
      }else{
            return(e1071::naiveBayes(class~.,dataTrain))
      }
}

