# library(RWeka)
# library(caret) #to use 'createFolds'
# library(nnet) #which.is.max

#'@title Saturation Filters
#'
#'@description Data complexity based filters for removing label noise from a dataset as a
#' preprocessing step of classification. For more information, see 'Details' and
#' 'References' sections.
#'
#'@param formula A formula describing the classification variable and the attributes to be used.
#'@param data,x Data frame containing the tranining dataset to be filtered.
#'@param noiseThreshold The threshold for removing noisy instances in the saturation filter.
#'Authors recommend values between 0.25 and 2. If it is set to \code{NULL}, the
#'threshold is appropriately chosen according to the number of training instances.
#'@param nfolds For \code{consensusSF} and \code{classifSF}, number of folds
#'in which the dataset is split.
#'@param consensusLevel For \code{consensusSF}, it sets the (minimum) number of 'noisy
#'votes' an instance must get in order to be removed. By default, the \code{nfolds-1} filters built
#'over each instance must label it as noise.
#'@param classColumn Positive integer indicating the column which contains the (factor of) classes.
#'By default, the last column is considered.
#'@param ... Optional parameters to be passed to other methods.
#'
#'@details
#'Based on theoretical studies about data complexity (Gamberger & Lavrac, 1997),
#'\code{saturationFilter} removes those
#'instances which most enable to reduce the CLCH (Complexity of the Least Complex Hypotheses)
#'of the training dataset. The full method can be looked up in (Gamberger et al., 1999), and
#'the previous step of \emph{literals} extraction is detailed in (Gamberger et al., 1996).
#'
#'\code{consensusSF} splits the dataset in \code{nfolds} folds, and applies
#'\code{saturationFilter} to every combination of \code{nfolds-1} folds. Those instances
#'with (at least) \code{consensusLevel} 'noisy votes' are removed.
#'
#'\code{classifSF} combines \code{saturationFilter} with a \code{nfolds}-folds cross validation
#'scheme (the latter in the spirit of filters such as \code{\link{EF}}, \code{\link{CVCF}}).
#'Namely, the dataset is split in \code{nfolds} folds and, for every combination
#'of \code{nfolds-1} folds, \code{saturationFilter} is applied and a classifier
#'(we implement a standard C4.5 tree) is built. Instances
#'from the excluded fold are removed according to this classifier.
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
#'Gamberger D., Lavrac N., Groselj C. (1999, June): Experiments with noise
#'filtering in a medical domain. In \emph{ICML} (pp. 143-151).
#'
#'Gamberger D., Lavrac N., Dzeroski S. (1996, January): Noise elimination in
#'inductive concept learning: A case study in medical diagnosis. In
#'\emph{Algorithmic Learning Theory} (pp. 199-212). Springer Berlin Heidelberg.
#'
#'Gamberger D., Lavrac N. (1997): Conditions for Occam's razor applicability
#'and noise elimination (pp. 108-123). Springer Berlin Heidelberg.
#'
#'@examples
#'# Next example is not run because saturation procedure is time-consuming.
#'\dontrun{
#'data(iris)
#'out1 <- saturationFilter(Species~., data = iris)
#'out2 <- consensusSF(Species~., data = iris)
#'out3 <- classifSF(Species~., data = iris)
#'print(out1)
#'print(out2)
#'print(out3)
#'}
#'@aliases consensusSF classifSF
#'@name saturationFilter
NULL

#' @export
saturationFilter <- function(x, ...)
{
  UseMethod("saturationFilter")
}

#' @rdname saturationFilter
#' @export
saturationFilter.formula <- function(formula,
  data,
  ...)
{
  if(!is.data.frame(data)){
    stop("data argument must be a data.frame")
  }
  modFrame <- model.frame(formula,data) # modFrame is a data.frame built from 'data' using the variables indicated in 'formula'. The first column of 'modFrame' is the response variable, thus we will indicate 'classColumn=1' when calling the HARF.default method in next line.
  attr(modFrame,"terms") <- NULL
  
  ret <- saturationFilter.default(x=modFrame,...,classColumn = 1)
  ret$call <- match.call(expand.dots = TRUE)
  ret$call[[1]] <- as.name("saturationFilter")
  # Next, we reconstruct the 'cleanData' from the removed and repaired indexes. Otherwise, the 'cleanData' would only contain those columns passed to the default method (for example imagine when running HARF(Species~Petal.Width+Sepal.Length,iris)).
  cleanData <- data
  if(!is.null(ret$repIdx)){
    cleanData[ret$repIdx,which(colnames(cleanData)==colnames(modFrame)[1])] <- ret$repLab  # This is not necessary in HARF because it only removes instances, it does not relabel. However, it must be used when the algorithm relabels instances (in our part there are some of them).
  }
  ret$cleanData <- cleanData[setdiff(1:nrow(cleanData),ret$remIdx),]
  return(ret)
}

#' @rdname saturationFilter
#' @export
saturationFilter.default <- function(x,
  noiseThreshold = NULL,
  classColumn = ncol(x),
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
  
  dataOrig <- x
  
  rownames(x) <- 1:nrow(x)
  inconsistencies <- (duplicated(x) != duplicated(x[-classColumn])) | 
    (duplicated(x, fromLast = TRUE) != duplicated(x[-classColumn], fromLast = TRUE))
  
  if(sum(inconsistencies) > 0){
    splitData <- split(x, as.list(x[,-classColumn]), drop = TRUE)
    splitData <- lapply(splitData, function(df){
      df[,classColumn] <- factor(levels(df[,classColumn])[nnet::which.is.max(table(df[,classColumn]))], levels = levels(df[,classColumn]))
      df
    })
    x <- do.call("rbind", splitData)
    rownames(x) <- unlist(lapply(splitData, rownames))
    x <- x[order(as.numeric(rownames(x))),]
    
    warning("There were inconsistencies (same attributes and different class) in the dataset to be filtered. \n To avoid stopping the process, majority class was internally assigned to all instances with the same attributes. Their real class was recovered at the end of the process")
  }
  
  names(x)[classColumn] <- "class"
  if(any(names(x)[-classColumn]=="class")){
    v <- names(x)[-classColumn]
    v[v=="class"] <- paste("classss",1:sum(v=="class"),sep="")
    names(x)[-classColumn] <- v
  }
  
  x$class <- factor(x$class) #unused levels are dropped
  numClasses <- nlevels(x$class)
  
  if(numClasses<2){
    stop("all the instances belong to the same class")
  }
  if(numClasses==2){
    indexesToRemove <- twoClassSaturationFilter(x,noiseThreshold)
    ##### Building the 'filter' object ###########
    cleanData <- dataOrig[setdiff(1:nrow(x),indexesToRemove),]
    remIdx  <- indexesToRemove
    repIdx <- NULL
    repLab <- NULL
    parameters <- list(noiseThreshold = noiseThreshold)
    call <- match.call()
    call[[1]] <- as.name("saturationFilter")
    
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
  indexesToRemove <- integer(0)
  dataCopy <- x
  dataCopy$class <- factor(rep("yes",nrow(x)),levels = c("yes","no"))
  for(i in 1:numClasses){
    dataCopy[which(x$class==levels(x$class)[i]),"class"] <- "yes"
    dataCopy[which(x$class!=levels(x$class)[i]),"class"] <- "no"
    indexesToRemove <- union(indexesToRemove,twoClassSaturationFilter(dataCopy,noiseThreshold))
  }
  ##### Building the 'filter' object ###########
  cleanData <- dataOrig[setdiff(1:nrow(x),indexesToRemove),]
  remIdx  <- indexesToRemove
  repIdx <- NULL
  repLab <- NULL
  parameters <- list(noiseThreshold = noiseThreshold)
  call <- match.call()
  call[[1]] <- as.name("saturationFilter")
  
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

#' @export
consensusSF <- function(x, ...)
{
  UseMethod("consensusSF")
}

#' @rdname saturationFilter
#' @export
consensusSF.formula <- function(formula,
  data,
  ...)
{
  if(!is.data.frame(data)){
    stop("data argument must be a data.frame")
  }
  modFrame <- model.frame(formula,data) # modFrame is a data.frame built from 'data' using the variables indicated in 'formula'. The first column of 'modFrame' is the response variable, thus we will indicate 'classColumn=1' when calling the HARF.default method in next line.
  attr(modFrame,"terms") <- NULL
  
  ret <- consensusSF.default(x=modFrame,...,classColumn = 1)
  ret$call <- match.call(expand.dots = TRUE)
  ret$call[[1]] <- as.name("consensusSF")
  # Next, we reconstruct the 'cleanData' from the removed and repaired indexes. Otherwise, the 'cleanData' would only contain those columns passed to the default method (for example imagine when running HARF(Species~Petal.Width+Sepal.Length,iris)).
  cleanData <- data
  if(!is.null(ret$repIdx)){
    cleanData[ret$repIdx,which(colnames(cleanData)==colnames(modFrame)[1])] <- ret$repLab  # This is not necessary in HARF because it only removes instances, it does not relabel. However, it must be used when the algorithm relabels instances (in our part there are some of them).
  }
  ret$cleanData <- cleanData[setdiff(1:nrow(cleanData),ret$remIdx),]
  return(ret)
}

#' @rdname saturationFilter
#' @export
consensusSF.default <- function(x,
  nfolds=10,
  consensusLevel = nfolds-1,
  noiseThreshold = NULL,
  classColumn = ncol(x),
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
  
  folds <- caret::createFolds(x[,classColumn],nfolds)
  votes <- integer(nrow(x))
  for(i in 1:nfolds){
    indNoise <- saturationFilter(x[-folds[[i]],], noiseThreshold, classColumn)$remIdx
    indTested <- unlist(folds[-i])
    votes[indTested[indNoise]] <- votes[indTested[indNoise]] + 1
  }
  ##### Building the 'filter' object ###########
  cleanData <- x[votes<consensusLevel,]
  remIdx  <- which(votes>=consensusLevel)
  repIdx <- NULL
  repLab <- NULL
  parameters <- list(nfolds=nfolds,
    consensusLevel = consensusLevel,
    noiseThreshold = noiseThreshold)
  call <- match.call()
  call[[1]] <- as.name("consensusSF")
  
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


#' @export
classifSF <- function(x, ...)
{
  UseMethod("classifSF")
}

#' @rdname saturationFilter
#' @export
classifSF.formula <- function(formula,
  data,
  ...)
{
  if(!is.data.frame(data)){
    stop("data argument must be a data.frame")
  }
  modFrame <- model.frame(formula,data) # modFrame is a data.frame built from 'data' using the variables indicated in 'formula'. The first column of 'modFrame' is the response variable, thus we will indicate 'classColumn=1' when calling the HARF.default method in next line.
  attr(modFrame,"terms") <- NULL
  
  ret <- classifSF.default(x=modFrame,...,classColumn = 1)
  ret$call <- match.call(expand.dots = TRUE)
  ret$call[[1]] <- as.name("classifSF")
  # Next, we reconstruct the 'cleanData' from the removed and repaired indexes. Otherwise, the 'cleanData' would only contain those columns passed to the default method (for example imagine when running HARF(Species~Petal.Width+Sepal.Length,iris)).
  cleanData <- data
  if(!is.null(ret$repIdx)){
    cleanData[ret$repIdx,which(colnames(cleanData)==colnames(modFrame)[1])] <- ret$repLab  # This is not necessary in HARF because it only removes instances, it does not relabel. However, it must be used when the algorithm relabels instances (in our part there are some of them).
  }
  ret$cleanData <- cleanData[setdiff(1:nrow(cleanData),ret$remIdx),]
  return(ret)
}

#' @rdname saturationFilter
#' @export
classifSF.default <- function(x,
  nfolds=10,
  noiseThreshold = NULL,
  classColumn = ncol(x),
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
  names(x)[classColumn] <- "class"
  if(any(names(x)[-classColumn]=="class")){
    v <- names(x)[-classColumn]
    v[v=="class"] <- paste("classss",1:sum(v=="class"),sep="")
    names(x)[-classColumn] <- v
  }
  
  folds <- caret::createFolds(x$class,nfolds)
  isNoise <- logical(nrow(x))
  for(i in 1:nfolds){
    isNoise[folds[[i]]] <- predict(RWeka::J48(class~.,saturationFilter(x[-folds[[i]],], noiseThreshold, classColumn)$cleanData),
      newdata=x[folds[[i]],])!=x$class[folds[[i]]]
  }
  ##### Building the 'filter' object ###########
  names(x) <- namesOrig
  cleanData <- x[!isNoise,]
  remIdx  <- which(isNoise)
  repIdx <- NULL
  repLab <- NULL
  parameters <- list(nfolds=nfolds,
    noiseThreshold = noiseThreshold)
  call <- match.call()
  call[[1]] <- as.name("classifSF")
  
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

