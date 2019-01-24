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


twoClassSaturationFilter <- function(data, noiseThreshold = NULL){
  initialInstances <- nrow(data)
  yes <- levels(data$class)[1]
  no <- levels(data$class)[2]
  
  if(is.null(noiseThreshold)){
    noiseThreshold <- getThreshold(min( sum(data$class==yes),sum(data$class==no) ))
  }
  
  cove <- coverage3array(data)
  
  indexesToRemove <- integer(0)
  KeepOn <- TRUE
  count <- 1
  while(KeepOn){
    #cat("Round:",count,"\n")
    if(dim(cove)[1]==0 || dim(cove)[2]==0){
      stop("a whole class was filtered!")
    }
    weightsPos <- numeric(dim(cove)[1])
    weightsNeg <- numeric(dim(cove)[2])
    literals <- optimLiterals(cove)
    for(k in literals){
      positive <- integer(0)
      negative <- integer(0)
      indexes <- arrayInd(which(cove[,,k] & apply(cove[,,literals],c(1,2),sum)==1),dim(cove)[c(1,2)])
      weightsPos[unique(indexes[,1])] <- weightsPos[unique(indexes[,1])]+1/length(unique(indexes[,1]))
      weightsNeg[unique(indexes[,2])] <- weightsNeg[unique(indexes[,2])]+1/length(unique(indexes[,2]))
    }
    if(max(c(weightsNeg,weightsPos))>noiseThreshold){
      count <- count+1
      if(max(weightsNeg)>max(weightsPos)||(max(weightsNeg)==max(weightsPos) && dim(cove)[1]<dim(cove)[2])){
        indexesToRemove <- c(indexesToRemove,which(data$class[setdiff(1:initialInstances,indexesToRemove)]==no)[which.max(weightsNeg)])
        cove <- cove[,- which.max(weightsNeg),]
      }
      else{
        indexesToRemove <- c(indexesToRemove,which(data$class[setdiff(1:initialInstances,indexesToRemove)]==yes)[which.max(weightsPos)])
        cove <- cove[- which.max(weightsPos),,]
      }
    }else{
      KeepOn <- FALSE
    }
  }
  indexesToRemove
}

coverage3array <- function(data){
  yes <- levels(data$class)[1]
  no <- levels(data$class)[2]
  posIndexes <- which(data$class==yes)
  negIndexes <- which(data$class==no)
  numLiterals <- 0
  out <- logical(0)
  for(col in 1:ncol(data)){
    if(names(data)[col]=="class"){
      next
    }
    if(class(data[,col])=="factor" || class(data[,col])=="integer"){
      #literals created from 'positive' instances
      literals <- unique(data[posIndexes,col])
      numLiterals <- numLiterals+length(literals)
      out <- c(out,sapply(literals,function(l){sapply(negIndexes,function(iNeg){sapply(posIndexes,function(iPos){data[iPos,col]==l && data[iNeg,col]!=l})})}))
      #literals created from 'negative' instances
      literals <- unique(data[negIndexes,col])
      numLiterals <- numLiterals+length(literals)
      out <- c(out,sapply(literals,function(l){sapply(negIndexes,function(iNeg){sapply(posIndexes,function(iPos){data[iPos,col]!=l && data[iNeg,col]==l})})}))
    }
    if(class(data[,col])=="numeric" || class(data[,col])=="integer"){
      literalsLQ <- numeric(0)
      literalsG <- numeric(0)
      v <- unique(data[posIndexes,col])
      names(v) <- rep(yes,length(v))
      w <- unique(data[negIndexes,col])
      names(w) <- rep(no,length(w))
      values <- ordering(c(v,w))
      currentClass <- names(values)[1]
      for(i in 1:(length(values)-1)){
        if(names(values)[i+1]!=currentClass){
          if(currentClass==yes){
            literalsLQ <- c(literalsLQ,mean(values[c(i,i+1)]))
          }
          else{
            literalsG <- c(literalsG,mean(values[c(i,i+1)]))
          }
          currentClass <- names(values)[i+1]
        }
      }
      #                   if(data[1,"class"]==yes){
      #                         cat("Column",col,":\n")
      #                         cat("Values:\n")
      #                         print(values)
      #                         cat("LiteralsLQ:\n",literalsLQ,"\n")
      #                         cat("LiteralsG:\n",literalsG,"\n")
      #                   }
      numLiterals <- numLiterals+length(literalsLQ)+length(literalsG)
      out <- c(out,sapply(literalsLQ,function(l){sapply(negIndexes,function(iNeg){sapply(posIndexes,function(iPos){data[iPos,col]<=l && data[iNeg,col]>=l && data[iPos,col]!=data[iNeg,col]})})}))
      out <- c(out,sapply(literalsG,function(l){sapply(negIndexes,function(iNeg){sapply(posIndexes,function(iPos){data[iPos,col]>=l && data[iNeg,col]<=l && data[iPos,col]!=data[iNeg,col]})})}))
    }
  }
  if(class(out)=="list"){
    out <- unlist(out)
  }
  dim(out) <- c(length(posIndexes),length(negIndexes),numLiterals)
  out
}

optimLiterals <- function(cove){
  coveCopy <- cove
  out <- integer(0)
  subround <- 1
  while(!all(is.na(cove))){
    #cat("Subround:",subround,"\n")
    #print(str(apply(cove,c(1,2),sum)))
    #print(all(is.na(apply(cove,c(1,2),sum))))
    #summing <- apply(cove,c(1,2),sum)
    #print(str(summing[which(!is.na(summing))]))
    weights <- 1/apply(cove,c(1,2),sum)
    #PRINT CLASS WEIGHTS TO MAKE SURE IT IS A MATRIX
    #print(str(weights))
    #print(all(is.na(weights)))
    #print(summary(as.vector(weights)))
    indMax <- arrayInd(which.max(weights),dim(cove)[c(1,2)])
    #print(weights[indMax])
    litToTest <- which(cove[indMax[1],indMax[2],])
    #print(length(litToTest))
    weightsLit <- numeric(length(litToTest))
    for(k in seq_along(litToTest)){
      arrIndexes <- arrayInd(which(!is.na(cove[,,litToTest[k]]) & cove[,,litToTest[k]]),dim(cove)[c(1,2)])
      #print(summary(apply(arrIndexes,2,function(ind){weights[ind]})))
      #print(class(arrIndexes))
      #apply(arrIndexes,2,function(ind){weights[ind]})
      weightsLit[k] <- sum(apply(arrIndexes,1,function(ind){weights[ind[1],ind[2]]}))
    }
    #print(weightsLit)
    litToKeep <- litToTest[which.max(weightsLit)]
    #print(weightsLit[which.max(weightsLit)])
    #print(str(weightsLit))
    out <- c(out,litToKeep)
    covered <- arrayInd(which(!is.na(cove[,,litToKeep]) & cove[,,litToKeep]),dim(cove)[c(1,2)])
    #print(length(which(!is.na(cove[,,litToKeep]) & cove[,,litToKeep])))
    for(i in 1:nrow(covered)){
      cove[covered[i,1],covered[i,2],] <- NA
    }
    cove <- cove[,,-litToKeep]
    #print(any(apply(cove,c(1,2),function(v){any(is.na(v)) && !all(is.na(v))})))
  }
  for(k in out){
    if(!any(coveCopy[,,k] & apply(coveCopy[,,out],c(1,2),sum)==1)){
      out <- setdiff(out,k)
    }
  }
  out
}

getThreshold <- function(x){
  if(x<=50){
    out <- 1.5
  }else if(50<x && x<=100){
    out <- 1
  }else if(100<x && x<=200){
    out <- 0.75
  }else{
    out <- 0.5
  }
  out
}

ordering <- function(values){
  out <- values[order(values)]
  for(i in 2:(length(out)-1)){
    if(out[i]==out[i+1] && names(out)[i+1]==names(out)[i-1]){
      names(out)[i+1] <- names(out)[i]
      names(out)[i] <- names(out)[i-1]
    }
  }
  #       if(out[1]==out[2] && names(out)[1]==names(out)[3]){
  #             names(out)[1] <- names(out)[2]
  #             names(out)[2] <- names(out)[3]
  #       }
  out
}
