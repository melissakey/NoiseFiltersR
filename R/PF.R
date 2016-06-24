# library(RWeka) # PART
# library(caret) # to use 'createFolds'
# library(rJava) # Dealing with Java objects in the step of "good rules selection".

#'@title Partitioning Filter
#'
#'@description Ensemble-based filter for removing label noise from a dataset as a
#' preprocessing step of classification. For more information, see 'Details' and
#' 'References' sections.
#'
#'@param formula A formula describing the classification variable and the attributes to be used.
#'@param data,x Data frame containing the tranining dataset to be filtered.
#'@param nfolds Number of partitions in each iteration.
#'@param p Real number between 0 and 1. It sets the minimum proportion of original
#'instances which must be tagged as noisy in order to go for another iteration.
#'@param s Positive integer setting the stop criterion together with \code{p}. The filter stops
#'after \code{s} iterations with not enough noisy instances removed (according to the proportion \code{p}).
#'@param y Real number between 0 and 1. It sets the proportion of good instances which
#'must be stored in each iteration.
#'@param theta Real number between 0 and 1. It sets the proportion of 'good rules' to be selected (see also
#''Details' section).
#'@param consensus Logical. If FALSE, majority voting scheme is used. If TRUE, consensus
#'voting scheme is applied.
#'@param classColumn Positive integer indicating the column which contains the (factor of) classes.
#'By default, the last column is considered.
#'@param ... Optional parameters to be passed to other methods.
#'
#'@details The full description of the method can be looked up in the provided references.
#'A PART rules set (from RWeka) is built in each of the \code{nfolds} partitions of \code{data}. After a
#''good rules selection' process based on the accuracy of each rule, the subsequent good rules sets are
#'tested in the whole dataset, and the removal of noisy instances is decided via consensus or
#'majority voting schemes. Finally, a proportion of good instances (i.e. those whose label agrees
#'with all the base classifiers) is stored and not considered in subsequent iterations. The process stops
#'after \code{s} iterations with not enough (according to the proportion \code{p}) noisy
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
#'Zhu X., Wu X., Chen Q. (2003, August): Eliminating class noise in large
#'datasets. \emph{International Conference in Machine Learning} (Vol. 3, pp. 920-927).
#'
#'Zhu X., Wu X., Chen Q. (2006): Bridging local and global data cleansing: Identifying class noise in
#'large, distributed data datasets. \emph{Data mining and Knowledge discovery}, 12(2-3), 275-308.
#'@note
#'The base rule classifier used is PART instead of C4.5rules used in the references.
#'
#'For the 'good rules selection' step, we implement the 'Best-L rules' scheme since, according to
#'the authors, it usually outperforms the others 'Adaptive Threshold' and 'Fixed Threshold' schemes.
#'
#'By means of a message, the number of noisy instances removed
#'in each iteration is displayed in the console.
#'@examples
#'# Next example is not run in order to save time
#'\dontrun{
#' data(iris)
#' # We fix a seed since there exists a random partition for the ensemble
#' set.seed(1)
#' out <- PF(Species~., data = iris, s = 1, nfolds = 3)
#' print(out)
#' identical(out$cleanData, iris[setdiff(1:nrow(iris),out$remIdx),])
#' }
#'@name PF
NULL

#' @export
PF <- function(x, ...)
{
      UseMethod("PF")
}

#' @rdname PF
#' @export
PF.formula <- function(formula,
                         data,
                         ...)
{
      if(!is.data.frame(data)){
            stop("data argument must be a data.frame")
      }
      modFrame <- model.frame(formula,data) # modFrame is a data.frame built from 'data' using the variables indicated in 'formula'. The first column of 'modFrame' is the response variable, thus we will indicate 'classColumn=1' when calling the HARF.default method in next line.
      attr(modFrame,"terms") <- NULL

      ret <- PF.default(x=modFrame,...,classColumn = 1)
      ret$call <- match.call(expand.dots = TRUE)
      ret$call[[1]] <- as.name("PF")
      # Next, we reconstruct the 'cleanData' from the removed and repaired indexes. Otherwise, the 'cleanData' would only contain those columns passed to the default method (for example imagine when running HARF(Species~Petal.Width+Sepal.Length,iris)).
      cleanData <- data
      if(!is.null(ret$repIdx)){
            cleanData[ret$repIdx,which(colnames(cleanData)==colnames(modFrame)[1])] <- ret$repLab  # This is not necessary in HARF because it only removes instances, it does not relabel. However, it must be used when the algorithm relabels instances (in our part there are some of them).
      }
      ret$cleanData <- cleanData[setdiff(1:nrow(cleanData),ret$remIdx),]
      return(ret)
}

#' @rdname PF
#' @export
PF.default <- function(x,
                       nfolds=5,
                       consensus=FALSE,
                       p=0.01,
                       s=3,
                       y=0.5,
                       theta=0.7,
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
      ##Disabling Java output###########
      # sys <- rJava::.jnew("java/lang/System")
      # ops <- sys$out
      # oes <- sys$err
      # nps <- rJava::.jnew("java/io/PrintStream","NUL:")
      # nes <- rJava::.jnew("java/io/PrintStream","NUL:")
      # sys$setOut(nps)
      # sys$setErr(nes)
      ##################################

      origSize <- nrow(x)
      namesOrig <- names(x)
      rownamesOrig <- attr(x,"row.names")
      names(x)[classColumn] <- "class"
      if(any(names(x)[-classColumn]=="class")){
            v <- names(x)[-classColumn]
            v[v=="class"] <- paste("classss",1:sum(v=="class"),sep="")
            names(x)[-classColumn] <- v
      }
      #setting some thresholds and auxiliary variables
      if(consensus){
            majThreshold <- nfolds-1 #the classifier built from the holding fold does not take part in the voting
      }else{
            majThreshold <- floor((nfolds-1)/2)+1
      }
      stopThreshold <- floor(origSize*p)
      KeepOn <- TRUE #will control the while loop
      counter <- 0 #will control how many consecutive times we filter few noisy data
      Dg <- data.frame() #Dg will store good instances
      countIter <- 0 # counts the total iterations

      while(KeepOn){
            countIter <- countIter+1
            VotesGeneral <- vector("integer",nrow(x)) #Counter for total votes
            VotesLocal <- vector("integer",nrow(x)) #Counter for votes from the base learner built using the particular instance

            folds <- caret::createFolds(x$class,nfolds)

            listClassif <- lapply(folds,function(ind){RWeka::PART(class~.,x[ind,])})
            numRules <- sum(sapply(listClassif,function(cl){cl$classifier$measureNumRules()}))
            L <- round(theta*numRules/nfolds) #parameter 'L' for the 'Best-L' scheme

            for(i in 1:nfolds){
                  goodRulesIndexes <- goodRulesSelection(cloningPARTobject(listClassif[[i]]), x[folds[[i]],], L)
                  votes <- votes(cloningPARTobject(listClassif[[i]]), goodRulesIndexes, x, folds[[i]])
                  VotesGeneral <- VotesGeneral+votes$general
                  VotesLocal <- VotesLocal+votes$local
            }
            NoisyIndexes <- which((VotesGeneral >= majThreshold) & (VotesLocal==1))
            GoodIndexes <- which(VotesGeneral==0 & VotesLocal==0)
            AmountGoodInst <- min(length(NoisyIndexes),round(y*length(GoodIndexes))) #good instances removed must be less than noisy (original paper requirement)
            GoodIndexesToKeep <- stratif(GoodIndexes,x$class[GoodIndexes],AmountGoodInst) #stratified sampling so that removed good instances belong to all the classes
            #removing good and noisy instances
            Dg <- rbind(Dg,x[GoodIndexesToKeep,])
            if(length(NoisyIndexes)>0){
                  x <- x[-c(NoisyIndexes,GoodIndexesToKeep),]
            }
            #refreshing stopping conditions
            if( length(NoisyIndexes) <= stopThreshold & counter+1==s) KeepOn <- FALSE
            if( length(NoisyIndexes) <= stopThreshold & counter+1<s) counter <- counter+1
            if(length(NoisyIndexes) > stopThreshold) counter <- 0

            message("Iteration ", countIter,": ", length(NoisyIndexes), " noisy instances removed")
      }
      finalData <- rbind(x,Dg)
      remIdx <- setdiff(1:origSize,as.integer(row.names(finalData)))
      cleanData <- finalData[order(as.integer(row.names(finalData))),]
      names(cleanData) <- namesOrig
      row.names(cleanData) <- rownamesOrig[as.integer(row.names(cleanData))]
      repIdx <- NULL
      repLab <- NULL
      parameters <- list(nfolds=nfolds,
                         consensus=consensus,
                         p=p,
                         s=s,
                         y=y,
                         theta=theta)
      call <- match.call()
      call[[1]] <- as.name("PF")

      ret <- list(cleanData = cleanData,
                  remIdx = remIdx,
                  repIdx=repIdx,
                  repLab=repLab,
                  parameters=parameters,
                  call = call,
                  extraInf = NULL
      )
      class(ret) <- "filter"
      ##Restoring Java output##########
      # sys$setOut(ops)
      # sys$setErr(oes)
      ################################
      return(ret)
}

cloningPARTobject <- function(original){
      cloned <- original
      cloned$classifier <- original$classifier$makeCopy(original$classifier)
      cloned
}

goodRulesSelection <- function(model,trainingData,numGoodRules){
      numRules <- model$classifier$measureNumRules()
      if(numRules <= numGoodRules)
            return(1:numRules)
      accuracyVector <- vector("numeric",numRules)
      coverageVector <- vector("integer",numRules)
      ClassPred2 <- predict(model,trainingData) #initializing coveniently the variable for the right performance of the 'for' loop
      for(i in numRules:1){
            classPred1 <- ClassPred2
            RemoveLastRule(model)
            classPred2 <- predict(model,trainingData)
            coveredIndexes <- which(!is.na(classPred1) & is.na(classPred2))
            correctClassif <- sum(classPred1[coveredIndexes]==trainingData$class[coveredIndexes])
            coverageVector[i] <- length(coveredIndexes)
            accuracyVector[i] <- correctClassif/coverageVector[i]
      }
      order(accuracyVector, coverageVector, decreasing = TRUE)[1:numGoodRules]
}

votes <- function(model,goodRulesIndexes,data,trainingIndexes){
      out <- list(general=vector("integer",nrow(data)),local=vector("integer",nrow(data)))
      numRules <- model$classifier$measureNumRules()
      for(i in numRules:1){
            if(!i%in%goodRulesIndexes){
                  RemoveLastRule(model)
            }else{
                  classPred1 <- predict(model,data)
                  RemoveLastRule(model)
                  classPred2 <- predict(model,data)
                  coveredIndexes <- which(!is.na(classPred1) & is.na(classPred2))
                  localIndexes <- intersect(coveredIndexes,trainingIndexes)
                  generalIndexes <- setdiff(coveredIndexes,localIndexes)
                  out$local[localIndexes] <- out$local[localIndexes]+(classPred1[localIndexes]!=data$class[localIndexes])
                  out$general[generalIndexes] <- out$general[generalIndexes]+(classPred1[generalIndexes]!=data$class[generalIndexes])
            }
      }
      out
}

RemoveLastRule <- function(model){
      modelPART <- model$classifier
      f <- modelPART$getClass()$getDeclaredField("m_root")
      f$setAccessible(TRUE)
      m_root <- f$get(modelPART)
      ff <- m_root$getClass()$getDeclaredField("theRules")
      ff$setAccessible(TRUE)
      theRules <- ff$get(m_root)
      theRules$remove(theRules$size()-1L)
      return(0)
}
