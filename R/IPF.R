# formu DONE
# library(RWeka)#to use J48 algorithm
# library(caret)#to use 'createFolds' function

#'@title Iterative Partitioning Filter
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
#'after \code{s} iterations with not enough noisy instances removed (according to the proportion \code{p}, see the 'Details' ).
#'@param y Real number between 0 and 1. It sets the proportion of good instances which
#'must be stored in each iteration.
#'@param consensus Logical. If FALSE, majority voting scheme is used. If TRUE, consensus
#'voting scheme is applied.
#'@param classColumn Positive integer indicating the column which contains the (factor of) classes.
#'By default, the last column is considered.
#'@param ... Optional parameters to be passed to other methods.
#'
#'@details The full description of the method can be looked up in the provided references.
#'A base classifier is built in each of the \code{nfolds} partitions of \code{data}. Then, they are
#'tested in the whole dataset, and the removal of noisy instances is decided via consensus or
#'majority voting schemes. Finally, a proportion of good instances (i.e. those whose label agrees
#'with all the base classifiers) is stored and removed for the next iteration. The process stops
#'after \code{s} iterations with not enough (according to the proportion \code{p}) noisy
#'instances removed. In this implementation, the base classifier used is C4.5.
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
#'Khoshgoftaar T. M., Rebours P. (2007): Improving software quality prediction by
#'noise filtering techniques. \emph{Journal of Computer Science and Technology}, 22(3), 387-396.
#'
#'Zhu X., Wu X., Chen Q. (2003, August): Eliminating class noise in large
#'datasets. \emph{International Conference in Machine Learning} (Vol. 3, pp. 920-927).
#'@note
#'By means of a message, the number of noisy instances removed
#'in each iteration is displayed in the console.
#'@examples
#'# Next example is not run in order to save time
#'\dontrun{
#'data(iris)
#'# We fix a seed since there exists a random folds partition for the ensemble
#'set.seed(1)
#'out <- IPF(Species~., data = iris, s = 2)
#'summary(out, explicit = TRUE)
#'identical(out$cleanData, iris[setdiff(1:nrow(iris),out$remIdx),])
#'}
#'@name IPF
NULL

#' @export
IPF <- function(x, ...){
      UseMethod("IPF")
}

#' @rdname IPF
#' @export
IPF.formula <- function(formula,
                         data,
                         ...){
      if(!is.data.frame(data)){
            stop("data argument must be a data.frame")
      }
      modFrame <- model.frame(formula,data) # modFrame is a data.frame built from 'data' using the variables indicated in 'formula'. The first column of 'modFrame' is the response variable, thus we will indicate 'classColumn=1' when calling the HARF.default method in next line.
      attr(modFrame,"terms") <- NULL

      ret <- IPF.default(x=modFrame,...,classColumn = 1)
      ret$call <- match.call(expand.dots = TRUE)
      ret$call[[1]] <- as.name("IPF")
      # Next, we reconstruct the 'cleanData' from the removed and repaired indexes. Otherwise, the 'cleanData' would only contain those columns passed to the default method (for example imagine when running HARF(Species~Petal.Width+Sepal.Length,iris)).
      cleanData <- data
      if(!is.null(ret$repIdx)){
            cleanData[ret$repIdx,which(colnames(cleanData)==colnames(modFrame)[1])] <- ret$repLab  # This is not necessary in HARF because it only removes instances, it does not relabel. However, it must be used when the algorithm relabels instances (in our part there are some of them).
      }
      ret$cleanData <- cleanData[setdiff(1:nrow(cleanData),ret$remIdx),]
      return(ret)
}

#' @rdname IPF
#' @export
IPF.default <- function(x,
                        nfolds=5,
                        consensus=FALSE,
                        p=0.01,
                        s=3,
                        y=0.5,
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

      dataOrig <- x
      origSize <- nrow(x)
      formu <- as.formula(paste(names(x)[classColumn],"~.",sep = ""))
      row.names(x) <- 1:nrow(x)

      #setting some thresholds and auxiliary variables
      if(consensus)
            majThreshold <- nfolds
      else
            majThreshold <- floor(nfolds/2)+1
      stopThreshold <- floor(nrow(x)*p)
      KeepOn <- TRUE #will control the while loop
      counter <- 0 #will control how many consecutive times we filter few nosiy data
      countIter <- 0 # counts the total iterations
      Dg <- data.frame() #Dg will store good instances

      while(KeepOn){
            countIter <- countIter+1
            VotesGeneral <- vector("integer",nrow(x)) #Counter for total votes
            VotesLocal <- vector("integer",nrow(x)) #Counter for votes from the base learner built using the particular instance
            folds <- caret::createFolds(x[,classColumn],nfolds)
            for(i in 1:nfolds){
                  classPreds <- predict(RWeka::J48(formu,x[folds[[i]],]),newdata=x)  #the prediction is made for every instance in data
                  VotesGeneral <- VotesGeneral+(classPreds!=x[,classColumn])
                  VotesLocal[folds[[i]]] <- VotesLocal[folds[[i]]]+(classPreds[folds[[i]]]!=x[,classColumn][folds[[i]]])
            }

            #predictions <- sapply(folds,function(ind){predict(J48(formu,x[ind,]),newdata=x)!=x[,classColumn]})
            #VotesGeneral <- rowSums(predictions)
            #VotesLocal <- rowSums(sapply(1:nfolds,function(i){predictions[-folds[[i]],i] <- FALSE;return(predictions[,i])}))

            #print(VotesGeneral)
            #print(VotesGeneral2)
            #cat("Iteration",counter,":\n",all(VotesLocal2==VotesLocal),"\n",all(VotesGeneral2==VotesGeneral),"\n")

            NoisyIndexes <- which((VotesGeneral >= majThreshold) & (VotesLocal==1))
            GoodIndexes <- which(VotesGeneral==0)
            AmountGoodInst <- min(length(NoisyIndexes),round(y*length(GoodIndexes))) #good instances removed must be less than noisy (original paper requirement)
            GoodIndexesToKeep <- stratif(GoodIndexes,x[,classColumn][GoodIndexes],AmountGoodInst) #stratified sampling so that removed good instances belong to all the classes
            #removing good and noisy instances
            Dg <- rbind(Dg,x[GoodIndexesToKeep,])
            if(length(NoisyIndexes)>0){
                  x <- x[-c(NoisyIndexes,GoodIndexesToKeep),]
            }
            #refreshing stopping conditions
            if( length(NoisyIndexes) <= stopThreshold & counter+1==s) KeepOn <- FALSE
            if( length(NoisyIndexes) <= stopThreshold & counter+1<s) counter <- counter+1
            if(length(NoisyIndexes) > stopThreshold) counter <- 0

            message("Iteration ", countIter,": ", length(NoisyIndexes), " noisy instances removed \n")
      }
      finalData <- rbind(x,Dg)
      ##### Building the 'filter' object ###########
      remIdx  <- setdiff(1:origSize,as.integer(row.names(finalData)))
      cleanData <- dataOrig[setdiff(1:origSize,remIdx),]
      repIdx <- NULL
      repLab <- NULL
      parameters <- list(nfolds=nfolds,
                         consensus=consensus,
                         p=p,
                         s=s,
                         y=y)
      call <- match.call()
      call[[1]] <- as.name("IPF")

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

#next auxiliary function randomly samples 'size' (integer) instances out of the 'values' vector in a stratified manner, i.e. taking
#into account the associated 'stratus' so that all stratus are equally present in the sample.
stratif <- function(values,stratus,size){
      if(size==0){
            return(vector("integer",0))
      }
      stratus <- as.factor(stratus)
      capacity <- summary(stratus) #we store the capacity of each stratus, i.e. the number of instances belonging to it

      subsizes <- splitting(size,capacity) #we obtain the amount of instances that must be sampled within each stratus

      groups <- split(values,stratus)
      unlist(mapply(sample,groups,subsizes))
}

#In next function 'size' is an integer and 'capacity' is a vector of integers. It returns how to equally distribute
#"size" elements in length(capacity) groups, taking into account the maximum capacity of each group (which is provided by the 'capacity' vector).
splitting <- function(size,capacity){
      if(size>sum(capacity)){
            stop("avoided an infinite while loop")
      }
      n <- length(capacity)
      out <- vector("integer",n)
      counter <- 1
      for(i in 1:size){
            KeepOn <- TRUE
            while(KeepOn){
                  if(out[modulus(counter,n)]<capacity[modulus(counter,n)]){
                        out[modulus(counter,n)] <- out[modulus(counter,n)]+1
                        KeepOn <- FALSE
                  }
            counter <- counter+1
            }
      }
      out
}

modulus <- function(D,d){
      out <- D%%d
      if(out==0)
            out <- d
      out
}
