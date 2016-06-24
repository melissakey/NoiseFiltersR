# formu DONE
# library(nnet) #Neural Network, which.is,max (ties at random)
# library(rpart) #CART decision tree
# library(e1071) #SVM
# library(kknn) #KNN

#'@title Hybrid Repair-Remove Filter
#'
#'@description Ensemble-based filter for removing or repairing label noise from a dataset as a
#' preprocessing step of classification. For more information, see 'Details' and
#' 'References' sections.
#'
#'@param formula A formula describing the classification variable and the attributes to be used.
#'@param data,x Data frame containing the tranining dataset to be processed.
#'@param consensus If set to \code{TRUE}, consensus voting scheme is applied to identify noisy instances. Otherwise (default),
#'majority approach is used.
#'@param noiseAction Character which can be set to "remove", "repair" or "hybrid". The filter accordingly decides
#'what to do with the identified noise (see Details).
#'@param classColumn Positive integer indicating the column which contains the (factor of) classes.
#'By default, the last column is considered.
#'@param ... Optional parameters to be passed to other methods.
#'
#'@details
#'As presented in (Miranda et al., 2009), \code{hybridRepairFilter} builds on the dataset an ensemble of four
#'classifiers: SVM, Neural Network, CART, KNN (combining k=1,3,5). According to their predictions and
#'majority or consensus voting schemes, a
#'subset of instances are labeled as noise. These are removed if \code{noiseAction} equals "remove", their class
#'is changed into the most voted among the ensemble if \code{noiseAction} equals "repair", and when the latter
#'is set to "hybrid", the vote of KNN decides whether remove or repair.
#'
#'All this procedure is repeated while the accuracy (over the original dataset) of the ensemble
#'trained with the processed dataset increases.
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
#'Miranda A. L., Garcia L. P. F., Carvalho A. C., Lorena A. C. (2009): Use of
#'classification algorithms in noise detection and elimination. In \emph{Hybrid Artificial
#'Intelligence Systems} (pp. 417-424). Springer Berlin Heidelberg.
#'@examples
#'# Next example is not run in order to save time
#'\dontrun{
#'data(iris)
#'out <- hybridRepairFilter(iris, noiseAction = "hybrid")
#'summary(out, explicit = TRUE)
#'}
#'@name hybridRepairFilter
NULL

#' @export
hybridRepairFilter <- function(x, ...)
{
      UseMethod("hybridRepairFilter")
}

#' @rdname hybridRepairFilter
#' @export
hybridRepairFilter.formula <- function(formula,
                                       data,
                                       ...)
{
      if(!is.data.frame(data)){
            stop("data argument must be a data.frame")
      }
      modFrame <- model.frame(formula,data) # modFrame is a data.frame built from 'data' using the variables indicated in 'formula'. The first column of 'modFrame' is the response variable, thus we will indicate 'classColumn=1' when calling the HARF.default method in next line.
      attr(modFrame,"terms") <- NULL

      ret <- hybridRepairFilter.default(x=modFrame,...,classColumn = 1)
      ret$call <- match.call(expand.dots = TRUE)
      ret$call[[1]] <- as.name("hybridRepairFilter")
      # Next, we reconstruct the 'cleanData' from the removed and repaired indexes. Otherwise, the 'cleanData' would only contain those columns passed to the default method (for example imagine when running HARF(Species~Petal.Width+Sepal.Length,iris)).
      cleanData <- data
      if(!is.null(ret$repIdx)){
            cleanData[ret$repIdx,which(colnames(cleanData)==colnames(modFrame)[1])] <- ret$repLab  # This is not necessary in HARF because it only removes instances, it does not relabel. However, it must be used when the algorithm relabels instances (in our part there are some of them).
      }
      ret$cleanData <- cleanData[setdiff(1:nrow(cleanData),ret$remIdx),]
      return(ret)
}

#' @rdname hybridRepairFilter
#' @export
hybridRepairFilter.default <- function(x,
                                       consensus = FALSE,
                                       noiseAction = "remove",
                                       classColumn = ncol(x),
                                       ...)
{
      if(!noiseAction %in% c("remove","repair","hybrid")){
            stop("parameter 'noiseAction' must exactly match one of 'remove', 'repair' or 'hybrid'.")
      }
      if(!is.data.frame(x)){
            stop("data argument must be a data.frame")
      }
      if(!classColumn%in%(1:ncol(x))){
            stop("class column out of range")
      }
      if(!is.factor(x[,classColumn])){
            stop("class column of data must be a factor")
      }

      originalData <- x
      formu <- as.formula(paste(names(x)[classColumn],"~.",sep = ""))
      row.names(x) <- 1:nrow(x)

      threshold <- ifelse(consensus,4,2)
      KeepOn <- TRUE
      correctClassifs <- 0
      counter <- 1
      while(KeepOn){
            preds <- predictions(trainn=x,test=x,formu)
            votesNoise <- rowSums(apply(preds,2,function(v){v!=x[,classColumn]}))
            noiseInds <- which(votesNoise>=threshold)

            if(noiseAction=="remove"){
                  dataTemp <- x[setdiff(1:nrow(x),noiseInds),]
            }
            if(noiseAction=="repair"){
                  dataTemp <- x
                  dataTemp[noiseInds,classColumn] <- apply(preds[noiseInds,,drop=FALSE],1,function(v){names(table(v))[nnet::which.is.max(table(v))]})
            }
            if(noiseAction=="hybrid"){
                  toRemoveInds <- noiseInds[which(preds[noiseInds,4]==x[,classColumn][noiseInds])]
                  toRepairInds <- setdiff(noiseInds,toRemoveInds)
                  dataTemp <- x
                  dataTemp[toRepairInds,classColumn] <- apply(preds[toRepairInds,,drop=FALSE],1,function(v){names(table(v))[nnet::which.is.max(table(v))]})
                  dataTemp <- dataTemp[setdiff(1:nrow(x),toRemoveInds),]
            }

            predsOriginal <- predictions(trainn=dataTemp,test=originalData,formu)
            correctClassifsNew <- sum(apply(predsOriginal,2,function(v){v==originalData[,classColumn]}))
            if(correctClassifsNew<=correctClassifs){
                  KeepOn <- FALSE
            }else{
                  correctClassifs <- correctClassifsNew
                  counter <- counter+1
                  x <- dataTemp
            }
      }
      ##### Building the 'filter' object ###########
      cleanData <- x
      names(cleanData) <- names(originalData)
      row.names(cleanData) <- attr(originalData,"row.names")[as.integer(row.names(x))]
      remIdx <- setdiff(1:nrow(originalData),as.integer(row.names(x)))
      isRepaired <- x[,classColumn]!=originalData[as.integer(row.names(x)),classColumn]
      repIdx <- as.integer(row.names(x))[isRepaired]
      repLab <- x[,classColumn][isRepaired]
      if(length(remIdx)==0){
            remIdx <- NULL
      }
      if(length(repIdx)==0){
            repIdx <- NULL
            repLab <- NULL
      }
      extraInf <- paste("The number of iterations was",counter,collapse = " ")
      parameters <- list(consensus = consensus,
                         noiseAction = noiseAction)
      call <- match.call()
      call[[1]] <- as.name("hybridRepairFilter")

      ret <- list(cleanData = cleanData,
                  remIdx = remIdx,
                  repIdx=repIdx,
                  repLab=repLab,
                  parameters=parameters,
                  call = call,
                  extraInf=extraInf
      )
      class(ret) <- "filter"
      return(ret)
}

predictions <- function(trainn,test,formu){
      out <- matrix(nrow = nrow(test), ncol = 4)
      out[,1] <- as.character(predict(e1071::svm(formu,trainn),test))
      rubbish <- utils::capture.output(model <- nnet::nnet(formu,trainn,size=10))
      out[,2] <- as.character(predict(model,test,type="class"))
      out[,3] <- as.character(predict(rpart::rpart(formu,trainn),test,type="class"))

      knnPreds <- matrix(nrow = nrow(test),ncol = 3)
      knnPreds[,c(1,2,3)] <- sapply(c(1,3,5),function(k){
            if(dim(trainn)==dim(test) && all(trainn==test)){
                  sapply(1:nrow(trainn),function(i){
                        kknn::kknn(formula = formu,
                                   train = trainn[-i,],
                                   test = test[i,],
                                   k = k,
                                   kernel = "rectangular")$fitted.values
                        })
                  }
            else{
                  kknn::kknn(formula = formu,
                             train = trainn,
                             test = test,
                             k = k,
                             kernel = "rectangular")$fitted.values
            }
            })

      out[,4] <- apply(knnPreds,1,function(v){names(table(v))[nnet::which.is.max(table(v))]})
      out
}
