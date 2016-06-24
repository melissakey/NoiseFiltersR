#formu DONE
#library(e1071) #SVM and Naive Bayes
#library(caret) #Multilayer Perceptron Neural Network
#library(RWeka) #C4.5 tree
#library(rpart) #CART tree
#library(randomForest) #Random Forest
#library(kknn) #KNN

#' @title Dynamic Classification Filter
#'
#' @description Ensemble-based filter for removing label noise from a dataset as a
#' preprocessing step of classification. For more information, see 'Details' and
#' 'References' sections.
#'
#' @param formula A formula describing the classification variable and the attributes to be used.
#' @param data,x Data frame containing the tranining dataset to be filtered.
#' @param nfolds Number of folds for the cross voting scheme.
#' @param m Number of classifiers to make up the ensemble. It must range between 1 and 9.
#' @param consensus If set to \code{TRUE}, consensus voting scheme is applied. Otherwise (default),
#' majority scheme is used.
#' @param classColumn Positive integer indicating the column which contains the (factor of) classes.
#' By default, the last column is considered.
#' @param ... Optional parameters to be passed to other methods.
#'
#' @details
#' \code{dynamicCF} (Garcia et al., 2012) follows the same approach as \code{\link{EF}}, but the ensemble of classifiers
#' is not fixed beforehand. Namely, \code{dynamicCF} trains 9 well-known classifiers in the
#' dataset to be filtered, and selects for the ensemble those with the \code{m} best predictions.
#' Then, a \code{nfolds}-folds cross voting scheme is applied, with consensus or majority strategies
#' depending on parameter \code{consensus}.
#'
#' The nine (standard) classifiers handled by \code{dynamicCF} are SVM, 3-KNN, 5-KNN, 9-KNN, CART, C4.5,
#' Random Forest, Naive Bayes and Multilayer Perceptron Neural Network.
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
#' Garcia L. P. F., Lorena A. C., Carvalho A. C. (2012, October): A study on class
#' noise detection and elimination. In \emph{Brazilian Symposium on Neural Networks (SBRN)}, pp. 13-18, IEEE.
#' @examples
#' # Next example is not run in order to save time
#' \dontrun{
#' data(iris)
#' trainData <- iris[c(1:20,51:70,101:120),]
#' # We fix a seed since there exists a random partition for the ensemble
#' set.seed(1)
#' out <- dynamicCF(Species~Petal.Length + Sepal.Length, data = trainData, nfolds = 5, m = 3)
#' summary(out, explicit = TRUE)
#' identical(out$cleanData, trainData[setdiff(1:nrow(trainData),out$remIdx),])
#' }
#' @name dynamicCF
NULL

#' @export
dynamicCF <- function(x, ...)
{
      UseMethod("dynamicCF")
}

#' @rdname dynamicCF
#' @export
dynamicCF.formula <- function(formula,
                              data,
                              ...)
{
      if(!is.data.frame(data)){
            stop("data argument must be a data.frame")
      }
      modFrame <- model.frame(formula,data) # modFrame is a data.frame built from 'data' using the variables indicated in 'formula'. The first column of 'modFrame' is the response variable, thus we will indicate 'classColumn=1' when calling the HARF.default method in next line.
      attr(modFrame,"terms") <- NULL

      ret <- dynamicCF.default(x=modFrame,...,classColumn = 1)
      ret$call <- match.call(expand.dots = TRUE)
      ret$call[[1]] <- as.name("dynamicCF")
      # Next, we reconstruct the 'cleanData' from the removed and repaired indexes. Otherwise, the 'cleanData' would only contain those columns passed to the default method (for example imagine when running HARF(Species~Petal.Width+Sepal.Length,iris)).
      cleanData <- data
      if(!is.null(ret$repIdx)){
            cleanData[ret$repIdx,which(colnames(cleanData)==colnames(modFrame)[1])] <- ret$repLab  # This is not necessary in HARF because it only removes instances, it does not relabel. However, it must be used when the algorithm relabels instances (in our part there are some of them).
      }
      ret$cleanData <- cleanData[setdiff(1:nrow(cleanData),ret$remIdx),]
      return(ret)
}

#' @rdname dynamicCF
#' @export
dynamicCF.default <- function(x,
                              nfolds = 10,
                              consensus = FALSE,
                              m = 3,
                              classColumn = ncol(x),
                              ...)
{
      if(m<1 || m>9){
            stop("parameter m must be set between 1 and 9")
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

      if(consensus){
            threshold <- m
      }else{
            threshold <- floor(m/2)+1
      }

      predictions <- matrix(nrow=nrow(x), ncol=9)
      predictions[,1:9] <- sapply(1:9,function(key){testing(key,x,x,classColumn)})
      ensembleInds <- order(colSums(predictions),decreasing=TRUE)[1:m]
      v <- c("SVM","KNN3","KNN5","KNN9","CART","C4.5","RandomForest","NaiveBayes","MultilayerPerceptron")
      extraInf <- paste(m,"selected classifiers:",paste(v[ensembleInds],collapse=" "))

      folds <- caret::createFolds(x[,classColumn],nfolds)
      votes <- integer(nrow(x))
      for(i in 1:nfolds){
           votes[folds[[i]]] <- m-rowSums(sapply(ensembleInds,function(ind){testing(ind,x[-folds[[i]],],x[folds[[i]],],classColumn)}))
      }
      cleanIndexes <- which(votes < threshold)
      ##### Building the 'filter' object ###########
      cleanData <- x[cleanIndexes,]
      remIdx  <- which(votes>=threshold)
      repIdx <- NULL
      repLab <- NULL
      parameters <- list(nfolds = nfolds,
                         consensus = consensus,
                         m = m)
      call <- match.call()
      call[[1]] <- as.name("dynamicCF")

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

testing <- function(key,trainn,test,classColumn){
      formu <- as.formula(paste(names(trainn)[classColumn],"~.",sep = ""))
      if(key==1){
            out <- predict(e1071::svm(formu,trainn),test)==test[,classColumn]
      }
      if(key%in%c(2,3,4) && !identical(trainn,test)){
            out <- kknn::kknn(formu,train = trainn,test = test,k = getK(key),kernel = "rectangular")$fitted.values==test[,classColumn]
      }
      if(key%in%c(2,3,4) && identical(trainn,test)){
            out <- sapply(1:nrow(trainn),function(i){kknn::kknn(formu,trainn[-i,],test[i,],k=getK(key),kernel = "rectangular")$fitted.values})==test[,classColumn]
      }
      if(key==5){
            out <- predict(rpart::rpart(formu,trainn),test,type="class")==test[,classColumn]
      }
      if(key==6){
            out <- predict(RWeka::J48(formu,trainn),test)==test[,classColumn]
      }
      if(key==7){
            out <- predict(randomForest::randomForest(formu,trainn),test)==test[,classColumn]
      }
      if(key==8){
            out <- predict(e1071::naiveBayes(formu,trainn),test)==test[,classColumn]
      }
      if(key==9){
            rubbish <- utils::capture.output(model <- nnet::nnet(formu,trainn,size=10))
            preds <- predict(model,test,type="class")
            out <- preds==test[,classColumn]
      }
      out
}

getK <- function(k){
      ifelse(k==2,3,ifelse(k==3,5,9))
}
