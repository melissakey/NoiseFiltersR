#' Mode Filter
#'
#' Similarity-based filter for removing or repairing label noise from a dataset as a
#' preprocessing step of classification. For more information, see 'Details' and
#' 'References' sections.
#'
#' \code{ModeFilter} estimates the most appropriate class for each instance based on the similarity metric
#' and the provided label. This can be addressed in three different ways (argument 'type'):
#'
#' In the classical approach, all labels are tried for all instances, and the one maximizing a metric
#' based on similarity is chosen. In the iterative approach, the same scheme is repeated until the proportion
#' of modified instances is less than \emph{epsilon} or the maximum number of iterations \emph{maxIter}
#' is reached. The weighted approach extends the classical one by assigning a weight for each instance, which
#' quantifies the reliability on its label. This weights is utilized in the computation of the metric to be
#' maximized.
#'
#' @param formula A formula describing the classification variable and the attributes to be used.
#' @param data,x Data frame containing the tranining dataset to be filtered.
#' @param type Character indicating the scheme to be used. It can be 'classical', 'iterative' or 'weighted'.
#' @param noiseAction Character indicating what to do with noisy instances. It can be either 'remove' or 'repair'.
#' @param epsilon If 'iterative' type is used, the loop will be stopped if the proportion of modified instances
#' is less or equal than this threshold.
#' @param maxIter Maximum number of iterations in 'iterative' type.
#' @param alpha Parameter used in the computation of the similarity between two instances.
#' @param beta It regulates the influence of the similarity metric in the estimation
#' of a new label for an instance.
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
#'
#' @references
#' Du W., Urahama K. (2010, November): Error-correcting semi-supervised pattern
#' recognition with mode filter on graphs.
#' In \emph{Aware Computing (ISAC), 2010 2nd International Symposium on} (pp. 6-11). IEEE.
#' @examples
#' # Next example is not run because in some cases it can be rather slow
#' \dontrun{
#'    data(iris)
#'    out <- ModeFilter(Species~., data = iris, type = "classical", noiseAction = "remove")
#'    print(out)
#'    identical(out$cleanData, iris[setdiff(1:nrow(iris),out$remIdx),])
#' }
#' @name ModeFilter
NULL

#' @export
ModeFilter <- function(x, ...)
{
      UseMethod("ModeFilter")
}

#' @rdname ModeFilter
#' @export
ModeFilter.formula <- function(formula,
                               data,
                               ...)
{
      if(!is.data.frame(data)){
            stop("data argument must be a data.frame")
      }
      modFrame <- model.frame(formula,data) # modFrame is a data.frame built from 'data' using the variables indicated in 'formula'. The first column of 'modFrame' is the response variable, thus we will indicate 'classColumn=1' when calling the HARF.default method in next line.
      attr(modFrame,"terms") <- NULL

      ret <- ModeFilter.default(x=modFrame,...,classColumn = 1)
      ret$call <- match.call(expand.dots = TRUE)
      ret$call[[1]] <- as.name("ModeFilter")
      # Next, we reconstruct the 'cleanData' from the removed and repaired indexes. Otherwise, the 'cleanData' would only contain those columns passed to the default method (for example imagine when running HARF(Species~Petal.Width+Sepal.Length,iris)).
      cleanData <- data
      if(!is.null(ret$repIdx)){
            cleanData[ret$repIdx,which(colnames(cleanData)==colnames(modFrame)[1])] <- ret$repLab  # This is not necessary in HARF because it only removes instances, it does not relabel. However, it must be used when the algorithm relabels instances (in our part there are some of them).
      }
      ret$cleanData <- cleanData[setdiff(1:nrow(cleanData),ret$remIdx),]
      return(ret)
}

#' @rdname ModeFilter
#' @export
ModeFilter.default <- function(x,
                               type = "classical",
                               noiseAction = "repair",
                               epsilon = 0.05,
                               maxIter = 100,
                               alpha = 1,
                               beta = 1,
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
      if(!type%in%c("classical","iterative","weighted")){
            stop("the argument 'type' must be set to 'classical', 'iterative' or 'weighted'")
      }
      if(!noiseAction%in%c("repair","remove")){
            stop("the argument 'noiseAction' must be set to 'repair' or 'remove'")
      }
      if(epsilon>1 | epsilon<0){
            stop("argument 'epsilon' must range between 0 and 1")
      }

  
      similarity <- matrix(NA, ncol = nrow(x), nrow = nrow(x))
      similarity <- sapply(1:nrow(x),function(i){
        c(rep(NA,i-1),
          sapply(i:nrow(x),
            function(j){exp(-alpha*distt(x[i,-classColumn],x[j,-classColumn])^2)}))})
      similarity <- similarity + t(similarity)
      labels <- levels(x[,classColumn])


      if(type=="classical"){
            newClass <- sapply(1:nrow(x),function(i){
                  sumsPerClass <- sapply(labels,function(label){
                        sum(similarity[x[,classColumn]==label,i]) * exp(-beta)*sum(similarity[x[,classColumn]!=label,i])
                  })
                  labels[nnet::which.is.max(sumsPerClass)]
            })
      }
      if(type=="iterative"){
            currentClass <- as.character(x[,classColumn])
            k <- 1
            convergenceRate <- 1.1
            while(k <= maxIter & convergenceRate > epsilon){
                  oldClass <- currentClass
                  currentClass <- sapply(1:nrow(x),function(i){
                        sumsPerClass <- sapply(labels,function(label){
                              sum(similarity[currentClass==label,i]) * exp(-beta)*sum(similarity[currentClass!=label,i])
                        })
                        labels[nnet::which.is.max(sumsPerClass)]
                  })
                  k <- k+1
                  convergenceRate <- sum(currentClass!=oldClass)/length(currentClass)
            }
            newClass <- currentClass
      }
      if(type=="weighted"){
            weights <- sapply(1:nrow(x),function(i){
                  (sum(similarity[i,x[,classColumn]==x[i,classColumn]]) * exp(-beta)*sum(similarity[i,x[,classColumn]!=x[i,classColumn]]))/sum(similarity[i,])
            })
            newClass <- sapply(1:nrow(x),function(i){
                  sumsPerClass <- sapply(labels,function(label){
                        indexEq <- x[,classColumn]==label
                        sum(similarity[indexEq,i]*weights[indexEq])+exp(-beta)*sum(similarity[!indexEq,i]*weights[!indexEq])
                  })
                  labels[nnet::which.is.max(sumsPerClass)]
            })
      }

      if(noiseAction=="remove"){
            remIdx <- which(newClass!=x[,classColumn])
            repIdx <- NULL
            repLab <- NULL
            cleanData <- x[setdiff(1:nrow(x),remIdx),]
      }
      if(noiseAction=="repair"){
            remIdx <- NULL
            repIdx <- which(newClass!=x[,classColumn])
            repLab <- factor(newClass[repIdx], levels = labels)
            cleanData <- x
            cleanData[,classColumn] <- newClass
      }

      ##### Building the 'filter' object ###########
      parameters <- list(type = type,
                         noiseAction = noiseAction,
                         maxIter = maxIter,
                         alpha = alpha,
                         beta = beta)
      call <- match.call()
      call[[1]] <- as.name("ModeFilter")

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
