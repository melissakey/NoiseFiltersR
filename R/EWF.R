#' Edge Weight Filter
#'
#' Similarity-based filter for removing or repairing label noise from a dataset as a
#' preprocessing step of classification. For more information, see 'Details' and
#' 'References' sections.
#'
#' \code{EWF} builds up a Relative Neighborhood Graph (RNG) from the dataset. Then, it identifies
#' as 'suspicious' those instances with a significant value of its\emph{local cut edge weight statistic}, which
#' intuitively means that they are surrounded by examples from a different class.
#'
#' Namely, the aforementioned statistic is the sum of the weights of edges joining
#' the instance (in the RNG graph) with instances from a different class.
#' Under the null hypothesis of the class label being independent of
#' the event 'being neighbors in the RNG graph', the distribution of this statistic can be approximated by a
#' gaussian one. Then, the p-value for the observed value is computed and contrasted with the
#' provided \code{threshold}.
#'
#' To handle 'suspicious' instances there are two approaches ('remove' or 'hybrid'), and the argument
#' 'noiseAction' determines which one to use. With 'remove', every suspect is removed from the dataset.
#' With the 'hybrid' approach, an instance is removed if it does not have \emph{good} (i.e. non-suspicious)
#' RNG-neighbors. Otherwise, it is relabelled with the majority class among its \emph{good} RNG-neighbors.
#'
#' @param formula A formula describing the classification variable and the attributes to be used.
#' @param data,x Data frame containing the tranining dataset to be filtered.
#' @param threshold Real number between 0 and 1. It sets the limit between good and suspicious instances. Its
#' default value is 0.25.
#' @param noiseAction Character being either 'remove' or 'hybrid'. It determines what to do with noisy
#' instances. By default, it is set to 'remove'.
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
#' @references
#' Muhlenbach F., Lallich S., Zighed D. A. (2004): Identifying and handling mislabelled
#' instances. \emph{Journal of Intelligent Information Systems}, 22(1), 89-109.
#' @examples
#' # Next example is not run because EWF is time-consuming
#' \dontrun{
#'    data(iris)
#'    trainData <- iris[c(1:20,51:70,101:120),]
#'    out <- EWF(Species~Petal.Length+Sepal.Length, data = trainData, noiseAction = "hybrid")
#'    print(out)
#' }
#' @name EWF
NULL

#' @export
EWF <- function(x, ...)
{
      UseMethod("EWF")
}

#' @rdname EWF
#' @export
EWF.formula <- function(formula,
                        data,
                        ...)
{
      if(!is.data.frame(data)){
            stop("data argument must be a data.frame")
      }
      modFrame <- model.frame(formula,data) # modFrame is a data.frame built from 'data' using the variables indicated in 'formula'. The first column of 'modFrame' is the response variable, thus we will indicate 'classColumn=1' when calling the HARF.default method in next line.
      attr(modFrame,"terms") <- NULL

      ret <- EWF.default(x=modFrame,...,classColumn = 1)
      ret$call <- match.call(expand.dots = TRUE)
      ret$call[[1]] <- as.name("EWF")
      # Next, we reconstruct the 'cleanData' from the removed and repaired indexes. Otherwise, the 'cleanData' would only contain those columns passed to the default method (for example imagine when running HARF(Species~Petal.Width+Sepal.Length,iris)).
      cleanData <- data
      if(!is.null(ret$repIdx)){
            cleanData[ret$repIdx,which(colnames(cleanData)==colnames(modFrame)[1])] <- ret$repLab  # This is not necessary in HARF because it only removes instances, it does not relabel. However, it must be used when the algorithm relabels instances (in our part there are some of them).
      }
      ret$cleanData <- cleanData[setdiff(1:nrow(cleanData),ret$remIdx),]
      return(ret)
}

#' @rdname EWF
#' @export
EWF.default <- function(x,
                        threshold = 0.25,
                        noiseAction = "remove",
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
      if(any(sapply(x,class)[-classColumn]=="factor")){
            stop("all attributes must be numerical or logical (0-1). Categorical variables should be
                  numerically coded in a manner suited to the problem")
      }
      if(!noiseAction%in%c("remove","hybrid")){
            stop("noiseAction argument must be either 'remove' or 'hybrid'")
      }


      #Standardize variables
      dataS <- data.frame(lapply(1:ncol(x),function(i){
            if(i==classColumn){
                  x[,i]
            }
            else{
                  x <- x[,i]
                  (x-mean(x))/stats::sd(x)
            }
      }))
      names(dataS) <- names(x)

      
      # Creating Relative Neighborhood Graph
      distMat <- as.matrix(
        dist(
          as.matrix(dataS[setdiff(which(sapply(dataS, class) == 'numeric'), classColumn)])
        )
      )
      graph <- isNeighbor(distMat, pmax)
      

      # Using cut edge weight statistic to indentify suspects
      isSuspect <- sapply(1:nrow(dataS),function(i){
            propClass <- sum(dataS[,classColumn]==dataS[i,classColumn])/nrow(dataS)
            neigh <- which(graph[i,])
            weights <- sapply(neigh,function(j){1/(1+distt(dataS[i,-classColumn],dataS[j,-classColumn]))})
            statisticValue <- sum(weights[dataS[neigh,classColumn]!=dataS[i,classColumn]])
            statisticMean <- (1-propClass)*sum(weights)
            statisticSD <- (1-propClass)*propClass*sum(weights^2)
            pvalue <- stats::pnorm(statisticValue,statisticMean,statisticMean)
            !(pvalue<threshold)
      })

      # Handling suspects
      if(noiseAction=="remove"){
            remIdx <- which(isSuspect)
            repIdx <- NULL
            repLab <- NULL
      }else{
            newLabel <- sapply(which(isSuspect),function(i){
                  goodNeigh <- !isSuspect & graph[i,]
                  if(any(goodNeigh)){
                        labels <- table(x[goodNeigh,classColumn])
                        return(names(labels)[nnet::which.is.max(labels)])
                  }
                  else{
                        return(NA)
                  }
            })
            remIdx <- which(isSuspect)[is.na(newLabel)]
            repIdx <- which(isSuspect)[!is.na(newLabel) & newLabel!=dataS[which(isSuspect),classColumn]]
            repLab <- newLabel[!is.na(newLabel) & newLabel!=dataS[which(isSuspect),classColumn]]
            repLab <- factor(repLab, levels = levels(x[,classColumn]))
      }

      ##### Building the 'filter' object ###########
      cleanData <- x
      if(!is.null(repIdx)){
            cleanData[repIdx,classColumn] <- repLab
      }
      cleanData <- cleanData[setdiff(1:nrow(x),remIdx),]
      parameters <- list(threshold = threshold,
                         noiseAction = noiseAction)
      call <- match.call()
      call[[1]] <- as.name("EWF")

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
