#' Editing with Neighbor Graphs
#'
#' Similarity-based filter for removing label noise from a dataset as a
#' preprocessing step of classification. For more information, see 'Details' and
#' 'References' sections.
#'
#' \code{ENG} builds a neighborhood graph which can be either \emph{Gabriel Graph (GG)} or
#' \emph{Relative Neighborhood Graph (RNG)} [S\'{a}nchez et al., 1997]. Then, an
#' instance is considered as 'potentially noisy' if most of its neighbors have a different class. To decide whether such
#' an instance 'X' is removed, let S be the subset given by 'X' together with its neighbors from the same class.
#' Compute the majority class 'C' among the neighbors of examples in S, and remove 'X' if its class is not 'C'.
#'
#' @param formula A formula describing the classification variable and the attributes to be used.
#' @param data,x Data frame containing the tranining dataset to be filtered.
#' @param graph Character indicating the type of graph to be constructed. It can be chosen between 'GG'
#' (Gabriel Graph) and 'RNG' (Relative Neighborhood Graph). See 'References' for more details on both graphs.
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
#' S\'{a}nchez J. S., Pla F., Ferri F. J. (1997): Prototype selection for the nearest
#' neighbour rule through proximity graphs. \emph{Pattern Recognition Letters}, 18(6), 507-513.
#' @examples
#' # The example is not run because the graph construction is quite time-consuming.
#' \dontrun{
#'    data(iris)
#'    trainData <- iris[c(1:20,51:70,101:120),]
#'    out <- ENG(Species~Petal.Length + Petal.Width, data = trainData, graph = "RNG")
#'    print(out)
#'    identical(out$cleanData,trainData[setdiff(1:nrow(trainData),out$remIdx),])
#' }
#' @name ENG
NULL

#' @export
ENG <- function(x, ...)
{
  UseMethod("ENG")
}

#' @rdname ENG
#' @export
ENG.formula <- function(formula,
  data,
  ...)
{
  if(!is.data.frame(data)){
    stop("data argument must be a data.frame")
  }
  modFrame <- model.frame(formula,data) # modFrame is a data.frame built from 'data' using the variables indicated in 'formula'. The first column of 'modFrame' is the response variable, thus we will indicate 'classColumn=1' when calling the HARF.default method in next line.
  attr(modFrame,"terms") <- NULL
  
  ret <- ENG.default(x=modFrame,...,classColumn = 1)
  ret$call <- match.call(expand.dots = TRUE)
  ret$call[[1]] <- as.name("ENG")
  # Next, we reconstruct the 'cleanData' from the removed and repaired indexes. Otherwise, the 'cleanData' would only contain those columns passed to the default method (for example imagine when running HARF(Species~Petal.Width+Sepal.Length,iris)).
  cleanData <- data
  if(!is.null(ret$repIdx)){
    cleanData[ret$repIdx,which(colnames(cleanData)==colnames(modFrame)[1])] <- ret$repLab  # This is not necessary in HARF because it only removes instances, it does not relabel. However, it must be used when the algorithm relabels instances (in our part there are some of them).
  }
  ret$cleanData <- cleanData[setdiff(1:nrow(cleanData),ret$remIdx),]
  return(ret)
}

#' @rdname ENG
#' @export
ENG.default <- function(x,
  graph = "RNG",
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
  if(!graph%in%c("RNG","GG")){
    stop("the 'graph' argument must be either 'GG' (Gabriel Graph) or 'RNG' (Relative Neighborhood Graph)")
  }
  
  distMat <- as.matrix(
    dist(
      as.matrix(x[setdiff(which(sapply(x, class) == 'numeric'), classColumn)])
    )
  )
  fun <- switch(graph,
    RNG = pmax,
    GG = function(x, y) x^2 + y^2
  )
  PG <- isNeighbor(distMat, fun)
  
  #First order graph
  isMisclassified <- sapply(1:nrow(x),function(i){
    classes <- table(x[PG[i,],classColumn])
    if(names(classes)[nnet::which.is.max(classes)]==x[i,classColumn]){
      out <- FALSE
    }else{
      out <- TRUE
    }
    return(out)
  })
  
  #Second order graph
  toRemove <- sapply(which(isMisclassified),function(i){
    sameClassNeigh <- c(i,which(PG[i,] & x[,classColumn]==x[i,classColumn]))
    classes <- character(0)
    for(j in sameClassNeigh){
      classes <- c(classes,as.character(x[PG[j,],classColumn]))
    }
    tableClasses <- table(classes)
    if(names(tableClasses)[nnet::which.is.max(tableClasses)]==x[i,classColumn]){
      out <- FALSE
    }
    else{
      out <- TRUE
    }
    return(out)
  })
  
  ##### Building the 'filter' object ###########
  remIdx  <- which(isMisclassified)[toRemove]
  cleanData <- x[setdiff(1:nrow(x),remIdx),]
  repIdx <- NULL
  repLab <- NULL
  parameters <- list(graph=graph)
  call <- match.call()
  call[[1]] <- as.name("ENG")
  
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