#' TomekLinks
#'
#' Similarity-based filter for removing label noise from a dataset as a
#' preprocessing step of classification. For more information, see 'Details' and
#' 'References' sections.
#'
#' The function \code{TomekLinks} removes "TomekLink points" from the dataset. These are introduced
#' in [Tomek, 1976], and are expected to lie on the border between classes.
#' Removing such points is a typical procedure for cleaning noise [Lorena, 2002].
#'
#' Since the computation of mean points is necessary for TomekLinks, only numeric attributes are allowed.
#' Moreover, only two different classes are allowed to detect TomekLinks.
#'
#' @param formula A formula describing the classification variable and the attributes to be used.
#' @param data,x Data frame containing the tranining dataset to be filtered.
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
#' Tomek I. (Nov. 1976): Two modifications of CNN, \emph{IEEE Trans. Syst., Man, Cybern.}, vol. 6, no. 11, pp. 769-772.
#'
#' Lorena A. C., Batista G. E. A. P. A., de Carvalho A. C. P. L. F., Monard M. C. (Nov. 2002): The influence of noisy patterns in the performance of learning methods in the splice junction recognition problem, in \emph{Proc. 7th Brazilian Symp. Neural Netw.}, Recife, Brazil, pp. 31-37.
#' @examples
#' # Next code fails since TomekLinks method is designed for two-class problems.
#' # Some decomposition strategy like OVO or OVA could be used to overcome this.
#' \dontrun{
#' data(iris)
#' out <- TomekLinks(Species~., data = iris)
#' }
#' @name TomekLinks
NULL

#' @export
TomekLinks <- function(x, ...)
{
      UseMethod("TomekLinks")
}

#' @rdname TomekLinks
#' @export
TomekLinks.formula <- function(formula,
                               data,
                               ...)
{
      if(!is.data.frame(data)){
            stop("data argument must be a data.frame")
      }
      modFrame <- model.frame(formula,data) # modFrame is a data.frame built from 'data' using the variables indicated in 'formula'. The first column of 'modFrame' is the response variable, thus we will indicate 'classColumn=1' when calling the HARF.default method in next line.
      attr(modFrame,"terms") <- NULL

      ret <- TomekLinks.default(x=modFrame,...,classColumn = 1)
      ret$call <- match.call(expand.dots = TRUE)
      ret$call[[1]] <- as.name("TomekLinks")
      # Next, we reconstruct the 'cleanData' from the removed and repaired indexes. Otherwise, the 'cleanData' would only contain those columns passed to the default method (for example imagine when running HARF(Species~Petal.Width+Sepal.Length,iris)).
      cleanData <- data
      if(!is.null(ret$repIdx)){
            cleanData[ret$repIdx,which(colnames(cleanData)==colnames(modFrame)[1])] <- ret$repLab  # This is not necessary in HARF because it only removes instances, it does not relabel. However, it must be used when the algorithm relabels instances (in our part there are some of them).
      }
      ret$cleanData <- cleanData[setdiff(1:nrow(cleanData),ret$remIdx),]
      return(ret)
}

#' @rdname TomekLinks
#' @export
TomekLinks.default <- function(x,
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
      if(nlevels(x[,classColumn])!=2){
            stop("TomekLinks method is designed for two-class problems")
      }
      if(!all(sapply(x,class)[-classColumn] %in% c("integer","numeric"))){
            stop("Numeric attributes are required, since TomekLinks involves computing mean points")
      }

      attributes <- x[,-classColumn]
      attributes <- as.matrix(attributes)

      class1Idx <- which(x[,classColumn]==levels(x[,classColumn])[1])
      class2Idx <- which(x[,classColumn]==levels(x[,classColumn])[2])
      tomekMatrix <- sapply(class1Idx,function(i){
            sapply(class2Idx,function(j){
                  meanPoint <- (attributes[i,]+attributes[j,])/2
                  dist1 <- apply(attributes[setdiff(class1Idx,i),],1,function(v){sum(abs(v-meanPoint))})
                  if(any(dist1<=sum(abs(attributes[i,]-meanPoint)))){
                        return(FALSE)
                  }
                  dist2 <- apply(attributes[setdiff(class2Idx,j),],1,function(v){sum(abs(v-meanPoint))})
                  if(any(dist2<=sum(abs(attributes[j,]-meanPoint)))){
                        return(FALSE)
                  }
                  return(TRUE)
            })
      })
      toRemove <- c(class1Idx[colSums(tomekMatrix)>0],class2Idx[rowSums(tomekMatrix)>0])

      ##### Building the 'filter' object ###########
      remIdx  <- sort(toRemove)
      cleanData <- x[setdiff(1:nrow(x),remIdx),]
      repIdx <- NULL
      repLab <- NULL
      parameters <- NULL
      call <- match.call()
      call[[1]] <- as.name("TomekLinks")

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
