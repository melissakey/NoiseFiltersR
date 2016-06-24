#' @export
print.filter <-
      function(x, ...)
      {
            cat("Call:", deparse(x$call, 0.8 * getOption("width")),sep = "\n")
            if(!is.null(x$parameters)){
                  cat("\nParameters:\n")
                  for (i in 1:length(x$parameters)) {
                        cat(attr(x$parameters,"names")[i],": ",x$parameters[[i]],"\n",sep =
                                  "")
                  }
            }
            cat("\nResults:\n")
            cat("Number of removed instances: ",length(x$remIdx)," (",length(x$remIdx)*100/(length(x$remIdx)+nrow(x$cleanData))," %)","\n",sep="")
            cat("Number of repaired instances: ",length(x$repIdx)," (",length(x$repIdx)*100/(length(x$remIdx)+nrow(x$cleanData))," %)","\n",sep="")
      }


#' Summary method for class filter
#'
#' This methods allows for appropriately displaying the most important
#' information about a filtered dataset, contained in the S3 class \code{filter}.
#'
#' The information offered is the following:
#' \itemize{
#'    \item Names of the dataset and the filter.
#'    \item Original call to the filter.
#'    \item Specific parameters used for the filter.
#'    \item Results: number of removed and repaired instances (absolute
#'    number and percentage).
#'    \item Additional information (if available, it depends on the filter).
#'    \item Optionally, if \code{explicit=TRUE}, the indexes for removed and
#'    repaired instances, as well as the new labels.
#' }
#'
#' @param object Object of class \code{filter}.
#' @param explicit If set to \code{TRUE}, the indexes for removed and repaired
#' instances (as well as new labels for the latters) are displayed. It defaults
#' to \code{FALSE}.
#' @param ... Additional arguments affecting the summary produced.
#' @examples
#' # Next example is not run in order to save time
#' \dontrun{
#' # Example of filter with additional information available.
#' data(iris)
#' out <- edgeBoostFilter(Species~., data = iris)
#' class(out)
#' summary(out)
#' summary(out, explicit = TRUE)
#' }
#' @export
summary.filter <- function(object, ..., explicit = FALSE)
      {
            if (explicit) {
                  structure(object,class = "summary1filter")
            }else{
                  structure(object,class = "summary2filter")
            }
      }

#' @export
print.summary2filter <- function(x, ...)
      {
            cat("Filter",x$call[[1]],"applied to dataset",x$call$data,"\n\n")
            print.filter(x)
            if(!is.null(x$extraInf)){
                  cat("\nAdditional information:\n")
                  cat(x$extraInf,"\n")
            }
      }

#' @export
print.summary1filter <- function(x, ...)
      {
            cat("Filter",x$call[[1]],"applied to dataset",x$call$data,"\n\n")
            print.filter(x)
            if(!is.null(x$extraInf)){
                  cat("\nAdditional information:\n")
                  cat(x$extraInf,"\n")
            }
            cat("\nExplicit indexes for removed instances:\n")
            cat(x$remIdx,"\n\n")
            if (length(x$repIdx) > 0) {
                  cat("Explicit indexes for repaired instances:\n")
                  cat(x$repIdx,"\n\n")
                  cat("New labels for repaired instances:\n")
                  cat(as.character(x$repLab),"\n\n")
            }
      }
