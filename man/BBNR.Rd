% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BBNR.R
\name{BBNR}
\alias{BBNR}
\alias{BBNR.default}
\alias{BBNR.formula}
\title{Blame Based Noise Reduction}
\usage{
\method{BBNR}{formula}(formula, data, ...)

\method{BBNR}{default}(x, k = 3, classColumn = ncol(x), ...)
}
\arguments{
\item{formula}{A formula describing the classification variable and the attributes to be used.}

\item{data, x}{Data frame containing the tranining dataset to be filtered.}

\item{...}{Optional parameters to be passed to other methods.}

\item{k}{Number of nearest neighbors to be used.}

\item{classColumn}{positive integer indicating the column which contains the
(factor of) classes. By default, the last column is considered.}
}
\value{
An object of class \code{filter}, which is a list with seven components:
\itemize{
   \item \code{cleanData} is a data frame containing the filtered dataset.
   \item \code{remIdx} is a vector of integers indicating the indexes for
   removed instances (i.e. their row number with respect to the original data frame).
   \item \code{repIdx} is a vector of integers indicating the indexes for
   repaired/relabelled instances (i.e. their row number with respect to the original data frame).
   \item \code{repLab} is a factor containing the new labels for repaired instances.
   \item \code{parameters} is a list containing the argument values.
   \item \code{call} contains the original call to the filter.
   \item \code{extraInf} is a character that includes additional interesting
   information not covered by previous items.
}
}
\description{
Similarity-based filter for removing label noise from a dataset as a
preprocessing step of classification. For more information, see 'Details' and
'References' sections.
}
\details{
\code{BBNR} removes an instance 'X' if: (i) it participates in the misclassification of other instance
(i.e. 'X' is among the \code{k} nearest neighbors of a misclassified instance and has a different class);
and (ii) its removal does not produce a misclassification in instances that, initially, were correctly
classified by 'X' (i.e. 'X' was initially among the \code{k} nearest neighbors and had the same class).
}
\examples{
# Next example is not run in order to save time
\dontrun{
data(iris)
out <- BBNR(iris, k = 5)
print(out)
identical(out$cleanData, iris[setdiff(1:nrow(iris),out$remIdx),])
}
}
\references{
Delany S. J., Cunningham P. (2004): An analysis of case-base editing in a spam filtering system.
In \emph{Advances in Case-Based Reasoning} (pp. 128-141). Springer Berlin Heidelberg.
}

