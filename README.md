## Description

*NoiseFiltersR* contains an extensive implementation of state-of-the-art and classical label noise
preprocessing algorithms for classification problems. Such a collection was missing for R statistical software.

Namely, *NoiseFiltersR* includes 30 label noise filters. All of them are appropriately
documented, with a general explanation of the method and the exact reference where it was first published.
Moreover, they can be called in a R-user-friendly manner, and their results are unified by means of the `filter` class, which also benefits from adapted `print` and `summary` methods.


## Installation

Use `install.packages` to install *NoiseFiltersR* and its dependencies from CRAN:
```R
install.packages("NoiseFiltersR")
```

Once installed, use the command `library` to attach the package:
```R
library("NoiseFiltersR")
```

## Example of use

Once the package is installed and attached, the user can apply any of the implemented algorithms.

Next instruction shows how to use the well-known *Iterative Partitioning Filter (IPF)* (Khoshgoftaar & Rebours, 2007) to
filter out class noise from the dataset `iris`. The *formula* allows us to indicate the classification
variable. Default parameters for the algorithm are considered:

```R
out <- IPF(Species~., data = iris)
```
Then, the variable `out` is an object of class `filter`. This is a list with seven elements:

* `cleanData`: a data frame containing the filtered dataset.
* `remIdx`: a vector of integers indicating the indexes for
removed instances (i.e. their row number with respect to the original data frame).
* `repIdx`: a vector of integers indicating the indexes for
repaired/relabelled instances (i.e. their row number with respect to the original data frame).
* `repLab`: a factor containing the new labels for repaired instances.
* `parameters`: a list containing the tuning parameters used for the filter.
* `call`: an expression that contains the original call to the filter.
* `extraInf`: a character that includes additional relevant
information not covered by previous items.

To appropriately display the information contained in a `filter` object, general functions
`print` and `summary` can be used (more details about their output can be found in the package vignette):

```R
print(out)
summary(out)
```

Finally, all the implemented algorithms can also be used without a *formula* argument, just
indicating the dataset to be preprocessed and the column that contains the classification variable
(last column is assumed by default):

```R
out <- IPF(iris, classColumn = 5)
```
For more specific information on how to use each filter, please refer to the
functions documentation page and the examples contained therein. For a general overview of the *NoiseFiltersR* package, please look up the associated vignette.

