
# prsCompaR

<!-- badges: start -->
<!-- badges: end -->

🚨 **This package is in under active development with no stable releases available yet** 🚨

prsCompaR is an R data package that contains performance metrics for polygenic risk score (PRS) development methods measured across five European biobanks...

## Installation

You can install the development version of prsCompaR using [devtools](https://devtools.r-lib.org):

``` r
devtools::install_github("intervene-EU-H2020/prsCompaR")
```

### Development dependencies

The Suggests: field in the DESCRIPTION file contains packages needed to run the data preprocessing scripts in data-raw (these make the rda files in data/).

The simplest way to install development dependencies is by running:

``` bash
$ git clone https://github.com/intervene-EU-H2020/prsCompaR.git
$ cd prsCompaR
$ R
R version 4.3.1 (2023-06-16) -- "Beagle Scouts"
...
> renv::restore()
```

You may need to install [renv](https://rstudio.github.io/renv/articles/renv.html) first.

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(prsCompaR)
data(metrics)
?metrics
```

## License

The data are licensed with [CC-BY-4.0](https://creativecommons.org/licenses/by/4.0/).

If you reuse data from this package in published work please cite:

> Remo's fantastic paper
