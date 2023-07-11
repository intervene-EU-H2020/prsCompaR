
# prsCompaR

<!-- badges: start -->
<!-- badges: end -->

🚨 **This package is in under active development with no stable releases available yet** 🚨

prsCompaR is an [R data package](https://r-pkgs.org/data.html) that contains performance metrics for polygenic risk score (PRS) development methods measured across five European biobanks.

## Installation

The fastest way to install the development version of prsCompaR using [devtools](https://devtools.r-lib.org):

``` r
devtools::install_github("intervene-EU-H2020/prsCompaR")
```

This data package only depends on base R. You can download the built release and install it locally using `install.packages()` also.

### Development dependencies

Development dependencies are required to run the scripts in `data-raw/` that process the raw data and save the `rda` files in `data/`.

The simplest way to install the development dependencies is to use `renv` and restore the development profile:

``` bash
$ git clone https://github.com/intervene-EU-H2020/prsCompaR.git
$ cd prsCompaR
$ R
R version 4.3.1 (2023-06-16) -- "Beagle Scouts"
...
> renv::activate(profile="dev")
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
