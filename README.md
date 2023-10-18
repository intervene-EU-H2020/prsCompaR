
# pgsCompaR

<!-- badges: start -->
<!-- badges: end -->

pgsCompaR is an [R data package](https://r-pkgs.org/data.html) that contains performance metrics for polygenic risk score (PGS) development methods measured across five European biobanks.

This data package doesn't provide any helpful functions for comparing PGS. It only contains processed experimental data and documentation. The raw experimental data [are also permissively licensed and publicly available](https://zenodo.org/records/10012996), but are more difficult to work with.

## Installation

The fastest way to install the development version of pgsCompaR using [devtools](https://devtools.r-lib.org):

``` r
devtools::install_github("intervene-EU-H2020/pgsCompaR")
```

This data package only depends on base R. You can download the built release and install it locally using `install.packages()` also.

### Development dependencies

Development dependencies are required to run the scripts in `data-raw/` that process the raw data and save the `rda` files in `data/`.

The simplest way to install the development dependencies is to use `renv` and restore the development profile:

``` bash
$ git clone https://github.com/intervene-EU-H2020/pgsCompaR.git
$ cd pgsCompaR
$ R
R version 4.3.1 (2023-06-16) -- "Beagle Scouts"
...
> renv::activate(profile="dev")
> renv::restore()
```

You may need to install [renv](https://rstudio.github.io/renv/articles/renv.html) first.

## Example

There are four datasets exported by this package:

| Dataset    | About                                                              |
|------------|--------------------------------------------------------------------|
| `metrics`  | Polygenic risk score performance metrics table for single biobanks |
| `meta_res` | Equivalent to `metrics`, but meta-analysed                         |
| `dst`      | Pairwise comparison of polygenic risk score development methods    |
| `pv_mrg`   | Equivalent to `dst`, but meta-analysed                             |

Dataset documentation can be viewed in R the normal way:

``` r
library(pgsCompaR)
data(metrics)
?metrics
```

## License

The data are licensed with [CC-BY-4.0](https://creativecommons.org/licenses/by/4.0/).

If you reuse data from this package in published work please cite:

> Remo's fantastic paper
