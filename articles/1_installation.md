# 1. Installation

### Development version using `devtools`

The simplest way to install the in-development version is via the
`devtools` package. You just have to indicate that you want to install
from the in-development `main` branch.

``` r
devtools::install_github("adamkocsis/orange@main")
```

### Manual install

You can also the install the package manually with the following steps:

1.  Clone the repository to your local hard drive.  
2.  Make sure that you have the package dependencies.  
3.  Open a terminal and navigate to the directory where you cloned. The
    `orange` directory should be visible from there.  
4.  Run this line in the terminal

&nbsp;

    R CMD INSTALL orange

- *If you see an error suggesting that `R` is not found, you have to add
  it to your `PATH` environmental variable.*  
- *If the R packages that `orange` depend on are not installed, you have
  to install them manually, or you will get an error.*
