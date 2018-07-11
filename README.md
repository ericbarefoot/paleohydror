# paleohydror - Reconstructing paleohydraulics from your field data

A package containing data and functions implemented in R for paleohydraulic analysis from field data.

Please let me know if you use it. It's all inefficient, dubiously documented, and poorly optimized, but maybe you'll find it useful for something, and I am interested to see what that is.

## Install

You will have to have a working `R` installation to use these routines. To get `R` -- which is free -- download and install from [the R website](https://www.r-project.org/)

For easy install, just open up your favorite terminal emulator, and start `R` first get `devtools` from CRAN. This package was written by Hadley Wickham. It's very good.
```r
install.packages('devtools')
```
Next, load the package.
```r
library(devtools)
```
Now you can install directly from *anyone's* github account. For this one, you would next enter:
```r
install_github('ericbarefoot/paleohydror')
```
And that's it! This package will now be in your `installed packages` library, and you can call it with `library` just like any other package. To check for updates or new documentation just install again, and you'll get the newest stuff.

## Usage

To use this package in your script, simply place this at the top:

```r
library(paleohydror)
```

There are many functions included in this package, including several different routines for calculating paleoslope. As a demonstration, a commonly used relationship currently was published by Trampush et al. (2014). This relationship depends on the median bedload grain size (D<sub>50</sub>), and flow depth (H).

This relationship is implemented as

```r
trampush_slp(Dbed, H, perc = 50)
```

The user supplies a vector of paired `Dbed` = D<sub>50</sub> and `H` = H values and the function returns a vector of estimates of slope. This relationship contains uncertainty in the empirical parameters, so to obtain estimates of slope with uncertainty, the slope can be calculated at several different percentiles (perc) of the parameter distributions.

To obtain 10<sup>th</sup> and 90<sup>th</sup> percentile estimates, simply pass 10 and 90 respectively to the `perc` variable.
