# udiff
An R package for estimating a unidiff model at the individual level

## Installation

``` r
remotes::install_github("wudongjie/udiff")
```

## Usage

Prepare a data with three variables.

For example, it could be the family of origin `y`, the family of destination `x` and the layer variable `z`.

```
#x,y,z should be categorical variables presented in numeric forms.
formula2 <- y ~ x + z
result <- udiff(formula2, data)
print(summary(result))
``` 