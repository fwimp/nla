# Perform nonlinear averaging of a temperature-dependent variable across a range of temperatures from empirically-derived data

Perform nonlinear averaging of a temperature-dependent variable across a
range of temperatures from empirically-derived data

## Usage

``` r
nla_tpc_empirical(tpcs, dailytemps, tempaxis)
```

## Arguments

- tpcs:

  A matrix of values measured across a temperature axis (represented in
  `tempaxis`).

- dailytemps:

  A matrix of temperatures across a 24h period, equating to each entry
  in `tempaxis`.

- tempaxis:

  The basis temperature axis.

## Value

A matrix of averaged values, where columns are temperatures, and rows
are individual tpc runs.

## Note

The precision of `tempaxis` MUST be greater than or equal to the values
of `dailytemps`.

The width of the `tpcs` matrix MUST be the same length as tempaxis.

## Input dimensionality

This function relies on carefully and appropriately-formatted inputs.

`tpcs` should be a 2D matrix where the rows correspond to individual
tpcs or bootstrap iterations of a tpc, and columns correspond to the
temperature axis in `tempaxis`.

`dailytemps` should be a 2D matrix where the rows correspond to sampling
times across a 24h window (so if your temperatures are hourly, this
should have 24 rows), and columns correspond to mean temperatures. So if
column 5 represents 10°C, `mean(dailytemps[,5]) == 10`

## Examples

``` r
if (FALSE) { # interactive()
nla_tpc_empirical()
}
```
