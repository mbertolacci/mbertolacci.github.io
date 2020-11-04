---
layout: page
title: Software
permalink: /software/
---

# R packages

## WoodburyMatrix ([github](https://github.com/mbertolacci/WoodburyMatrix) | [CRAN](https://cran.r-project.org/package=WoodburyMatrix))

A hierarchy of classes and methods for manipulating matrices formed implicitly from the sums of the inverses of other matrices, a situation commonly encountered in spatial statistics and related fields. Enables easy use of the Woodbury matrix identity and the matrix determinant lemma to allow computation (e.g., solving linear systems) without having to form the actual matrix.

## armspp ([github](https://github.com/mbertolacci/armspp) | [CRAN](https://CRAN.R-project.org/package=armspp))

An efficient Rcpp implementation of the Adaptive Rejection Metropolis Sampling (ARMS) algorithm proposed by [Gilks, W. R., Best, N. G. and Tan, K. K. C. (1995)](https://doi.org/10.2307/2986138). This allows for sampling from a univariate target probability distribution specified by its (potentially unnormalised) log density. This was my first CRAN package!

## bomdata ([github](https://github.com/mbertolacci/bomdata/))

An R package for interfacing with the [Australian Bureau of Meteorology's](http://www.bom.gov.au/) [climate data service](http://www.bom.gov.au/climate/data). It retrieves site metadata and daily precipitation values, with facilities for efficiently bulk downloading many sites.

## climatedata ([github](https://github.com/mbertolacci/climatedata/))

An R package that simplifies downloading a few climate indices.

## savepointr ([github](https://github.com/mbertolacci/savepointr))

An R package containing functions to create and manage 'savepoints', a method of saving intermediate states of long running processes in a fault tolerant manner.

# Other

## lorem-rss ([github](https://github.com/mbertolacci/lorem-rss))

A webservice that generates Lorem Ipsum RSS at specified intervals, available at [http://lorem-rss.herokuapp.com](http://lorem-rss.herokuapp.com).
