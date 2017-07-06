---
layout: post
title: New R package on Github - climatedata
excerpt_separator: <!--more-->
tags:
- R
---

I've recently put an R package on Github, [climatedata](https://github.com/mbertolacci/climatedata/), that I've put together
as part of my PhD work. At present, the idea is to help the package user
download climate index data, either directly or as calculated from source data.

At this point, there are three indices available to download:

- the Indian Ocean Dipole, via the Dipole Mode Index provided by [JAMSTEC](http://www.jamstec.go.jp/frsgc/research/d1/iod/iod/dipole_mode_index.html);
- the Southern Oscillation Index, as provided by the [Australian Bureau of Meteorology](http://www.bom.gov.au/climate/current/soi2.shtml); and
- the Southern Annular Mode, either the [Marshall index](https://legacy.bas.ac.uk/met/gjma/sam.html), or calculated from [HadSLP2](http://www.metoffice.gov.uk/hadobs/hadslp2/).

I expect to add more functionality in dribs and drabs as time passes. The package is not yet submitted to CRAN, but it can be installed from Github using devtools:

    devtools::install_github('mbertolacci/climatedata')

Feel free to reach out to me if this package is useful to you or you have anything you'd like to add to it—or, even better, send a pull request!
