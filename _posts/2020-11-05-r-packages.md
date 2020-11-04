---
layout: post
title: "Some R packages: armspp and WoodburyMatrix"
excerpt_separator: <!--more-->
tags:
- R
---

Over the past year I've published two R packages that I'd like to highlight in this blog post.

<!--more-->

The first is armspp ([github](https://github.com/mbertolacci/armspp), [CRAN](https://CRAN.R-project.org/package=armspp)), which provides an efficient Rcpp implementation of the Adaptive Rejection Metropolis Sampling (ARMS) algorithm. The algorithm can be called from R with a user-specified target log density, or it can be called from C++ directly through inclusion of the header-only implementation.

The second package is WoodburyMatrix ([github](https://github.com/mbertolacci/WoodburyMatrix), [CRAN](https://cran.r-project.org/package=WoodburyMatrix)). It provides a hierarchy of classes and methods for manipulating matrices formed implicitly from the sums of the inverses of other matrices, a situation commonly encountered in spatial statistics and related fields. It makes it easy to use  the Woodbury matrix identity and the matrix determinant lemma to allow computation (e.g., solving linear systems) without having to form the actual matrix.

There is a more general idea here, which is to define S4 classes to encapsulate implicitly formed matrices. For example, one could make a matrix-product class to encapsulate the product of two matrices \\(AB\\) without explicitly calculating the product. Then the class could provide operations like `%*%` and `solve`, just as in the WoodburyMatrix package. This could be advantageous in a few situations, for example when \\(A\\) is dense and \\(B\\) is sparse (or vice-versa), and the product may be dense. Or possibly the matrices could themselves be implicit. I don't think such a package exists yet, so if I need it, I'll write it and release it. But if you do so first, please tell me about it!
