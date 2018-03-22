library(magick)
library(Matrix)
library(raster)
library(spdep)

acar_M_inverse <- function(input_nb) {
  Diagonal(x = spdep::card(input_nb))
}

acar_W <- function(input_nb, canonical = FALSE) {
  sqrt_n_neighbours <- sqrt(spdep::card(input_nb))

  adjacencies <- do.call(rbind, lapply(1 : length(input_nb), function(index) {
    neighbours <- input_nb[[index]]
    # Takes advantage of symmetric = TRUE in the sparseMatrix call below (and
    # also has the effect of removing 0 for islands)
    neighbours <- neighbours[neighbours > index]

    weights <- sqrt_n_neighbours[neighbours] * sqrt_n_neighbours[index]

    cbind(
      rep(index, length(neighbours)),
      neighbours,
      weights
    )
  }))

  sparseMatrix(
    i = adjacencies[, 1],
    j = adjacencies[, 2],
    x = adjacencies[, 3],
    symmetric = TRUE
  )
}

acar_Q <- function(input_nb, phi) {
  acar_M_inverse(input_nb) - phi * acar_W(input_nb)
}

set.seed(314159)

grid_raster <- raster(nrows = 10, ncols = 10, xmn = 0, xmx = 1, ymn = 0, ymx = 1)
grid_poly <- rasterToPolygons(grid_raster)
grid_nb <- poly2nb(grid_poly)

permutation <- sample.int(100)
Q_100 <- acar_Q(grid_nb, 0.12)[permutation, permutation]
str(Q_100)

plot_image <- image_graph(width = 960, height = 960, res = 200)
image(Q_100)
dev.off()
plot_image %>%
  image_crop('+0+80') %>%
  image_write(path = '../../assets/2018-03-22-sparse-gmrf-Q.png')

chol_Q_100 <- t(chol(Q_100))

str(chol_Q_100)

plot_image <- image_graph(width = 960, height = 960, res = 200)
image(chol_Q_100)
dev.off()
plot_image %>%
  image_crop('+0+80') %>%
  image_write(path = '../../assets/2018-03-22-sparse-gmrf-Q-chol.png')

chol_Q_100_permuted <- Cholesky(Q_100, LDL = FALSE, perm = TRUE)

str(chol_Q_100_permuted)

P <- as(chol_Q_100_permuted, 'pMatrix')

plot_image <- image_graph(width = 740, height = 740, res = 180)
image(P %*% Q_100 %*% t(P), main = 'Q (permuted)')
image(chol_Q_100_permuted, main = 'Cholesky')
dev.off()

plot_image %>%
  image_append() %>%
  image_write(path = '../../assets/2018-03-22-sparse-gmrf-Q-permuted.png')

z <- rnorm(nrow(Q_100))
y <- as.vector(solve(chol_Q_100_permuted, solve(chol_Q_100_permuted,
  z,
  system = 'Lt'
), system = 'Pt'))
print(y)

values(grid_raster) <- y[invPerm(permutation)]

plot_image <- image_graph(width = 960, height = 960, res = 200)
plot(grid_raster)
dev.off()

plot_image %>%
  image_crop('960x690+0+150') %>%
  image_write(path = '../../assets/2018-03-22-sparse-gmrf-Q-sample.png')
