#' Squared distance matrix
#'
#' @description Computes the squared distance matrix of a set of covariance matrices
#'
#' @param Sm a (p x p x m) array of covariance matrices,
#' where p is the number of variables and m the number of groups.
#' @param dist. "Riemannian" or "Euclidean"
#' @param method an integer for the method of matrix inversion
#' @param pa an integer for the parameter of matrix inversion
#'
#' @return The matrix of squared Riemannian or Euclidean distances
#'
#' @seealso See \code{\link{minv}} for the method and the parameter used for the matrix inversion
#'
#' @examples
#'
#' # Data matrix of 2D landmark coordinates
#' data("Tropheus.IK.coord")
#' coords <- which(names(Tropheus.IK.coord) == "X1"):which(names(Tropheus.IK.coord) == "Y19")
#' proc.coord <- as.matrix(Tropheus.IK.coord[coords])
#'
#' # Data reduction
#' phen.pca <- prcomp(proc.coord, rank. = 5, tol = sqrt(.Machine$double.eps))
#' pc.scores <- phen.pca$x
#'
#' # Covariance matrix of each population
#' S.phen.pop <- cov.group(pc.scores, groups = Tropheus.IK.coord$POP.ID)
#'
#' # Squared Riemannian distance matrix of the covariance matrices of all populations
#' eigen.phen.r <- mat.sq.dist(S.phen.pop, dist. = "Riemannian")
#'
#' # Squared Euclidean distance matrix of the covariance matrices of all populations
#' eigen.phen.e <- mat.sq.dist(S.phen.pop, dist. = "Euclidean")
#'
#' @export
mat.sq.dist <-
  function (Sm, dist. = "Riemannian", method = 0, pa = 0) {

    k <- dim(Sm)[[3]]
    tol <- .Machine$double.eps * k  # Machine tolerance
    V <- matrix(0, nrow = k, ncol = k)
    for (l in 1:k) {
      for (m in 1:k) {
        if (m != l) {
          if (dist. == "Euclidean") {
            E_lm <- euclidean.dist(Sm[, , l], Sm[, , m])
          }
          if (dist. == "Riemannian") {
            E_lm <- relative.eigen(Sm[, , l], Sm[, , m], method, pa)$distCov
          }
          if ((E_lm)^2 > tol) {
            V[l, m] <- (E_lm) ^ 2
          }
        }
      }
    }

    rownames(V) <- dimnames(Sm)[[3]]
    colnames(V) <- dimnames(Sm)[[3]]

    return(V)

  }
