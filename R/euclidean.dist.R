#' Euclidean distance between two covariance matrices
#'
#' @description Computes the Euclidean distance (Frobenius norm) between two variance-covariance matrices of same dimensions
#'
#' @param S1 a variance-covariance matrix
#' @param S2 a variance-covariance matrix
#'
#' @return Euclidean distance between S1 and S2 following Dryden et al. (2009).
#'
#' @references Dryden IL, Koloydenko A, Zhou D (2009)
#' Non-Euclidean statistics for covariance matrices, with applications to diffusion tensor imaging.
#' \emph{The Annals of Applied Statistics 3}:1102-1123.
#' \url{https://projecteuclid.org/euclid.aoas/1254773280}
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
#' # Euclidean distance between the covariance matrices of 2 populations
#' # (IKA1 relative to IKS5)
#' dist.a1s5 <- euclidean.dist(S.phen.pop[, , "IKA1"], S.phen.pop[, , "IKS5"])
#'
#' @export
euclidean.dist <-
  function (S1, S2) {

    if (is.data.frame(S2))
      S2 <- as.matrix(S2)
    if (is.data.frame(S1))
      S1 <- as.matrix(S1)
    if (is.null(S1) | is.null(S2))
      stop("supply both 'S1' and 'S2'")
    if (!is.matrix(S1) | !is.matrix(S2))
      stop("'S1' and 'S2' must be matrices or data frames")
    if (!all.equal(dim(S1), dim(S2)))
      stop("'S1' and 'S2' must be square matrices of the same dimensions")

    S12 <- S1 - S2
    M <- t(S12) %*% S12
    trM <- sum(diag(M))
    distEucl <- sqrt(trM)

    return(distEucl)

  }
