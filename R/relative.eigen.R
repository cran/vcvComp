#' Relative eigenanalysis
#'
#' @description Computes the Riemanian distance between two variance-covariance matrices of same dimensions and the relative eigenvectors and eigenvalues of S1 with respect to S2
#'
#' @param S1 a variance-covariance matrix
#' @param S2 a variance-covariance matrix
#' @param method an integer for the method of matrix inversion (see function 'minv')
#' @param pa an integer for the parameter of matrix inversion (see function 'minv')
#'
#' @return
#' A list containing the following named components:
#' \item{relValues}{the vector of relative eigenvalues}
#' \item{relVectors}{the matrix of relative eigenvectors}
#' \item{distCov}{the distance between the two covariance matrices}
#' \item{relGV}{the product of the nonzero relative eigenvalues = the ratio of the generalized variances.
#' The generalized variance corresponds to the determinant of the covariance matrix.}
#' \item{logGV}{the log ratio of the generalized variances}
#' \item{q}{the number of nonzero eigenvalues}
#'
#' @seealso See \code{\link{minv}} for the method and the parameter used for the matrix inversion
#'
#' @references Bookstein F, Mitteroecker P (2014)
#' Comparing covariance matrices by relative eigenanalysis, with applications to organismal biology.
#' \emph{Evolutionary Biology 41}: 336-350.
#' \url{https://doi.org/10.1007/s11692-013-9260-5}
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
#' # Relative PCA = relative eigenanalysis between 2 covariance matrices
#' # (population IKA1 relative to IKS5)
#' relEigen.a1s5 <- relative.eigen(S.phen.pop[, , "IKA1"], S.phen.pop[, , "IKS5"])
#'
#' @export
relative.eigen <-
  function (S1, S2, method = 0, pa = 0) {

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

    # Computation of S2^-1 S1
    invS2 <- minv(S2, method, pa)
    M <- invS2 %*% S1

    # Eigenvectors and eigenvalues
    eigenM <- eigen(M)
    relVectors <- Re(eigenM$vectors)
    rownames(relVectors) <- colnames(S1)

    D <- Re(eigenM$values)
    D0 <- rep(0, length(D))
    # Keeps only the nonzero eigenvalues (above tol or below -tol)
    tol <- .Machine$double.eps * max(dim(M)) * max(D)  # Machine tolerance value
    for (i in 1:length(D)) {
      if (abs(D[i]) > tol) {
        D0[i] <- D[i]
        q <- i
      }
    }
    relValues <- D0

    # Riemannian distance between S1 and S2
    sqLogVal <- (log(relValues[1:q])) ^ 2
    distCov <- sqrt(sum(sqLogVal))

    # Product of the nonzero relative eigenvalues = ratio of the generalized variances
    relGV <- prod(relValues[1:q])
    logGV <- log(relGV)

    relEigen <- list("relValues" = relValues,
                     "relVectors" = relVectors,
                     "distCov" = distCov,
                     "relGV" = relGV,
                     "logGV" = logGV,
                     "q" = q)
    return(relEigen)

  }
