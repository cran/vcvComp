#' Scaling factor between two matrices
#'
#' @description Computes the maximum-likelihood estimate
#' of the scaling factor between two proportional covariance matrices.
#' Note that the scaling factor between the two matrices
#' is equal to the arithmetic mean of their relative eigenvalues.
#'
#' @param S1 a variance-covariance matrix
#' @param S2 a variance-covariance matrix
#' @param method an integer for the method of matrix inversion (see function 'minv')
#' @param pa an integer for the parameter of matrix inversion (see function 'minv')
#'
#' @return The scaling factor between the two matrices.
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
#' # Between-group (B) and within-group (W) covariance matrices for all populations
#' B <- cov.B(proc.coord, groups = Tropheus.IK.coord$POP.ID, sex = Tropheus.IK.coord$Sex)
#' W <- cov.W(proc.coord, groups = Tropheus.IK.coord$POP.ID, sex = Tropheus.IK.coord$Sex)
#'
#' # ML estimate of the scaling factor between B and W
#' sc <- scaling.BW(B, W)
#'
#' # Scaling of B to W
#' Bsc <- B / sc
#'
#' @export
scaling.BW <-
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

    p <- dim(S1)[1]
    M <- minv(S2, method = 0, pa = 0) %*% S1
    trM <- sum(diag(M))
    k <- trM / p
    return(k)

  }
