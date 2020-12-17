#' Matrix pseudoinverse
#'
#' @description Computes the inverse or the pseudoinverse of a matrix
#'
#' @param M a numeric matrix (square matrix)
#' @param method an integer for the method of inversion.
#' If method = 0, only the nonzero eigenvalues are kept;
#' if method = 1, only the eigenvalues above a threshold are kept;
#' if method = 2, only the several first eigenvalues are kept;
#' if method = 3, a Tikhonov regularization (= ridge regression) is performed.
#' @param pa an integer for the parameter of inversion.
#' If method = 1, pa is the threshold below which the eigenvalues are not kept;
#' if method = 2, pa is an positive integer number corresponding to number of eigenvalues that are kept;
#' if method = 3, pa is the scaling factor for the identity matrix
#'
#' @return A numeric matrix corresponding to the pseudoinverse of M
#'
#' @examples
#'
#' # Data matrix of 2D landmark coordinates
#' data("Tropheus.IK.coord")
#' coords <- which(names(Tropheus.IK.coord) == "X1"):which(names(Tropheus.IK.coord) == "Y19")
#' proc.coord <- as.matrix(Tropheus.IK.coord[coords])
#'
#' # Covariance matrix of each population
#' S.phen.pop <- cov.group(proc.coord, groups = Tropheus.IK.coord$POP.ID)
#'
#' # Pseudo-inversion of a square matrix (covariance matrix of the population IKS5)
#' S2 <- S.phen.pop[, , "IKS5"]
#' invS2 <- minv(S2, method = 0, pa = 0)  # Pseudoinverse keeping non-zero eigenvalues
#' invS2 <- minv(S2, method = 1, pa = 10^-8)  # Pseudoinverse keeping eigenvalues above 10^-8
#' invS2 <- minv(S2, method = 2, pa = 5)  # Pseudoinverse keeping the first five eigenvalues
#' invS2 <- minv(S2, method = 3, pa = 0.5)  # Ridge regression with Tikhonov factor of 0.5
#'
#' @export
minv <-
  function (M, method = 0, pa = 0) {

    # Checkings
    if (length(dim(M)) > 2L || !is.numeric(M))
      stop("'M' must be a numeric matrix")
    if (!is.matrix(M))
      M <- as.matrix(M)

    # Pseudoinverse of M
    if (method == 0 | method == 1 | method == 2) {
      E <- eigen(M)
      D <- E$values
      D0 <- rep(0, length(D))
      # Keeps only the values above 'pa'
      if (method == 0 | method == 1) {
        if (method == 0) { pa <- .Machine$double.eps * max(dim(M)) * max(D) }  # Machine tolerance value
        for (i in 1:length(D0)) {
          if (abs(D[i]) > pa) { D0[i] <- 1 / (D[i]) }
        }
      }
      # Keeps only the 'pa' first eigenvalues
      if (method == 2) {
        if (pa > length(D)) { pa <- length(D) }
        for (i in 1:pa) { D0[i] <- 1 / (D[i]) }
      }
    }

    # Ridge regression = Tikhonov regularization
    if (method == 3) {
      ID <- diag(nrow(M))  # Identity matrix
      M0 <- pa * ID + M
      E <- eigen(M0)
      D <- E$values
      D0 <- 1 / D
    }

    V <- E$vectors
    MF <- V %*% diag(D0) %*% t(V)

    invisible(MF)

  }
