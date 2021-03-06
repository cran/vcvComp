#' Difference test for successive relative eigenvalues
#'
#' @description Tests the difference between two successive relative eigenvalues
#'
#' @param n the sample size(s), given as a number or a vector of length 2
#' @param relValues a vector of relative eigenvalues
#'
#' @return The P-values for the test of difference between successive eigenvalues
#'
#' @seealso \code{\link{relative.eigen}} for the computation of relative eigenvalues,
#' @seealso \code{\link[stats:Chisquare]{pchisq}} for Chi-squared distribution
#'
#' @references Mardia KV, Kent JT, Bibby JM (1979)
#' \emph{Multivariate analysis}. Academic Press, London.
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
#' # Test of the difference between 2 successives eigenvalues
#' # of the covariance matrix of IKA1 relative to IKS5
#' n_ika1 <- length(which(Tropheus.IK.coord$POP.ID == "IKA1"))  # sample size for IKA1
#' n_iks5 <- length(which(Tropheus.IK.coord$POP.ID == "IKS5"))  # sample size for IKS5
#' eigen.test(n = c(n_ika1, n_iks5), relValues = relEigen.a1s5$relValues)
#'
#' @export
eigen.test <-
  function (n, relValues) {

    if (is.null(n))
      stop("supply the sample size 'n'")
    if (!is.vector(n) | !is.numeric(n))
      stop("supply the sample size 'n' as a number or a numeric vector")
    if (length(n) < 1 | length(n) > 2)
      stop("supply the sample size 'n' as a single number or a vector of length 2")
    if (length(n) == 2)
      n <- 2 / (1 / n[1] + 1 / n[2])  # harmonic mean

    if (!is.vector(relValues) | !is.numeric(relValues))
      stop("supply the relative eigenvalues 'relValues' as a numeric vector")
    if (length(relValues) < 2)
      stop("supply at least 2 numbers in 'relValues'")

    p <- length(relValues)  # number of relative eigenvalues

    pValues <- rep(1, (p - 1))
    for (i in 1:(p - 1)) {
      val <- 2 * n * log((relValues[i] + relValues[i + 1]) / (2 * (relValues[i] * relValues[i + 1]) ^ 0.5))
      pValues[i] <- pchisq(q = val, df = 2, lower.tail = FALSE)
    }

    return(pValues)

  }
