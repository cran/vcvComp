#' Proportionality test of two variance-covariance matrices
#'
#' @description Tests the proportionality of two variance-covariance matrices
#'
#' @param n the sample size(s), given as a number or a vector of length 2
#' @param S1 a variance-covariance matrix
#' @param S2 a variance-covariance matrix
#' @param method an integer for the method of matrix inversion (see function 'minv')
#' @param pa an integer for the parameter of matrix inversion (see function 'minv')
#'
#' @return The P-value for the test of proportionality between two variance-covariance matrices
#'
#' @importFrom stats pchisq
#'
#' @seealso \code{\link{relative.eigen}} for the computation of relative eigenvalues,
#' @seealso \code{\link{minv}} for the method and the parameter used for the matrix inversion,
#' @seealso \code{\link[stats]{pchisq}} for Chi-squared distribution
#'
#' @references Mardia KV, Kent JT, Bibby JM (1979)
#' \emph{Multivariate analysis}. Academic Press, London.
#'
#' @examples
#'
#' # Data matrix of 2D landmark coordinates
#' data("Tropheus")
#' PHEN <- as.matrix(Tropheus[which(names(Tropheus) == "X1"):which(names(Tropheus) == "Y19")])
#'
#' # Procrustes superimposition
#' library("geomorph")
#' PHEN_array <- arrayspecs(PHEN, p = 19, k = 2)
#' phen.gpa <- gpagen(PHEN_array, print.progress = FALSE)
#' proc.coord <- two.d.array(phen.gpa$coords)
#'
#' # Data reduction
#' phen.pca <- prcomp(proc.coord, rank. = 5, tol = sqrt(.Machine$double.eps))
#' pc.scores <- phen.pca$x
#'
#' # Covariance matrix of each population
#' S.phen.pop <- cov.group(pc.scores, groups = Tropheus$POP.ID)
#'
#' # Maximum likelihood test of proportionality between 2 covariance matrices
#' # (IKA1 relative to IKS5) - 71 and 75 are the sample sizes
#' prop.vcv.test(n = c(71, 75), S.phen.pop[,,"IKA1"], S.phen.pop[,,"IKS5"])
#'
#' @export
prop.vcv.test <-
  function (n, S1, S2, method = 0, pa = 0) {

    if (is.null(n))
      stop("supply the sample size 'n'")
    if (!is.vector(n) | !is.numeric(n))
      stop("supply the sample size 'n' as a number or a numeric vector")
    if (length(n) < 1 | length(n) > 2)
      stop("supply the sample size 'n' as a single number or a vector of length 2")
    if (length(n) == 2)
      n <- 2 / (1 / n[1] + 1 / n[2])  # harmonic mean

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

    relValues <- relative.eigen(S1, S2, method = 0, pa = 0)$relValues  # relative eigenvalues
    p <- length(relValues)  # number of relative eigenvalues
    if (n < 10 * p) {
      warning("The sample size is not very large compared to the number of relative eigenvalues.")
    }

    a <- mean(relValues)  # arithmetic mean
    g <- prod(relValues) ^ (1 / p)  # geometric mean

    val <- 0.5 * n * p * log(a / g)
    ddl <- (p - 1) * (p + 2) / 2

    pValue <- pchisq(q = val, df = ddl, lower.tail = FALSE)

    return(pValue)

  }
