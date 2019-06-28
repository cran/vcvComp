#' Principal coordinates ordination
#'
#' @description Performs a principal coordinates analysis of a distance matrix
#'
#' @param V a square distance matrix
#'
#' @return
#' A list containing the following named components:
#' \item{k}{the number of groups (value)}
#' \item{vectors}{the eigenvectors of the centered inner product matrix (matrix)}
#' \item{values}{the eigenvalues of the centered inner product matrix (vector)}
#' \item{PCoords}{the principal coordinates = scaled eigenvectors (matrix)}
#' \item{Variance}{a dataframe containing the following named variables:
#' \describe{
#' \item{eigenvalues}{eigenvalues of the centered inner product matrix}
#' \item{variance}{variance of each principal coordinate}
#' \item{exVar}{proportion of the total variation accounted by each principal coordinate}
#' \item{cumVar}{cumulative proportion of the total variation accounted by principal coordinate}
#' }
#' }
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
#' # Squared distance matrix of the covariance matrices of all populations
#' eigen.phen.pop <- mat.sq.dist(S.phen.pop, dist. = "Riemannian")  # Riemannian distances
#'
#' # Ordination of the squared distance matrix
#' prcoa.pop <- pr.coord(eigen.phen.pop)
#'
#' # Visualization
#' plot(prcoa.pop$PCoords[, 1], prcoa.pop$PCoords[, 2])
#' abline(h = 0) ; abline(v = 0)
#' text(prcoa.pop$PCoords[, 1], prcoa.pop$PCoords[, 1], labels = rownames(prcoa.pop$PCoords))
#'
#' @export
pr.coord <-
  function (V) {

    if (is.data.frame(V))
      V <- as.matrix(V)
    else if (!is.matrix(V))
      stop("'V' must be a matrix or a data frame")
    if (!all(is.finite(V)))
      stop("'V' must contain finite values only")
    if (dim(V)[1] != dim(V)[2])
      stop("'V' must be a square matrix")

    # Centered inner product matrix
    k <- dim(V)[1]
    H <- diag(k) - matrix((1 / k), nrow = k, ncol = k)  # centering matrix
    D <- - 0.5 * H %*% V %*% H

    # Number of principal coordinates
    max_pc <- k - 1

    # Eigenanalysis
    E <- eigen(D)
    vectors <- E$vectors[, 1:max_pc]
    rownames(vectors) <- rownames(V)
    colnames(vectors) <- paste("PCo", 1:max_pc, sep = "")

    L <- E$values[1:max_pc]
    L0 <- rep(0, length(L))
    # Keeps only the nonzero eigenvalues (above tol or below -tol)
    tol <- .Machine$double.eps * max(dim(D)) * max(L)  # Machine tolerance value
    for (i in 1:length(L)) {
      if (abs(L[i]) > tol) {
        L0[i] <- L[i]
      }
    }
    values <- L0
    PCoords <- vectors %*% diag(sqrt(values))
    colnames(PCoords) <- paste("PCo", 1:max_pc, sep = "")

    variance <- values / max_pc
    exVar <- values / sum(values)
    cumVar <- exVar
    for (i in 2:max_pc) {
      cumVar[i] <- cumVar[i - 1] + exVar[i]
    }
    Variance <- data.frame("eigenvalues" = values, "variance" = variance, "exVar" = exVar, "cumVar" = cumVar)

    prCoord <- list("k" = k, "vectors" = vectors, "values" = values, "PCoords" = PCoords, "Variance" = Variance)
    return(prCoord)

  }
