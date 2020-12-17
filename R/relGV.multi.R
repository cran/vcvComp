#' Ratio of generalized variances
#'
#' @description Computes the (log-transformed) ratios of the generalized variances
#' of a set of covariance matrices
#'
#' @param Sm a (p x p x m) array of covariance matrices,
#' where p is the number of variables and m the number of groups.
#' @param logGV a logical argument to indicate if the ratios should be log-transformed
#'
#' @return The matrix of the (log-transformed) ratios of the generalized variances.
#' For each row, the ratio corrresponds to the group of the row
#' relative to the group of a column.
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
#' # Ratio of the generalized variances of 2 populations (IKA1 and IKS5)
#' relGV.multi(S.phen.pop[, , c("IKA1", "IKS5")], logGV = FALSE)
#'
#' @export
relGV.multi <-
  function (Sm, logGV = TRUE) {

    k <- dim(Sm)[[3]]
    V <- matrix(1, nrow = k, ncol = k)
    for (l in 1:k) {
      for (m in 1:k) {
        if (m != l) {
          V[l, m] <- det(Sm[, , l]) / det(Sm[, , m])
        }
      }
    }

    if (logGV == TRUE) {
      V <- log(V)
    }

    rownames(V) <- dimnames(Sm)[[3]]
    colnames(V) <- dimnames(Sm)[[3]]

    return(V)

  }
