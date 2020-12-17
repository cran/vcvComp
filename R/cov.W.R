#' Within-group covariance matrix
#'
#' @description Computes the pooled within-group covariance matrix.
#' The effect of sexual dimorphism can be removed by using, for each group,
#' the average of the covariance matrix of males and the covariance matrix of females.
#'
#' @param X a data matrix with variables in columns and group names as row names
#' @param groups a character / factor vector containing grouping variable
#' @param sex NULL (default). A character / factor vector containing sex variable,
#' to remove sexual dimorphism by averaging males and females in each group
#' @param weighted logical. If FALSE (default), the average of all the within-group covariance matrices is used.
#' If TRUE, the within-group covariance matrices are weighted by their sample size.
#'
#' @return The pooled within-group covariance matrix
#'
#' @importFrom stats cov
#'
#' @seealso \code{\link[stats:cor]{cov}}
#'
#' @examples
#'
#' # Data matrix of 2D landmark coordinates
#' data("Tropheus.IK.coord")
#' coords <- which(names(Tropheus.IK.coord) == "X1"):which(names(Tropheus.IK.coord) == "Y19")
#' proc.coord <- as.matrix(Tropheus.IK.coord[coords])
#'
#' # Pooled within-group covariance matrix for all populations (weighted by sample size)
#' W <- cov.W(proc.coord, groups = Tropheus.IK.coord$POP.ID, weighted = TRUE)
#'
#' # Pooled within-group covariance matrix for all populations (unweighted)
#' W <- cov.W(proc.coord, groups = Tropheus.IK.coord$POP.ID)
#'
#' # Within-group covariance matrix for all populations, pooled by sex
#' W.mf <- cov.W(proc.coord, groups = Tropheus.IK.coord$POP.ID, sex = Tropheus.IK.coord$Sex)
#'
#' @export
cov.W <-
  function (X, groups, sex = NULL, weighted = FALSE) {

    if (is.data.frame(X))
      X <- as.matrix(X)
    else if (!is.matrix(X))
      stop("'X' must be a matrix or a data frame")
    if (!all(is.finite(X)))
      stop("'X' must contain finite values only")

    # Groups
    groups <- factor(groups)
    glev <- levels(groups)
    nlev <- length(glev)
    gsizes <- as.vector(table(groups))
    if (1 %in% gsizes) {
      warning("group with one entry found")
    }

    # Weighting
    wt <- gsizes
    if (weighted == FALSE) {
      wt <- rep(2, nlev)
    }

    # Sex
    slev <- 0
    if (!is.null(sex)) {
      sex <- factor(sex)
      slev <- levels(sex)
      if (length(slev) != 2) {
        warning("The number of sex categories is different from two.
                Sexual dimorphism will not be removed.")
      }
    }

    p <- ncol(X)  # number of variables
    Gvcv <- array(NA, dim = c(p, p, nlev))
    for (i in 1:nlev) {
      # No correction for sexual dimorphism
      if (is.null(sex) || length(slev) != 2) {
        Xi <- X[which(groups == glev[i]), ]
        Gvcv[, , i] <- cov(Xi) * (wt[i] - 1)
      }
      # Correction for sexual dimorphism: mean males / females
      if (!is.null(sex) & length(slev) == 2) {
        X1 <- X[which(groups == glev[i] & sex == slev[1]), ]
        X2 <- X[which(groups == glev[i] & sex == slev[2]), ]
        Gsex <- array(c(cov(X1), cov(X2)), dim = c(p, p, 2))
        Gvcv[, , i] <- apply(Gsex, c(1, 2), mean) * (wt[i] - 1)
      }
    }
    dimnames(Gvcv) <- list(colnames(p), colnames(p), glev)
    W <- apply(Gvcv, c(1, 2), sum)
    W <- W / (sum(wt) - nlev)

    return(W)

  }
