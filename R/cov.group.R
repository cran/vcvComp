#' Group covariance matrices
#'
#' @description Computes the covariance matrix of each group.
#' The effect of sexual dimorphism can be removed by using, for each group,
#' the average of the covariance matrix of males and the covariance matrix of females.
#'
#' @param X a data matrix with variables in columns and group names as row names
#' @param groups a character / factor vector containing grouping variable
#' @param sex NULL (default). A character / factor vector containing sex variable,
#' to remove sexual dimorphism by averaging males and females in each group
#' @param use an optional character string giving a method for computing covariances in the presence of missing values.
#' This must be (an abbreviation of) one of the strings "everything", "all.obs", "complete.obs", "na.or.complete", or "pairwise.complete.obs".
#'
#' @return A (p x p x m) array of covariance matrices,
#' where p is the number of variables and m the number of groups.
#'
#' @importFrom stats cov
#'
#' @seealso \code{\link[stats:cor]{cov}} and \code{\link[base]{scale}}
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
#' # Covariance matrix of each population, pooled by sex
#' S.phen.pooled <- cov.group(proc.coord,
#' groups = Tropheus.IK.coord$POP.ID, sex = Tropheus.IK.coord$Sex)
#'
#' @export
cov.group <-
  function (X, groups, sex = NULL, use = "everything") {

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
        Gvcv[, , i] <- cov(Xi, use = use)
      }
      # Correction for sexual dimorphism: mean males / females
      if (!is.null(sex) & length(slev) == 2) {
        X1 <- X[which(groups == glev[i] & sex == slev[1]), ]
        X2 <- X[which(groups == glev[i] & sex == slev[2]), ]
        Gsex <- array(c(cov(X1, use = use), cov(X2, use = use)), dim = c(p, p, 2))
        Gvcv[, , i] <- apply(Gsex, c(1, 2), mean)
      }
    }

    dimnames(Gvcv) <- list(colnames(X), colnames(X), glev)
    return(Gvcv)

  }
