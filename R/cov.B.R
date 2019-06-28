#' Between-group covariance matrix
#'
#' @description Computes the between-group covariance matrix.
#' The effect of sexual dimorphism can be removed by using, for each group,
#' the average of the mean of males and the mean of females.
#'
#' @param X a data matrix with variables in columns and group names as row names
#' @param groups a character / factor vector containing grouping variable
#' @param sex NULL (default). A character / factor vector containing sex variable,
#' to remove sexual dimorphism by averaging males and females in each group
#' @param center either a logical value or a numeric vector of length equal to the number of columns of X
#' @param weighted logical. Should the between-group covariance matrix be weighted?
#'
#' @return The between-group covariance matrix
#'
#' @importFrom stats cov cov.wt
#'
#' @seealso \code{\link[stats]{cov}}, \code{\link[stats]{cov.wt}}
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
#' # Between-group covariance matrix for all populations
#' B <- cov.B(proc.coord, groups = Tropheus$POP.ID)
#'
#' # Between-group covariance matrix for all populations, pooled by sex
#' B.mf <- cov.B(proc.coord, groups = Tropheus$POP.ID, sex = Tropheus$Sex)
#'
#' @export
cov.B <-
  function (X, groups, sex = NULL, center = FALSE, weighted = FALSE) {

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

    p <- ncol(X)
    Gmeans <- matrix(NA, nrow = nlev, ncol = p, dimnames = list(glev, colnames(X)))
    for (i in 1:nlev) {
      # No correction for sexual dimorphism
      if (is.null(sex) || length(slev) != 2) {
        Gmeans[i, ] <- apply(X[which(groups == glev[i]), ], 2, mean)
      }
      # Correction for sexual dimorphism: mean males / females
      if (!is.null(sex) & length(slev) == 2) {
        Gsex1 <- apply(X[which(groups == glev[i] & sex == slev[1]), ], 2, mean)
        Gsex2 <- apply(X[which(groups == glev[i] & sex == slev[2]), ], 2, mean)
        Gsex <- rbind(Gsex1, Gsex2)
        Gmeans[i, ] <- apply(Gsex, 2, mean)
      }
    }

    if (weighted == TRUE) {
      wt <- gsizes / sum(gsizes)
      wcov <- cov.wt(Gmeans, wt, cor = FALSE, center)
      B <- wcov$cov
    }

    if (weighted == FALSE) {
      B <- cov(Gmeans)
    }

    return(B)

  }
