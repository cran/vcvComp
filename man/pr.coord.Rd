% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pr.coord.R
\name{pr.coord}
\alias{pr.coord}
\title{Principal coordinates ordination}
\usage{
pr.coord(V)
}
\arguments{
\item{V}{a square distance matrix}
}
\value{
A list containing the following named components:
\item{k}{the number of groups (value)}
\item{vectors}{the eigenvectors of the centered inner product matrix (matrix)}
\item{values}{the eigenvalues of the centered inner product matrix (vector)}
\item{PCoords}{the principal coordinates = scaled eigenvectors (matrix)}
\item{Variance}{a dataframe containing the following named variables:
\describe{
\item{eigenvalues}{eigenvalues of the centered inner product matrix}
\item{variance}{variance of each principal coordinate}
\item{exVar}{proportion of the total variation accounted by each principal coordinate}
\item{cumVar}{cumulative proportion of the total variation accounted by principal coordinate}
}
}
}
\description{
Performs a principal coordinates analysis of a distance matrix
}
\examples{

# Data matrix of 2D landmark coordinates
data("Tropheus.IK.coord")
coords <- which(names(Tropheus.IK.coord) == "X1"):which(names(Tropheus.IK.coord) == "Y19")
proc.coord <- as.matrix(Tropheus.IK.coord[coords])

# Data reduction
phen.pca <- prcomp(proc.coord, rank. = 5, tol = sqrt(.Machine$double.eps))
pc.scores <- phen.pca$x

# Covariance matrix of each population
S.phen.pop <- cov.group(pc.scores, groups = Tropheus.IK.coord$POP.ID)

# Squared distance matrix of the covariance matrices of all populations
eigen.phen.pop <- mat.sq.dist(S.phen.pop, dist. = "Riemannian")  # Riemannian distances

# Ordination of the squared distance matrix
prcoa.pop <- pr.coord(eigen.phen.pop)

# Visualization
plot(prcoa.pop$PCoords[, 1], prcoa.pop$PCoords[, 2])
abline(h = 0) ; abline(v = 0)
text(prcoa.pop$PCoords[, 1], prcoa.pop$PCoords[, 1], labels = rownames(prcoa.pop$PCoords))

}
