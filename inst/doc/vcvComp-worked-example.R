## ----load, echo = TRUE---------------------------------------------------
library("vcvComp")
data("Tropheus")

## ----var_SexPop----------------------------------------------------------
outliers <- c(18, 56, 155, 351, 624)
Tropheus.IK <- Tropheus[- outliers, ]
# Sample reduced to six populations
Tropheus.IK <- subset(Tropheus.IK, subset = POP.ID %in% levels(POP.ID)[1:6])
Tropheus.IK$POP.ID <- factor(Tropheus.IK$POP.ID)
# New variable combining population and sex
Tropheus.IK$SexPop <- paste(Tropheus.IK$POP.ID, Tropheus.IK$Sex, sep = "_")
Tropheus.IK$SexPop <- as.factor(Tropheus.IK$SexPop)

## ----LM------------------------------------------------------------------
PHEN <- as.matrix(Tropheus.IK[which(names(Tropheus.IK) == "X1"):
                                which(names(Tropheus.IK) == "Y19")])
rownames(PHEN) <- Tropheus.IK$List_TropheusData_ID

## ----GPA-----------------------------------------------------------------
library("geomorph")  # load packages geomorph, rgl and RRPP
# conversion matrix -> array (19 landmarks, 2 dimensions)
PHEN_array <- arrayspecs(PHEN, p = 19, k = 2)
# Procrustes superimposition
phen.gpa <- gpagen(PHEN_array, print.progress = FALSE)
# conversion array -> matrix of Procrustes coordinates
proc.coord <- two.d.array(phen.gpa$coords)
colnames(proc.coord) <- colnames(PHEN)

## ----PCA-----------------------------------------------------------------
phen.pca <- prcomp(proc.coord, rank. = 5, tol = sqrt(.Machine$double.eps))
pc.scores <- phen.pca$x

## ----Pop vcv-------------------------------------------------------------
S.phen.pooled <- cov.group(pc.scores, groups = Tropheus.IK$POP.ID, sex = Tropheus.IK$Sex)

## ----Pop_PCoA------------------------------------------------------------
eigen.phen <- mat.sq.dist(S.phen.pooled, dist. = "Riemannian")  # Riemannian distances
prcoa <- pr.coord(eigen.phen)  # ordination
prcoa$Variance  # variance explained

## ----Pop_PCoA_screePlot, fig.height = 3, fig.width = 4, fig.cap = paste("**Figure 1**", "Fraction of variance explained by each principal coordinate in the ordination of the six *Tropheus moorii* populations.")----
# Visualization of the variance explained by each dimension (fig. 1)
barplot(prcoa$Variance$exVar, las = 1, col = "darkblue",
        names.arg = 1:nrow(prcoa$Variance), cex.axis = 0.8, cex  = 0.8,
        xlab = "Dimensions", ylab = "Variance explained")

## ----Pop_PCoA_ordination, fig.height = 3.8, fig.width = 4.5, fig.cap = paste("**Figure 2**", "Scatterplot of the first three principal coordinates (PCoord) in the ordination of the six *Tropheus moorii* populations. Populations living in sympatry are shown in dark blue, allopatric populations in light blue.")----
# Visualization of PCo1, PCo2 and PCo3 (fig. 2)
coul.pop <- c(rep("blue", 3), rep("darkblue", 3))  # colors
pch.pop <- c(rep(19, 3), rep(15, 3))  # symbols
pco3d <- c(1, 2, 3)  # dimensions
xyzlab <- c(paste("PCoord", pco3d[1]), 
            paste("PCoord", pco3d[2]), 
            paste("PCoord", pco3d[3]))
s3d <- scatterplot3d:::scatterplot3d(prcoa$PCoords[, pco3d[1:3]],
              xlab = xyzlab[1], ylab = xyzlab[2], zlab = xyzlab[3],
              color = coul.pop, pch = 19, angle = 55,
              type = "h", lty.hplot = 3, 
              cex.symbols = 1, cex.axis = 0.8)
s3d.coords <- s3d$xyz.convert(prcoa$PCoords[, pco3d[1:3]])
text(s3d.coords$x, s3d.coords$y, 
     labels = row.names(prcoa$PCoords), 
     pos = 4, cex = 0.7, col = coul.pop)

## ----IKA1-IKS5_ML_test---------------------------------------------------
table(Tropheus.IK$POP.ID)  # sample sizes
prop.vcv.test(n = c(69,75), S.phen.pooled[,,"IKA1"], S.phen.pooled[,,"IKS5"])  # ML test

## ----IKA1-IKS5_relPCA, fig.height = 3, fig.width = 4, fig.cap = paste("**Figure 3**", "Relative eigenvalues (on a log scale) of the population IKA1 relative to IKS5.")----
# Ratio of generalized variances of IKA1 and IKS5
relGV.multi(S.phen.pooled[, , c("IKA1", "IKS5")], logGV = FALSE)
# Relative PCA = relative eigenanalysis
relEigen.a1s5 <- relative.eigen(S.phen.pooled[, , "IKA1"], S.phen.pooled[, , "IKS5"])
relEigen.a1s5$relValues  # relative eigenvalues
# Visualization of the relative eigenvalues (fig. 3)
plot(relEigen.a1s5$relValues[1:relEigen.a1s5$q], 
     log = "y",  las = 1, col = "blue", type = "b", 
     main = "IKA1 relative to IKS5", cex = 0.8, 
     cex.main = 1, cex.axis = 0.8, cex.sub = 0.7, 
     sub = paste("Relative generalized variance =", relEigen.a1s5$relGV), 
     xlab = NA, ylab = "Relative eigenvalues")
abline(h = 1)

## ----IKA1-IKS5_vcvPattern, fig.height = 3, fig.width = 7, fig.cap = paste("**Figure 4**", "Visualization as thin-plate spline (TPS) deformation grids (Bookstein, 1991) of the shape pattern corresponding to the first relative PC, which has the maximal excess of variance in IKA1 relative to IKS5.")----
# Shape patterns corresponding to the relative eigenvectors
a1s5 <- c(which(Tropheus.IK$POP.ID %in% "IKA1"), 
               which(Tropheus.IK$POP.ID %in% "IKS5"))  # specimens
REF.A1S5 <- mshape(phen.gpa$coords[, , a1s5])  # average shape
A.A1S5 <- arrayspecs(t(phen.pca$rotation %*% relEigen.a1s5$relVectors), 
                     p = 19, k = 2)  # loadings
# Graphical parameters
WF <- cbind(c(1, 1, 2, 2, 3, 4, 5, 6, 7, 8, 9, 1, 12, 14, 14), 
            c(19, 18, 18, 3, 4, 5, 6, 7, 8, 9, 10, 10, 11, 11, 15))
gp3 <- gridPar(grid.col = "grey", tar.link.col = "blue", 
               tar.pt.size = 0.7, tar.pt.bg = "blue")
# Visualization of the first dimension (fig. 4)
par(new = FALSE, mfrow = c(1, 2), mar = c(0.5, 0.5, 0.5, 0.5))
plotRefToTarget(REF.A1S5, (REF.A1S5 - 0.01 * A.A1S5[, , 1]), 
                mag = 7, method = "TPS", gridPar = gp3, links = WF)
title(main = "-", line = -1)
plotRefToTarget(REF.A1S5, (REF.A1S5 + 0.01 * A.A1S5[, , 1]), 
                mag = 7, method = "TPS", gridPar = gp3, links = WF)
title(main = "+", line = -1)
title("First relative eigenvector", outer = TRUE, line = - 1)

## ----Sex_Pop_PCoA--------------------------------------------------------
S.phen.mf <- cov.group(pc.scores, groups = Tropheus.IK$SexPop)  # covariance matrices
eigen.phen.mf <- mat.sq.dist(S.phen.mf, dist. = "Riemannian")  # Riemannian distances
prcoa.mf <- pr.coord(eigen.phen.mf)  # ordination
prcoa.mf$Variance  # variance explained

## ----Sex_Pop_PCoA_screePlot, fig.height = 3, fig.width = 4, fig.cap = paste("**Figure 5**", "Fraction of variance explained by each principal coordinate of the 12 sex-specific covariance matrices.")----
# Visualization of the variance explained by each dimension (fig. 5)
barplot(prcoa.mf$Variance$exVar, las = 1, col = "darkblue", 
        names.arg = 1:nrow(prcoa.mf$Variance), cex.axis = 0.8, cex  = 0.8,
        xlab = "Dimensions", ylab = "Variance explained")

## ----Sex_Pop_PCoA_ordination, fig.height = 3.5, fig.width = 6, fig.cap = paste("**Figure 6**", "Principal coordinates ordination of the 12 sex-specific covariance matrices. Males in blue, females in red. Populations living in sympatry with *Tropheus polli* in dark colors.")----
# Visualization of PCo1 and PCo2 (fig. 6)
coul.mf <- c(rep(c("red", "blue"), 3), rep(c("darkred", "darkblue"), 3))  # colors
pch.mf <- c(rep(19, 6), rep(15, 6))  # symbols
pco <- c(1, 2)  # dimensions
plot(prcoa.mf$PCoords[, pco[1]], prcoa.mf$PCoords[, pco[2]], 
     xlab = paste("Principal coordinate", pco[1]), 
     ylab = paste("Principal coordinate", pco[2]), 
     asp = 1,  las = 1, pch = pch.mf, col = coul.mf, cex.axis = 0.8)
abline(h = 0) ; abline(v = 0)
text(prcoa.mf$PCoords[, pco[1]], prcoa.mf$PCoords[, pco[2]], 
     labels = rownames(prcoa.mf$PCoords), 
     adj = 1.5, cex = 0.6, col = coul.mf)

## ----Sex_IKA1_relPCA-----------------------------------------------------
pop.ika1 <- grep("IKA1", levels(Tropheus.IK$SexPop)) 
relEigen.ika1 <- relative.eigen(S.phen.mf[, , pop.ika1[1]], S.phen.mf[, , pop.ika1[2]])
relEigen.ika1$relGV  # ratio of generalized variances
relEigen.ika1$relValues  # relative eigenvalues

## ----Sex_IKA1_relVal, fig.height = 3, fig.width = 4, fig.cap = paste("**Figure 7**", "Relative eigenvalues (maximal ratios of variance) of females relative to males for the population IKA1.")----
# Visualization of the relative eigenvalues (fig. 7)
plot(relEigen.ika1$relValues[1:relEigen.ika1$q], 
     log = "y",  las = 1, col = "blue", type = "b", 
     main = "IKA1: females / males", cex.main = 1, cex.axis = 0.8, 
     xlab = NA, ylab = "Relative eigenvalues")
abline(h = 1)

## ----Sex_IKA1_vcvPattern, fig.height = 5, fig.width = 7, fig.cap = paste("**Figure 8**", "Visualization as thin-plate spline (TPS) deformation grids (Bookstein, 1991) of the shape patterns corresponding to the first relative PC, which has the maximal excess of variance in females relative to males for the population IKA1 (*top*), and the last relative PC, which has the maximal excess of variance in males relative to females (*bottom*).")----
# Population IKA1: average shape and loadings
ika1 <- which(Tropheus.IK$POP.ID %in% "IKA1")  # specimens
REF.IKA1 <- mshape(phen.gpa$coords[, , ika1])  # average shape
A.IKA1 <- arrayspecs(t(phen.pca$rotation %*% relEigen.ika1$relVectors), 
                     p = 19, k = 2)  # loadings
# Graphical parameters
WF <- cbind(c(1, 1, 2, 2, 3, 4, 5, 6, 7, 8, 9, 1, 12, 14, 14), 
            c(19, 18, 18, 3, 4, 5, 6, 7, 8, 9, 10, 10, 11, 11, 15))
gp3 <- gridPar(grid.col = "grey", tar.link.col = "blue", 
               tar.pt.size = 0.7, tar.pt.bg = "blue")
par(new = FALSE, mfrow = c(2, 2), mar = c(0.5, 0.5, 1, 0.5))
# Visualization of the first dimension (fig. 8: top)
plotRefToTarget(REF.IKA1, (REF.IKA1 - 0.01 * A.IKA1[, , 1]), 
                mag = 7, method = "TPS", gridPar = gp3, links = WF)
title(main = "-", line = - 1)
plotRefToTarget(REF.IKA1, (REF.IKA1 + 0.01 * A.IKA1[, , 1]), 
                mag = 7, method = "TPS", gridPar = gp3, links = WF)
title(main = "+", line = - 1)
title("First relative eigenvector", outer = TRUE, line = - 1)
# Visualization of the last dimension (fig. 8: bottom)
plotRefToTarget(REF.IKA1, (REF.IKA1 - 0.01 * A.IKA1[, , 5]), 
                mag = 7, method = "TPS", gridPar = gp3, links = WF)
title(main = "-", line = - 1)
plotRefToTarget(REF.IKA1, (REF.IKA1 + 0.01 * A.IKA1[, , 5]), 
                mag = 7, method = "TPS", gridPar = gp3, links = WF)
title(main = "+", line = - 1)
title("Last relative eigenvector", outer = TRUE, line = -16)

## ----BW_ML_test----------------------------------------------------------
# Computation of B and W (pooled by sex)
B <- cov.B(pc.scores, groups = Tropheus.IK$POP.ID, sex = Tropheus.IK$Sex)
W <- cov.W(pc.scores, groups = Tropheus.IK$POP.ID, sex = Tropheus.IK$Sex)
# Proportionality test between B and W = ML test
prop.vcv.test(n = c(6, 511), B, W)  # 6 groups, 511 specimens

## ----BW_Pop_PCoA---------------------------------------------------------
Bsc <- B / scaling.BW(B, W)  # scale B to W
# Create an array of group covariance matrices, B and W
S.bw <- array(c(S.phen.pooled, W, Bsc), 
              dim = c(dim(S.phen.pooled)[[1]], 
                      dim(S.phen.pooled)[[2]], 
                      dim(S.phen.pooled)[[3]] + 2))
dimnames(S.bw) <- list(dimnames(S.phen.pooled)[[1]], 
                       dimnames(S.phen.pooled)[[2]], 
                       c(dimnames(S.phen.pooled)[[3]], "W", "B"))
# Ordination
eigen.phen.bw <- mat.sq.dist(S.bw, dist. = "Riemannian")
prcoa.bw <- pr.coord(eigen.phen.bw)

## ----BW_Pop_PCoA_ordination, fig.height = 3.9, fig.width = 5, fig.cap = paste("**Figure 9**", "Principal coordinates ordination of the six populations (males and females pooled), along with their between-group (**B**) and their within-group (**W**) covariance matrices.")----
# Visualization of PCo1, PCo2 and PCo3 (fig. 9)
coul.bw <- c(coul.pop, rep("darkgreen", 2))  # colors
pco3d <- c(1, 2, 3)  # dimensions
xyzlab <- c(paste("PCoord", pco3d[1]), 
            paste("PCoord", pco3d[2]), 
            paste("PCoord", pco3d[3]))
s3d <- scatterplot3d:::scatterplot3d(prcoa.bw$PCoords[, pco3d[1:3]],
              xlab = xyzlab[1], ylab = xyzlab[2], zlab = xyzlab[3],
              color = coul.bw, pch = 19, angle = 55,
              type = "h", lty.hplot = 3, 
              cex.symbols = 1, cex.axis = 0.8)
s3d.coords <- s3d$xyz.convert(prcoa.bw$PCoords[, pco3d[1:3]])
text(s3d.coords$x, s3d.coords$y, labels = row.names(prcoa.bw$PCoords), 
     pos = 4, cex = 0.7, col = coul.bw)

## ----BW_relPCA, fig.height = 3, fig.width = 4, fig.cap = paste("**Figure 10**", "Relative eigenvalues of the between-group covariance matrix versus the within-group covariance matrix for the six *Tropheus* populations.")----
# Relative PCA of B with respect to to W
relEigenBW <- relative.eigen(B, W)
relEigenBW$relValues  # relative eigenvalues
# Test differences between two successive relative eigenvalues
eigen.test(n = c(6, 511), relValues = relEigenBW$relValues)
# Visualization of the relative eigenvalues (fig. 10)
plot(relEigenBW$relValues[1:relEigenBW$q], 
     log = "y",  las = 1, col = "blue", type = "b",
     main = "B relative to W", cex.main = 1, cex.axis = 0.8, 
     xlab = NA, ylab = "Relative eigenvalues")
abline(h = 1)

## ----BW_vcvPattern, fig.height = 5, fig.width = 7, fig.cap = paste("**Figure 11**", "Visualization of the shape patterns corresponding to the variance within populations (*top*), and the last relative PC, which has the maximal excess of variance within populations relative to that between populations (*bottom*).")----
# Shape patterns corresponding to the relative eigenvectors
REF <- mshape(phen.gpa$coords)  # average shape
A <- arrayspecs(t(phen.pca$rotation %*% relEigenBW$relVectors), p=19, k=2)  # loadings
# Graphical parameters
WF <- cbind(c(1, 1, 2, 2, 3, 4, 5, 6, 7, 8, 9, 1, 12, 14, 14), 
            c(19, 18, 18, 3, 4, 5, 6, 7, 8, 9, 10, 10, 11, 11, 15))
gp3 <- gridPar(grid.col = "grey", tar.link.col = "blue", 
               tar.pt.size = 0.7, tar.pt.bg = "blue")
par(new = FALSE, mfrow = c(2, 2), mar = c(0.5, 0.5, 1, 0.5))
# Visualization of the first dimension (fig. 11: top)
plotRefToTarget(REF, (REF - 0.01 * A[, , 1]), 
                mag = 7, method = "TPS", gridPar = gp3, links = WF)
title(main = "-", line = - 1)
plotRefToTarget(REF, (REF + 0.01 * A[, , 1]), 
                mag = 7, method = "TPS", gridPar = gp3, links = WF)
title(main = "+", line = - 1)
title("First relative eigenvector", outer = TRUE, line = - 1)
# Visualization of the last dimension (fig. 11: bottom)
plotRefToTarget(REF, (REF - 0.01 * A[, , 5]), 
                mag = 7, method = "TPS", gridPar = gp3, links = WF)
title(main = "-", line = - 1)
plotRefToTarget(REF, (REF + 0.01 * A[, , 5]), 
                mag = 7, method = "TPS", gridPar = gp3, links = WF)
title(main = "+", line = - 1)
title("Last relative eigenvector", outer = TRUE, line = - 16)

