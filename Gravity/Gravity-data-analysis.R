library(dplyr) 
library(gravity)
library(igraph)
library(R.matlab)

data("gravity_no_zeros")

#### Part 1: Fit the gravity model and obtain p-values

Res <- gravity_no_zeros
fit <- ppml(
dependent_variable = "flow",
distance = "distw",
additional_regressors = c("rta", "comlang_off", "comcur", "contig", "iso_o", "iso_d"), data = Res
)
Res$fitted <- fit$fitted.values
Res$pvals1 <- ppois(Res$flow, Res$fitted, lower.tail=FALSE)    #small p-value: observed >> fitted
Res$pvals2 <- ppois(Res$flow, Res$fitted, lower.tail=TRUE)   #small p-value: observed << fitted 

# Create matrices of P-values
Pvals1 <- as_adjacency_matrix(graph_from_data_frame(Res[,c(1:2,12)]), attr="pvals1", sparse=FALSE)
Pvals2 <- as_adjacency_matrix(graph_from_data_frame(Res[,c(1:2,13)]), attr="pvals2", sparse=FALSE)
countries <- rownames(Pvals1)

# Extract GDP
n <- length(countries)
GDP <- numeric(length=n)
for (i in 1:length(countries)) {
	temp <- Res$gdp_o[Res$iso_o==countries[i]]
	GDP[i] <- temp[1]
}
sortTemp <- sort.int(GDP, decreasing=TRUE, index.return=TRUE)
topInd <- sortTemp$ix[1:15]
otherInd <- sortTemp$ix[16:n]


#### Part 2: Network mixed membership analysis

library(RSpectra)
A1 <- 1*as.matrix( (Pvals1!=0 & Pvals1 < 0.05) | (t(Pvals1)!=0 & t(Pvals1) < 0.05))
A2 <- 1*as.matrix( (Pvals2!=0 & Pvals2 < 0.05) | (t(Pvals2)!=0 & t(Pvals2) < 0.05))
diag(A1) <- 0
diag(A2) <- 0

source("MixedSCORE.R")
source("MixedSCORE-add.R")

## First network
K <- 3
eig.out <- RSpectra::eigs(A1, k=K)
ev <- eig.out$vectors[,1:K]
eig.out$values
R <- ev[,2:K] / matrix(rep(ev[,1], K-1), ncol = K-1, byrow = F)
R <- as.matrix(R, nrow = nrow(ev))
centers <- kmeans(R, 40, iter.max = 100, nstart = 100)$centers 
vh.out <- vertexSearch(centers, K = K)
vertices <- vh.out$vertices

pdf(file="NetworkPos.pdf")
plot(R[,1],R[,2], pch=16, cex=0.5, xlab="", ylab="")
text(R[otherInd,1], R[otherInd,2]-0.05, countries[otherInd], col="blue", cex=0.7)
text(R[topInd,1],R[topInd,2]-0.05, countries[topInd],col="orange",cex=0.7)
points(vertices[,1], vertices[,2], col="red", pch=16, cex=1.5)
lines(vertices[c(1,2),1], vertices[c(1,2),2], col="red", lty=2, lwd=2)
lines(vertices[c(2,3),1], vertices[c(2,3),2], col="red", lty=2, lwd=2)
lines(vertices[c(1,3),1], vertices[c(1,3),2], col="red", lty=2, lwd=2)
dev.off()

mscoreFit1 <- Estimation(R, vertices, K, eig.out$values, eig.out$vectors)
newdyadic <- mscoreFit1$memberships %*% mscoreFit1$P %*% t(mscoreFit1$memberships)
lognewdyadic <- log(newdyadic - min(newdyadic)+0.001)

R1 <- R
vertices1 <- vertices
lognewdyadic1 <- lognewdyadic



## Second network
K <- 3
eig.out <- RSpectra::eigs(A2, k=K)
ev <- eig.out$vectors[,1:K]
eig.out$values
R <- ev[,2:K] / matrix(rep(ev[,1], K-1), ncol = K-1, byrow = F)
R <- as.matrix(R, nrow = nrow(ev))
centers <- kmeans(R, 40, iter.max = 100, nstart = 100)$centers 
vh.out <- vertexSearch(centers, K = K)
vertices <- vh.out$vertices



pdf(file="NetworkNeg.pdf")
plot(R[,1],R[,2], pch=16, cex=0.5, xlab="", ylab="")
text(R[otherInd,1], R[otherInd,2]-0.05, countries[otherInd], col="blue", cex=0.7)
text(R[topInd,1],R[topInd,2]-0.05, countries[topInd],col="orange",cex=0.7)
points(vertices[,1], vertices[,2], col="red", pch=16, cex=1.5)
lines(vertices[c(1,2),1], vertices[c(1,2),2], col="red", lty=2, lwd=2)
lines(vertices[c(2,3),1], vertices[c(2,3),2], col="red", lty=2, lwd=2)
lines(vertices[c(1,3),1], vertices[c(1,3),2], col="red", lty=2, lwd=2)
dev.off()


mscoreFit2 <- Estimation(R, vertices, K, eig.out$values, eig.out$vectors)
newdyadic <- mscoreFit2$memberships %*% mscoreFit2$P %*% t(mscoreFit2$memberships)
lognewdyadic <- log(newdyadic - min(newdyadic)+0.001)

R2 <- R
vertices2 <- vertices
lognewdyadic2 <- lognewdyadic


N <- dim(Res)[1]
Res$lognewdyadic1 <- numeric(length=N)
Res$lognewdyadic2 <- numeric(length=N)
Res$logdeg1_o <- numeric(length=N)
Res$logdeg1_d <- numeric(length=N)
Res$logdeg2_o <- numeric(length=N)
Res$logdeg2_d <- numeric(length=N)
Res$logGDP_o <- log(Res$gdp_o)
Res$logGDP_d <- log(Res$gdp_d)
for (i in 1:N) {
    id_o <- which(countries==Res$iso_o[i])
    id_d <- which(countries==Res$iso_d[i])
	Res$lognewdyadic1[i] <- lognewdyadic1[id_o,id_d]
	Res$lognewdyadic2[i] <- lognewdyadic2[id_o,id_d]
	Res$logdeg1_o[i] <- log(mscoreFit1$degrees[id_o])
	Res$logdeg1_d[i] <- log(mscoreFit1$degrees[id_d])
	Res$logdeg2_o[i] <- log(mscoreFit2$degrees[id_o])
	Res$logdeg2_d[i] <- log(mscoreFit2$degrees[id_d])
}
Res$logdistw <- log(Res$distw)


fit_new <- ppml(
dependent_variable = "flow",
distance = "distw",
additional_regressors = c("rta", "comlang_off", "comcur", "contig", "lognewdyadic1", "lognewdyadic2", "iso_o", "iso_d"), data = Res
)


save(Res, fit, fit_new,Pvals1, Pvals2, countries, A1, A2, R1, vertices1, mscoreFit1, R2, vertices2, mscoreFit2, file="GravityResults-Revision-final.RData")





