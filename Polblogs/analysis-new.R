library(igraph)
polblogs <- read.graph("polblogs.gml", format = "gml")
components <- decompose.graph(polblogs)
giantIndx <- which.max(sapply(components, vcount))
pbGiant <- components[giantIndx]
pbGiant <- pbGiant[[1]]

label <- V(pbGiant)$value
url <- V(pbGiant)$label
A <- get.adjacency(pbGiant)
A <- as.matrix(A)
A <- pmax(A, t(A))
n <- dim(A)[1]
degs <- rowSums(A)

library(RSpectra)
EigsResults <- eigs_sym(A,2)
V <- EigsResults$vectors
if (V[1,1]<0) {
	V[,1] <- -V[,1]
}
R <- V[,2]/V[,1]
Lambda <- diag(EigsResults$values)



temp = kmeans(R, 2)
centers = temp$centers
upperVX = max(centers)
lowerVX = min(centers)

w = (R - lowerVX)/(upperVX-lowerVX)
w = pmax(w,0)
w = pmin(w,1)
Q = matrix(c(1, 1, upperVX, lowerVX), ncol=2, nrow=2)
P0 = Q %*% Lambda %*% t(Q)
b1 = 1/sqrt(diag(P0))
P = diag(b1) %*% P0 %*% diag(b1)
pi_star = cbind(w, 1-w) %*% diag(1/b1)
pi = pi_star[,1]/rowSums(pi_star)
theta = V[,1]/(pi*b1[1]+(1-pi)*b1[2])

  


colors <- vector("character", length=n)
colors[label==0] <- "blue"
colors[label==1] <- "red"
plot(R, pch=20, cex=2, col=colors)
abline(h=upperVX, lty=2, lwd=2)
abline(h=lowerVX, lty=2, lwd=2)


plot(degs, pi, pch=20, cex=2, col=colors)
plot(pi, pch=20, cex=2, col=colors)

TempColors <- matrix(c(1,0,0, 0, 0, 1), nrow=2,ncol=3, byrow=TRUE)
newcolor <- cbind(pi, 1-pi) %*% TempColors

summary <- data.frame(R=R, label=label,pi=pi, theta=theta, url=url,col=newcolor)
temp <- sort(summary$theta, decreasing=TRUE, index.return=TRUE)
summary <- summary[temp$ix, ]
plot(summary$pi, summary$theta, pch=20, cex=2, col = rgb(summary$col.1, summary$col.2, summary$col.3), xlab="political orientation", ylab="influence", cex.lab=1.5, xaxt='n')

text(summary$pi[1], summary$theta[1]-0.05, labels=summary$url[1], cex=1.2)
text(0.1, summary$theta[2], labels=summary$url[2], cex=1.2, adj=0)
text(-0.02, summary$theta[3]-0.06, labels=summary$url[3], cex=1.2, adj=0)
text(summary$pi[4]+0.01, summary$theta[4]+0.02, labels=summary$url[4], cex=1.2, adj=0)
text(0.2, summary$theta[5]-0.05, labels=summary$url[5], cex=1.2, adj=0)
text(summary$pi[6]-0.02, summary$theta[6], labels=summary$url[6], cex=1.2, adj=1)
text(summary$pi[7]-0.02, summary$theta[7]+0.02, labels=summary$url[7], cex=1.2, adj=1)
text(summary$pi[8]-0.02, summary$theta[8]-0.03, labels=summary$url[8], cex=1.2, adj=1)
text(summary$pi[37], summary$theta[37]-0.05, labels=summary$url[37], cex=1.2)
text(summary$pi[156], summary$theta[156]-0.05, labels=summary$url[156], cex=1.2)
text(summary$pi[118], summary$theta[118]-0.05, labels=summary$url[118], cex=1.2)

liberal <- summary[label==0 & theta>0.2, ]
temp <- sort(liberal$pi, decreasing=TRUE, index.return=TRUE)
liberal <- liberal[temp$ix,]
conservative <- summary[label==1 & theta>0.2,]
temp <- sort(conservative$pi, decreasing=FALSE, index.return=TRUE)
conservative <- conservative[temp$ix,]


#Change the x-axis
axis(1,at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), labels =c(-1, -0.6, -0.2, 0.2, 0.6, 1))
