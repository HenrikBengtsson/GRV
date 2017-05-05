## Source: https://wwwf.imperial.ac.uk/~gmontana/software/grv/example.R
library(GRV)

# small GRV example:
# ------------------

# codes required for example:

"rmnorm" <- function (n = 1, mean = rep(0, d), varcov) # returns random multivariate Normal observation
{
    d <- if (is.matrix(varcov))
        ncol(varcov)
    else 1
    z <- matrix(rnorm(n * d), n, d) %*% chol(varcov)
    y <- t(mean + t(z))
    return(y)
}

"rwishart" <- function(df, p) # returns random Wishart matrix
{
    Z <- matrix(0, p, p)
    diag(Z) <- sqrt(rchisq(p, df:(df-p+1)))
    if(p > 1)
    {
        pseq <- 1:(p-1)
        Z[rep(p*pseq, pseq) + unlist(lapply(pseq, seq))] <- rnorm(p*(p-1)/2)
    }
    return(crossprod(Z))
}

"gen.snp.data" <- function(N,Q,rmafs) # returns N samples of Q snps with different minor allele frequencies (mafs)
{
	return(mapply(function(i) sample(c(0,1,2),N,prob=c((1-rmafs[i])^2,2*rmafs[i]*(1-rmafs[i]),rmafs[i]^2),replace=TRUE), 1:Q))
}


"vector.dists" <- function(X,dist) # a selection of distance measures for vectorial observations
{
	"norm.mutual.info.dist" <- function(x1,x2) 
	{
		n <- length(x1)
		nbins <- floor(sqrt(n)) # #choice of nbins as sqrt(n) comes from paper - 'evaluation of gene expression clustering via mutual information distance measure' (2007)
		jnt.density <- hist2d(x1,x2,nbins,show="FALSE")$counts/n
		marginal.density.x1 <- apply(jnt.density,1,sum)
		marginal.density.x2 <- apply(jnt.density,2,sum)
		entropy.x1 <- -sum(ifelse(marginal.density.x1 > 0, marginal.density.x1*log(marginal.density.x1), 0) )
		entropy.x2 <- -sum(ifelse(marginal.density.x2 > 0, marginal.density.x2*log(marginal.density.x2), 0) )
		jnt.density1 <- as.vector(jnt.density)
		jnt.entropy <- -sum(ifelse(jnt.density1 > 0, jnt.density1*log(jnt.density1), 0))
		mi <- entropy.x1 + entropy.x2 -  jnt.entropy
		return(1 - mi/max(entropy.x1,entropy.x2))
	}
	if(dist=="euclidean"){
		n <- nrow(X)
		XXt <- X%*%t(X)
		d <- diag(XXt)
		ones <- matrix(1,nrow=n,ncol=1)
		DX <- sqrt(d%*%t(ones) + ones%*%t(d) - 2*XXt)
	}
	if(dist=="mahalanobis"){
		DX <- sqrt(as.matrix(distance(X,method="mahalanobis")))
	}
	if(dist=="manhattan"){
		DX <- as.matrix(distance(X,method="manhattan"))
	}
	if(dist=="maximum"){
		DX <- as.matrix(dist(X,method="maximum"))
	}
	if(dist=="bray.curtis"){
		DX <- as.matrix(distance(X,method="bray-curtis"))
	}
	if(dist=="correlation.pearson"){
		DX <- 1 - cor(t(X),use="pairwise.complete.obs",method="pearson")
	}
	if(dist=="correlation.spearman"){
		DX <- 1 - cor(t(X),use="pairwise.complete.obs",method="spearman")
	}
	if(dist=="cosine.angle"){
		n <- nrow(X)
		DX <- matrix(0,nrow=n,ncol=n) 
		for(i in 1:(n-1)){
			for(j in 2:n){
				DX[i,j] <- DX[j,i] <-  1- X[i,]%*%X[j,]/(sqrt(sum(X[i,]^2)*sum(X[j,]^2)))
			}
		}
	}
	if(dist=="norm.mutual.information"){ 
		n <- nrow(X)
		DX <- matrix(0,nrow=n,ncol=n)
		for(i in 1:(n-1)){
			for(j in 2:n){
				DX[i,j] <- DX[j,i] <- norm.mutual.info.dist(X[i,],X[j,])	
			}
		}
	}
	return(DX)
}

"snp.dists" <- function(data,method) # snp distances taken from 'Similarity Measures for Clustering SNP Data' (2005)
{
	n <- nrow(data)
	vec.pairs <- combn(c(1:n),2)
	no.vec.pairs <- ncol(vec.pairs)
	dmat <- matrix(0,nrow=n,ncol=n)
	if(method=="IBS"){
		sim.vals <- sapply(1:no.vec.pairs, function(i) 0.5*mean(2 - abs(data[vec.pairs[1,i],] - data[vec.pairs[2,i],])))      
		dmat[lower.tri(dmat)] <- 1-sim.vals
		dmat <- dmat+t(dmat)
	}
	if(method=="simple.match"){
		sim.vals <- sapply(1:no.vec.pairs, function(i) mean(data[vec.pairs[1,i],] - data[vec.pairs[2,i],]==0))      
		dmat[lower.tri(dmat)] <- 1-sim.vals
		dmat <- dmat+t(dmat)
	}
	if(method=="sokal.sneath"){
		p <- ncol(data)
		sim.vals <- sapply(1:no.vec.pairs, function(i) {m.plus <- sum(data[vec.pairs[1,i],] - data[vec.pairs[2,i],]==0)
			 							m.minus <- p - m.plus
										m.plus/(m.plus+0.5*m.minus)})
		dmat[lower.tri(dmat)] <- 1-sim.vals
		dmat <- dmat+t(dmat)
	}
	if(method=="rogers.tanimoto.1"){
		p <- ncol(data)
		sim.vals <- sapply(1:no.vec.pairs, function(i) {m.plus <- sum(data[vec.pairs[1,i],] - data[vec.pairs[2,i],]==0)
										m.minus <- p - m.plus
										m.plus/(m.plus+2*m.minus)})
		dmat[lower.tri(dmat)] <- 1-sim.vals
		dmat <- dmat+t(dmat)
	}
	if(method=="hamman.1"){
		p <- ncol(data)
		sim.vals <- sapply(1:no.vec.pairs, function(i) {m.plus <- sum(data[vec.pairs[1,i],] - data[vec.pairs[2,i],]==0)
										m.minus <- p - m.plus
										(m.plus-m.minus)/p})
		sim.vals2 <- sim.vals + abs(min(sim.vals))
		dmat[lower.tri(dmat)] <- 1-(sim.vals2/max(sim.vals2))
		dmat <- dmat+t(dmat)
	}
	return(dmat)
}

library(MASS) 
library(graphics)
library(stats)
library(Matrix)
library(combinat)
library(ecodist)
library(gplots)

set.seed(101)

# First we demonstrate how to use the GRV test:
# generate (1) X = N x P multivariate data matrix
#          (2) Y = N x Q snp data matrix
#	     (3) obtain corresponding N x N distance matrices DX and DY and apply GRV test

# simulate X:

N <- 50
P <- 10
X <- rmnorm(N, mean=runif(P,-6,6), varcov=rwishart(df=P,p=P))	

# simulate Y:

Q <- 5
Y <- gen.snp.data(N,Q,rmafs=runif(Q,0.1,0.5))

# apply mahalanobis distance to X, and IBS distance to Y, and apply GRV test:

DX.mah <- vector.dists(X,dist="mahalanobis")
DY.ibs <- snp.dists(Y,method="IBS")

GRV.test(DX.mah,DY.ibs)
#grv.statistic   grv.p.value 
#   0.14236262    0.03878647 

# plot approximate probability density function (PDF) of GRV against sampling distribution using 10^4 Monte Carlo permutations:

GX <- create.G(N,DX.mah)
GY <- create.G(N,DY.ibs)

ys <- seq(0,0.2,len=1000)
pdf.GRV.vals <- pdf.GRV(DX.mah,DY.ibs,ys)
GRV.perms <- sapply(1:10000, function(i) {perms <- sample(1:N,replace=FALSE)
							GRV(GX[perms,perms],GY)})
plot(ys,pdf.GRV.vals,type="l",ylab="PDF(y)",xlab="y",main="PDF of GRV")
hist(GRV.perms,prob=TRUE,add=TRUE)

