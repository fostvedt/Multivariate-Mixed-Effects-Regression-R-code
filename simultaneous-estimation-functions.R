require(MASS)
require(reshape)
require(magic)
require(Matrix)

############################################## 
# Functions for the EM algorithm for the simulataneous
# modeling approach
############################################## 

Egamma.i <- function(D,Y,Z,P,X,Beta){ return(D %*% t(Z) %*% P %*% (Y - X %*% Beta) )  }
Vgamma.i <- function(D,Z,P){ return(D - D %*% t(Z) %*% P %*% Z %*% D ) }
Eerror.i <- function(Sigma,Y,X,Beta,P,ni){ return( (Sigma %x% diag(ni) ) %*% P %*% (Y- X %*% Beta) )}
Verror.i <- function(Sigma,P,ni){return( (Sigma %x% diag(ni)) - (Sigma %x% diag(ni)) %*% P %*% (Sigma %x% diag(ni)))}


FisherInf.2level <- function(X,D,Z,Sigma){
	term1 <- matrix(0,nparx,nparx)
	for (i in 1:N){
	ind <- which(ID==i)
	X1 <- X[ind,-1]
	Z1 <- Z[ind,-1]	
	ni <- dim(X1)[1]/nvar
	V <- Z1 %*% D %*% t(Z1) + Sigma %x% diag(ni)
	P <- chol2inv( chol(V) )
	term1 <- term1 + t(X1) %*% P %*% X1
	}
	return(solve(term1))	
}

FisherInf.3level <- function(X,D1,D2,Z1,Z2,Sigma,ID){
	term1 <- matrix(0,nparx,nparx)
	for (i in 1:N){
	ind <- which(ID==i)
	X10 <- X[ind,-1]
	Z10 <- Z1[ind,-1]
	Z20 <- Z2[ind,-1]		
	ni <- dim(X10)[1]/nvar
	V <-  Z10 %*% D1 %*% t(Z10) + Z20 %*% D2 %*% t(Z20)  + Sigma %x% diag(ni) 
	P <- chol2inv(chol(V))
	term1 <- term1 + t(X10) %*% P %*% X10
	}
	return(solve(term1))	
}

################################################
# Two-Level Model Estimation functions
################################################




makedata.2level <- function(data1,data2,m1,m2,slope=TRUE){
X01 <- model.matrix(m1) 
X02 <- model.matrix(m2)
Y <-  matrix(c(data1$Score,data2$Score),ncol=1) 
L2ID <- matrix(c(data1$L2ID,data2$L2ID),ncol=1) 
X <- cbind(L2ID,adiag(X01,X02))
if(slope==TRUE){
Z <- cbind(L2ID,adiag(X01[,c(1,2)], X02[,c(1,2)]))	
}
else{
Z01 <- as.matrix(X01[,1],ncol=1)
Z02 <- as.matrix(X02[,1],ncol=1)
Z <- cbind(L2ID,adiag(Z01,Z02))	
}
return(list(X=X,Z=Z,Y=Y,L2ID=L2ID))
}


startvals.2level <- function(X,Z,Y,ID,nYvar){
# Get start val for Sigma
p <- q <- nYvar
n <- dim(Y)[1]/nYvar
Sigma <- matrix(0,nYvar,nYvar)
for (k in 1:q){
for (j in 1:p){
Sigma[k,j]  <- cov(Y[((k-1)*n+1):(k*n),],Y[((j-1)*n+1):(j*n),])		}
}
# Get start val for D
coef <- NULL
for(i in unique(ID)){
ind <- which(ID==i)
coef <- rbind(coef,as.numeric(lm(Y[ind,]~-1+Z[ind,])$coefficients))
}
m <- apply(coef,2,mean,na.rm=T)
D <- matrix(0,length(m),length(m))
for(i in 1:dim(coef)[1]){
	D <- D+(coef[i,]-m)%*%t(coef[i,]-m)
}
D <- D/dim(coef)[1]
# Get start val for Beta
nparx <- dim(X)[2]
term1 <- matrix(0,nparx,nparx)
term2 <- matrix(0,nparx,1)
N <- length(unique(ID))	
for (i in 1:N){
X1 <- X[which(ID==i),]
Z1 <- Z[which(ID==i),]
Y1 <- Y[which(ID==i)]
ni <- dim(X1)[1]/nYvar
V <- Z1 %*% D %*% t(Z1) + Sigma %x% diag(ni) 
P <- chol2inv( chol(V ))
term1 <- term1 + t(X1) %*% P %*% X1
term2 <- term2 + t(X1) %*% P %*% Y1
}
Beta <- chol2inv(chol(term1)) %*% term2	
return(list(Beta=Beta, D=D, Sigma=Sigma))
}




estimatepars.2lev <- function(D,Sigma,Beta,X,Z,Y,ID,nYvar,tol){
if(tol<0) tol=0.1
nparx <- dim(X)[2]
ni <- as.vector(table(ID))[1]/nYvar
p <- q <- nYvar
nparz <- dim(Z)[2]
ll <- matrix(0,nrow=1,ncol=length(unique(ID)))
repeat{
loglik <- NULL
	D.new <- matrix(0,nparz,nparz)
	Sigma1 <- Sigma.new <- matrix(0,p,p)
	term1 <- matrix(0,nparx,nparx)
    term2 <- matrix(0,nparx,1)
	N <- length(unique(ID))
	for (i in 1:N){ 
	ind <- which(ID==i)
	X10 <- X[ind,]
	Z10 <- Z[ind,]
	Y1 <- Y[ind,]
	ni <- dim(X10)[1]/nYvar
	V <- Z10 %*% D %*% t(Z10) + Sigma %x% diag(ni)
	P <- chol2inv( chol(V) )
	term1 <- term1 + t(X10) %*% P %*% X10
	term2 <- term2 + t(X10) %*% P %*% Y1	
	EE1 <- Egamma.i(D=D,Y=Y1,Z=Z10,P,X=X10,Beta)
	VV1 <- Vgamma.i(D,Z10,P)
	D.new <- D.new + EE1 %*% t(EE1) + VV1	
	
	EE <- Eerror.i(Sigma,Y=Y1,X=X10,Beta,P,ni=ni)
	VE <- Verror.i(Sigma,P,ni)
if(ni==1){
for (k in 1:q){
	for (j in 1:p){
		Sigma.new[k,j] <- Sigma.new[j,k] <- t(EE[(ni*(j-1) + 1):((j*ni))]) %*%
                       EE[(ni*(k-1) + 1):((k*ni))] + 
                   VE[(ni*(j-1) + 1):((j*ni)),(ni*(k-1) + 1):((k*ni))]
		}
	Sigma1 <- Sigma1 + Sigma.new	
	}
}
else{
	for (k in 1:q){
	for (j in 1:p){
    Sigma.new[k,j] <- Sigma.new[j,k] <- t(EE[(ni*(j-1) + 1):((j*ni))]) %*%
                       EE[(ni*(k-1) + 1):((k*ni))] + 
                    sum(diag(VE[(ni*(j-1) + 1):((j*ni)),(ni*(k-1) + 1):((k*ni))]))
		}
	Sigma1 <- Sigma1 + Sigma.new	
	}}
	llik <- -0.5*sum(log(eigen(V, symmetric = TRUE, only.values = TRUE)$values)) - 0.5* (t(Y1 - X10%*%Beta)%*%P%*%(Y1 - X10%*%Beta))
	loglik <- c(loglik, llik)
	}
Beta <- chol2inv(chol(term1)) %*% term2;# Beta
SE <- solve(term1)
D <- D.new/N #;D
Sigma <- Sigma1/dim(X)[1] #;Sigma
ll <- rbind(ll,loglik)
cond <- ll[dim(ll)[1],] -  ll[dim(ll)[1]-1,]
if(sum(abs(cond)) < tol) break
}
return(list(Beta=Beta,Sigma=Sigma,D=D,SE=SE,iterations=dim(ll)[1],Logliklihood=sum(ll[dim(ll)[1],])))
}



################################################
# Three-Level Model Estimation functions
################################################



makedata.3level <- function(data1,data2,m1,m2){
z1 <- names(ranef(m1))
z2 <- names(ranef(m2))

X01 <- model.matrix(m1) 
X02 <- model.matrix(m2)
Y <-  matrix(c(data1$Score,data2$Score),ncol=1) 
#L1ID <- matrix(c(data1$L1ID,data2$L1ID),ncol=1) 
L2ID <- matrix(c(data1$L2ID,data2$L2ID),ncol=1) 
L3ID <- matrix(c(data1$L3ID,data2$L3ID),ncol=1) 

X <- cbind(L2ID,L3ID,adiag(X01,X02))
indz11 <- match(names(ranef(m1)[[z1[1]]]),names(fixef(m1)))
indz21 <- match(names(ranef(m1)[[z1[2]]]),names(fixef(m1)))
indz12 <- match(names(ranef(m2)[[z2[1]]]),names(fixef(m2)))
indz22 <- match(names(ranef(m2)[[z2[2]]]),names(fixef(m2)))


mat11 <- as.matrix(X01[,indz11],ncol=length(indz11))
mat12 <- as.matrix(X02[,indz12],ncol=length(indz12))
mat21 <- as.matrix(X01[,indz21],ncol=length(indz21))
mat22 <- as.matrix(X02[,indz22],ncol=length(indz22))

Z1 <- cbind(L3ID,adiag(mat11 , mat12))
Z2 <- cbind(L2ID,adiag(mat21 , mat22))	
return(list(X=X,Z1=Z1,Z2=Z2,Y=Y,L2ID=L2ID,L3ID=L3ID))	
}



startvals.3level <- function(X1,Z1,Z2,Y,L2,L3,nYvar){
# Get start val for Sigma
p <- q <- nYvar
n <- dim(Y)[1]/nYvar
Sigma <- matrix(0,nYvar,nYvar)
for (k in 1:q){
for (j in 1:p){
Sigma[k,j]  <- cov(Y[((k-1)*n+1):(k*n),],Y[((j-1)*n+1):(j*n),])		}
}
# Get start val for D
coef1 <- NULL
for(i in unique(L2)){
ind <- which(L2==i)
coef1 <- rbind(coef1,as.numeric(lm(Y[ind,]~-1+Z1[ind,])$coefficients))
}
m1 <- apply(coef1,2,mean)
D1 <- matrix(0,ncol=length(m1),nrow=length(m1))
for(i in 1:dim(coef1)[1]){
	D1 <- D1+(coef1[i,]-m1)%*%t(coef1[i,]-m1)
}
D1 <- D1/dim(coef1)[1]
# Get start val for D
coef2 <- NULL
for(i in unique(L3)){
ind <- which(L3==i)
coef2 <- rbind(coef2,as.numeric(lm(Y[ind,]~-1+Z2[ind,])$coefficients))
}
m2 <- apply(coef2,2,mean)
D2 <- matrix(0,ncol=length(m2),nrow=length(m2))
for(i in 1:dim(coef2)[1]){
	D2 <- D2+(coef2[i,]-m2)%*%t(coef2[i,]-m2)
}
D2 <- D2/dim(coef2)[1]

# Get start val for Beta
nparx <- dim(X1)[2]
term1 <- matrix(0,nparx,nparx)
term2 <- matrix(0,nparx,1)
N <- length(unique(L3))	
for (i in 1:N){
ind <- which(L3==i)
X10 <- X1[ind,]
Z10 <- Z1[ind,]
Z20 <- Z2[ind,]
Y1 <- Y[ind]
ni <- dim(X10)[1]/nYvar
V <- Z10 %*% D1 %*% t(Z10) + Z20 %*% D2 %*% t(Z20)  + Sigma %x% diag(ni)
P <- chol2inv( chol(V) )
term1 <- term1 + t(X10) %*% P %*% X10
term2 <- term2 + t(X10) %*% P %*% Y1
}
Beta <- chol2inv(chol(term1)) %*% term2	
return(list(Beta=Beta,D1=D1,D2=D2,Sigma=Sigma))
}


#D1=b$D1 
#D2=b$D2 
#Sigma=b$Sigma 
#Beta=b$Beta 
#X=a$X[,-c(1,2)] 
#Z1=a$Z1[,-1]
#Z2=a$Z2[,-1]
#Y=a$Y
#ID=a$L3ID
#nYvar=2
#tol=0.1

estimatepars.3lev <- function(D1,D2,Sigma,Beta,X,Z1,Z2,Y,ID,nYvar,tol){
if(tol<0) tol=0.1
nparx <- dim(X)[2]
ni <- as.vector(table(ID))[1]/nYvar
p <- q <- nYvar
nparz1 <- dim(Z1)[2]
nparz2 <- dim(Z2)[2]
ll <- matrix(0,nrow=1,ncol=length(unique(ID)))
repeat{
loglik <- NULL
	D1.new <- matrix(0,nparz1,nparz1)
	D2.new <- matrix(0,nparz2,nparz2)
	Sigma1 <- Sigma.new <- matrix(0,p,p)
	term1 <- matrix(0,nparx,nparx)
    term2 <- matrix(0,nparx,1)
	N <- length(unique(ID))
	for (i in 1:N){ 
	ind <- which(ID==i)
	X10 <- X[ind,]
	Z10 <- Z1[ind,]
	Z20 <- Z2[ind,]
	Y1 <- Y[ind,]
	ni <- dim(X10)[1]/nYvar
	V <- Z10 %*% D1 %*% t(Z10) + Z20 %*% D2 %*% t(Z20)  + Sigma %x% diag(ni) 
	P <- chol2inv( chol(V) )
	term1 <- term1 + t(X10) %*% P %*% X10
	term2 <- term2 + t(X10) %*% P %*% Y1	
	EE1 <- Egamma.i(D=D1,Y=Y1,Z=Z10,P,X=X10,Beta)
	VV1 <- Vgamma.i(D1,Z10,P)
	D1.new <- D1.new + EE1 %*% t(EE1) + VV1	
	
	EE2 <- Egamma.i(D2,Y1,Z20,P,X10,Beta)
	VV2 <- Vgamma.i(D2,Z20,P)
	D2.new <- D2.new + EE2 %*% t(EE2) + VV2	
	EE <- Eerror.i(Sigma,Y=Y1,X=X10,Beta,P,ni=ni)
	VE <- Verror.i(Sigma,P,ni)
if(ni==1){
for (k in 1:q){
	for (j in 1:p){
		Sigma.new[k,j] <- Sigma.new[j,k] <- t(EE[(ni*(j-1) + 1):((j*ni))]) %*%
                       EE[(ni*(k-1) + 1):((k*ni))] + 
                   VE[(ni*(j-1) + 1):((j*ni)),(ni*(k-1) + 1):((k*ni))]
		}
	Sigma1 <- Sigma1 + Sigma.new	
	}
}
else{
	for (k in 1:q){
	for (j in 1:p){
    Sigma.new[k,j] <- Sigma.new[j,k] <- t(EE[(ni*(j-1) + 1):((j*ni))]) %*%
                       EE[(ni*(k-1) + 1):((k*ni))] + 
                    sum(diag(VE[(ni*(j-1) + 1):((j*ni)),(ni*(k-1) + 1):((k*ni))]))
		}
	Sigma1 <- Sigma1 + Sigma.new	
	}}
	llik <- -0.5*sum(log(eigen(V, symmetric = TRUE, only.values = TRUE)$values)) - 0.5* (t(Y1 - X10%*%Beta)%*%P%*%(Y1 - X10%*%Beta))
	loglik <- c(loglik, llik)
	}
Beta <- chol2inv(chol(term1)) %*% term2;# Beta
SE <- solve(term1)
D1 <- D1.new/N #;D
D2 <- D2.new/N
Sigma <- Sigma1/dim(X)[1] #;Sigma
ll <- rbind(ll,loglik)
cond <- ll[dim(ll)[1],] -  ll[dim(ll)[1]-1,]
if(abs(sum(cond)) < tol) break
}
return(list(Beta=Beta,Sigma=Sigma,D1=D1,D2=D2,SE=SE,iterations=dim(ll)[1],Logliklihood=sum(ll[dim(ll)[1],])))
}




