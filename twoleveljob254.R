setwd("/Users/fostvedt/Jobs/")
#setwd("//Users//lukefostvedt//Documents//Simultaneous Inferences//")
source("simultaneous-estimation-functions.R")
source("Two-level-simulating-functions.R")

rep <- 100
tol <- 0.00001
Beta.stor.long2541 <- matrix(0,nrow=rep,ncol=10)
D.stor.long2541 <- array(0,dim=c(4,4,rep))
Sigma.stor.long2541 <- array(0,dim=c(2,2,rep))
SE.stor.long2541 <- array(0,dim=c(10,10,rep))
iterations.stor.long2541 <- NULL
Logliklihood.stor.long2541 <- NULL

se.ind.long2541 <- matrix(0,nrow=rep,ncol=10)
fix.ind.long2541 <- matrix(0,nrow=rep,ncol=10)
D.ind.long2541 <- array(0,dim=c(4,4,rep))
Sigma.ind.long2541 <- array(0,dim=c(2,2,rep))

for(k in 1:rep){
dat <- sim.2p.2l.slope.long(n=100,r=4,error1=16,error2=16, BY1=c(150,15,7,2), BY2=c(150,15,4,6),varB=c(240,12,250,10),TRT=c(4,4), corD=c(0.5,0.5,0.5,0.5,0.5,0.5),cory=0.5,C=1)
d1 <- dat$data[which(dat$data$Yvar==1),]
d2 <- dat$data[which(dat$data$Yvar==2),]
m1 <- lmer(Score ~ Grade+TRT+x3+x4+(Grade|L2ID),data=d1,REML=F)
m2 <- lmer(Score ~ Grade+TRT+x3+x4+(Grade|L2ID),data=d2,REML=F)
a <- makedata.2level(d1,d2,m1,m2)
b <- startvals.2level(X=a$X[,-1],Z=a$Z[,-1],Y=a$Y,ID=a$L2ID,nYvar=2)
c <- estimatepars.2lev(D=b$D, Sigma=b$Sigma, Beta=b$Beta, X=a$X[,-1], Z=a$Z[,-1], Y=a$Y, ID=a$L2ID, nYvar=2, tol=tol) 


Dm1 <- diag(attr(VarCorr(m1)$L2ID,"stddev"))%*%attr(VarCorr(m1)$L2ID,"correlation")%*%diag(attr(VarCorr(m1)$L2ID,"stddev"))
Dm2 <- diag(attr(VarCorr(m2)$L2ID,"stddev"))%*%attr(VarCorr(m2)$L2ID,"correlation")%*%diag(attr(VarCorr(m2)$L2ID,"stddev"))
se.ind.long2541[k,] <- c(sqrt(diag(vcov(m1))),sqrt(diag(vcov(m2))))
fix.ind.long2541[k,] <- c(fixef(m1),fixef(m2))
D.ind.long2541[,,k] <- adiag(Dm1,Dm2)
Sigma.ind.long2541[,,k] <- adiag(attr(VarCorr(m1),"sc")^2 ,attr(VarCorr(m2),"sc")^2 )


Beta.stor.long2541[k,] <- t(c$Beta)
SE.stor.long2541[,,k] <- c$SE
D.stor.long2541[,,k] <- c$D

Sigma.stor.long2541[,,k] <- c$Sigma
iterations.stor.long2541[k] <- c$iterations
Logliklihood.stor.long2541[k] <- c$Logliklihood


print(c(k,c$iterations))
save(se.ind.long2541,fix.ind.long2541,D.ind.long2541,Sigma.ind.long2541,Beta.stor.long2541,SE.stor.long2541, Sigma.stor.long2541, D.stor.long2541, iterations.stor.long2541, Logliklihood.stor.long2541, file="long2541.results")
	
}
save(se.ind.long2541,fix.ind.long2541,D.ind.long2541,Sigma.ind.long2541,Beta.stor.long2541,SE.stor.long2541, Sigma.stor.long2541, D.stor.long2541, iterations.stor.long2541, Logliklihood.stor.long2541, file="long2541.results")



