setwd("/Users/fostvedt/Jobs/")
#setwd("//Users//lukefostvedt//Documents//Simultaneous Inferences//")
source("simultaneous-estimation-functions.R")
source("Two-level-simulating-functions.R")

rep <- 100
tol <- 0.00001
Beta.stor.h2545 <- matrix(0,nrow=rep,ncol=10)
D.stor.h2545 <- array(0,dim=c(2,2,rep))
Sigma.stor.h2545 <- array(0,dim=c(2,2,rep))
SE.stor.h2545 <- array(0,dim=c(10,10,rep))
iterations.stor.h2545 <- NULL
Logliklihood.stor.h2545 <- NULL

se.ind.h2545 <- matrix(0,nrow=rep,ncol=10)
fix.ind.h2545 <- matrix(0,nrow=rep,ncol=10)
D.ind.h2545 <- array(0,dim=c(2,2,rep))
Sigma.ind.h2545 <- array(0,dim=c(2,2,rep))

for(k in 1:rep){
dat <- sim.2p.2l.hier(teacher=50, student=15, error1=16, error2=16, BY1=c(150,15,4,6),BY2=c(150,15,4,6),TRT=c(4,4), varB=c(20,20), cor1=0.5, cory=0.5)
d1 <- dat$data[which(dat$data$Yvar==1),]
d2 <- dat$data[which(dat$data$Yvar==2),]
m1 <- lmer(Score ~ Grade+TRT+x3+x4+(1|L2ID),data=d1,REML=F)
m2 <- lmer(Score ~ Grade+TRT+x3+x4+(1|L2ID),data=d2,REML=F)
a <- makedata.2level(d1,d2,m1,m2,slope=F)
b <- startvals.2level(X=a$X[,-1],Z=a$Z[,-1],Y=a$Y,ID=a$L2ID,nYvar=2)
c <- estimatepars.2lev(D=b$D, Sigma=b$Sigma, Beta=b$Beta, X=a$X[,-1], Z=a$Z[,-1], Y=a$Y, ID=a$L2ID, nYvar=2, tol=tol) 


Dm1 <- attr(VarCorr(m1)$L2ID,"stddev")^2
Dm2 <- attr(VarCorr(m2)$L2ID,"stddev")^2
se.ind.h2545[k,] <- c(sqrt(diag(vcov(m1))),sqrt(diag(vcov(m2))))
fix.ind.h2545[k,] <- c(fixef(m1),fixef(m2))
D.ind.h2545[,,k] <- adiag(Dm1,Dm2)
Sigma.ind.h2545[,,k] <- adiag(attr(VarCorr(m1),"sc")^2 ,attr(VarCorr(m2),"sc")^2 )


Beta.stor.h2545[k,] <- t(c$Beta)
SE.stor.h2545[,,k] <- c$SE
D.stor.h2545[,,k] <- c$D

Sigma.stor.h2545[,,k] <- c$Sigma
iterations.stor.h2545[k] <- c$iterations
Logliklihood.stor.h2545[k] <- c$Logliklihood


print(k)
save(se.ind.h2545,fix.ind.h2545,D.ind.h2545,Sigma.ind.h2545,Beta.stor.h2545,SE.stor.h2545, Sigma.stor.h2545, D.stor.h2545, iterations.stor.h2545, Logliklihood.stor.h2545, file="h2545.results")
	
}
save(se.ind.h2545,fix.ind.h2545,D.ind.h2545,Sigma.ind.h2545,Beta.stor.h2545,SE.stor.h2545, Sigma.stor.h2545, D.stor.h2545, iterations.stor.h2545, Logliklihood.stor.h2545, file="h2545.results")


