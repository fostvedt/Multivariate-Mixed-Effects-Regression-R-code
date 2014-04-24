require(MASS) 
require(plyr) 
require(lme4)  
require(magic)
require(reshape)
require(Matrix)
#############################################
#Three level simulation functions
#############################################


###########################################
# Simulating longitudinal data 
############### INPUTS #################### 
sim.p.long.lmer <- function(s,n,r,error1,error2, BY1, BY2,varB,varS,TRT, corD,cors,cory,C){
varb0 <- varB[1];varb1 <- varB[2];varb2 <- varB[3];varb3 <- varB[4]
vars1 <- varS[1];vars2 <- varS[2]
B0 <- BY1[1]; B1<- BY1[2];  B3<- BY1[3]; B4<- BY1[4]
B5 <- BY2[1]; B6<- BY2[2];  B8<- BY2[3]; B9<- BY2[4]
B2<- TRT[1];
B7<- TRT[2];
cor12<-corD[1];cor13<-corD[2];cor14<-corD[3]
cor23<-corD[4];cor24<-corD[5];cor34<-corD[6]

cov12 <- cor12*sqrt(varb0*varb1) #correlation between random effects
cov13 <- cor13*sqrt(varb0*varb2)	
cov14 <- cor14*sqrt(varb0*varb3)	
cov23 <- cor23*sqrt(varb1*varb2)	
cov24 <- cor24*sqrt(varb1*varb3)	
cov34 <- cor34*sqrt(varb2*varb3)	

d1 <- matrix(c(varb0, cov12, cov13, cov14,
               cov12, varb1, cov23, cov24,
               cov13, cov23, varb2, cov34,
               cov14, cov24, cov34, varb3),4,4)#var-cov matrix of random effects for y1
ind1 <- mvrnorm(n,c(0,0,0,0),d1) #generate bivariate random effects for y1
b0 <- ind1[,1] # individual intercepts' deviation from fixed intercept 
b1 <- ind1[,2] # individual slopes' deviation from fixed slope 
b5 <- ind1[,3] # individual intercepts' deviation from fixed intercept 
b6 <- ind1[,4] # individual slopes' deviation from fixed slope 
#zero mean 
ind.slo1 <- B1+b1
ind.int1 <- B0+b0
ind.int2 <- B5+b5
ind.slo2 <- B6+b6
rand.eff1 <- cbind(ind.int1,ind.slo1,ind.int2,ind.slo2) 

#generate school random effect
covs <- cors*sqrt(vars1*vars1)	
d2 <- matrix(c(vars1,covs,covs,vars2),2,2) #var-cov matrix of random effects for y2
inds <- mvrnorm(s,c(0,0),d2) #generate bivariate random effects for y2
bs1 <- inds[,1] # individual intercepts' deviation from fixed intercept 
bs2 <- inds[,2] # individual slopes' deviation from fixed slope 
#zero mean 
s1.int <- bs1
s2.int <- bs2
rand.eff2 <- cbind(s1.int,s2.int) 

# generating multivariate normal error terms with zero mean for both variables
covy <- cory*error1*error2
sigma <- matrix(c(error1^2,covy,covy, error2^2),2,2)
d <- sigma %x% diag(r) # var-cov matrix of error terms at time points 
err1 <- mvrnorm(n,rep(0,2*r),d) # generate multivariate normal error terms with
dim(err1)
############################################## 
#r <- 4; s <- 30;
sa <- c(sample(c(1:s),s),sample(c(1:s),n-s,replace=TRUE))
if(C==1){
Trt <- matrix(0,s,r)
Ts <- sample(1:s,floor(s/2),replace=FALSE)
Trt[Ts,] <- rep(0,r)
Trt[-Ts,] <- matrix(rep(c(rep(0,1),rep(1,r-1)),s-length(Ts)), nrow=s-length(Ts),ncol=r,byrow=TRUE)
}

else {
samp <- sample(c(1:r),s,replace=TRUE)
Trt <- matrix(0,s,r)
for(tt in 1:s) Trt[tt,] <- c(rep(0,samp[tt]),rep(1,r-samp[tt]))
#Ts <- t(apply(T,1,cumsum))
}

x3 <- rbinom(n,1,0.3)
x4 <- rbinom(n,1,0.5)

############################################## 
data <- matrix(nrow=2*n,ncol=r) 

for(ii in 1:n){
for(kk in 1:r){
data[ii,kk] = B0+ bs1[sa[ii]] + B1*(kk-1) + B2*Trt[sa[ii],kk] + (b0[ii] + b1[ii]*Trt[sa[ii],kk])  +  B3*x3[ii]  + B4*x4[ii] 
data[n+ii,kk] = B5 + bs2[sa[ii]] + B6*(kk-1) + B7*Trt[sa[ii],kk] + (b5[ii] + b6[ii]*Trt[sa[ii],kk]) + B8*x3[ii]  + B9*x4[ii] 
}}
err <- rbind(err1[,1:r],err1[,(r+1):(2*r)])
data2 <- data + err 

############################################## 
data2 <- as.data.frame(data2)

names <- c()
for(i in 1:r) names[i]=paste("Score",i,sep="") 
colnames(data2) <- names ### Add column names to data 
mynames <- colnames(data2) 
data2$L3ID <- rep(sa,2) 
data2$L2ID <- rep(1:n,2) 
data2$Yvar <- rep(1:2,each=n)
d <- reshape(data2,varying=mynames,idvar=c("Yvar","L2ID"), v.names="Score",timevar="Grade",times=1:r,direction="long") 
d <- d[order(d$Yvar,d$L3ID,d$L2ID),]
TRT <- NULL
for(tt in 1:s) {
times <- as.numeric(table(d$L3ID)/r)[tt]
TRT <- c(TRT,rep(Trt[tt,],times))
}
d$TRT <- TRT
d$x3 <- rep(rep(x3,2),each=r)
d$x4 <- rep(rep(x4,2),each=r)
list(data=d, student.ran.eff = rand.eff1, school.rand.eff=rand.eff2 )
} 
# # # # # # # # # # End of Function



############################################
# three-level hierarchical model with random intercepts. 
############################################



sim.2p.3l.hier <- function(school,teacher, student, error1, error2, Grades, BY1,BY2, TRT, varB, varT, corD,cort,cory){

varb0 <- varB[1];varb1 <- varB[2];varb2 <- varB[3];varb3 <- varB[4]
vart1 <- varT[1];vart2 <- varT[2]
T0 <- TRT[1]; T1 <- TRT[2]
B0 <- BY1[1]; B1<- BY1[2]; B2<- BY1[3]; B3<- BY1[4]; 
B4 <- BY2[1]; B5<- BY2[2]; B6<- BY2[3]; B7<- BY2[4]; 
cor12<-corD[1];cor13<-corD[2];cor14<-corD[3]
cor23<-corD[4];cor24<-corD[5];cor34<-corD[6]

cov12 <- cor12*sqrt(varb0*varb1) #correlation between random effects
cov13 <- cor13*sqrt(varb0*varb2)	
cov14 <- cor14*sqrt(varb0*varb3)	
cov23 <- cor23*sqrt(varb1*varb2)	
cov24 <- cor24*sqrt(varb1*varb3)	
cov34 <- cor34*sqrt(varb2*varb3)	

#school level effects since trt is at school leve
Sig <- matrix(c(varb0, cov12, cov13, cov14,
               cov12, varb1, cov23, cov24,
               cov13, cov23, varb2, cov34,
               cov14, cov24, cov34, varb3),4,4)#var-cov matrix of random effects for y1
d1 <- Sig
ind1 <- mvrnorm(school,rep(0,dim(Sig)[1]),d1) #generate bivariate random effects for y1
b0 <- ind1[,1] # individual intercepts' deviation from fixed intercept 
b1 <- ind1[,2] # individual slopes' deviation from fixed slope 
b4 <- ind1[,3] # individual intercepts' deviation from fixed intercept 
b5 <- ind1[,4] # individual slopes' deviation from fixed slope 
#zero mean 
ind.slo1 <- B1+b1
ind.int1 <- B0+b0
ind.int2 <- B4+b4
ind.slo2 <- B5+b5
rand.eff1 <- cbind(ind.int1,ind.slo1,ind.int2,ind.slo2) 	
	
#teacher level random effects since teachers are random
covt1 <- cort*sqrt(vart2*vart1) #correlation between random effects
d1t <- matrix(c(vart1,covt1,covt1,vart2),2,2)
ind1t <- mvrnorm(teacher*school,c(0,0),d1t) 
b0t <- ind1t[,1] # teacher intercepts' deviation from fixed intercept 
b1t <- ind1t[,2] 
ind.int2.t <- B1+b1t
ind.int1.t <- B0+b0t
rand.eff2 <- cbind(ind.int1.t,ind.int2.t) 
terr1 <- matrix(rep(ind1t[,1],each=5),nrow=student*teacher,ncol=school)
terr2 <- matrix(rep(ind1t[,2],each=5),nrow=student*teacher,ncol=school)
terr <- rbind(terr1,terr2)
# generating multivariate normal error terms with zero mean for both variables
covy <- cory*error1*error2
sigma <- matrix(c(error1^2,covy,covy,error2^2),2,2)
d <- sigma %x% diag(school) # var-cov matrix of error terms at time points 
err1 <- mvrnorm(student*teacher,rep(0,2*school),d) # generate multivariate normal error terms with
############################################## 
# Y1 = int + trt + demographics + continuous covariate + error + level1 error
# Y2 = int + trt + demographics + continuous covariate + error + level1 error
##############################################

t <- matrix(0,nrow=teacher*student, ncol=school)
sa <- sample(1:school,floor(school/2))
t[,sa] <- 1
G <- matrix(0,nrow=teacher*student,ncol=school)
for(g in 1:school) {
l1 <- min(teacher,length(Grades))
s1 <- sample(Grades,l1,replace=FALSE)
s2 <- sample(Grades,teacher-l1,replace=TRUE)
G[,g] <- rep(c(s1,s2),each=student)
}
#else G <- matrix( rep( sample( Grades, teacher*school, replace=TRUE), each=student), student*teacher, school, byrow=FALSE)
x3 <- matrix( rbinom(teacher*student*school,1,0.3), teacher*student, school)
x4 <- matrix( rbinom(teacher*student*school,1,0.5), teacher*student, school)
x5 <- matrix( rbinom(teacher*student*school,1,0.2), teacher*student, school)

data <- matrix(nrow=2*student*teacher,ncol=school) 
for(i in 1:(teacher*student)) { 
for(k in 1:school) { 
data[i,k] = B0 + B1*G[i,k]+ T0*t[i,k] + B2*x3[i,k] + B3*x4[i,k] + b0[k] + b1[k]*G[i,k]
data[teacher*student+i,k] = B4 +B5*G[i,k]+ T1*t[i,k] + B6*x3[i,k] + B7*x4[i,k] + b4[k] + b5[k]*G[i,k]
}} 
err <- rbind(err1[,1:school],err1[,(school+1):(2*school)])
data2 <- data + err +terr

############################################## 
data2 <- as.data.frame(data2)
names <- c()
for(i in 1:school) names[i]=paste("school",i,sep="") 
colnames(data2) <- names ### Add column names to data 
mynames <- colnames(data2) 
d <- melt(data2,idvar=mynames)
tt0 <- matrix(rep(1:(school*teacher),each=student), nrow=teacher*student,ncol=school)
d$L2ID <- as.vector(rbind(tt0,tt0))
st0 <- matrix(1:(school*teacher*student), nrow=teacher*student,ncol=school)
d$L1ID <- as.vector(rbind(st0,st0))
d$Yvar <- rep(rep(1:2,each=student*teacher),school)
d$Grade = as.vector(rbind(G,G))
d$Score <- d$value
d$L3ID <- as.numeric(as.factor(d$variable))
d$x3 <- as.vector(rbind(x3,x3))
d$x4 <- as.vector(rbind(x4,x4))
d$x5 <- as.vector(rbind(x5,x5))
d$TRT <- as.vector(rbind(t,t))
#rand.eff <- cbind(rand.eff1,)
#list(data=d,rand.eff=rand.eff)
return(data=d[,-c(1,2)])
} 


