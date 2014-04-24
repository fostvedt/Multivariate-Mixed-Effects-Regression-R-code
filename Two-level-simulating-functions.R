require(MASS) 
require(plyr) 
require(lme4)  
require(magic)
require(reshape)
require(Matrix)
#############################################
#Two level simulation functions
#############################################



############################################
# two-level hierarchical model with random intercepts. 
############################################


sim.2p.2l.hier <- function(teacher, student, error1, error2, BY1,BY2,TRT, varB, cor1,cory){
#Fixed Effects 
B0 <- BY1[1]; B1<- BY1[2];  B2<- BY1[3]; B3<- BY1[4]
B4 <- BY2[1]; B5<- BY2[2];  B6<- BY2[3]; B7<- BY2[4]
T0<- TRT[1];
T1<- TRT[2];
#Variance Terms
varb0 <- varB[1];varb1 <- varB[2]
#correlation between random effects
cov1 <- cor1*sqrt(varb0*varb1) 	
d1 <- matrix(c(varb0,cov1,cov1,varb1),2,2)#var-cov matrix of random effects for y1
ind1 <- mvrnorm(teacher,c(0,0),d1) #generate bivariate random effects for y1
b0 <- ind1[,1] # individual intercepts' deviation from fixed intercept 
b1 <- ind1[,2] # individual slopes' deviation from fixed slope 
ind.int2 <- B1+b1
ind.int1 <- B0+b0
rand.eff1 <- cbind(ind.int1,ind.int2) 
# generating multivariate normal error terms with zero mean for both variables
covy <- cory*error1*error2
sigma <- matrix(c(error1^2,covy,covy,error2^2),2,2)
d <- sigma %x% diag(student) # var-cov matrix of error terms at time points 
err1 <- mvrnorm(teacher,rep(0,2*student),d) # generate multivariate normal error terms with
############################################## 
# Y1 = int + trt + demographics + continuous covariate + error + level1 error
# Y2 = int + trt + demographics + continuous covariate + error + level1 error
##############################################

t <-  matrix(rep(sample(rep(c(0,1),each=teacher/2)),each=student),teacher, student,byrow=T)
G <-  matrix( rep( sample( c(3,4,5,6), teacher, replace=T), each=student), teacher, student, byrow=T)
x3 <- matrix(rbinom(teacher*student,1,0.3),teacher, student)
x4 <- matrix(rbinom(teacher*student,1,0.5),teacher, student)
x5 <- matrix(rbinom(teacher*student,1,0.2),teacher, student)

data <- matrix(nrow=2*teacher,ncol=student) 
for(i in 1:teacher) { 
for(k in 1:student) { 
data[i,k] = B0 + B1*G[i,k]+ T0*t[i,k] + B2*x3[i,k] + B3*x4[i,k] + b0[i]
data[teacher+i,k] = B4 +B5*G[i,k]+ T1*t[i,k] + B6*x3[i,k] + B7*x4[i,k] + b1[i]
}} 
err <- rbind(err1[,1:student],err1[,(student+1):(2*student)])
data2 <- data + err 

############################################## 
data2 <- as.data.frame(data2)
names <- c()
for(i in 1:student) names[i]=paste("Student",i,sep="") 
colnames(data2) <- names ### Add column names to data 
mynames <- colnames(data2) 
data2$L2ID <- rep(1:teacher,2) 
data2$Yvar <- rep(1:2,each=teacher)
d <- reshape(data2,varying=mynames,idvar=c("L2ID","Yvar"), v.names="Score",timevar="L1ID",times=1:student,direction="long") 
d$Grade = rep(G[,1],2*student)
d$x3 <- as.vector(rbind(x3,x3))
d$x4 <- as.vector(rbind(x4,x4))
d$TRT <- as.vector(rbind(t,t))

#d
d <- d[order(d$Yvar,d$L2ID),]
d$L1ID <- rep(1:(teacher*student),2)
rand.eff <- cbind(rand.eff1)
list(data=d,rand.eff=rand.eff)
} 



###############
# Simulating 2 level longitudinal data with random slope
###############


sim.2p.2l.slope.long <- function(n,r,error1,error2, BY1, BY2,varB,TRT, corD,cory,C){
#Fixed Effects 
B0 <- BY1[1]; B1<- BY1[2];  B3<- BY1[3]; B4<- BY1[4]
B5 <- BY2[1]; B6<- BY2[2];  B8<- BY2[3]; B9<- BY2[4]
B2<- TRT[1];
B7<- TRT[2];
#Variance Terms
varb0 <- varB[1];varb1 <- varB[2];varb2 <- varB[3];varb3 <- varB[4]
cor12<-corD[1];cor13<-corD[2];cor14<-corD[3]
cor23<-corD[4];cor24<-corD[5];cor34<-corD[6]

#correlation between random effects
cov12 <- cor12*sqrt(varb0*varb1) 
cov13 <- cor13*sqrt(varb0*varb2)	
cov14 <- cor14*sqrt(varb0*varb3)	
cov23 <- cor23*sqrt(varb1*varb2)	
cov24 <- cor24*sqrt(varb1*varb3)	
cov34 <- cor34*sqrt(varb2*varb3)	

#var-cov matrix of random effects for y1
d1 <- matrix(c(varb0, cov12, cov13, cov14,
               cov12, varb1, cov23, cov24,
               cov13, cov23, varb2, cov34,
               cov14, cov24, cov34, varb3),4,4)
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

# generating multivariate normal error terms with zero mean for both variables
covy <- cory*error1*error2
sigma <- matrix(c(error1^2,covy,covy, error2^2),2,2)
d <- sigma %x% diag(r) # var-cov matrix of error terms at time points 
err1 <- mvrnorm(n,rep(0,2*r),d) # generate multivariate normal error terms with
dim(err1)
############################################## 
if(C==1){
Trt <- matrix(0,n,r)
Ts <- sample(1:n,floor(n/2),replace=FALSE)
Trt[Ts,] <- rep(0,r)
Trt[-Ts,] <- matrix(rep(c(rep(0,1),rep(1,r-1)),n-length(Ts)), nrow=n-length(Ts),ncol=r,byrow=TRUE)
}

else {
samp <- sample(c(1:r),n,replace=TRUE)
Trt <- matrix(0,n,r)
for(tt in 1:n) Trt[tt,] <- c(rep(0,samp[tt]),rep(1,r-samp[tt]))
}

x3 <- rbinom(n,1,0.3)
x4 <- rbinom(n,1,0.5)

############################################## 
data <- matrix(nrow=2*n,ncol=r) 

for(ii in 1:n){
for(kk in 1:r){
data[ii,kk] = B0 + b0[ii] + (B1+b1[ii])*(kk-1) + B2*Trt[ii,kk]  + B3*x3[ii]  + B4*x4[ii] 
data[n+ii,kk] = B5 + b5[ii] + (B6+b6[ii])*(kk-1) + B7*Trt[ii,kk] +  B8*x3[ii]  + B9*x4[ii] 
}}
err <- rbind(err1[,1:r],err1[,(r+1):(2*r)])
data2 <- data + err 

############################################## 
data2 <- as.data.frame(data2)

names <- c()
for(i in 1:r) names[i]=paste("Score",i,sep="") 
colnames(data2) <- names ### Add column names to data 
mynames <- colnames(data2) 
data2$L2ID <- rep(1:n,2) 
data2$Yvar <- rep(1:2,each=n)
d <- reshape(data2,varying=mynames,idvar=c("Yvar","L2ID"), v.names="Score",timevar="Grade",times=1:r,direction="long") 
d <- d[order(d$Yvar,d$L2ID),]
TRT <- NULL
for(tt in 1:n) TRT <- c(TRT,Trt[tt,])
d$TRT <- rep(TRT,2)
d$x3 <- rep(rep(x3,2), each=r)
d$x4 <- rep(rep(x4,2), each=r)
list(data=d, student.ran.eff = rand.eff1)
} 


###############
# Simulating 2 level longitudinal data 
# with only random intercepts
###############

sim.2p.2l.noint.long <- function(n,r,error1,error2, BY1, BY2,varB, TRT,corD,cory,C){
#Fixed Effects 
B0 <- BY1[1]; B1<- BY1[2];  B3<- BY1[3]; B4<- BY1[4]
B5 <- BY2[1]; B6<- BY2[2];  B8<- BY2[3]; B9<- BY2[4]
B2<- TRT[1];
B7<- TRT[2];
#Variance Terms
varb0 <- varB[1];varb1 <- varB[2]
cor12<-corD

#correlation between random effects
cov12 <- cor12*sqrt(varb0*varb1) 
	

#var-cov matrix of random effects for y1
d1 <- matrix(c(varb0, cov12, cov12, varb1),2,2)
ind1 <- mvrnorm(n,c(0,0),d1) #generate bivariate random effects for y1
b0 <- ind1[,1] # individual intercepts' deviation from fixed intercept 
b5 <- ind1[,2] # individual intercepts' deviation from fixed intercept 
#zero mean 
ind.int1 <- B0+b0
ind.int2 <- B5+b5
rand.eff1 <- cbind(ind.int1,ind.int2) 

# generating multivariate normal error terms with zero mean for both variables
covy <- cory*error1*error2
sigma <- matrix(c(error1^2,covy,covy, error2^2),2,2)
d <- sigma %x% diag(r) # var-cov matrix of error terms at time points 
err1 <- mvrnorm(n,rep(0,2*r),d) # generate multivariate normal error terms with
dim(err1)
############################################## 
if(C==1){
Trt <- matrix(0,n,r)
Ts <- sample(1:n,floor(n/2),replace=FALSE)
Trt[Ts,] <- rep(0,r)
Trt[-Ts,] <- matrix(rep(c(rep(0,1),rep(1,r-1)),n-length(Ts)), nrow=n-length(Ts),ncol=r,byrow=TRUE)
}

else {
samp <- sample(c(1:r),n,replace=TRUE)
Trt <- matrix(0,n,r)
for(tt in 1:n) Trt[tt,] <- c(rep(0,samp[tt]),rep(1,r-samp[tt]))
}

x3 <- rbinom(n,1,0.3)
x4 <- rbinom(n,1,0.5)

############################################## 
data <- matrix(nrow=2*n,ncol=r) 

for(ii in 1:n){
for(kk in 1:r){
data[ii,kk] = B0 + b0[ii]+ B1*(kk-1) + B2*Trt[ii,kk] +  B3*x3[ii]  + B4*x4[ii] 
data[n+ii,kk] = B5 +b5[ii] + B6*(kk-1) + B7*Trt[ii,kk] +  B8*x3[ii]  + B9*x4[ii] 
}}
err <- rbind(err1[,1:r],err1[,(r+1):(2*r)])
data2 <- data + err 

############################################## 
data2 <- as.data.frame(data2)

names <- c()
for(i in 1:r) names[i]=paste("Score",i,sep="") 
colnames(data2) <- names ### Add column names to data 
mynames <- colnames(data2) 
data2$L2ID <- rep(1:n,2) 
data2$Yvar <- rep(1:2,each=n)
d <- reshape(data2,varying=mynames,idvar=c("Yvar","L2ID"), v.names="Score",timevar="Grade",times=1:r,direction="long") 
d <- d[order(d$Yvar,d$L2ID),]
TRT <- NULL
for(tt in 1:n) TRT <- c(TRT,Trt[tt,])
d$TRT <- rep(TRT,2)
d$x3 <- rep(rep(x3,2),each=r)
d$x4 <- rep(rep(x4,2),each=r)
list(data=d, student.ran.eff = rand.eff1)
} 



