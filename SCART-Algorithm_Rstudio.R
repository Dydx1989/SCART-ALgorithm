library(survival)
library(cmprsk)
library(rpart)
library(mvpart)
library(party)
library("evtree")
library(timereg)
library(randomForest)
library(randomForestSRC)
library(earth)
library(ggplot2)
Pbc=read.csv("pbc2.csv")
attach(Pbc)
#$View(Pbc)
Pbc_dat=data.frame(Time,status,Age,sex,ascites,hepato,spiders,edema,bili,chol,albumin,copper,alk.phos,ast,trig,platelet,protime,stage)
table(Pbc$status)
table(Pbc$trt,Pbc$status)
#### Cumulative Incidence Function######
dis=factor(trt, levels=c(1,2), labels=c("D-penicillamine","Placebo"))
tapply(Time,list(dis,status),median)
############# CIF ################################
##################################################
ci <- cuminc(ftime=Time, fstatus=status)
ci
plot(ci, curvlab = c("Transplant", "Death from other Diseases"),
     xlab = "Time in Days",col=c(2,3),lwd=3,ylab="Cumulative Incidence Function (CIF)",cex=1.8,cex.lab=1.4)
ci2 <- cuminc(ftime=Time, fstatus=status, group = trt)
plot(ci2, lty = c(1,2, 3,4), col = c("darkred", "blue","orange", "darkolivegreen"),
curvlab = c("Transplant + D-penicillamine", "Transplant + Placebo","Death from other Diseases + D-penicillamine","Death other Diseases + Placebo"), 
xlab = "Time in Days",cex.lab=1.5,cex=1.7,lwd=3,
ylab="Cumulative Incidence Function (CIF)")
calcna = function(time,event) {
  na.fit = survfit(coxph(Surv(time,event)~1), type="aalen")
  jumps = c(0, na.fit$time, max(time))
  # need to be careful at the beginning and end
  surv = c(1, na.fit$surv, na.fit$surv[length(na.fit$surv)])
  
  # apply appropriate transformation
  neglogsurv = -log(surv)
  
  # create placeholder of correct length
  naest = numeric(length(time))
  for (i in 2:length(jumps)) {
    naest[which(time>=jumps[i-1] & time<=jumps[i])] =
      neglogsurv[i-1]   # snag the appropriate value
  }
  return(naest)
}
Hazard1 = calcna(Time,(status==1)^2)
hazard1=cbind(Hazard1)
Hazard2=calcna(Time,(status==2)^2)
hazard2=cbind(Hazard2)
#### Joint Hazard Function #######
HA=hazard1+hazard2
################ Joint Survival Function #########################
##################################################################  
ssurv=exp(-HA)
##### Joint Cumulative Hazard Function ####
CIF=1-ssurv
plot(sort(Time),sort(CIF),type="l",lwd=3,col="blue",xlab=" Time in Days ",ylab=" Marginal Cumulative Incidence Function ", 
     main="Marginal Cumulative Incidence Function For PBC Data", cex.lab=1.5, cex.axis=1.0)   
#qplot(sort(Time),sort(CIF), geom=c("point","smooth"),xlab="Failure Time in Days ",ylab="Joint Cumulative Incidence Function ", 
#main="Joint Cumulative Incidence Function For PBC Data")   
## Data frame containing CIF
dat=data.frame(CIF,Age,as.factor(sex),as.factor(ascites),as.factor(hepato),as.factor(spiders),
as.factor(edema),bili,chol,albumin,copper,alk.phos,ast,trig,platelet,protime,as.factor(stage))
dat <- as.data.frame(dat)
names(dat)[1] <- "JointCIF"
head(dat)
set.seed(34567)
P=0.80
trainIndex<-sample.int(n = nrow(dat), size = floor(P*nrow(dat)), replace = F)
trainData <- dat[trainIndex,]
testData<- dat[-trainIndex,]
trainX <-trainData[,-1] # Create training feature data frame
testX <- testData[,-1] # Create test feature data frame
trainY=trainData[,1]
datatrain=data.frame(trainY,trainX)
head(datatrain)
############ Decision Tree #############
## Sum of Square 
ModelJCIF_SS<-rpart(trainY~. ,method="anova",data = datatrain, maxdepth=20,minbucket=20)
#printcp(ModelJCIF_SS)
#plotcp(ModelJCIF_SS) # visualize cross-validation results
summary(ModelJCIF_SS)
JCIFSS <- as.party(ModelJCIF_SS)
JCIFSS
plot(JCIFSS)
#### performance Using Root Mean Square Error 
Predc_JointCIFSS <- predict(ModelJCIF_SS,testData)
MSESS_JCIF<-mean((Predc_JointCIFSS-testData$JointCIF)^2)
round(MSESS_JCIF,4)
RMSESS_JCIF<- sqrt(mean((Predc_JointCIFSS-testData$JointCIF)^2))
round(RMSESS_JCIF,4)
##### Absolute Sum 
ModelJCIF_ABS<-rpart(trainY~. ,method="anova", data = datatrain,maxdepth=20,cp=0,minbucket=20,dissim="man")
#printcp(ModelJCIF_ABS)
#plotcp(ModelJCIF_ABS) # visualize cross-validation results
summary(ModelJCIF_ABS)
JCIFABS <- as.party(ModelJCIF_ABS)
JCIFABS
plot(JCIFABS)
#### performance Using Root Mean Square Error 
Predc_JointCIFABS <- predict(ModelJCIF_ABS,testData)
MSEABS_JCIF<-mean((Predc_JointCIFABS-testData$JointCIF)^2)
round(MSEABS_JCIF,4)
RMSEABS_JCIF<- sqrt(mean((Predc_JointCIFABS-testData$JointCIF)^2))
round(RMSEABS_JCIF,4)
##################### Joint Cumulative Probability ###########
########################################################
trainData <- dat[trainIndex,]
trainX <-trainData[,-1] # Create training feature data frame
trainY=trainData[,1]
datatrain=data.frame(trainY,trainX)
## Node Subset ######
Node3<- datatrain$trainY[which(datatrain$bili>=2.35&datatrain$edema>=0.25)]
Node5<- datatrain$trainY[which(datatrain$bili>=2.35&datatrain$edema<0.25&datatrain$stage>=3.5)]
Node6<- datatrain$trainY[which(datatrain$bili>=2.35&datatrain$edema<0.25&datatrain$stage<3.5)]
Node8<- datatrain$trainY[which(datatrain$bili<2.35&datatrain$albumin<3.325)]
Node12<- datatrain$trainY[which(datatrain$bili<2.35&datatrain$albumin>=3.325&datatrain$alk.phos<1871.5&
                                  datatrain$copper>=37.5&datatrain$chol>=309.75)]
Node13<- datatrain$trainY[which(datatrain$bili<2.35&datatrain$albumin>=3.325&datatrain$alk.phos<1871.5&
                                  datatrain$copper>=37.5&datatrain$chol<309.75)]
Node14<- datatrain$trainY[which(datatrain$bili<2.35&datatrain$albumin>=3.325&datatrain$alk.phos<1871.5&
                                  datatrain$copper<37.5)]
Node15<- datatrain$trainY[which(datatrain$bili<2.35&datatrain$albumin>=3.325&datatrain$alk.phos>=1871.5)]
### Time Subset
Time13<- datatrain$Time[which(datatrain$bili<2.35&datatrain$albumin>=3.325&datatrain$alk.phos<1871.5&datatrain$copper>=37.5&
                                datatrain$chol<309.75)]
plot(sort(Time13),sort(Node13),type="l") # returns 2    
plot(sort(Time13),sort(Node13),col=3,type="l",lwd=3,main="Power of the Five Selected Test Through Cause 1",
     xlab="Sample Sizes",ylab="Power of the Test",cex.lab=1.3,ylim=c(min(Node13),max(Node13)))
Time3<- datatrain$Time[which(datatrain$bili>=2.35&datatrain$edema>=0.25)]
lines(sort(Time3),sort(Node3),type="l") # returns 2 
Time5<- datatrain$Time[which(datatrain$bili>=2.35&datatrain$edema<0.25&datatrain$stage>=3.5)]
lines(sort(Time5),sort(Node5),type="l") # returns 2 
Time6<- datatrain$Time[which(datatrain$bili>=2.35&datatrain$edema<0.25&datatrain$stage<3.5)]
lines(sort(Time6),sort(Node6),type="l") # returns 2 
Time8<- datatrain$Time[which(datatrain$bili<2.35&datatrain$albumin<3.325)]
lines(sort(Time8),sort(Node8),type="l") # returns 2 
Time12<- datatrain$Time[which(datatrain$bili<2.35&datatrain$albumin>=3.325&datatrain$alk.phos<1871.5&
                                datatrain$copper>=37.5&datatrain$chol>=309.75)]
lines(sort(Time12),sort(Node12),type="l") # returns 2 
Time14<- datatrain$Time[which(datatrain$bili<2.35&datatrain$albumin>=3.325&datatrain$alk.phos<1871.5&
                                datatrain$copper<37.5)]
lines(sort(Time14),sort(Node14),type="l") # returns 2
Time15<- datatrain$Time[which(datatrain$bili<2.35&datatrain$albumin>=3.325&datatrain$alk.phos>=1871.5)]
lines(sort(Time15),sort(Node15),type="l") # returns 2  
################ Conditional CIF  ###################################################
#####################################################################################  
Pbc_dat2=data.frame(Time,status,trt,Age,sex,ascites,hepato,spiders,
                    edema,bili,chol,albumin,copper,alk.phos,ast,trig,platelet,protime,stage)
##### Joint Cumulative Hazard Function ####   
CumInc<- cuminc(Time,  # failure time variable
                status)  # variable with distinct codes for different causes of failure
timepoints <- function(w,times) {
  # w=list (see cuminc or km), times= times you want estimates for.
  # output is a list with components est giving the estimates, var giving
  # the variances,
  if (!is.null(w$Tests)) w <- w[names(w) != 'Tests']
  ng <- length(w)
  times <- sort(times)
  nt <- length(times)
  storage.mode(times) <- "double"
  storage.mode(nt) <- "integer"
  ind <- matrix(0,ncol=nt,nrow=ng)
  oute <- matrix(NA,ncol=nt,nrow=ng)
  outv <- oute
  storage.mode(ind) <- "integer"
  slct <- rep(TRUE,ng)
  for (i in 1:ng) {
    if (is.null((w[[i]])$est)) { slct[i] <- FALSE} else {
      z <- .Fortran("tpoi",as.double(w[[i]][[1]]),
                    as.integer(length(w[[i]][[1]])),ind[i,],times,nt,PACKAGE = "cmprsk")
      ind[i,] <- z[[3]]
      oute[i,ind[i,]>0] <- w[[i]][[2]][z[[3]]]
      if (length(w[[i]])>2) outv[i,ind[i,]>0] <- w[[i]][[3]][z[[3]]]
    }
  }
  dimnames(oute) <- list(names(w)[1:ng],as.character(times))
  dimnames(outv) <- dimnames(oute)
  list(est=oute[slct,,drop=FALSE],var=outv[slct,,drop=FALSE])
}                 
CumInc2<-timepoints(CumInc, times =Pbc$Time)
CumInc2=CumInc2$est
IF1= CumInc2[1,]
IF2= CumInc2[2,]
# Conditional Incidence Function pepe
CIF=(IF1)/ (1-IF2)
write.csv(CIF,"CIFPBC.csv")
CIF=read.csv("CIFPBC.csv")
plot(sort(Pbc$Tim),CIF[,2],type="l",ylab=" Conditional Incidence Function (CIF)",cex.lab=1.5, cex.axis=1.3,
     xlab="Time in Days", col="blue", lwd=3, main="Conditional Incidence Function (CIF) of PBC Data ") 
#CPF.mean=mean(CPF,na.rm=TRUE)
#CPF<- ifelse(is.na(CPF), CPF.mean, CPF)
#Ttime=unique(Time)
#plot(Ttime,CPF,type="l",ylim=c(0,1)) 
## Data frame containing CIF
CPF=(IF21+IF22)/ (IF21+IF22+IF11+IF12)
CPF_1=write.csv(CPF,"CPF_1.csv")
CPF_2=read.csv("CPF_1.csv")
Conditional_IF=as.matrix(CPF_2)
CPF.mean=mean(Conditional_IF,na.rm=TRUE)
#Conditional_IF=data.frame(CPF)
Conditional_IF=ifelse(is.na(Conditional_IF==T), CPF.mean, Conditional_IF)
#Conditional_IF=data.frame(Conditional_IF)
datCPF=data.frame(Conditional_IF,Age,sex,ascites,hepato,spiders,edema,bili,chol,albumin,copper,alk.phos,ast,trig,platelet,protime,stage)
datCPF <- as.data.frame(datCPF)
names(datCPF)[1] <- "ConditionaIF"
head(datCPF)
set.seed(22345)
P=0.85
trainIndex<-sample.int(n = nrow(datCPF), size = floor(P*nrow(datCPF)), replace = F)
traindatCPF <- datCPF[trainIndex,]
testdatCPF<- datCPF[-trainIndex,]
trainX_CPF <-traindatCPF[,-1] # Create training feature data frame
testX_CPF <- testdatCPF[,-1] # Create test feature data frame
trainY_CPF=traindatCPF[,1]
datatrain=data.frame(trainY_CPF,trainX_CPF)
############ Decision Tree #############
## Sum of Square 
ModelCPF_SS<-rpart(trainY_CPF~. ,method="anova",data = datatrain, maxdepth=20,minbucket=20)
printcp(ModelCPF_SS)
plotcp(ModelCPF_SS) # visualize cross-validation results
summary(ModelCPF_SS)
CPFSS <- as.party(ModelCPF_SS)
CPFSS
plot(CPFSS)
#### performance Using Root Mean Square Error 
Predc_ConditionalIFSS <- predict(ModelCPF_SS,testdatCPF)
MSESS_CPF<-mean((Predc_ConditionalIFSS-testdatCPF$ConditionaIF)^2)
round(MSESS_CPF,4)
RMSESS_CPF<- sqrt(mean((Predc_ConditionalIFSS-testdatCPF$ConditionaIF)^2))
round(RMSESS_CPF,4)
##### Absolute Sum 
ModelCPF_ABS<-rpart(trainY_CPF~.,method="anova", data = datatrain,maxdepth=20,cp=0,minbucket=20,dissim="man")
printcp(ModelCPF_ABS)
plotcp(ModelCPF_ABS) # visualize cross-validation results
summary(ModelCPF_ABS)
CPFABS <- as.party(ModelCPF_ABS)
CPFABS
plot(CPFABS)
#### performance Using Root Mean Square Error 
Predc_ConditionalIFABS <- predict(ModelCPF_ABS,testdatCPF)
MSEABS_CPF<-mean((Predc_ConditionalIFABS-testdatCPF$ConditionaIF)^2)
round(MSEABS_CPF,4)
RMSEABS_CPF<- sqrt(mean((Predc_ConditionalIFABS-testdatCPF$ConditionaIF)^2))
round(RMSEABS_CPF,4)
################ Simultaneous Log-Odds Residual Based Tree ##########################
#####################################################################################  
############# Log-Odd ################################
##################################################
calcna = function(time,event) {
  na.fit = survfit(coxph(Surv(time,event)~1), type="aalen")
  jumps = c(0, na.fit$time, max(time))
  # need to be careful at the beginning and end
  surv = c(1, na.fit$surv, na.fit$surv[length(na.fit$surv)])
  
  # apply appropriate transformation
  neglogsurv = -log(surv)
  
  # create placeholder of correct length
  naest = numeric(length(time))
  for (i in 2:length(jumps)) {
    naest[which(time>=jumps[i-1] & time<=jumps[i])] =
      neglogsurv[i-1]   # snag the appropriate value
  }
  return(naest)
}
Hazard1 = calcna(Time,(status==1)^2)
hazard1=cbind(Hazard1)
Hazard2=calcna(Time,(status==2)^2)
hazard2=cbind(Hazard2)
#### Joint Hazard Function #######
HA=hazard1+hazard2
################ Joint Survival Function #########################
##################################################################  
ssurv=exp(-HA)
##### Joint Cumulative Hazard Function ####
CIF=1-ssurv
Log_Odd=log((1-CIF)/CIF)
### Empirical Distribution of Log-Odds Residual 

########### Residual Plots #####################
###### Log-Odds Residual 
par(mfrow=c(2,2))
plot(Time,sort(Log_Odd) ,
     xlab="Time in Days)", ylab="Log-Odds Residuals", col="blue",
     main=' Log-Odds Residual for PBC Data',cex.lab=1.5, cex.axis=1.3)

lines(lowess(Time, Log_Odd ),col='red',lwd=3 )
qplot(Time,sort(Log_Odd), geom=c("point","smooth"),xlab="Waiting Time to Transplant Days",ylab="Log-Odds Residuals",
      main="Plot of Log-Odds Residuals against Patient Waiting Time") 
##### Martingale Residual
fit2=coxph(Surv(Time,(status==2)^2)~1,method='breslow')
rr=resid(fit2,type='martingale')
plot(Time, sort(rr),xlab="Waiting Time to Death(Days)", ylab="Martingale Residuals",
     main='Martingale Residual for PBC Data',cex.lab=1.5, cex.axis=1.3)
lines(lowess(Time, resid(fit2)),col='red',lwd=3)
qplot(Time,sort(rr), geom=c("point","smooth"),xlab="Waiting Time to Transplant Days",ylab="Martingale Residuals", 
      main="Plot of martingale residuals against Patient Waiting Time") 
###
plot(Time, sort(rr),xlab="Waiting Time to Death(Days)", ylab="Martingale Residuals",
     main='Martingale Residual for PBC Data',cex.lab=1.5, cex.axis=1.3,col=ifelse(rr>=mean(rr), "blue", "red"),
     pch=ifelse(rr>=mean(rr), 19, 17),cex =1.0)
lines(lowess(Time, resid(fit2)),col='brown',lwd=3)
#Add Legend to graph.  You can change the size of the box by changing cex = 0.75  Large # makes it larger.
legend("topright", c("Smoothing Curve", "Cenored", "Uncensored"), col = c("brown", "blue","red"), cex = 1.0,
       text.col = "black", lty = c(1 ,-1, -1), pch = c(-1, 19, 17),
       merge = TRUE, bg = 'gray90',lwd=3)
plot(Year, EMD, col=ifelse(D_EMD, "black", "red"),ylab = "EMD mg/L", pch=ifelse(D_EMD, 19, 17), cex = 0.7)

# Apply loess smoothing using the default span value of 0.8.  You can change the curve by changing the span value.
y.loess <- loess(y ~ x, span=0.8, data.frame(x=Year, y=EMD))
######  Residual plots against Covariates ###
library(ggplot2)
par(mfrow=c(1,3))
qplot(Age,rr, geom=c("point","smooth"),xlab="Patient Age",ylab="Martingale Residuals", main="Figure 7: Plot of martingale residuals against Age")   
qplot(Age,Log_Odd, geom=c("point","smooth"),xlab="Patient Age",ylab="Log-Odds Residuals", main="Figure 9: Plot of Log-Odds  residuals against Age") 
######## Data frame containing Log-Odss Residual
dat_LO=data.frame(Log_Odd,Age,sex,ascites,hepato,spiders,edema,bili,chol,albumin,copper,alk.phos,ast,trig,platelet,protime,stage)
dat_LO <- as.data.frame(dat_LO)
names(dat_LO)[1] <- "LOresiduals"
head(dat_LO)
set.seed(1245)
P=0.85
trainIndex<-sample.int(n = nrow(dat_LO), size = floor(P*nrow(dat_LO)), replace = F)
trainData_LO <- dat_LO[trainIndex,]
testData_LO<- dat_LO[-trainIndex,]
trainX_LO <-trainData_LO[,-1] # Create training feature data frame
testX_LO <- testData_LO[,-1] # Create test feature data frame
trainY_LO=trainData_LO[,1]
datatrain_LO=data.frame(trainY_LO,trainX_LO)
############ Decision Tree #############
## Sum of Square 
ModelLO_SS<-rpart(trainY_LO~. ,method="anova",data = datatrain_LO, maxdepth=20,minbucket=20)
printcp(ModelLO_SS)
plotcp(ModelLO_SS) # visualize cross-validation results
summary(ModelLO_SS)
LOSS <- as.party(ModelLO_SS)
LOSS
plot(LOSS)
#### performance Using Root Mean Square Error 
Predc_LOSS <- predict(ModelLO_SS,testData_LO)
MSESS_LO<-mean((Predc_LOSS-testData_LO$LOresiduals)^2)
round(MSESS_LO,4)
RMSESS_LO<- sqrt(mean((Predc_LOSS-testData_LO$LOresiduals)^2))
round(RMSESS_LO,4)
##### Absolute Sum 
ModelLO_ABS<-rpart(trainY_LO~. ,method="anova", data = datatrain_LO,maxdepth=20,cp=0,minbucket=20,dissim="man")
printcp(ModelLO_ABS)
plotcp(ModelLO_ABS) # visualize cross-validation results
summary(ModelLO_ABS)
LOABS <- as.party(ModelLO_ABS)
LOABS
plot(LOABS)
#### performance Using Root Mean Square Error 
Predc_LOABS <- predict(ModelLO_ABS,testData_LO)
MSEABS_LO<-mean((Predc_LOABS-testData_LO$LOresiduals)^2)
round(MSEABS_LO,4)
RMSEABS_LO<- sqrt(mean((Predc_LOABS-testData_LO$LOresiduals)^2))
round(RMSEABS_LO,4)
############## Classical Method ################################
################################################################
#library(pec)
Status=Status
Pbc_dat=data.frame(Time,status,Age,sex,ascites,hepato,spiders,edema,bili,chol,albumin,copper,alk.phos,ast,trig,platelet,protime,stage)
##### Mean Square Error ########
P=0.75
trainIndex<-sample.int(n = nrow(Pbc_dat), size = floor(P*nrow(Pbc_dat)), replace = F)
trainData <- Pbc_dat[trainIndex,]
testData<- Pbc_dat[-trainIndex,]
trainX <-trainData[,-c(1,2)] # Create training feature data frame
testX <- testData[,-c(1,2)] # Create test feature data frame
testY=testData[,1]
traintime= trainData[,1]
trainstatus= trainData[,2]
testtime= testData[,1]
teststaus= (testData[,2]==2)^2
datatrain=data.frame(traintime,trainstatus,trainX)
datatest=data.frame(testtime,teststaus,testX)
head(datatrain)
head(datatest)
train.fit <- coxph(Surv(traintime, (trainstatus==2)^2) ~ Age+sex+ascites+hepato+spiders+edema+bili+
chol+albumin+copper+alk.phos+ast+trig+platelet+protime+stage, data=datatrain,x=TRUE, y=TRUE, method="breslow")
### lpnew <- predict(train.fit, type="expected")
### lpnew <- predict(train.fit, type="expected")
#Predictedhazard <- predict(train.fit, newdata=testData)
lp.pred <- predict(train.fit,type="lp",newdata=testData)
# Baseline Function = 1
Pred.valhazard <- exp(lp.pred)
Pred.valsurv=exp(-Pred.valhazard) 
ClassicalMSE<-mean((testData[,2]-Pred.valsurv)^2)
round(ClassicalMSE,4)
######## Prediction Accuracy #################
Coxmodel1=coxph(Surv(traintime,(trainstatus==1)^2)~Age+sex+ascites+hepato+spiders+edema+bili+
                  chol+albumin+copper+alk.phos+ast+trig+platelet+protime+stage, data=datatrain)
Coxmodel2=coxph(Surv(traintime,(trainstatus==2)^2)~Age+sex+ascites+hepato+spiders+edema+bili+
                  chol+albumin+copper+alk.phos+ast+trig+platelet+protime+stage, data=datatrain)
library(prodlim)
library(survival)
library(pec)
# Cause 1
perror1=pec(list(Cox=Coxmodel1),Hist(traintime,(trainstatus==1)^2)~Age+sex+ascites+hepato+spiders+edema+bili+
             chol+albumin+copper+alk.phos+ast+trig+platelet+protime+stage, data=datatrain)
##  prediction error
ibs(perror1)
# Cause 2 
perror2=pec(list(Cox=Coxmodel2),Hist(traintime,(trainstatus==2)^2)~Age+sex+ascites+hepato+spiders+edema+bili+
             chol+albumin+copper+alk.phos+ast+trig+platelet+protime+stage, data=datatrain)
##  prediction error
ibs(perror2)
####### AUC #############
library(risksetROC)
trainstatus=(trainstatus==1)^2
coxfit <- coxph(Surv(traintime, (trainstatus==1)^2) ~ Age+sex+ascites+hepato+spiders+edema+bili+
                  chol+albumin+copper+alk.phos+ast+trig+platelet+protime+stage, data=datatrain)
eta= coxfit$linear.predictors
nobs=length(traintime[trainstatus==1])
span=1.0*(nobs^(-0.2))
ROCcox=risksetROC(Stime=traintime,status=trainstatus,marker=eta,predict.time=median(traintime),method="Cox",
                  main=" Classical Model ROC ", col="red", xlab="1-Specificity (FPR)", ylab="Sensitivity (TPR)")
######## JCIF ################### ################
CIF[median(CIF)<=CIF]<- 1
CIF[median(CIF)>CIF]<- 0
CIF=floor(CIF)
CIF=ifelse(CIF==1,"Lowrisk","Highrisk")
dat=data.frame(CIF,Age,sex,ascites,hepato,spiders,edema,bili,chol,albumin,copper,alk.phos,ast,trig,platelet,protime,stage)
dat <- as.data.frame(dat)
names(dat)[1] <- "CIF"
head(dat)                                                                                                                            
set.seed(34567)
P=0.75
trainIndex<-sample.int(n = nrow(dat), size = floor(P*nrow(dat)), replace = F)
trainData <- dat[trainIndex,]
testData<- dat[-trainIndex,]
trainX <-trainData[,-1] # Create training feature data frame
testX <- testData[,-1] # Create test feature data frame
trainY=trainData[,1]
testY=testData[,1]
datatrain=data.frame(trainY,trainX)
head(datatrain) 
############ AUC Decision Tree #############
## Sum of Square        randomForest
ModelJCIF_SS<-rpart(trainY~.,data = datatrain)
printcp(ModelJCIF_SS)
## randomForest
ModelJCIF_SS<-randomForest(trainY~.,data = datatrain)
#plotcp(ModelJCIF_SS) # visualize cross-validation results
library(caret)
rpartPred <- predict(ModelJCIF_SS, testData,type="class")
#rpartPred=ifelse(rpartPred==2,"Lowrisk","Highrisk")
confusionMatrix(rpartPred, testY) # requires 2 factor vectors
library(pROC)
library(ROCR)
#ModelJCIF_SS1<-rpart(trainY~. ,data = datatrain)
rpartProbs <- predict(ModelJCIF_SS, testData, type = "prob")
head(rpartProbs)
lev <- c("Highrisk","Lowrisk")
rpartROC <- roc(testY, rpartProbs[, "Highrisk"], levels = rev(lev))
plot(rpartROC, type = "S", print.thres = .5, col="blue",main="Marginal CIF ROC", 
     xlab="1-Specificity (FPR)", ylab="Sensitivity (TPR)")
rpartROC
###### Conditional ##############
CPF_2=read.csv("CPF_1.csv")
Conditional_IF=as.matrix(CPF_2)
CPF.mean=mean(Conditional_IF,na.rm=TRUE)
#Conditional_IF=data.frame(CPF)
Conditional_IF=ifelse(is.na(Conditional_IF==T), CPF.mean, Conditional_IF)
Conditional_IF[median(Conditional_IF)<=Conditional_IF]<- 1
Conditional_IF[median(Conditional_IF)>Conditional_IF]<- 0
Conditional_IF=floor(Conditional_IF)
Conditional_IF=ifelse(Conditional_IF==1,"Lowrisk","Highrisk")
datCPF=data.frame(Conditional_IF,Age,sex,ascites,hepato,spiders,edema,bili,chol,albumin,copper,alk.phos,ast,trig,platelet,protime,stage)
datCPF <- as.data.frame(datCPF)
names(datCPF)[1] <- "ConditionaIF"
head(datCPF)
set.seed(22345)
P=0.85
trainIndex<-sample.int(n = nrow(datCPF), size = floor(P*nrow(datCPF)), replace = F)
traindatCPF <- datCPF[trainIndex,]
testdatCPF<- datCPF[-trainIndex,]
trainX_CPF <-traindatCPF[,-1] # Create training feature data frame
testX_CPF <- testdatCPF[,-1] # Create test feature data frame
trainY_CPF=traindatCPF[,1]
testY_CPF=testdatCPF[,1]
datatrain=data.frame(trainY_CPF,trainX_CPF)
Model_CIF_SS<-rpart(trainY_CPF~.,data = datatrain)
printcp(Model_CIF_SS)
#plotcp(ModelJCIF_SS) # visualize cross-validation results
library(caret)
rpartPred_2 <- predict(Model_CIF_SS, testdatCPF,type="class")
#rpartPred=ifelse(rpartPred==2,"Lowrisk","Highrisk")
confusionMatrix(rpartPred_2, testY_CPF) # requires 2 factor vectors
rpartProbs_2 <- predict(Model_CIF_SS, testdatCPF, type = "prob")
#head(rpartProbs_2)
lev <- c("Highrisk","Lowrisk")
rpartROC_2 <- roc(testY_CPF, rpartProbs_2[, "Highrisk"], levels = rev(lev))
plot(rpartROC_2, type = "S", print.thres = 0.5, col="blue",main="Conditional CIF ROC", xlab="1-Specificity (FPR)", ylab="Sensitivity (TPR)")
rpartROC_2 
#################### Log-Odds Residual ##############
Log_Odd=log((1-CIF)/CIF)
Log_Odd[median(Log_Odd)<=Log_Odd]<- 1
Log_Odd[median(Log_Odd)>Log_Odd]<- 0
Log_Odd=floor(Log_Odd)
Log_Odd=ifelse(Log_Odd==1,"Lowrisk","Highrisk")
dat_LO=data.frame(Log_Odd,Age,sex,ascites,hepato,spiders,edema,bili,chol,albumin,copper,alk.phos,ast,trig,platelet,protime,stage)
dat_LO <- as.data.frame(dat_LO)
names(dat_LO)[1] <- "Log_Odd"
head(dat_LO)
set.seed(1245)
P=0.85
trainIndex<-sample.int(n = nrow(dat_LO), size = floor(P*nrow(dat_LO)), replace = F)
trainData_LO <- dat_LO[trainIndex,]
testData_LO<- dat_LO[-trainIndex,]
trainX_LO <-trainData_LO[,-1] # Create training feature data frame
testX_LO <- testData_LO[,-1] # Create test feature data frame
trainY_LO=trainData_LO[,1]
testY_LO=testData_LO[,1]
datatrain_LO=data.frame(trainY_LO,trainX_LO)
Model_Log<-rpart(trainY_LO~.,data = datatrain_LO)
printcp(Model_Log)
rpartPred_3 <- predict(Model_Log, testData_LO,type="class")
#rpartPred=ifelse(rpartPred==2,"Lowrisk","Highrisk")
confusionMatrix(rpartPred_3, testY_LO) # requires 2 factor vectors
rpartProbs_3 <- predict(Model_Log, testData_LO, type = "prob")
#head(rpartProbs_2)
lev <- c("Highrisk","Lowrisk")
rpartROC_3 <- roc(testY_LO, rpartProbs_3[, "Highrisk"], levels = rev(lev))
plot(rpartROC_3, type = "S", print.thres = 0.5, col="blue",main="Log-Odds CIF ROC",
     xlab="1-Specificity (FPR)", ylab="Sensitivity (TPR)")
rpartROC_3 
############ Tuning Model ##############
########################################
library(ipred)
library(e1071)
library(caret)
### Three repeats 10-fold Cross Validation ##
#cvCtrl <- trainControl(method = "repeatedcv", repeats = 3)
#train(Class ~ ., data = training, method = "rpart",
#tuneLength = 30,trControl = cvCtrl)
###########################################################
cvCtrl <- trainControl(method = "repeatedcv", repeats = 3,
                       summaryFunction = twoClassSummary,classProbs = TRUE)
set.seed(1234)
rpartTune <- train(Class ~ ., data = training, method = "rpart",
                   tuneLength = 30,metric = "ROC",trControl = cvCtrl)
######## JCIF ################### 
CIF[median(CIF)<=CIF]<- 1
CIF[median(CIF)>CIF]<- 0
CIF=floor(CIF)
CIF=ifelse(CIF==1,"Lowrisk","Highrisk")
dat=data.frame(CIF,Age,sex,ascites,hepato,spiders,edema,bili,chol,albumin,copper,alk.phos,ast,trig,platelet,protime,stage)
data_scaling=data.frame(Age,sex,ascites,hepato,spiders,edema,bili,chol,albumin,copper,alk.phos,ast,trig,platelet,protime,stage)
# Assuming goal class is column 10
#preObj <- preProcess(data_scaling, method=c("center", "scale"))BoxCox
preObj <- preProcess(data_scaling, method=c("center", "scale","BoxCox"))
newData <- predict(preObj,data_scaling ) 
dat=data.frame(CIF,newData) 
dat <- as.data.frame(dat)
names(dat)[1] <- "CIF"
head(dat)                                                                                                                         
set.seed(34567)
P=0.75
trainIndex<-sample.int(n = nrow(dat), size = floor(P*nrow(dat)), replace = F)
trainData <- dat[trainIndex,]
testData<- dat[-trainIndex,]
trainX <-trainData[,-1] # Create training feature data frame
testX <- testData[,-1] # Create test feature data frame
trainY=trainData[,1]
testY=testData[,1]
datatrain=data.frame(trainY,trainX)
head(datatrain) 
#cvCtrl <- trainControl(method = "repeatedcv", repeats = 3,summaryFunction = twoClassSummary,classProbs = TRUE)
cvCtrl <- trainControl(method = "LGOCV",summaryFunction = twoClassSummary,classProbs = TRUE,
                       savePredictions = TRUE,number = 1000)
set.seed(1234)
rpartTune_JCIF <- train(trainY ~ ., data = datatrain, method = "rf",
                        tuneLength = 80,metric = "ROC",trControl = cvCtrl)
## plot
plot(rpartTune_JCIF, scales = list(x = list(log = 10)))
### Prediction
#predict(rpartTune_JCIF, newdata)
rpartPred_JCIF<-predict(rpartTune_JCIF,testData)
confusionMatrix(rpartPred_JCIF,testY)
#### Prediction Probability 
rpartProbs_CIF<-predict(rpartTune_JCIF,testData,type="prob")
#head(rpartProbs)
### ROC
library(pROC)
lev <- c("Highrisk","Lowrisk")
rpartROC_CIF<-roc(testY,rpartProbs_CIF[,"Highrisk"],levels=rev(lev))
plot(rpartROC_CIF,type="S",main="Marginal CIF ROC", 
     xlab="1-Specificity (FPR)", ylab="Sensitivity (TPR)")
rpartROC_CIF
### Boosting 
dat_cif=data.frame(CIF,Age,sex,ascites,hepato,spiders,edema,bili,chol,albumin,copper,alk.phos,ast,trig,platelet,protime,stage)
P=0.75
trainIndex<-sample.int(n = nrow(dat_cif), size = floor(P*nrow(dat_cif)), replace = F)
trainData <- dat_cif[trainIndex,]
testData<- dat_cif[-trainIndex,]
trainX <-trainData[,-1] # Create training feature data frame
testX <- testData[,-1] # Create test feature data frame
trainY=trainData[,1]
testY=testData[,1]
datatrain=data.frame(trainY,trainX)
X=dat_cif[,2:16]
Y= dat_cif[,1]
grid <- expand.grid(.model = "tree",.trials = c(1:100),.winnow = FALSE)
c5Tune <- train(X, Y ,method = "C5.0",preProcess = c("center", "scale"),metric = "Accuracy",tuneGrid = grid,trControl = cvCtrl)
c5Pred<-predict(c5Tune,dat_cif)
confusionMatrix(c5Pred,Y)
###### Conditional ##############
CPF_2=read.csv("CPF_1.csv")
Conditional_IF=as.matrix(CPF_2)
CPF.mean=mean(Conditional_IF,na.rm=TRUE)
#Conditional_IF=data.frame(CPF)
Conditional_IF=ifelse(is.na(Conditional_IF==T), CPF.mean, Conditional_IF)
Conditional_IF[median(Conditional_IF)<=Conditional_IF]<- 1
Conditional_IF[median(Conditional_IF)>Conditional_IF]<- 0
Conditional_IF=floor(Conditional_IF)
Conditional_IF=ifelse(Conditional_IF==1,"Lowrisk","Highrisk")
datCPF=data.frame(Conditional_IF,Age,sex,ascites,hepato,spiders,edema,bili,chol,albumin,copper,alk.phos,ast,trig,platelet,protime,stage)
data_scaling_CPF=data.frame(Age,sex,ascites,hepato,spiders,edema,bili,chol,albumin,copper,alk.phos,ast,trig,platelet,protime,stage)
# Assuming goal class is column 10
#preObj <- preProcess(data_scaling, method=c("center", "scale"))BoxCox
preObj_CPF <- preProcess(data_scaling_CPF, method=c("center", "scale","BoxCox"))
newDataCPF <- predict(preObj_CPF ,data_scaling_CPF) 
datCPF=data.frame(Conditional_IF,newDataCPF) 
datCPF <- as.data.frame(dat)                 
names(datCPF)[1] <- "ConditionaIF"
head(datCPF)
set.seed(22345)
P=0.75
trainIndex<-sample.int(n = nrow(datCPF), size = floor(P*nrow(datCPF)), replace = F)
traindatCPF <- datCPF[trainIndex,]
testdatCPF<- datCPF[-trainIndex,]
trainX_CPF <-traindatCPF[,-1] # Create training feature data frame
testX_CPF <- testdatCPF[,-1] # Create test feature data frame
trainY_CPF=traindatCPF[,1]
testY_CPF=testdatCPF[,1]
datatrainCPF=data.frame(trainY_CPF,trainX_CPF)
rpartTune_CPF <- train(trainY_CPF ~ ., data = datatrainCPF, method = "rf",
                       tuneLength = 80,metric = "ROC",trControl = cvCtrl)
## plot
plot(rpartTune_CPF, scales = list(x = list(log = 10)))
### Prediction
#predict(rpartTune_JCIF, newdata)
rpartPred_CPF<-predict(rpartTune_CPF,testdatCPF)
confusionMatrix(rpartPred_CPF,testY_CPF)
#### Prediction Probability 
rpartProbs_CPF<-predict(rpartTune_CPF,testdatCPF,type="prob")
#head(rpartProbs)
### ROC
library(pROC)
lev <- c("Highrisk","Lowrisk")
rpartROC_CPF<-roc(testY_CPF,rpartProbs_CPF[,"Highrisk"],levels=rev(lev))
plot(rpartROC_CPF,type="S",main="Conditional CIF ROC", xlab="1-Specificity (FPR)", ylab="Sensitivity (TPR)")
rpartROC_CPF
######### Log-Odds Residual #########
####################################
Log_Odd=log((1-CIF)/CIF)
Log_Odd[median(Log_Odd)<=Log_Odd]<- 1
Log_Odd[median(Log_Odd)>Log_Odd]<- 0
Log_Odd=floor(Log_Odd)
Log_Odd=ifelse(Log_Odd==1,"Lowrisk","Highrisk")
dat_LO=data.frame(Log_Odd,Age,sex,ascites,hepato,spiders,edema,bili,chol,albumin,copper,alk.phos,ast,trig,platelet,protime,stage)
dat_LO <- as.data.frame(dat_LO)
names(dat_LO)[1] <- "Log_Odd"
head(dat_LO)
set.seed(1245)
P=0.75
trainIndex<-sample.int(n = nrow(dat_LO), size = floor(P*nrow(dat_LO)), replace = F)
trainData_LO <- dat_LO[trainIndex,]
testData_LO<- dat_LO[-trainIndex,]
trainX_LO <-trainData_LO[,-1] # Create training feature data frame
testX_LO <- testData_LO[,-1] # Create test feature data frame
trainY_LO=trainData_LO[,1]
testY_LO=testData_LO[,1]
datatrain_LO=data.frame(trainY_LO,trainX_LO) 
rpartTune_LO <- train(trainY_LO ~ ., data = datatrain_LO, method = "rf",
                      tuneLength = 80,metric = "ROC",trControl = cvCtrl)
## plot
#plot(rpartTune_CPF, scales = list(x = list(log = 10)))
### Prediction
#predict(rpartTune_JCIF, newdata)
rpartPred_LO<-predict(rpartTune_LO,testData_LO)
confusionMatrix(rpartPred_LO,testY_LO)
#### Prediction Probability 
rpartProbs_LO<-predict(rpartTune_LO,testData_LO,type="prob")
#head(rpartProbs)
### ROC
library(pROC)
lev <- c("Highrisk","Lowrisk")
rpartROC_LO<-roc(testY_LO,rpartProbs_LO[,"Highrisk"],levels=rev(lev))
plot(rpartROC_LO,type="S",main="Log-Odds CIF ROC",
     xlab="1-Specificity (FPR)", ylab="Sensitivity (TPR)")
rpartROC_LO

######### Total plot 
plot(rpartROC_CIF,type="S", xlab="1-Specificity (FPR)", ylab="Sensitivity (TPR)", col="blue",lty=1,lwd=3 )
lines(rpartROC_CPF,type="S",main="Conditional CIF ROC", 
      xlab="1-Specificity (FPR)", ylab="Sensitivity (TPR)", col="black",lty=2, lwd=3)
lines(rpartROC_LO,type="S",main="Log-Odds CIF ROC",
      xlab="1-Specificity (FPR)", ylab="Sensitivity (TPR)", col="red",lty=3, lwd=3)
legend("bottomright",c("Log-Odds CIF","Conditional CIF ","Marginal CIF ROC"),lwd=3,lty=1:3,col=c("blue","black","red"))