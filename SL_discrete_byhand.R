#################################################################################
#################################################################################
#Super learner for time-to-event outcomes tutorial
#
#Method: discrete-time superLearner
#Polley and van der Laan 2011
#
#Implemented by hand 
#################################################################################
#################################################################################

source("packages.R")
source("rotterdam_setup.R")
source("functions.R")
source("censoring_weights.R")


#---
#divide training data into 5 folds
n<-dim(dta_train)[1]
n.fold<-5
set.seed(1430)
dta_train$cv.group<-cut(runif(n,0,n.fold),breaks=seq(0,n.fold,1),include.lowest = T,labels = F)

#---
#put data into counting process format using periods of length 0.1 years
dta_train_split<-survSplit(Surv(time,status)~.,data=dta_train,cut = seq(0,max(dta_train$time),0.1),start="tstart",end="tstop",event="status")

#generate time^2 and time^3 variables to be used as covariates in the SL
dta_train_split$tstart2<-dta_train_split$tstart^2
dta_train_split$tstart3<-dta_train_split$tstart^3

#---
#data that we will use to predict hazards at all times
times<-data.frame(tstart=seq(0,max(dta_train$time),0.1))
dta_train_haz<-merge(dta_train,times)
dta_train_haz<-dta_train_haz[order(dta_train_haz$pid,dta_train_haz$tstart),]
dta_train_haz$tstart2<-dta_train_haz$tstart^2
dta_train_haz$tstart3<-dta_train_haz$tstart^3

#---
#storing predicted hazards and survival probabilities (cross-fitted) from each method
n.methods=4
haz.pred<-matrix(nrow=dim(dta_train_split)[1],ncol=n.methods)
haz.pred.all<-matrix(nrow=dim(dta_train_haz)[1],ncol=n.methods)
surv.pred<-matrix(nrow=n,ncol=n.methods)

#storing expected loss from each method and each fold
loss.L2<-matrix(nrow=n.fold,ncol=n.methods)
loss.loglik<-matrix(nrow=n.fold,ncol=n.methods)
loss.ipcw<-matrix(nrow=n.fold,ncol=n.methods)

#storing censoring weights
# cens.wt<-rep(NA,n.fold)
cens.wt<-rep(NA,n)

#---
#start of cross-validation loop
for(i in 1:n.fold){
  print(i)
  
  cvtrain<-dta_train_split[!dta_train_split$cv.group==i,] #data used to fit models
  cvtest<-dta_train_haz[dta_train_haz$cv.group==i,] #test data - note this has 1 row for each time for each person (including times after the observed event/cens time)
  cvtest.rows<-which(dta_train_haz$cv.group==i) #used for storing results
  
  cvtest.1row<-dta_train[dta_train$cv.group==i,]#test data with 1 row per person, used for estimating expected loss
  cvtest.split<-dta_train_split[dta_train_split$cv.group==i,] #test data with multiple rows per person, used for estimating expected loss
  cvtest.rows.1row<-which(dta_train$cv.group==i) #used for storing results
  
  #-------------------------------------------------------
  #Fit candidate models and obtain cross-fitted estimates
  #for individuals in the excluded fold we obtain: (a) the hazard at each time, (b)the survival probability at time 3000
  #-------------------------------------------------------
  
  #---
  #method 1: glm with no covariates
  method.glm0<-glm(status~1,data=cvtrain,family = "binomial")
  haz.pred.all[cvtest.rows,1]<-predict(method.glm0,newdata=cvtest,type="response")
  haz.pred[dta_train_split$cv.group==i,1]<-haz.pred.all[cvtest.rows,1][cvtest$tstart<cvtest$time]
  surv.pred[dta_train$cv.group==i,1]<-ave(1-haz.pred.all[cvtest.rows,1],cvtest$pid,FUN=cumprod)[cvtest$tstart==10]
  
  #---
  #method 2: glm including time plus covariates
  method.glm1<-glm(status~tstart+tstart2+tstart3+year1+year2+age+meno+size1+size2+grade+nodes+pgr+er+hormon+chemo,
                   data=cvtrain,family = "binomial")
  haz.pred.all[cvtest.rows,2]<-predict(method.glm1,newdata=cvtest,type="response")
  haz.pred[dta_train_split$cv.group==i,2]<-haz.pred.all[cvtest.rows,2][cvtest$tstart<cvtest$time]
  surv.pred[dta_train$cv.group==i,2]<-ave(1-haz.pred.all[cvtest.rows,2],cvtest$pid,FUN=cumprod)[cvtest$tstart==10]

  #---
  #method 3: lasso
  xmat.cvtrain<-model.matrix(method.glm1)[,-1]
  xmat.cvtest<-model.matrix(method.glm1,data=cvtest)[,-1]
  method.lasso <- cv.glmnet(x=xmat.cvtrain,y=cvtrain$status,family = "binomial",alpha=1)
  haz.pred.all[cvtest.rows,3]<-predict(method.lasso, newx = xmat.cvtest, s = "lambda.min",type="response")
  haz.pred[dta_train_split$cv.group==i,3]<-haz.pred.all[cvtest.rows,3][cvtest$tstart<cvtest$time]
  surv.pred[dta_train$cv.group==i,3]<-ave(1-haz.pred.all[cvtest.rows,3],cvtest$pid,FUN=cumprod)[cvtest$tstart==10]
  
  #---
  #method 4: random forest
  method.rf <- ranger(status ~ tstart+year1+year2+age+meno+size1+size2+grade+nodes+pgr+er+hormon+chemo,
                      data=cvtrain,oob.error = FALSE)
  haz.pred.all[cvtest.rows,4] <- predict(method.rf, data = cvtest)$predictions
  haz.pred[dta_train_split$cv.group==i,4]<-haz.pred.all[cvtest.rows,4][cvtest$tstart<cvtest$time]
  surv.pred[dta_train$cv.group==i,4]<-ave(1-haz.pred.all[cvtest.rows,4],cvtest$pid,FUN=cumprod)[cvtest$tstart==10]
  
  #-------------------------------------------------------
  #Calculate expected loss for each method
  #Using loss functions: L2 version, loglik version, IPCW version
  #-------------------------------------------------------
  
  #estimate censoring survival distribution (assuming independent censoring), which is needed as a weight in the ipcw loss function
  cens.km<-survfit(Surv(time,1-status)~1,data=cvtest.1row)
  # cens.wt[i]<-1/summary(cens.km,times = 10)$surv
  
  cens.km.step<-stepfun(cens.km$time,c(1,cens.km$surv))
  cens.wt[cvtest.rows.1row]<-1/cens.km.step(pmin(cvtest.1row$time,10))
  
  for(j in 1:n.methods){
    #L2 loss
    loss.L2[i,j]<-mean((cvtest.split$status-haz.pred[dta_train_split$cv.group==i,j])^2)
    
    #loglik loss
    loss.loglik[i,j]<-mean((log(haz.pred[dta_train_split$cv.group==i,j])^cvtest.split$status)*
                        (log(1-haz.pred[dta_train_split$cv.group==i,j])^(1-cvtest.split$status)))
    
    # #ipcw loss
    # loss.ipcw[i,j]<-mean(cens.wt[i]*I(cvtest.1row$status==1|(cvtest.1row$status==0 & cvtest.1row$time>10))*
    #                       (I(cvtest.1row$time>10)-surv.pred[dta_train$cv.group==i,j])^2)
    
    #ipcw loss
    loss.ipcw[i,j]<-mean(cens.wt[cvtest.rows.1row]*I(cvtest.1row$status==1|(cvtest.1row$status==0 & cvtest.1row$time>10))*
                           (I(cvtest.1row$time>10)-surv.pred[dta_train$cv.group==i,j])^2)
  }
}

#-------------------------------------------------------
#For discrete SL: take a look at the mean expected loss for each method
#-------------------------------------------------------

meanloss.L2<-colMeans(loss.L2)
meanloss.loglik<-colMeans(loss.loglik)
meanloss.ipcw<-colMeans(loss.ipcw)

discrete.SL<-cbind(meanloss.L2,meanloss.loglik,meanloss.ipcw)
rownames(discrete.SL)<-c("1. GLM: mean only",
                         "2. GLM: time+covariates",
                         "3. GLM lasso: time+covariates",
                         "4. Random forest: time+covariates")
discrete.SL

#-------------------------------------------------------
#ensemble step
#Using loss functions: L2 version, loglik version, IPCW version
#-------------------------------------------------------

#---
#ensemble based on L2 loss function, using non-negative least squares
comb.mod.L2<-nnls(b=dta_train_split$status,A=haz.pred)

#scale the coefficients to sum to 1
weights.L2<-comb.mod.L2$x/sum(comb.mod.L2$x)

#---
#ensemble based on loglik loss function

#obtain logit transformations of the estimated hazards
g.haz.loglik<-log(haz.pred/(1-haz.pred))

#fit combining model, without an intercept
comb.mod.loglik<-glm(dta_train_split$status~-1+g.haz.loglik,family="binomial")

#no need to scale the coefficients for this loss function
weights.loglik<-comb.mod.loglik$coefficients

#---
#ensemble based on ipcw loss function, using non-negative weighted least squares
#This is achieved by multiplying the outcome and explanatory variables by the square root of the weights. 
#https://stackoverflow.com/questions/47888996/weighted-nonnegative-least-squares-in-r
# cens.wt.dta<-data.frame(cens.wt=cens.wt,cv.group=1:n.fold)
# dta_train2<-merge(dta_train,cens.wt.dta)

dta_train2<-dta_train
dta_train2$cens.wt<-cens.wt

dta_train2<-dta_train2[order(dta_train2$pid),]
dta_train2$elig<-I(dta_train2$status==1|(dta_train2$status==0 & dta_train2$time>10))
comb.mod.ipcw<-nnls(b=I(dta_train2$status[dta_train2$elig==1]==0)*sqrt(dta_train2$cens.wt[dta_train2$elig==1]),
                    A=surv.pred[dta_train2$elig==1,]*sqrt(dta_train2$cens.wt[dta_train2$elig==1]))

#scale the coefficients to sum to 1
weights.ipcw<-comb.mod.ipcw$x/sum(comb.mod.ipcw$x)

#---
#Look at weights for ensembles based on different loss functions

ensemble.SL<-cbind(weights.L2,weights.loglik,weights.ipcw)
rownames(ensemble.SL)<-c("1. GLM: mean only",
                         "2. GLM: time+covariates",
                         "3. GLM lasso: time+covariates",
                         "4. Random forest: time+covariates")
ensemble.SL

#-------------------------------------------------------
#-------------------------------------------------------
#Obtaining predictions in the test data based on the ensemble
#-------------------------------------------------------
#-------------------------------------------------------

#---
#first fit each candidate model on the full training data, and obtain predictions for people in dta_test
#---

#---
#generate test data in the long format needed to obtain predictions
dta_test_haz<-merge(dta_test,times)
dta_test_haz<-dta_test_haz[order(dta_test_haz$pid,dta_test_haz$tstart),]
dta_test_haz$tstart2<-dta_test_haz$tstart^2
dta_test_haz$tstart3<-dta_test_haz$tstart^3

#---
#storing predicted hazards from each method
n.methods=4
haz.pred.all<-matrix(nrow=dim(dta_test_haz)[1],ncol=n.methods)

#storing predicted survival probabilities from each method
n.test<-dim(dta_test)[1]
surv.pred<-matrix(nrow=n.test,ncol=n.methods)

#---
#method 1: glm with no covariates
method.glm0<-glm(status~1,data=dta_train_split,family = "binomial")
haz.pred.all[,1]<-predict(method.glm0,newdata=dta_test_haz,type="response")
surv.pred[,1]<-ave(1-haz.pred.all[,1],dta_test_haz$pid,FUN=cumprod)[dta_test_haz$tstart==10]

#---
#method 2: glm including time plus covariates
method.glm1<-glm(status~tstart+tstart2+tstart3+year1+year2+age+meno+size1+size2+grade+nodes+pgr+er+hormon+chemo,
                 data=dta_train_split,family = "binomial")
haz.pred.all[,2]<-predict(method.glm1,newdata=dta_test_haz,type="response")
surv.pred[,2]<-ave(1-haz.pred.all[,2],dta_test_haz$pid,FUN=cumprod)[dta_test_haz$tstart==10]

#---
#method 3: lasso
xmat.train<-model.matrix(method.glm1)[,-1]
dta_test_haz$status=99#just a trick so that we can use model.matrix to get xmat.test
xmat.test<-model.matrix(method.glm1,data=dta_test_haz)[,-1]
method.lasso <- cv.glmnet(x=xmat.train,y=dta_train_split$status,family = "binomial",alpha=1)
haz.pred.all[,3]<-predict(method.lasso, newx = xmat.test, s = "lambda.min",type="response")
surv.pred[,3]<-ave(1-haz.pred.all[,3],dta_test_haz$pid,FUN=cumprod)[dta_test_haz$tstart==10]

#---
#method 4: random forest
method.rf <- ranger(status ~ tstart+year1+year2+age+meno+size1+size2+grade+nodes+pgr+er+hormon+chemo,
                    data=dta_train_split,oob.error = FALSE)
haz.pred.all[,4] <- predict(method.rf, data = dta_test_haz)$predictions
surv.pred[,4]<-ave(1-haz.pred.all[,4],dta_test_haz$pid,FUN=cumprod)[dta_test_haz$tstart==10]

#---
#now combine using the ensemble weights
#---

#---
#ensemble SL based on L2 loss function
haz.pred.SL.L2<-weights.L2[1]*haz.pred.all[,1]+weights.L2[2]*haz.pred.all[,2]+
  weights.L2[3]*haz.pred.all[,3]+weights.L2[4]*haz.pred.all[,4]

surv.pred.SL.L2<-ave(1-haz.pred.SL.L2,dta_test_haz$pid,FUN=cumprod)[dta_test_haz$tstart==10]

#---
#ensemble SL based on loglik loss function

#obtain logit transformations of the estimated hazards
logit<-function(x){log(x/(1-x))}
expit<-function(x){exp(x)/(1+exp(x))}

g.haz.pred.SL.loglik<-
  weights.loglik[1]*logit(haz.pred.all[,1])+
  weights.loglik[2]*logit(haz.pred.all[,2])+
  weights.loglik[3]*logit(haz.pred.all[,3])+
  weights.loglik[4]*logit(haz.pred.all[,4])

haz.pred.SL.loglik<-expit(g.haz.pred.SL.loglik)

surv.pred.SL.loglik<-ave(1-haz.pred.SL.loglik,dta_test_haz$pid,FUN=cumprod)[dta_test_haz$tstart==10]

#---
#ensemble SL based on ipcw

surv.pred.SL.ipcw<-weights.ipcw[1]*surv.pred[,1]+
  weights.ipcw[2]*surv.pred[,2]+
  weights.ipcw[3]*surv.pred[,3]+
  weights.ipcw[4]*surv.pred[,4]

#-------------------------------------------------------
#-------------------------------------------------------
#obtaining measures of predictive performance
#-------------------------------------------------------
#-------------------------------------------------------

#the line below can be altered depending on which SL we want to evaluate (surv.pred.SL.L2,surv.pred.SL.loglik,surv.pred.SL.ipcw)
risk.pred<-1-surv.pred.SL.ipcw

#---------------------------------
#Calibration plots
#---------------------------------

#---
#calibration plot
cut_points=c(0,quantile(risk.pred,probs=seq(0.1,0.9,0.1)),1)
risk_group=cut(risk.pred,breaks=cut_points,include.lowest = T,labels = F)
calib_risk_group<-sapply(FUN=function(x){mean(risk.pred[risk_group==x])},1:10)

km.grp=survfit(Surv(time,status)~strata(risk_group),data=dta_test)

risk_obs_grp<-rep(NA,10)
for(k in 1:10){
  group.exists<-paste0("strata(risk_group)=risk_group=",k)%in%summary(km.grp,cens=T)$strata
  if(group.exists){
    step.grp=stepfun(km.grp$time[summary(km.grp,cens=T)$strata==paste0("strata(risk_group)=risk_group=",k)],
                     c(1,km.grp$surv[summary(km.grp,cens=T)$strata==paste0("strata(risk_group)=risk_group=",k)]))
    risk_obs_grp[k]<-1-step.grp(10)
  }
}

plot(calib_risk_group,risk_obs_grp,type="both",xlab="Predicted risk",ylab="Estimated actual risk",xlim=c(0,1),ylim=c(0,1))
abline(0,1)

#---------------------------------
#Brier score, IPA, and integrated Brier score
#---------------------------------

#---
#Brier score and IPA

Brier(time=dta_test$time, status=dta_test$status, risk=risk.pred, 
      seq.time=10, weights=dta_test$cens.wt)

ipa(time=dta_test$time, status=dta_test$status, risk=risk.pred, 
    seq.time=10, weights=dta_test$cens.wt)

#---------------------------------
#C-index, AUC and AUCt
#---------------------------------

#---
#C-index - using our function
c_index_ties(time=dta_test$time,status=dta_test$status, risk=risk.pred, tau=10, weightmatrix = wt_matrix_eventsonly)

#c-index - using concordance
concordance(Surv(dta_test$time, dta_test$status) ~ risk.pred,
            newdata=dta_test,
            reverse = TRUE,
            timewt = "n/G2")$concordance

#---
#AUC - using timeROC
timeROC( T = dta_test$time,delta = dta_test$status,marker = risk.pred,
         cause = 1,weighting = "marginal",times = 10,iid = FALSE)

#---
#C/D AUCt - using our function
max.event.time<-max(dta_test$time[dta_test$status==1])
wCD_AUCt(time=dta_test$time,status=dta_test$status, risk=risk.pred, seq.time =max.event.time, weightmatrix = wt_matrix_eventsonly)





