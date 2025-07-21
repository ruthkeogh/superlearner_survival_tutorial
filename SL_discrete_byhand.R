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

#---
#divide training data into 5 folds
n<-dim(dta_train)[1]
n.fold<-5
set.seed(1430)
dta_train$cv.group<-cut(runif(n,0,n.fold),breaks=seq(0,n.fold,1),include.lowest = T,labels = F)

#---
#put data into counting process format using 30 day periods
dta.split<-survSplit(Surv(time,status)~.,data=dta_train,cut = seq(0,max(dta_train$time),0.1),start="tstart",end="tstop",event="status")

#generate time^2 and time^3 variables to be used as covariates in the SL
dta.split$tstart2<-dta.split$tstart^2
dta.split$tstart3<-dta.split$tstart^3

#---
#data that we will use to predict hazards at all times
times<-data.frame(tstart=seq(0,max(dta_train$time),0.1))
dta.haz<-merge(dta_train,times)
dta.haz<-dta.haz[order(dta.haz$pid,dta.haz$tstart),]
dta.haz$tstart2<-dta.haz$tstart^2
dta.haz$tstart3<-dta.haz$tstart^3

#---
#storing predicted hazards from each method
n.methods=4
haz.pred<-matrix(nrow=dim(dta.split)[1],ncol=n.methods)
haz.pred.all<-matrix(nrow=dim(dta.haz)[1],ncol=n.methods)

#storing predicted survival probabilities from each method
surv.pred<-matrix(nrow=n,ncol=n.methods)

#storing expected loss from each method and each fold
loss.L2<-matrix(nrow=n.fold,ncol=n.methods)
loss.loglik<-matrix(nrow=n.fold,ncol=n.methods)
loss.ipcw<-matrix(nrow=n.fold,ncol=n.methods)

#---
#start of cross-validation loop
for(i in 1:n.fold){
  print(i)
  
  train<-dta.split[!dta.split$cv.group==i,] #data used to fit models
  test<-dta.haz[dta.haz$cv.group==i,] #test data - note this has 1 row for each time for each person (including times after the observed event/cens time)
  test.rows<-which(dta.haz$cv.group==i) #used for storing results
  
  test.1row<-dta_train[dta_train$cv.group==i,]#test data with 1 row per person, used for estimating expected loss
  test.split<-dta.split[dta.split$cv.group==i,] #test data with multiple rows per person, used for estimating expected loss
  
  #-------------------------------------------------------
  #Fit candidate models and obtain cross-fitted estimates
  #for individuals in the excluded fold we obtain: (a) the hazard at each time, (b)the survival probability at time 3000
  #-------------------------------------------------------
  
  #---
  #method 1: glm with no covariates
  method.glm0<-glm(status~1,data=train,family = "binomial")
  haz.pred.all[test.rows,1]<-predict(method.glm0,newdata=test,type="response")
  haz.pred[dta.split$cv.group==i,1]<-haz.pred.all[test.rows,1][test$tstart<test$time]
  surv.pred[dta_train$cv.group==i,1]<-ave(1-haz.pred.all[test.rows,1],test$pid,FUN=cumprod)[test$tstart==10]
  
  #---
  #method 2: glm including time plus covariates
  method.glm1<-glm(status~tstart+tstart2+tstart3+year1+year2+age+meno+size1+size2+grade+nodes+pgr+er+hormon+chemo,
                   data=train,family = "binomial")
  haz.pred.all[test.rows,2]<-predict(method.glm1,newdata=test,type="response")
  haz.pred[dta.split$cv.group==i,2]<-haz.pred.all[test.rows,2][test$tstart<test$time]
  surv.pred[dta_train$cv.group==i,2]<-ave(1-haz.pred.all[test.rows,2],test$pid,FUN=cumprod)[test$tstart==10]

  #---
  #method 3: lasso
  xmat.train<-model.matrix(method.glm1)[,-1]
  test$status=99#just a trick so that we can use model.matrix to get xmat.test
  xmat.test<-model.matrix(method.glm1,data=test)[,-1]
  method.lasso <- cv.glmnet(x=xmat.train,y=train$status,family = "binomial",alpha=1)
  haz.pred.all[test.rows,3]<-predict(method.lasso, newx = xmat.test, s = "lambda.min",type="response")
  haz.pred[dta.split$cv.group==i,3]<-haz.pred.all[test.rows,3][test$tstart<test$time]
  surv.pred[dta_train$cv.group==i,3]<-ave(1-haz.pred.all[test.rows,3],test$pid,FUN=cumprod)[test$tstart==10]
  
  #---
  #method 4: random forest
  method.rf <- ranger(status ~ tstart+year1+year2+age+meno+size1+size2+grade+nodes+pgr+er+hormon+chemo,
                      data=train,oob.error = FALSE)
  haz.pred.all[test.rows,4] <- predict(method.rf, data = test)$predictions
  haz.pred[dta.split$cv.group==i,4]<-haz.pred.all[test.rows,4][test$tstart<test$time]
  surv.pred[dta_train$cv.group==i,4]<-ave(1-haz.pred.all[test.rows,4],test$pid,FUN=cumprod)[test$tstart==10]
  
  #-------------------------------------------------------
  #Calculate expected loss for each method
  #Using loss functions: L2 version, loglik version, IPCW version
  #-------------------------------------------------------
  
  #estimate censoring survival distribution (assuming independent censoring), which is needed as a weight in the ipcw loss function
  cens.km<-survfit(Surv(time,1-status)~1,data=test.1row)
  cens.wt<-1/summary(cens.km,times = 10)$surv
  
  for(j in 1:n.methods){
    #L2 loss
    loss.L2[i,j]<-sum((test.split$status-haz.pred[dta.split$cv.group==i,j])^2)
    
    #loglik loss
    loss.loglik[i,j]<-sum((log(haz.pred[dta.split$cv.group==i,j])^test.split$status)*
                        (log(1-haz.pred[dta.split$cv.group==i,j])^(1-test.split$status)))
    
    #ipcw loss
    loss.ipcw[i,j]<-sum(cens.wt*I(test.1row$status==1|(test.1row$status==0 & test.1row$time>10))*
                          (I(test.1row$time>10)-surv.pred[dta_train$cv.group==i,j])^2)
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
comb.mod.L2<-nnls(b=dta.split$status,A=haz.pred)

#scale the coefficients to sum to 1
weights.L2<-comb.mod.L2$x/sum(comb.mod.L2$x)

#---
#ensemble based on loglik loss function

#obtain logit transformations of the estimated hazards
g.haz.loglik<-log(haz.pred/(1-haz.pred))

#fit combining model, without an intercept
comb.mod.loglik<-glm(dta.split$status~-1+g.haz.loglik,family="binomial")

#no need to scale the coefficients for this loss function
weights.loglik<-comb.mod.loglik$coefficients

#---
#ensemble based on ipcw loss function, using non-negative weighted least squares
#This is achieved by multiplying the outcome and explnatory variables by the square root of the weights. 
#https://stackoverflow.com/questions/47888996/weighted-nonnegative-least-squares-in-r
cens.wt.dta<-data.frame(cens.wt=cens.wt,cv.group=1:10)
dta_train2<-merge(dta_train,cens.wt.dta)
dta_train2<-dta_train2[order(dta_train2$pid),]
dta_train2$elig<-I(dta_train2$status==1|(dta_train2$status==0 & dta_train2$time>10))
dta_train2$status<-ifelse(dta_train2$status==1 & dta_train2$time<=10,1,0)
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
dta.haz.test<-merge(dta_test,times)
dta.haz.test<-dta.haz.test[order(dta.haz.test$pid,dta.haz.test$tstart),]
dta.haz.test$tstart2<-dta.haz.test$tstart^2
dta.haz.test$tstart3<-dta.haz.test$tstart^3

#---
#storing predicted hazards from each method
n.methods=4
haz.pred.all<-matrix(nrow=dim(dta.haz.test)[1],ncol=n.methods)

#storing predicted survival probabilities from each method
n.test<-dim(dta_test)[1]
surv.pred<-matrix(nrow=n.test,ncol=n.methods)

#---
#method 1: glm with no covariates
method.glm0<-glm(status~1,data=dta.split,family = "binomial")
haz.pred.all[,1]<-predict(method.glm0,newdata=dta.haz.test,type="response")
surv.pred[,1]<-ave(1-haz.pred.all[,1],dta.haz.test$pid,FUN=cumprod)[dta.haz.test$tstart==10]

#---
#method 2: glm including time plus covariates
method.glm1<-glm(status~tstart+tstart2+tstart3+year1+year2+age+meno+size1+size2+grade+nodes+pgr+er+hormon+chemo,
                 data=dta.split,family = "binomial")
haz.pred.all[,2]<-predict(method.glm0,newdata=dta.haz.test,type="response")
surv.pred[,2]<-ave(1-haz.pred.all[,2],dta.haz.test$pid,FUN=cumprod)[dta.haz.test$tstart==10]

#---
#method 3: lasso
xmat.train<-model.matrix(method.glm1)[,-1]
dta.haz.test$status=99#just a trick so that we can use model.matrix to get xmat.test
xmat.test<-model.matrix(method.glm1,data=dta.haz.test)[,-1]
method.lasso <- cv.glmnet(x=xmat.train,y=dta.split$status,family = "binomial",alpha=1)
haz.pred.all[,3]<-predict(method.lasso, newx = xmat.test, s = "lambda.min",type="response")
surv.pred[,3]<-ave(1-haz.pred.all[,3],dta.haz.test$pid,FUN=cumprod)[dta.haz.test$tstart==10]

#---
#method 4: random forest
method.rf <- ranger(status ~ tstart+year1+year2+age+meno+size1+size2+grade+nodes+pgr+er+hormon+chemo,
                    data=dta.split,oob.error = FALSE)
haz.pred.all[,4] <- predict(method.rf, data = dta.haz.test)$predictions
surv.pred[,4]<-ave(1-haz.pred.all[,4],dta.haz.test$pid,FUN=cumprod)[dta.haz.test$tstart==10]

#---
#now combine using the ensemble weights
#---

#---
#ensemble SL based on L2 loss function
haz.pred.SL.L2<-weights.L2[1]*haz.pred.all[,1]+weights.L2[2]*haz.pred.all[,2]+
  weights.L2[3]*haz.pred.all[,3]+weights.L2[4]*haz.pred.all[,4]

surv.pred.SL.L2<-ave(1-haz.pred.SL.L2,dta.haz.test$pid,FUN=cumprod)[dta.haz.test$tstart==10]

#---
#ensemble SL based on loglik loss function
logit<-function(x){log(x/(1-x))}
expit<-function(x){exp(x)/(1+exp(x))}

g.haz.pred.SL.loglik<-
  weights.loglik[1]*logit(haz.pred.all[,1])+
  weights.loglik[2]*logit(haz.pred.all[,2])+
  weights.loglik[3]*logit(haz.pred.all[,3])+
  weights.loglik[4]*logit(haz.pred.all[,4])

haz.pred.SL.loglik<-expit(g.haz.pred.SL.loglik)

surv.pred.SL.loglik<-ave(1-haz.pred.SL.loglik,dta.haz.test$pid,FUN=cumprod)[dta.haz.test$tstart==10]

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
risk.pred<-1-surv.pred.SL.loglik

#---------------------------------
#Calibration plots
#---------------------------------

#---
#calibration plot - using riskRegression
#Not sure this is possible

#---
#calibration plot - using our function
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
#Brier score and IPA - using riskRegression
#Not sure this is possible

#---
#Brier score and IPA - using our function

Brier(time=dta_test$time, status=dta_test$status, risk=risk.pred, 
      seq.time=10, weights=dta_test$cens.wt)

ipa(time=dta_test$time, status=dta_test$status, risk=risk.pred, 
    seq.time=10, weights=dta_test$cens.wt)

#---
#Integrated Brier score - using riskRegression
#?

#---
#Integrated Brier score - using pec
#?

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
#AUC - using Score
#?

#---
#AUC - using timeROC
timeROC( T = dta_test$time,delta = dta_test$status,marker = risk.pred,
         cause = 1,weighting = "marginal",times = 10,iid = FALSE)

#---
#C/D AUCt - using our function
#it is very slightly lower using our function
wCD_AUCt(time=dta_test$time,status=dta_test$status, risk=risk.pred, seq.time = 10, weightmatrix = wt_matrix_eventsonly)




