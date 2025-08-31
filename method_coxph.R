#################################################################################
#################################################################################
#Super learner for time-to-event outcomes tutorial
#
#Method: standard Cox regression.
#################################################################################
#################################################################################

source("packages.R")
source("rotterdam_setup.R")
source("functions.R")
source("censoring_weights.R")

#---------------------------------
#---------------------------------
#fit Cox model
#---------------------------------
#---------------------------------

cox.mod<-coxph(Surv(time,status)~year1+year2+age+meno+size1+size2+grade+nodes+pgr+er+hormon+chemo,data=dta_train,x=TRUE)

#---------------------------------
#---------------------------------
#obtain predictions
#---------------------------------
#---------------------------------

cox.pred<-predictRisk(cox.mod,newdata=dta_test,times=10)

#---------------------------------
#---------------------------------
#Calibration plots
#---------------------------------
#---------------------------------

#---
#calibration plot - using riskRegression

xs=Score(list(Cox=cox.mod),Surv(time,status)~1,data=dta_test,
         plots="cal",times=10,metrics=NULL)
plotCalibration(x=xs,models="Cox",method="quantile",q=10,cens.method = "local")

#---
#calibration plot - using our function
#This gives same results to those obtained above using riskRegression

cut_points=c(0,quantile(cox.pred,probs=seq(0.1,0.9,0.1)),1)
risk_group=cut(cox.pred,breaks=cut_points,include.lowest = T,labels = F)
calib_risk_group<-sapply(FUN=function(x){mean(cox.pred[risk_group==x])},1:10)

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
#---------------------------------
#Brier score, IPA, and integrated Brier score
#---------------------------------
#---------------------------------

#---
#Brier score and IPA - using riskRegression
#IPA=1-brier(Cox model)/brier(null model)
brier<-Score(list(Cox=cox.mod),Surv(time,status)~1,data=dta_test,
         metrics="brier",times=10)

#---
#Brier score and IPA - using our function
#gives same results as above

Brier(time=dta_test$time, status=dta_test$status, risk=cox.pred, 
      seq.time=10, weights=dta_test$cens.wt)

ipa(time=dta_test$time, status=dta_test$status, risk=cox.pred, 
      seq.time=10, weights=dta_test$cens.wt)

#---
#Integrated Brier score - using riskRegression
Score(list(Cox=cox.mod),Surv(time,status)~1,data=dta_test,
               metrics="brier",times=seq(0,10,0.1),summary="ibs")

#---
#Integrated Brier score - using pec
#This also gives the integrated Brier score
#Results are similar to those obtained above using Score, provided we use a fine time grid for 'times' in Score
pec(list(Cox=cox.mod),data=dta_test,formula=Surv(time,status)~1)

#---------------------------------
#---------------------------------
#C-index, AUC and AUCt
#---------------------------------
#---------------------------------

#---
#C-index - using our function
c_index_ties(time=dta_test$time,status=dta_test$status, risk=cox.pred, tau=10, weightmatrix = wt_matrix_eventsonly)

#c-index - using concordance
#same result (to 3 decimal places)
concordance(Surv(dta_test$time, dta_test$status) ~ cox.pred,
            newdata=dta_test,
            reverse = TRUE,
            timewt = "n/G2")$concordance

#---
#AUC - using riskRegression
Score(list(Cox=cox.mod),Surv(time,status)~1,data=dta_test,
               metrics="auc",times=10)

#---
#AUC - using timeROC
#same result
timeROC( T = dta_test$time,delta = dta_test$status,marker = cox.pred,
  cause = 1,weighting = "marginal",times = 10,iid = FALSE)

#---
#C/D AUCt - using our function
max.event.time<-max(dta_test$time[dta_test$status==1])
wCD_AUCt(time=dta_test$time,status=dta_test$status, risk=risk.pred, seq.time =max.event.time, weightmatrix = wt_matrix_eventsonly)






