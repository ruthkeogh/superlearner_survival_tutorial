#################################################################################
#################################################################################
#Super learner for time-to-event outcomes tutorial
#
#Method: continuous-time superLearner
#Munch and Gerds 2024/25
#################################################################################
#################################################################################

source("packages.R")
source("rotterdam_setup.R")
source("functions.R")
source("censoring_weights_KM.R")#Kaplan-Meier estimates of censoring weights
# source("censoring_weights_continuoustime_SL_MunchGerds.R")#uncomment to use SL estimates of censoring weights instead

#---------------------------------
#---------------------------------
#superlearner
#---------------------------------
#---------------------------------

#make data into data.table, as required by the jossl code
dta_train2<-dta_train
# dta_train2$status<-dta_train2$status2
dta_train2<-as.data.table(dta_train2)

#list of learners to be used
learners <- list(
  N_Aa = list(model = "cox", x_form = ~1),
  cox = list(model = "cox", x_form = ~ year1+year2+age+meno+size1+size2+grade+nodes+pgr+er+hormon+chemo),
  cox_penalty = list(model = "GLMnet", x_form = ~year1+year2+age+meno+size1+size2+grade+nodes+pgr+er+hormon+chemo),
  rf = list(model = "rfsrc", x_form = ~year1+year2+age+meno+size1+size2+grade+nodes+pgr+er+hormon+chemo,n.time=0)
)

#apply jossl
set.seed(1)
sl = jossl(learners = list(cause1 = learners,
                                  censor = learners),
                  data = dta_train2,
                  time = 10,
                  integrate = TRUE,
                  split.method = "cv5",
                  B=1,
                  time_grid_length=100,
                  vars = NULL,
                  collapse = TRUE)

sl$fitted_winners 
sl$cv_fit
#For cause 1 the best fitting model is rfsrc

#---------------------------------
#---------------------------------
#obtain predictions
#---------------------------------
#---------------------------------

#refit the best fitting model on the full training data
set.seed(1)
rfsrc.fit <- rfsrc(Surv(time, status) ~ year1+year2+age+meno+size1+size2+grade+nodes+pgr+er+hormon+chemo, data = dta_train,ntime=0)

#obtain predictions in the test data
rfsrc.pred<-predict(rfsrc.fit, newdata=dta_test, importance='none')
rfsrc.surv<-rfsrc.pred$survival

approx.surv <- c(t(sapply(1:nrow(rfsrc.surv), function(i) {
  stats::approx(c(0,rfsrc.pred$time.interest), c(1,rfsrc.surv[i,]), method='constant', xout = 10, rule = 2)$y
})))

risk.pred<-1-approx.surv

#---------------------------------
#---------------------------------
#Calibration plots
#---------------------------------
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
#---------------------------------
#Brier score and Scaled Brier (IPA)
#---------------------------------
#---------------------------------

#---
#Brier score and IPA - using riskRegression
Score(list(rf=rfsrc.fit),Surv(time,status)~1,data=dta_test,
      metrics="brier",times=10)

#---
#Brier score and IPA - using our function
#gives same results as above

Brier(time=dta_test$time, status=dta_test$status, risk=risk.pred, 
      seq.time=10, weights=dta_test$cens.wt)

ipa(time=dta_test$time, status=dta_test$status, risk=risk.pred, 
    seq.time=10, weights=dta_test$cens.wt)

#---------------------------------
#---------------------------------
#C-index and AUC
#---------------------------------
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
Score(list(rf=rfsrc.fit),Surv(time,status)~1,data=dta_test,
      metrics="auc",times=10)

#---
#C/D AUCt - using our function
max.event.time<-max(dta_test$time[dta_test$status==1])
wCD_AUCt(time=dta_test$time,status=dta_test$status, risk=risk.pred, seq.time =max.event.time, weightmatrix = wt_matrix_eventsonly)


#---
#AUC - using timeROC
timeROC( T = dta_test$time,delta = dta_test$status,marker = risk.pred,
         cause = 1,weighting = "marginal",times = 10,iid = FALSE)
