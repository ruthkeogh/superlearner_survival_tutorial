#################################################################################
#################################################################################
#Super learner for time-to-event outcomes tutorial
#
#Method: continuous-time superLearner
#Westling et al. 2023
#
#As in method_continuoustime_SL_Westling.R but with a larger set of candidate learners,
#using the full set of learners built into survSuperLearner.
#See https://github.com/tedwestling/survSuperLearner
#################################################################################
#################################################################################

source("packages.R")
source("rotterdam_setup.R")
source("functions.R")
source("censoring_weights.R")

#---------------------------------
#---------------------------------
#superlearner
#---------------------------------
#---------------------------------

#Candidate learners
event.SL.library <- list("survSL.km", "survSL.coxph", "survSL.expreg", 
                         "survSL.weibreg", "survSL.loglogreg", "survSL.gam"
                         , "survSL.rfsrc", "survSL.pchreg")

cens.SL.library<-event.SL.library

covs_train<-dta_train[,c("year1","year2","age","meno","size1","size2","grade","nodes","pgr","er","hormon","chemo")]
covs_test<-dta_test[,c("year1","year2","age","meno","size1","size2","grade","nodes","pgr","er","hormon","chemo")]
set.seed(1)
sl <- survSuperLearner(time = dta_train$time, event = dta_train$status, X = covs_train, newX = covs_test, new.times = 10, 
                       event.SL.library = event.SL.library, cens.SL.library = cens.SL.library, verbose = TRUE,
                       cvControl=list(V = 5),
                       control=list(saveFitLibrary = TRUE))

#coefficients in the ensemble
sl$event.coef
sl$cens.coef

#---------------------------------
#---------------------------------
#obtain predictions
#---------------------------------
#---------------------------------

#risk predictions in the test data
#This appears to give predictions up to the last time point (here time 10)
risk.pred<-1-sl$event.SL.predict

#---------------------------------
#---------------------------------
#Calibration plots
#---------------------------------
#---------------------------------

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
#---------------------------------
#Brier score, IPA, and integrated Brier score
#---------------------------------
#---------------------------------

#---
#Brier score and IPA - using our function

Brier(time=dta_test$time, status=dta_test$status, risk=risk.pred, 
      seq.time=10, weights=dta_test$cens.wt)

ipa(time=dta_test$time, status=dta_test$status, risk=risk.pred, 
    seq.time=10, weights=dta_test$cens.wt)

#---------------------------------
#---------------------------------
#C-index, AUC and AUCt
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
#AUC - using timeROC
timeROC( T = dta_test$time,delta = dta_test$status,marker = risk.pred,
         cause = 1,weighting = "marginal",times = 10,iid = FALSE)

#---
#C/D AUCt - using our function
max.event.time<-max(dta_test$time[dta_test$status==1])
wCD_AUCt(time=dta_test$time,status=dta_test$status, risk=risk.pred, seq.time =max.event.time, weightmatrix = wt_matrix_eventsonly)






