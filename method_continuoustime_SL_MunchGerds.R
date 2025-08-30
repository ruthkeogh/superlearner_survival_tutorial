#################################################################################
#################################################################################
#Super learner for time-to-event outcomes tutorial
#
#Method: continuous-time superLearner
#Munch and Gerds 2024
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

#make data into data.table, as required by the statelearner code
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

#apply statelearner

set.seed(1)
sl = statelearner(learners = list(cause1 = learners,
                                  censor = learners),
                  data = dta_train2,
                  time = 10,
                  integrate = TRUE,
                  split.method = "cv5",
                  B=1,
                  time_grid_length=100,
                  vars = NULL,
                  collapse = TRUE)

sl$fitted_winners #For cause 1 the best fitting model is rfsrc
sl$cv_fit

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

#---
#Note that the statelearner code suggests we should be able to use predict to get predictions based on the best fitting model, as shown below.
#But this gives somewhat different results to just refitting the survival random forest, with poor predictive performance
#It is possible that this is getting the cumulative incidence for the event being observed in a world in which censoring occurs, i.e. the F function 

# dta_test2<-dta_test
# dta_test2<-dta_test2[,c("year1","year2","age","meno","size1","size2","grade","nodes","pgr","er","hormon","chemo")]
# 
# #this gives risk predictions at 150 times
# risk.pred.alltimes<-predict(object = sl$fitted_winners$cause1, newdata = dta_test2, 
#                             onlySL = TRUE)
# 
# #predictions at time 10
# risk.pred2<-risk.pred.alltimes$cif[,150,1]

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
#Brier score, IPA, and integrated Brier score
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

#---
#Integrated Brier score - using riskRegression
Score(list(rf=rfsrc.fit),Surv(time,status)~1,data=dta_test,
      metrics="brier",times=seq(0,10,0.1),summary="ibs")

#---
#Integrated Brier score - using pec
#This also gives the integrates Brier score
#Results are similar to those obtained above using Score, provided we use a fine time grid for 'times' in Score
pec(list(rf=rfsrc.fit),data=dta_test,formula=Surv(time,status)~1)

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
#AUC - using Score
Score(list(rf=rfsrc.fit),Surv(time,status)~1,data=dta_test,
      metrics="auc",times=10)

#---
#AUC - using timeROC
timeROC( T = dta_test$time,delta = dta_test$status,marker = risk.pred,
         cause = 1,weighting = "marginal",times = 10,iid = FALSE)

#---
#C/D AUCt - using our function
#it is very slightly lower using our function
wCD_AUCt(time=dta_test$time,status=dta_test$status, risk=risk.pred, seq.time = 10, weightmatrix = wt_matrix_eventsonly)
