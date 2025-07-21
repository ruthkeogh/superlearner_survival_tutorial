#################################################################################
#################################################################################
#Super learner for time-to-event outcomes tutorial
#
#Method: Cox-lasso
#################################################################################
#################################################################################

#---------------------------------
#---------------------------------
#fit Cox-lasso model
#---------------------------------
#---------------------------------

cox.mod<-coxph(Surv(time,status)~year1+year2+age+meno+size1+size2+grade+nodes+pgr+er+hormon+chemo,data=dta_train,x=TRUE)
xmat.train<-as.matrix(dta_train[,c("year1","year2","age","meno","size1","size2","grade","nodes","pgr","er","hormon","chemo")])
xmat.test<-as.matrix(dta_test[,c("year1","year2","age","meno","size1","size2","grade","nodes","pgr","er","hormon","chemo")])
set.seed(1)
cox.lasso <- cv.glmnet(x=xmat.train,y=Surv(dta_train$time,dta_train$status),family = "cox",alpha=1,type.measure = "C")

#take a look at which covariates have non-zero coefficients
predict(cox.lasso, newx = xmat.test, s = "lambda.min",type="coefficients")
predict(cox.lasso, newx = xmat.test, s = "lambda.min",type="nonzero")

#---------------------------------
#---------------------------------
#obtain predictions
#---------------------------------
#---------------------------------

#see https://cran.r-project.org/web/packages/glmnet/vignettes/Coxnet.pdf
#survfit:coxnet gives survival probabilities at all times, and we can then obtain the estimates of risk up to time 10
cox.pred.alltimes <- survfit(cox.lasso, s = "lambda.min", x = xmat.train, y = Surv(dta_train$time,dta_train$status), newx = xmat.test)
cox.pred <- 1-cox.pred.alltimes$surv[cox.pred.alltimes$time==10,]

#---------------------------------
#---------------------------------
#Calibration plots
#---------------------------------
#---------------------------------

#---
#calibration plot - using our function
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
#Brier score and IPA - using our function

Brier(time=dta_test$time, status=dta_test$status, risk=cox.pred, 
      seq.time=10, weights=dta_test$cens.wt)

ipa(time=dta_test$time, status=dta_test$status, risk=cox.pred, 
    seq.time=10, weights=dta_test$cens.wt)

#---
#Integrated Brier score - using riskRegression
#Not sure how to get this to work: code is from standard coxph
# Score(list(Cox=cox.mod),Surv(time,status)~1,data=dta_test,
#       metrics="brier",times=seq(0,10,0.1),summary="ibs")

#---
#Integrated Brier score - using pec
#Not working
# risk.matrix<-1-t(cox.pred.alltimes$surv)
# 
# pec(risk.matrix,times=cox.pred.alltimes$time,data=dta_test,formula=Surv(time,status)~1,exact=F)

#---------------------------------
#---------------------------------
#C-index, AUC and AUCt
#---------------------------------
#---------------------------------

#---
#C-index - using our function
c_index_ties(time=dta_test$time,status=dta_test$status, risk=cox.pred, tau=10, weightmatrix = wt_matrix_eventsonly)

#c-index - using concordance
concordance(Surv(dta_test$time, dta_test$status) ~ cox.pred,
            newdata=dta_test,
            reverse = TRUE,
            timewt = "n/G2")$concordance

#---
#AUC - using Score
#Not sure how to get this to work: code is from standard coxph
# Score(list(Cox=cox.mod),Surv(time,status)~1,data=dta_test,
#       metrics="auc",times=10)

#---
#AUC - using timeROC
timeROC( T = dta_test$time,delta = dta_test$status,marker = cox.pred,
         cause = 1,weighting = "marginal",times = 10,iid = FALSE)

#---
#C/D AUCt - using our function
#it is very slightly lower using our function
wCD_AUCt(time=dta_test$time,status=dta_test$status, risk=cox.pred, seq.time = 10, weightmatrix = wt_matrix_eventsonly)





