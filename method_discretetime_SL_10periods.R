#################################################################################
#################################################################################
#Super learner for time-to-event outcomes tutorial
#
#Method: discrete-time superLearner, using 10 time periods and indicators for time period as covariates in the SL
#Polley and van der Laan 2011
#################################################################################
#################################################################################

source("packages.R")
source("rotterdam_setup.R")
source("functions.R")
source("censoring_weights_KM.R")#Kaplan-Meier estimates of censoring weights
# source("censoring_weights_discretetime_SL.R")#uncomment to use SL estimates of censoring weights instead

#---------------------------------
#---------------------------------
#superlearner
#---------------------------------
#---------------------------------

#time grid to be used for discrete time analysis
times<-seq(0,9,1)

#split training data into discrete-time format
dta_train.split<-survSplit(Surv(time,status)~.,data=dta_train,cut = times,start="tstart",end="tstop",event="event")

#in the test data we need one row for all times in the grid 'times'
times.dta<-data.frame(tstart=times)
dta_test.pred<-merge(dta_test,times.dta)
dta_test.pred<-dta_test.pred[order(dta_test.pred$pid,dta_test.pred$tstart),]

#generate time period indicator variables to be used as covariates in the SL

dta_train.split$tstart0<-I(dta_train.split$tstart==0)
dta_train.split$tstart1<-I(dta_train.split$tstart==1)
dta_train.split$tstart2<-I(dta_train.split$tstart==2)
dta_train.split$tstart3<-I(dta_train.split$tstart==3)
dta_train.split$tstart4<-I(dta_train.split$tstart==4)
dta_train.split$tstart5<-I(dta_train.split$tstart==5)
dta_train.split$tstart6<-I(dta_train.split$tstart==6)
dta_train.split$tstart7<-I(dta_train.split$tstart==7)
dta_train.split$tstart8<-I(dta_train.split$tstart==8)
dta_train.split$tstart9<-I(dta_train.split$tstart==9)

dta_test.pred$tstart0<-I(dta_test.pred$tstart==0)
dta_test.pred$tstart1<-I(dta_test.pred$tstart==1)
dta_test.pred$tstart2<-I(dta_test.pred$tstart==2)
dta_test.pred$tstart3<-I(dta_test.pred$tstart==3)
dta_test.pred$tstart4<-I(dta_test.pred$tstart==4)
dta_test.pred$tstart5<-I(dta_test.pred$tstart==5)
dta_test.pred$tstart6<-I(dta_test.pred$tstart==6)
dta_test.pred$tstart7<-I(dta_test.pred$tstart==7)
dta_test.pred$tstart8<-I(dta_test.pred$tstart==8)
dta_test.pred$tstart9<-I(dta_test.pred$tstart==9)

#apply SuperLearner using NNLS
#change to method = "method.NNloglik" for loss function based on loglik
set.seed(1)
sl<-SuperLearner(Y=dta_train.split$event,
                    X=dta_train.split[,c("tstart0","tstart1","tstart2","tstart3","tstart4","tstart5",
                                         "tstart6","tstart7","tstart8","tstart9",
                                         "year1","year2","age","meno","size1","size2","grade","nodes","pgr",
                                         "er","hormon","chemo")] ,
                    newX=dta_test.pred[,c("tstart0","tstart1","tstart2","tstart3","tstart4","tstart5",
                                          "tstart6","tstart7","tstart8","tstart9",
                                          "year1","year2","age","meno","size1","size2","grade","nodes","pgr",
                                          "er","hormon","chemo")],
                    family = binomial(), 
                    SL.library= list("SL.mean", "SL.glm", "SL.glmnet", "SL.ranger"),
                    method = "method.NNLS", id = dta_train.split$pid, verbose = TRUE,
                    cvControl = list(V=5))

#---------------------------------
#---------------------------------
#obtain predictions
#---------------------------------
#---------------------------------

pred.dat<-cbind(dta_test.pred[,c("pid","tstart")],pred=sl$SL.predict)

pred.dat<-pred.dat%>%group_by(pid)%>%mutate(surv=cumprod(1-pred))

risk.pred<-1-pred.dat$surv[pred.dat$tstart==9]

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
#Brier score and Scaled Brier (IPA)
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
#C/D AUCt - using our function
max.event.time<-max(dta_test$time[dta_test$status==1])
wCD_AUCt(time=dta_test$time,status=dta_test$status, risk=risk.pred, seq.time =max.event.time, weightmatrix = wt_matrix_eventsonly)

#---
#AUC - using timeROC
timeROC( T = dta_test$time,delta = dta_test$status,marker = risk.pred,
         cause = 1,weighting = "marginal",times = 10,iid = FALSE)