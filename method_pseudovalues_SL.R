#################################################################################
#################################################################################
#Super learner for time-to-event outcomes tutorial
#
#Method: pseudo values super learner
#Sachs et al 2019
#
#Implemented using superLearner package, 
#and making use of Michael Sach's code from https://github.com/sachsmc/pseupersims
#################################################################################
#################################################################################

source("packages.R")
source("rotterdam_setup.R")
source("functions.R")
source("functions_pseudovalues.R")
source("censoring_weights_KM.R")#Kaplan-Meier estimates of censoring weights
# source("censoring_weights_discretetime_SL.R")#uncomment to use SL estimates of censoring weights instead

#extra packages
library(pseudo)
library(kernlab)

#---------------------------------
#---------------------------------
#calculate pseudo values at times 1:10
#this uses the add_pseudo_obs which is modified form the code of Sachs.
#Create a new version of the data which contains pseudo values at times 1:10 and a time variable (pvtime), 
#---------------------------------
#---------------------------------

add_pseudo_obs <- function(data, tme = 1:10) {
  
  psuo <- pseudoci(data$time, event = data$status, tmax = tme)
  data <- do.call(rbind, lapply(1:length(tme), function(x) cbind(data, cause1.pseudo = psuo$pseudo$cause1[, x], pvtime = tme[x])))
  data
}

#dta_train_pv is stacked across times 1:10
dta_train_pv<-add_pseudo_obs(data=dta_train,tme=1:10)

#---------------------------------
#---------------------------------
#Superlearner using pseudovalues at 10 time horizons
#---------------------------------
#---------------------------------

#add time of interest to the test data
dta_test$pvtime<-10

sl.full <- SuperLearner(Y = dta_train_pv$pv, 
                        X = dta_train_pv[,c("pvtime","year1","year2","age","meno","size1","size2","grade","nodes","pgr",
                                            "er","hormon","chemo")],
                        newX=dta_test[,c("pvtime","year1","year2","age","meno","size1","size2","grade","nodes","pgr",
                                         "er","hormon","chemo")],
                        SL.library = SL.library, 
                        id = dta_train_pv$id,
                        verbose = TRUE, 
                        method = "method.pseudoAUC", 
                        control = list(timedex = dta_train_pv$pvtime == 10),
                        cvControl = list(V=5))


#using 1 time horizon only

# dta_train_pv10<-dta_train_pv[dta_train_pv$pvtime==10,]
# 
# 
# sl.full10 <- SuperLearner(Y = dta_train_pv10$pv, 
#                           X = dta_train_pv10[,c("year1","year2","age","meno","size1","size2","grade","nodes","pgr",
#                                                 "er","hormon","chemo")],
#                           newX=dta_test[,c("year1","year2","age","meno","size1","size2","grade","nodes","pgr",
#                                            "er","hormon","chemo")],
#                           SL.library = list("SL.mean", "SL.glm", "SL.glmnet", "SL.ranger"), 
#                           id = dta_train_pv10$id,
#                           verbose = TRUE, 
#                           method = "method.pseudoAUC", control = list(timedex = dta_train_pv10$pvtime == 10))
# 
# sl.full10$SL.predict

#---------------------------------
#---------------------------------
#Superlearner using pseudovalues at 10 time horizons
#---------------------------------
#---------------------------------

#standardise continuous variables

mean.age<-mean(dta_train$age)
sd.age<-sd(dta_train$age)

mean.nodes<-mean(dta_train$nodes)
sd.nodes<-sd(dta_train$nodes)

mean.pgr<-mean(dta_train$pgr)
sd.pgr<-sd(dta_train$pgr)

mean.er<-mean(dta_train$er)
sd.er<-sd(dta_train$er)

dta_train_pv$age_std<-(dta_train_pv$age-mean.age)/sd.age
dta_train_pv$nodes_std<-(dta_train_pv$nodes-mean.nodes)/sd.nodes
dta_train_pv$pgr_std<-(dta_train_pv$pgr-mean.pgr)/sd.pgr
dta_train_pv$er_std<-(dta_train_pv$er-mean.er)/sd.er

dta_test$age_std<-(dta_test$age-mean.age)/sd.age
dta_test$nodes_std<-(dta_test$nodes-mean.nodes)/sd.nodes
dta_test$pgr_std<-(dta_test$pgr-mean.pgr)/sd.pgr
dta_test$er_std<-(dta_test$er-mean.er)/sd.er

#add time of interest to the test data
dta_test$pvtime<-10

#set up learners 
tune = list(ntrees = c(200),
            max_depth = 2,
            shrinkage = c(0.01, 0.1, .2))
learners = create.Learner("SL.xgboost", tune = tune, detailed_names = T, name_prefix = "xgb")

length(learners$names)
SL.library <-  c("SL.glm",  "SL.gam",  "SL.ksvm", "SL.ranger",
                 "SL.rpart", "SL.glmnet","SL.polymars", learners$names)

sl.full.std <- SuperLearner(Y = dta_train_pv$pv, 
                            X = dta_train_pv[,c("pvtime","year1","year2","age_std","meno","size1","size2","grade","nodes_std","pgr_std",
                                                "er_std","hormon","chemo")],
                            newX=dta_test[,c("pvtime","year1","year2","age_std","meno","size1","size2","grade","nodes_std","pgr_std",
                                             "er_std","hormon","chemo")],
                            SL.library = SL.library, 
                            id = dta_train_pv$id,
                            verbose = TRUE, 
                            method = "method.pseudoAUC", control = list(timedex = dta_train_pv$pvtime == 10))

list("SL.mean", "SL.glm", "SL.glmnet", "SL.ranger")


#---------------------------------
#---------------------------------
#obtain predictions
#---------------------------------
#---------------------------------

risk.pred<-sl.full.std$SL.predict

#truncate to range 0 to 1
risk.pred<-ifelse(risk.pred<0,0,risk.pred)
risk.pred<-ifelse(risk.pred>1,1,risk.pred)

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
#it is very slightly lower using our function
wCD_AUCt(time=dta_test$time,status=dta_test$status, risk=risk.pred, seq.time = 10, weightmatrix = wt_matrix_eventsonly)

#---
#AUC - using timeROC
timeROC( T = dta_test$time,delta = dta_test$status,marker = risk.pred,
         cause = 1,weighting = "marginal",times = 10,iid = FALSE)






