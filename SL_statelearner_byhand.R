#################################################################################
#################################################################################
#Super learner for time-to-event outcomes tutorial
#
#Method: continuous-time superLearner
#Munch and Gerds 2024
#
#Implemented by hand 
#################################################################################
#################################################################################

#divide data into 5 folds
n<-dim(dta_train)[1]
set.seed(1430)
dta_train$cv.group<-cut(runif(n,0,5),breaks=seq(0,5,1),include.lowest = T,labels = F)

#-------------------------------------------------------
#candidate learners: Nelson-Aalen, Cox, Cox-Lasso, random survival forest
#fit model excluding each fold in turn and then obtain predictions in the excluded fold
#-------------------------------------------------------

IBS<-matrix(nrow=5,ncol=16)

#cross-validation loop
for(i in 1:5){
  print(i)
  
  train<-dta_train[!dta_train$cv.group==i,]
  test<-dta_train[dta_train$cv.group==i,]
  n.test<-dim(test)[1]
  
  #unique event/censoring times in training data
  train.times<-sort(unique(train$time))
  n.tt<-length(train.times)
  
  #---
  #Nelson-Aalen
  
  #for event
  method.na.event<-coxph(Surv(time,status2==1)~1,data=train,x=TRUE)
  chaz.na.event<-predictCox(method.na.event,newdata=test,times=train.times,type="cumhazard")$cumhazard
  haz.na.event<-t(sapply(1:n.test,FUN=function(x){diff(c(0,chaz.na.event[x,]))}))
  
  #censoring
  method.na.cens<-coxph(Surv(time,status2==-1)~1,data=train,x=TRUE)
  chaz.na.cens<-predictCox(method.na.cens,newdata=test,times=train.times,type="cumhazard")$cumhazard
  haz.na.cens<-t(sapply(1:n.test,FUN=function(x){diff(c(0,chaz.na.cens[x,]))}))
  
  #---
  #Cox model
  
  #for event
  method.coxph.event<-coxph(Surv(time,status2==1)~year1+year2+age+meno+size+grade+nodes+pgr+er+hormon+chemo,data=train,x=TRUE)
  chaz.coxph.event<-predictCox(method.coxph.event,newdata=test,times=train.times,type="cumhazard")$cumhazard
  haz.coxph.event<-t(sapply(1:n.test,FUN=function(x){diff(c(0,chaz.coxph.event[x,]))}))
    
  #for censoring
  method.coxph.cens<-coxph(Surv(time,status2==-1)~year1+year2+age+meno+size+grade+nodes+pgr+er+hormon+chemo,data=train,x=TRUE)
  chaz.coxph.cens<-predictCox(method.coxph.cens,newdata=test,times=train.times,type="cumhazard")$cumhazard
  haz.coxph.cens<-t(sapply(1:n.test,FUN=function(x){diff(c(0,chaz.coxph.cens[x,]))}))
  
  #---
  #Cox with lasso
  
  #for event
  xmat.train<-as.matrix(train[,c("year1","year2","age","meno","size1","size2","grade","nodes","pgr","er","hormon","chemo")])
  xmat.test<-as.matrix(test[,c("year1","year2","age","meno","size1","size2","grade","nodes","pgr","er","hormon","chemo")])
  method.coxlasso.event <- cv.glmnet(x=xmat.train,y=Surv(train$time,train$status2==1),family = "cox",alpha=1,type.measure = "C")
  chaz.coxlasso.event<-t(survfit(method.coxlasso.event, s = "lambda.min", x = xmat.train, y = Surv(train$time,train$status2==1), newx = xmat.test)$cumhaz)
  haz.coxlasso.event<-t(sapply(1:n.test,FUN=function(x){diff(c(0,chaz.coxlasso.event[x,]))}))

  #for censoring
  method.coxlasso.cens <- cv.glmnet(x=xmat.train,y=Surv(train$time,train$status2==-1),family = "cox",alpha=1,type.measure = "C")
  chaz.coxlasso.cens<-t(survfit(method.coxlasso.cens, s = "lambda.min", x = xmat.train, y = Surv(train$time,train$status2==-1), newx = xmat.test)$cumhaz)
  haz.coxlasso.cens<-t(sapply(1:n.test,FUN=function(x){diff(c(0,chaz.coxlasso.cens[x,]))}))
  
  #---
  #survival random forest
  
  #for event
  method.rfsrc.event <- rfsrc(Surv(time, status) ~ year1+year2+age+meno+size1+size2+grade+nodes+pgr+er+hormon+chemo, data = train,ntime=0)
  pred.rfsrc.event<-predict(method.rfsrc.event, newdata=test, importance='none')
  chaz.rfsrc.event<-matrix(nrow=n.test,ncol=length(train.times))
  for(id in 1:n.test){
    step.rfsrc.event<-stepfun(pred.rfsrc.event$time.interest,c(0,pred.rfsrc.event$chf[id,]))
    chaz.rfsrc.event[id,]<-step.rfsrc.event(train.times)
  }
  haz.rfsrc.event<-t(sapply(1:n.test,FUN=function(x){diff(c(0,chaz.rfsrc.event[x,]))}))

  #for censoring
  train$status.cens<-ifelse(train$status2==-1,1,0)
  test$status.cens<-ifelse(test$status2==-1,1,0)
  method.rfsrc.cens <- rfsrc(Surv(time, status.cens) ~ year1+year2+age+meno+size1+size2+grade+nodes+pgr+er+hormon+chemo, data = train,ntime=0)
  pred.rfsrc.cens<-predict(method.rfsrc.cens, newdata=test, importance='none')
  chaz.rfsrc.cens<-matrix(nrow=n.test,ncol=length(train.times))
  for(id in 1:n.test){
    step.rfsrc.cens<-stepfun(pred.rfsrc.cens$time.interest,c(0,pred.rfsrc.cens$chf[id,]))
    chaz.rfsrc.cens[id,]<-step.rfsrc.cens(train.times)
  }
  haz.rfsrc.cens<-t(sapply(1:n.test,FUN=function(x){diff(c(0,chaz.rfsrc.cens[x,]))}))
  
  #---
  #F functions for each combination (equation (6) in Munch and Gerds)

  #arrays of cumulative hazards from each model
  chaz.all.event<-array(cbind(chaz.na.event,chaz.coxph.event,chaz.coxlasso.event,chaz.rfsrc.event),
                        dim=c(n.test,n.tt,4))
  chaz.all.cens<-array(cbind(chaz.na.cens,chaz.coxph.cens,chaz.coxlasso.cens,chaz.rfsrc.cens),
                        dim=c(n.test,n.tt,4))
  
  #arrays of hazards from each model
  haz.all.event<-array(cbind(haz.na.event,haz.coxph.event,haz.coxlasso.event,haz.rfsrc.event),
                        dim=c(n.test,n.tt,4))
  haz.all.cens<-array(cbind(haz.na.cens,haz.coxph.cens,haz.coxlasso.cens,haz.rfsrc.cens),
                       dim=c(n.test,n.tt,4))
  
  
  #overall event-free survival
  #This contains a matrix of cumulative incidences for each person, at each time,
  #and for each of the 16 combinations of models for event and censoring 
  F.none<-array(dim=c(n.test,n.tt,16))
  m<-1
  for(j in 1:4){
    for(k in 1:4){
      F.none[,,m]<-exp(-chaz.all.event[,,j]-chaz.all.cens[,,k]) 
      m<-m+1
    }
  }

  #lags of F.none, used in the next step
  
  chaz.all.event.lags<-array(cbind(cbind(rep(0,n.test),chaz.na.event[,-n.tt]),
                                   cbind(rep(0,n.test),chaz.coxph.event[,-n.tt]),
                                   cbind(rep(0,n.test),chaz.coxlasso.event[,-n.tt]),
                                   cbind(rep(0,n.test),chaz.rfsrc.event[,-n.tt])),
                        dim=c(n.test,length(train.times),4))
  chaz.all.cens.lags<-array(cbind(cbind(rep(0,n.test),chaz.na.cens[,-n.tt]),
                                   cbind(rep(0,n.test),chaz.coxph.cens[,-n.tt]),
                                   cbind(rep(0,n.test),chaz.coxlasso.cens[,-n.tt]),
                                   cbind(rep(0,n.test),chaz.rfsrc.cens[,-n.tt])),
                             dim=c(n.test,n.tt,4))
  lagF.none<-array(dim=c(n.test,n.tt,16))
  m<-1
  for(j in 1:4){
    for(k in 1:4){
      lagF.none[,,m]<-exp(-chaz.all.event.lags[,,j]-chaz.all.cens.lags[,,k]) 
      m<-m+1
    }
  }
  
  #1: event-na, cens-na
  #2: event-na, cens-cox
  #3: event-na, cens-coxlasso
  #4: event-na, cens-rfsrc
  #5: event-cox, cens-na
  #6: event-cox, cens-cox
  #7: event-cox, cens-coxlasso
  #8: event-cox, cens-rfsrc
  #9: event-coxlasso, cens-na
  #10: event-coxlasso, cens-cox
  #11: event-coxlasso, cens-coxlasso
  #12: event-coxlasso, cens-rfsrc
  #13: event-rfsrc, cens-na
  #14: event-rfsrc, cens-cox
  #15: event-rfsrc, cens-coxlasso
  #16: event-rfsrc, cens-rfsrc
  
  
  #cumulative incidence for event
  F.event<-array(dim=c(n.test,n.tt,16))
  for(j in 1:16){
    # print(j)
    event.model<-ceiling(j/4)
    F.event[,,j]<-t(sapply(1:n.test,
                    FUN=function(x){cumsum(lagF.none[x,,j]*haz.all.event[x,,event.model])}))
  }
  
  # above uses combination of models (1-16) -- event model (1-4)
  # 1,2,3,4 -- 1,1,1,1
  # 5,6,7,8 -- 2,2,2,2
  # 9,10,11,12 -- 3,3,3,3
  # 13,14,15,16 -- 4,4,4,4
  
  #cumulative incidence for censoring
  F.cens<-array(dim=c(n.test,n.tt,16))
  for(j in 1:16){
    # print(j)
    cens.model<-(j-1)%%4+1
    F.cens[,,j]<-t(sapply(1:n.test,
                           FUN=function(x){cumsum(lagF.none[x,,j]*haz.all.cens[x,,cens.model])}))
  }
  
  # above uses combination of models (1-16) -- cens model (1-4)
  # 1,2,3,4 -- 1,2,3,4
  # 5,6,7,8 -- 1,2,3,4
  # 9,10,11,12 -- 1,2,3,4
  # 13,14,15,16 -- 1,2,3,4
  
  #---
  #indicators of which state a person is in at each time
  #denoted eta(t) in Munch and Gerds
  
  eta.event<-t(sapply(1:n.test,FUN=function(x){I(test$time[x]<=train.times & test$status[x]==1)}))
  eta.cens<-t(sapply(1:n.test,FUN=function(x){I(test$time[x]<=train.times & test$status[x]==-1)}))
  eta.none<-t(sapply(1:n.test,FUN=function(x){I(test$time[x]<=train.times & test$status[x]==0)}))
  
  #---
  #Brier scores
  
  #Brier score over all three outcomes
  B<-array(dim=c(n.test,n.tt,16))
  for(j in 1:16){
    B[,,j]<-(F.event[,,j]-eta.event)^2+(F.cens[,,j]-eta.cens)^2+(F.none[,,j]-eta.none)^2
    #check Munch code - it might be as below
    # B[,,j]<-(F.event[,,j]-eta.event)^2+(F.cens[,,j]-eta.cens)^2
  }
  
  #integrated Brier score up to the last time point (10 years)
  Bsum<-matrix(nrow=n.test,ncol=16)
  for(j in 1:16){
    Bsum[,j]<-rowSums(B[,,j]*diff(c(0,train.times)))
  }

  #mean integrated Brier score
  IBS[i,]<-colMeans(Bsum)
 
}

IBS.cv<-matrix(colMeans(IBS),nrow=4,ncol=4,byrow=T)
rownames(IBS.cv)<-paste0("event-",c("na","cox","coxlasso","rfsrc"))
colnames(IBS.cv)<-paste0("cens-",c("na","cox","coxlasso","rfsrc"))

IBS.cv

#The lowest Brier score is for the combination (event-rfsrc, cens-na)
