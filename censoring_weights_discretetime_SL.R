#################################################################################
#################################################################################
#Super learner for time-to-event outcomes tutorial
#
#The file obtains censoring weights in the test data (dta_test)
#, which are required for estimating predictive performance measures using the functions in functions.R, and in some functions used from packages. 
#
#In this file the weights are obtained using the discrete-time super learner.
#
#For some measures the weights are required only at the person's observed event/censoring time. 
#For some measures (C-index and C/D AUCt) the weights are required for each person at every observed event time (i.e. a matrix of weights for each person at each time).
#################################################################################
#################################################################################

#---------------------------------
#---------------------------------
# time-fixed censoring weights at observed event or censoring times
# obtained using discrete-time SL
# This gives dta_test$cens.wt
#---------------------------------
#---------------------------------

#time grid to be used for discrete time analysis
times<-seq(0,max(dta_test$time),0.1)

#split data into discrete-time format
dta_test.split<-survSplit(Surv(time,1-status)~.,data=dta_test,cut = times,start="tstart",end="tstop",event="cens")

#generate time^2 and time^3 variables to be used as covariates in the SL
dta_test.split$tstart2<-dta_test.split$tstart^2
dta_test.split$tstart3<-dta_test.split$tstart^3

#we do not want to include administrative censoring at just after 10 years, as this is a different censoring process
#there only use censorings up to 10 years
dta_test.split$cens<-ifelse(dta_test.split$tstop==3651/365,0,dta_test.split$cens)

#apply SuperLearner using NNLS
set.seed(1)
sl.cens<-SuperLearner(Y=dta_test.split$cens,X=dta_test.split[,c("tstart","tstart2","tstart3","year1","year2","age","meno","size1","size2","grade","nodes","pgr","er","hormon","chemo")] ,
                      newX=dta_test.split[,c("tstart","tstart2","tstart3","year1","year2","age","meno","size1","size2","grade","nodes","pgr","er","hormon","chemo")],
                      family = binomial(), SL.library= list("SL.mean", "SL.glm", "SL.glmnet", "SL.ranger"),
                      method = "method.NNLS", id = dta_test.split$pid, verbose = TRUE,
                      cvControl = list(V=5))

#obtain predicted discrete-time hazards
pred.dat<-cbind(dta_test.split[,c("pid","tstart")],pred=sl.cens$SL.predict)

#obtain censoring probabilities based on the discrete-time hazards
#For each person up to their observed event/cens time
pred.dat<-pred.dat%>%group_by(pid)%>%mutate(cens.surv=cumprod(1-pred))

pred.dat<-pred.dat%>%group_by(pid)%>%mutate(rownum=row_number())%>%mutate(maxrow=max(rownum))

cens.prob.sl1<-pred.dat$cens.surv[pred.dat$rownum==pred.dat$maxrow]

#censoring weights
dta_test$cens.wt<-1/cens.prob.sl1

#---------------------------------
#---------------------------------
# time-dependent censoring weights
# obtained using discrete-time SL
# This gives wt_matrix and wt_matrix_eventsonly
#---------------------------------
#---------------------------------

# unique event time points
event.times <- sort(unique(dta_test$time[dta_test$status==1]))
# we add the evaluation time point t=10 as we need it for calculation of CD AUCt
event.times<-c(event.times,10)

#use pred.dat from above to obtain the censoring weight for each person at each event time
ids<-unique(pred.dat$pid)
wt_matrix<-matrix(nrow=dim(dta_test)[1],ncol=length(event.times))

for(i in 1:dim(dta_test)[1]){
  temp<-pred.dat[pred.dat$pid==ids[i],]
  step.temp<-stepfun(temp$tstart,c(1,temp$cens.surv))
  cens.surv.temp<-step.temp(event.times)
  wt_matrix[i,]<-1/ifelse(event.times>dta_test$time[i],NA,cens.surv.temp)
}

# remove last column from weightmatrix as this is not needed in the cindex calculation
wt_matrix_eventsonly <- wt_matrix[,-ncol(wt_matrix)] 

