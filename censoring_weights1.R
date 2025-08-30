#################################################################################
#################################################################################
#Super learner for time-to-event outcomes tutorial
#
#The file obtains censoring weights in the test data (dta_test)
#, which are required for estimating predictive performance measures using the functions in functions.R, and in some functions used from packages. 
#
#Two methods are used, resulting in different sets of weights:
#(1) Censoring assumed to be independent of covariates: weights obtained using Kaplan-Meier estimator.
#(2) Censoring allowed to depend on covariates: weights obtained using super learner. This is done separately using each of the three SL methods. 
#
#For some measures the weights are required only at the person's observed event/censoring time. 
#For some measures (C-index and C/D AUCt) the weights are required for each person at every observed event time (i.e. a matrix of weights for each person at each time).
#################################################################################
#################################################################################

#---------------------------------
#---------------------------------
# time-fixed censoring weights at observed event or censoring times
# obtained using Kaplan-Meier
# This gives dta_test$cens.wt
#---------------------------------
#---------------------------------

cens.km<-survfit(Surv(time,1-status)~1,data=dta_test)
cens.func<-stepfun(cens.km$time,c(1,cens.km$surv))

dta_test$cens.wt<-1/cens.func(dta_test$time)
dta_test$cens.wt[dta_test$time==3651/365] <- 1/cens.func(10)

#---------------------------------
#---------------------------------
# time-dependent censoring weights
# obtained using Kaplan-Meier
# This gives wt_matrix and wt_matrix_eventsonly
#---------------------------------
#---------------------------------

# unique event time points
event.times <- sort(unique(dta_test$time[dta_test$status==1]))
# we add the evaluation time point t=10 as we need it for calculation of CD AUCt
event.times<-c(event.times,10)

#split data on event.times
dta_test_split <- survSplit( Surv(time,status) ~., data = dta_test, cut = event.times)

# calculate weights for standard censoring at the end time points
dta_test_split$cens.wt<-1/cens.func(dta_test_split$time)
dta_test_split$cens.wt[dta_test_split$time==3651/365] <- 1/cens.func(10)

#--------------
# construct weights matrix needed for cindex / auct 
# select the weights at event time points + make sure subjects who are censored before the first event time point are kept in the dataset
# note that the weight in NA after the person's event/cens time
# put weights in wide format
dta_test_wide<- dta_test_split %>% 
  filter(time %in% event.times | (time < min(event.times))) %>% 
  select(c("pid","time","cens.wt")) %>% 
  spread(time, cens.wt)  

# subjects with censoring time before first event time should add a row but not a column to the weightsmatrix (only columns for event time points are needed)
n.cens.before.first.event <- sum(dta_test_split$time < min(event.times))
dta_test_wide <- as.matrix(dta_test_wide[,-1:-(1+n.cens.before.first.event)]) #one additional column is deleted ("id")

wt_matrix <- dta_test_wide

# remove last column from weightmatrix as this is not needed in the cindex calculation
wt_matrix_eventsonly <- wt_matrix[,-ncol(wt_matrix)] 

#---------------------------------
#---------------------------------
# time-fixed censoring weights at observed event or censoring times
# obtained using discrete-time SL
# This gives dta_test$cens.wt.sl1
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
dta_test$cens.wt.sl1<-1/cens.prob.sl1

#---------------------------------
#---------------------------------
# time-dependent censoring weights
# obtained using discrete-time SL
# This gives wt_matrix_sl1 and wt_matrix_sl1_eventsonly
#---------------------------------
#---------------------------------

# unique event time points
event.times <- sort(unique(dta_test$time[dta_test$status==1]))
# we add the evaluation time point t=10 as we need it for calculation of CD AUCt
event.times<-c(event.times,10)

#use pred.dat from above to obtain the censoring weight for each person at each event time
ids<-unique(pred.dat$pid)
wt_matrix_sl1<-matrix(nrow=dim(dta_test)[1],ncol=length(event.times))

for(i in 1:dim(dta_test)[1]){
  temp<-pred.dat[pred.dat$pid==ids[i],]
  step.temp<-stepfun(temp$tstart,c(1,temp$cens.surv))
  cens.surv.temp<-step.temp(event.times)
  wt_matrix_sl1[i,]<-1/ifelse(event.times>dta_test$time[i],NA,cens.surv.temp)
}

# remove last column from weightmatrix as this is not needed in the cindex calculation
wt_matrix_sl1_eventsonly <- wt_matrix_sl1[,-ncol(wt_matrix_sl1)] 

#---------------------------------
#---------------------------------
# time-fixed censoring weights at observed event or censoring times
# obtained using continuous-time SL, using survSuperLearner
# This gives dta_test$cens.wt.sl2
#---------------------------------
#---------------------------------

#times at which we want to obtain censoring survival probabilities
#this is all times except the administrative censoring time
new.times<-dta_test$time
new.times<-ifelse(new.times==3651/365,10,new.times)

#Screening algorithms
#Note: Not sure how to allow different screening options for different algorithms
event.SL.library <- cens.SL.library <- lapply(c("survSL.km", "survSL.coxph", "survSL.rfsrc"), 
                                              function(alg) {
                                                c(alg,"survscreen.glmnet","All")
                                              })

covs_test<-dta_test[,c("year1","year2","age","meno","size1","size2","grade","nodes","pgr","er","hormon","chemo")]
set.seed(1)
sl <- survSuperLearner(time = dta_test$time, event = dta_test$status, X = covs_test, newX = covs_test, new.times = new.times, 
                       event.SL.library = event.SL.library, cens.SL.library = cens.SL.library, verbose = TRUE,
                       cvControl=list(V = 5),
                       control=list(saveFitLibrary = TRUE))

#obtain censoring probabilities
#For each person up to their observed event/cens time
cens.prob.sl2<-sl$cens.SL.predict

#censoring weights - use time 10 for people administratively censored just after time 10
dta_test$cens.wt.sl2<-1/diag(cens.prob.sl2)

#---------------------------------
#---------------------------------
# time-dependent censoring weights
# obtained using continuous-time SL, using survSuperLearner
# This gives wt_matrix_sl2 and wt_matrix_sl2_eventsonly
#---------------------------------
#---------------------------------

# unique event time points
event.times <- sort(unique(dta_test$time[dta_test$status==1]))
# we add the evaluation time point t=10 as we need it for calculation of CD AUCt
event.times<-c(event.times,10)

#SL: is same as above, except we obtain the predictions at event.times
set.seed(1)
sl <- survSuperLearner(time = dta_test$time, event = dta_test$status, X = covs_test, newX = covs_test, new.times = event.times, 
                       event.SL.library = event.SL.library, cens.SL.library = cens.SL.library, verbose = TRUE,
                       cvControl=list(V = 5),
                       control=list(saveFitLibrary = TRUE))

#obtain censoring probabilities
#For each person up to each event time
cens.prob.sl2<-sl$cens.SL.predict

#censoring weight for each person at each event time
wt_matrix_sl2<-1/cens.prob.sl2

# remove last column from weightmatrix as this is not needed in the cindex calculation
wt_matrix_sl2_eventsonly <- wt_matrix_sl2[,-ncol(wt_matrix_sl2)] 

#---------------------------------
#---------------------------------
# time-fixed censoring weights at observed event or censoring times
# obtained using continuous-time SL, using statelearner
# This gives dta_test$cens.wt.sl3
#---------------------------------
#---------------------------------

#make data into data.table, as required by the statelearner code
dta_test2<-dta_test
dta_test2$status<-dta_test2$status2
dta_test2<-as.data.table(dta_test2)

#list of learners to be used
learners <- list(
  N_Aa = list(model = "cox", x_form = ~1),
  cox = list(model = "cox", x_form = ~ year1+year2+age+meno+size1+size2+grade+nodes+pgr+er+hormon+chemo),
  cox_penalty = list(model = "GLMnet", x_form = ~year1+year2+age+meno+size1+size2+grade+nodes+pgr+er+hormon+chemo),
  rf = list(model = "rfsrc", x_form = ~year1+year2+age+meno+size1+size2+grade+nodes+pgr+er+hormon+chemo, ntree = 50)
)

#apply statelearner
set.seed(1)
sl = statelearner(learners = list(cause1 = learners,
                                  censor = learners),
                  data = dta_test2,
                  time = 10,
                  integrate = TRUE,
                  split.method = "cv5",
                  B=1,
                  time_grid_length=100,
                  vars = NULL,
                  collapse = TRUE)

#the best fitting model for censoring is a Cox model, an we use this below to obtain the censoring weights
sl$fitted_winners

#fit Cox model for censoring: noting that we use the status2=-1 indicator (which denotes censoring before time 10, i.e. excludes administrative censoring.
cox.cens<-coxph(Surv(time,status2==-1)~year1+year2+age+meno+size1+size2+grade+nodes+pgr+er+hormon+chemo,data=dta_test,x=TRUE)

cens.prob.sl3<-survfit(cox.cens,newdata=dta_test)

#censoring weights
dta_test$cens.wt.sl3<-1/sapply(1:dim(dta_test)[1],FUN=function(x){cens.prob.sl3$surv[which(dta_test$time[x]==cens.prob.sl3$time),x]})

#---------------------------------
#---------------------------------
# time-dependent censoring weights
# obtained using continuous-time SL, using statelearner
# This gives wt_matrix_sl3 and wt_matrix_sl3_eventsonly
#---------------------------------
#---------------------------------

# unique event time points
event.times <- sort(unique(dta_test$time[dta_test$status==1]))
# we add the evaluation time point t=10 as we need it for calculation of CD AUCt
event.times<-c(event.times,10)

#use pred.dat from above to obtain the censoring weight for each person at each event time
ids<-unique(pred.dat$pid)
wt_matrix_sl3<-matrix(nrow=dim(dta_test)[1],ncol=length(event.times))

for(i in 1:dim(dta_test)[1]){
  step.temp<-stepfun(cens.prob.sl3$time,c(1,cens.prob.sl3$surv[,1]))
  cens.surv.temp<-step.temp(event.times)
  wt_matrix_sl3[i,]<-1/ifelse(event.times>dta_test$time[i],NA,cens.surv.temp)
}

# remove last column from weightmatrix as this is not needed in the cindex calculation
wt_matrix_sl3_eventsonly <- wt_matrix_sl3[,-ncol(wt_matrix_sl3)] 
