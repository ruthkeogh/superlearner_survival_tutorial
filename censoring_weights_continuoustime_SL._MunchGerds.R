#################################################################################
#################################################################################
#Super learner for time-to-event outcomes tutorial
#
#The file obtains censoring weights in the test data (dta_test)
#, which are required for estimating predictive performance measures using the functions in functions.R, and in some functions used from packages. 
#
#In this file the weights are obtained using the continuous-time super learner of Munch & Gerds.
#
#For some measures the weights are required only at the person's observed event/censoring time. 
#For some measures (C-index and C/D AUCt) the weights are required for each person at every observed event time (i.e. a matrix of weights for each person at each time).
#################################################################################
#################################################################################

#---------------------------------
#---------------------------------
# time-fixed censoring weights at observed event or censoring times
# obtained using continuous-time SL, using statelearner
# This gives dta_test$cens.wt
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
dta_test$cens.wt<-1/sapply(1:dim(dta_test)[1],FUN=function(x){cens.prob.sl3$surv[which(dta_test$time[x]==cens.prob.sl3$time),x]})

#---------------------------------
#---------------------------------
# time-dependent censoring weights
# obtained using continuous-time SL, using statelearner
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
  step.temp<-stepfun(cens.prob.sl3$time,c(1,cens.prob.sl3$surv[,1]))
  cens.surv.temp<-step.temp(event.times)
  wt_matrix[i,]<-1/ifelse(event.times>dta_test$time[i],NA,cens.surv.temp)
}

# remove last column from weightmatrix as this is not needed in the cindex calculation
wt_matrix_eventsonly <- wt_matrix[,-ncol(wt_matrix)] 
