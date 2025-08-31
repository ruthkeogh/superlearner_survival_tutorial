#################################################################################
#################################################################################
#Super learner for time-to-event outcomes tutorial
#
#The file obtains censoring weights in the test data (dta_test)
#, which are required for estimating predictive performance measures using the functions in functions.R, and in some functions used from packages. 
#
#In this file the weights are obtained using Kaplan-Meier estimator, 
#with censoring assumed to be independent of covariates
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

