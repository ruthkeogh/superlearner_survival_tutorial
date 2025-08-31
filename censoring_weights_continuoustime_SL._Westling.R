#################################################################################
#################################################################################
#Super learner for time-to-event outcomes tutorial
#
#The file obtains censoring weights in the test data (dta_test)
#, which are required for estimating predictive performance measures using the functions in functions.R, and in some functions used from packages. 
#
#In this file the weights are obtained using the continuous-time super learner of Westling et al.
#
#For some measures the weights are required only at the person's observed event/censoring time. 
#For some measures (C-index and C/D AUCt) the weights are required for each person at every observed event time (i.e. a matrix of weights for each person at each time).
#################################################################################
#################################################################################

#---------------------------------
#---------------------------------
# time-fixed censoring weights at observed event or censoring times
# obtained using continuous-time SL, using survSuperLearner
# This gives dta_test$cens.wt
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
dta_test$cens.wt<-1/diag(cens.prob.sl2)

#---------------------------------
#---------------------------------
# time-dependent censoring weights
# obtained using continuous-time SL, using survSuperLearner
# This gives wt_matrix and wt_matrix_eventsonly
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
wt_matrix<-1/cens.prob.sl2

# remove last column from weightmatrix as this is not needed in the cindex calculation
wt_matrix_eventsonly <- wt_matrix[,-ncol(wt_matrix)] 
