#################################################################################
#################################################################################
#Super learner for time-to-event outcomes tutorial
#
#This file sets up the rotterdam data and creates the training and test data sets
#################################################################################
#################################################################################

#---------------------------------
#---------------------------------
#data set-up
#---------------------------------
#---------------------------------

head(rotterdam)

dta<-rotterdam

n<-dim(dta)[1]

#implement administrative censoring  at 3651 days (just after 10 years=3650 days)
#status: 1=death, -1=censoring, 0=administrative censoring
#status2: 1=death, -1=censoring, 0=administrative censoring

dta$status<-dta$death
dta$status<-ifelse(dta$dtime>3650,0,dta$status)
dta$time<-ifelse(dta$dtime>3650,3651,dta$dtime)#we want the censoring time to be just after time 10 years, to avoid messing up various functions later

dta$status2<-ifelse(dta$death==1,1,-1)
dta$status2<-ifelse(dta$dtime>3650,0,dta$status2)

dta$time<-dta$time/365 #make time in years

#create dummy variables for non-binary categorical variable 'size'
dta$size1<-ifelse(dta$size=="20-50",1,0)
dta$size2<-ifelse(dta$size==">50",1,0)

#make 'grade' a 0/1 variable (it is originally coded as 2/3)
dta$grade<-ifelse(dta$grade==2,0,1)

#categorize year into three periods
dta$year0<-ifelse(dta$year<=1985,1,0)
dta$year1<-ifelse(dta$year>=1986 & dta$year<=1989,1,0)
dta$year2<-ifelse(dta$year>=1990 & dta$year<=1993,1,0)

#create training and test data sets
set.seed(1)
train.sample<-sort(sample(1:n,2087,replace = F))
dta_train<-dta[train.sample,]
dta_test<-dta[!(1:n)%in%train.sample,]

