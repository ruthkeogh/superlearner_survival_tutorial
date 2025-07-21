#################################################################################
#################################################################################
#Super learner for time-to-event outcomes tutorial
#
#packages needed
#################################################################################
#################################################################################

library(tidyverse)
library(gam)
library(nnls)
library(survival)
library(riskRegression)
library(pec)
library(timeROC)
library(ranger)
library(randomForest)
library(randomForestSRC)
library(data.table)
library(glmnet)

#---
#Discrete-time superLearner: Polley and van der Laan 2011
library(SuperLearner)

#---
#Continuous-time superLearner: Westling et al. 2023
# devtools::install_github("tedwestling/survSuperLearner")
library(survSuperLearner)

#---
#Continuous-time superLearner (statelearner): Munch & Gerds 2024
#Need to run the functions in this folder
#https://github.com/amnudn/joint-survival-super-learner/blob/main/R-code/functions


