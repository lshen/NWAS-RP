
####################################################################################
#   GWAS reprioritization
#   Reprioritize GWAS results using SVM, SVR and Ridge regression,
#   for next step module identification
# --------------------------------------------------------------------
#   Input:
#   - wm: genetics function interaction matrix
#   - wm.train: subset rows of wm, for training model
#   - y.train: subset of GWAS p-values corresponding to wm.train, for training model
#   - paras, parameters: C.svm, C.svr,epsilon, lambda (unknown, tuned from other functions)
# 
#   Output:
#   - y.pred: distance from SVM hyperplan, or regression response from SVR and Ridge
# ------------------------------------------
#   Author: Xiaohui Yao, yao2@umail.iu.edu
#   Date created: Oct-20-2016
#   Date updated: Nov-06-2016
#   @Indiana University.
####################################################################################

############
### SVM ####
############
# load required library
library(kernlab)
threshold = 0.01 # nominal significant or non-significant

NetWAS.SVM <- function(wm, y, wm.train, y.train, C.svm)
{
  # train svm model using training data, we use 'C-svc' and 'rbf' kernel in experiment
  svm.model <- ksvm(wm.train, y.train, type='C-svc', kernel='rbfdot', C=C.svm)
  
  # predict for all GWAS
  y.pred.svm <- predict(svm.model, wm, type='decision')
  return(y.pred.svm)
}

############
### SVR ####
############
# load required library
library(kernlab)

NetWAS.SVR <- function(wm, y, wm.train, y.train, C.svr, epsilon)
{
  # train SVR model using training data, we use 'eps-svr' and 'rbf' kernel in experiment
  svr.model <- ksvm(wm.train, y.train, type='eps-svr', kernel='rbfdot', C=C.svr, epsilon=epsilon)
  
  # predict for all GWAS
  y.pred.svr <- predict(svr.model, wm, type='response') 
  return(y.pred.svr)
}


##############
### Ridge ####
##############
# load required library
library(glmnet)

NetWAS.Ridge <- function(wm, y, wm.train, y.train, lambda)
{
  # train Ridge model using training data
  ridge.model <- glmnet(wm.train, y.train, family='gaussian', alpha=0, lambda=lambda)
  # predict for all GWAS
  y.pred.ridge <- predict(ridge.model, wm, type='response')
  return(y.pred.ridge)
}