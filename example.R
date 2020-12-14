
###################################################
# Run this example to see how to implement 
# machine learning methods to GWAS reprioritization
# ------------------------------------
# Author: Xiaohui Yao, yao2@umail.iu.edu
# @Indiana University
################################################### 

source('repriori.R')

# load data
load('example.RData') # wm.example; y.example.classification; y.example.regression; gene

# set parameters, should be tuned before running.
# As an example, we fix them.
C.svm = 1
C.svr = 1
epsilon = 0.01
lambda = 20

n <- nrow(wm.example)
# 1/3 for train model
idx.train <- sample(1:n)[1:(0.3*n)]
wm.train <- wm.example[idx.train,]
y.train.c <- y.example.classification[idx.train]
y.train.r <- y.example.regression[idx.train]

# run svm, svr and ridge method to re-prioritize GWAS results
### SVM
repri.svm <- NetWAS.SVM(wm.example, y.example, wm.train, y.train.c, C.svm)
result.svm <- cbind(gene.example, repri.svm)
result.svm <- result.svm[order(-result.svm$repri.svm),]

### SVR
repri.svr <- NetWAS.SVR(wm.example, y.example, wm.train, y.train.r, C.svr, epsilon)
result.svr <- cbind(gene.example, repri.svr)
result.svr <- result.svr[order(-result.svr$repri.svr),]

### Ridge
repri.ridge <- NetWAS.Ridge(wm.example, y.example, wm.train, y.train.r, lambda)
result.ridge <- cbind(gene.example, repri.ridge)
result.ridge <- result.ridge[order(-result.ridge$s0),]