# -*- coding: utf-8 -*-
"""
Created on Fri May 18 15:43:32 2018

@author: gfarhani
"""
import pandas as pd
import numpy as np
import numpy
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import LinearSVC
from sklearn import preprocessing
from sklearn.model_selection import cross_val_score
from sklearn.naive_bayes import GaussianNB

# creat classifiers
lr = LogisticRegression()
svc = LinearSVC(C = 1)
rfc = RandomForestClassifier(n_estimators=100)
gnb = GaussianNB()
dataset = numpy.loadtxt("Eureka_1.csv", delimiter=",")
# split into input (X) and output (Y) variables
X = dataset[:,0:-1]
Y = dataset[:,-1]
#X = preprocessing.normalize(X_n)
#X = preprocessing.scale(X)
from sklearn.model_selection import train_test_split
X_train,  X_test, Y_train, Y_test = train_test_split(X,Y, test_size = .4)

lag = lr.fit(X_train,Y_train)
tag = svc.fit(X_train,Y_train)
rag = rfc.fit(X_train,Y_train)
nag = gnb.fit(X_train,Y_train)

pred_lr = lr.predict(X_test)
pred_tag = svc.predict(X_test)
pred_rag = rfc.predict(X_test)
pred_nag = gnb.predict(X_test)

print(pd.crosstab(Y_test,pred_lr,rownames=['Actual Value'],colnames=['Predicted Value']))
print(pd.crosstab(Y_test,pred_tag,rownames=['Actual Value'],colnames=['Predicted Value']))
print(pd.crosstab(Y_test,pred_rag,rownames=['Actual Value'],colnames=['Predicted Value']))
print(pd.crosstab(Y_test,pred_nag,rownames=['Actual Value'],colnames=['Predicted Value']))

A= svc.score(X_test, Y_test)
print('The Score for svc:') 
print(A)

B= lr.score(X_test, Y_test)
print('The Score for lr:') 
print(B)

C = rfc.score(X_test, Y_test)
print('The Score for rfc:') 
print(C)


D = gnb.score(X_test, Y_test)
print('The Score for rfc:') 
print(D)

scores = cross_val_score(lag, X, Y, cv=5)
print(scores)
print(numpy.mean(scores))

scores1 = cross_val_score(tag, X, Y, cv=5)
print('The cross validation on the Score for svc:') 
print(scores1)
print(numpy.mean(scores1))

scores2 = cross_val_score(rag, X, Y, cv=5)
print('The cross validation on the Score for random forest:') 
print(scores2)
print(numpy.mean(scores2))

scores3 = cross_val_score(nag, X, Y, cv=5)
print('The cross validation on the Score for naive Bayessian:') 
print(scores3)
print(numpy.mean(scores3))

##########
from sklearn.decomposition import PCA
pca = PCA(n_components=10)
pca.fit(X)
X_reduced = pca.transform(X)
import pylab as pl
pl.scatter(X_reduced[:, 1], X_reduced[:, 2], X_reduced[:, 3], c=Y,
           cmap='RdYlBu')



from sklearn.cluster import KMeans
k_means = KMeans(n_clusters=3, random_state=0) # Fixing the RNG in kmeans


k_means.fit(X_reduced)
y_pred = k_means.predict(X_reduced)
import matplotlib.pyplot as pl
pl.scatter(X[:, 0], X[:, 1], X[:,3], X[:,4], c=y_pred, cmap='RdYlBu');
