wdbc=read.csv(file.choose(),sep=",",header =TRUE)
dim(wdbc)
#convert the features of the data: wdbs.data
wdbc.data=as.matrix(wdbc[,c(3:12)])
#set the row names of wdbc.data
row.names(wdbc.data)=wdbc$id
#create diagnosis vector
diagnosis=as.numeric(wdbc$diagnosis=="M")
head(diagnosis)
#summary of data
summary(wdbc.data)
str(wdbc.data)
# total no of observation of malignant diagnosis
table(wdbc$diagnosis)
# what is mean of each of the columns ?
round(colMeans(wdbc.data),2)
# what is sd of each of the columns ?
roundSD=function(x){
  round(sd(x),2)
}
apply(wdbc.data,2,roundSD)
# how the variables related to each other ?
library(corrplot)
corMatrix=wdbc[,c(3:12)]
# rename the columns ?
cNames=c( "rad_w","txt_w","per_w","are_w","smt_w","cmp_w", "con_w","ccp_w","sym_w","frd_w")
colnames(corMatrix)=cNames
# create the correlation matrix 
M=round(cor(corMatrix),2)

# create corrplot
corrplot(M,diag=FALSE,method="color",order="FPC",tl.srt=90)
# from the corrplot it is evident that there 
# are many variable that are highly correlated with each other


#Principle component Analysis
#  why PCA ? Due to the number of variables in the model 
# ,we can try using a dimensionality reduction technique to unevil 
# any patterns in the data. As mentioned in the Exploratory 
# Data Analysis section, there are thirty variables that when combined
# can be used to model each patient's diagnosis.Using PCA we can combine
# our many variables into different linear combinations that each
# explain a part of the variance of model . By proceeding with PCA
# we are assuming the linearity of the of our variables within
# dataset.By choosing only the linear combinations that provide a majority (>=85%)
# of the covariance, we can reduce the complexity of our model.
# We can then more easily see how the model works and provide meaningful 
# graphs and representations of our complext dataset.
# 
# The first step in doing a PCA , is to ask ourselves whether or not 
# tha data should be scaled to unit variance .That is ,to bring all the numeric
# variables to the same scale .In other words, we are trying to determine 
# whether we should use a correlation matrix or covariance matrix in our calculationns of eigen values 
# and eigen vectors .


# Running PCA using correlationn matrix:
# when the correlation matrix is used to calculate the eigen
# values and eigen vectors, we use the prcomp()funtion.
wdbc.pr=prcomp(wdbc.data,scale=TRUE,center=TRUE)
attributes(wdbc.pr)
summary(wdbc.pr)
library(psych)
pairs.panels(wdbc.pr$x[,(1:3)],gap=0,bg=c("green","red")[wdbc$diagnosis],pch=21,main="scatter plot of PC")


# Let's visualize this using a scree plot
#Calculate variablity of each component 
pr.var=wdbc.pr$sdev^2
pr.var

# Variance explained by each principal component :pve
pve=pr.var/sum(pr.var)
pve
# eigen values
round(pr.var,2)
# percent variation explained
round(pve,2)
# cumulative percent explained
round(cumsum(pve),2)
# create a plot of variance explained for each principal
# component .
plot(pve,xlab="principal component",ylab="Proportion of
     variance explained ",ylim=c(0,1),type="b",main="Scree plot")
# plot cumulative proportion of variance explained 
plot(cumsum(pve),xlab="principal component ",ylab ="
     cumulative var explained",
     ylim=c(0,1),type="b",main="scree plot (cumulative var)")

# 88% variation is explained by the first six PC's.Moreover,
# the eigen values associated with the first 6 PC's are greater 
# than 1 .We will use this criteria to decide on how many PC's 
# to include in the model building phase .

# Next ,lets create a scaater plot observations by principal components
# one and two :
plot(wdbc.pr$x[,c(1,2)],col=(diagnosis+1),xlab ="PC1",
     ylab="PC2",main="scatter plot of PC1 VS PC2")
legend(x="topleft",pch=1,col=c("red","black"),
       legend=c("B","M"))

# There is clear separation of diagnosis(M or B) that is
# evident in the PC1 VS PC2 plot.


# conclusion:By using PCA we took a complex model of 30 predictors the model
# down to six linear combinations of the various predictors.


# Linear Discriminant Analysis(LDA)
#   From the principal component's scatter plots it is evident 
#   that there is some clustering of benign and malignant point.
#   This suggests that we could build a linear discriminant 
#   function using these principal components.Now that we have 
#   our chosen principal components we can perform the linear 
#   discriminant analysis.

#         -----Model building and validation-------
# Here's the high level process followed:
#  * Build the model using training data
#  * Predict using the test data
#  * Evaluate model performance using ROC and AUC

# Our next tast is to use the first 6 PCs to build a Linear discrimainant
# function using the lda() function in R.
# From the wdbc.pr object ,we need to extract the first six PCs
# To do this lets first check available for this object .
ls(wdbc.pr)
# We are interested in the rotation (also called loadings) of the 
# first six PCs multiplied by the scaled data,which are called
# scores (basically pc  transformed data)
wdbc.pcs=wdbc.pr$x[,1:3]
head(wdbc.pcs,20)
# here the rownames help us see the how PC transformed data  
# looks like .
# Now ,we need to append the diagnosis column to this PC transformed 
# data frame wdbc.pcs.Lets call the new data frame as
# wdbc.pcst.
wdbc.pcst=wdbc.pcs
wdbc.pcst=cbind(wdbc.pcs,diagnosis)
head(wdbc.pcst,25)
# Here ,diagnosis==1 represents malignant
# and diagnosis==0 represents benign


#   ------split the dataset into training/test  data----
# using the training data we can build the LDA function .
# Next ,we use the test data to make predictions.
# calculate N
N=nrow(wdbc.pcst)
N
ind=sample(2,nrow(wdbc.pcst),replac=TRUE,prob=c(0.8,0.2))
wdbc.pcst.train=wdbc.pcst[ind==1,]
wdbc.pcst.test=wdbc.pcst[ind==2,]
nrow(wdbc.pcst.train)
nrow(wdbc.pcst.test)
# so 442 observations are in the training dataset
#  and 127 observations are in the test dataset.
#  We will use the training dataset to calculate the 
# linear discriminant function  by passing it to the lda()
# fucntion to the MASS PACAKAGE
library(MASS)

wdbc.pcst.train.df=wdbc.pcst.train
# convert matrix to a dataframe 
wdbc.pcst.train.df=as.data.frame(wdbc.pcst.train)
wdbc.pcst.test.df=as.data.frame(wdbc.pcst.test)
#       PERFORMANCE LDA ON DIAGNOSIS
wdbc.lda=lda(diagnosis~PC1+PC2+PC3,data=wdbc.pcst.train.df)
# lets summarize the LDA OUTPUT
attributes(wdbc.lda)
head(wdbc.lda$prior)
wdbc.lda$counts  
wdbc.lda$scaling
p=predict(wdbc.lda,wdbc.pcst.train.df)$class
# confusion matrix and accuracy ~ training data       
tab=table(predicted=p,actual=wdbc.pcst.train.df$diagnosis)
tab
table(wdbc.pcst.train.df$diagnosis)

# accuracy of training data
sum(diag(tab))/sum(tab)
# confusion matrix and accuracy ~ testting  data
p1=predict(wdbc.lda,wdbc.pcst.test.df)$class
tab1=table(predicted=p1,actual=wdbc.pcst.test.df$diagnosis)
tab1
table(wdbc.pcst.test.df$diagnosis)
# accuracy of testing data
sum(diag(tab1))/sum(tab1)
# model peroformance evaluation
library(ROCR)
pred=predict(wdbc.lda,wdbc.pcst.train.df,type='prob')     
hist(pred$posterior[,2])
pred=prediction(pred$posterior[,2],wdbc.pcst.train.df$diagnosis)
roc=performance(pred,"tpr","fpr")
plot(roc,colorize=T)
abline(a=0,b=1)
# area under the curve
auc=performance(pred,"auc")
auc=unlist(slot(auc,"y.values"))
auc
auc=round(auc,4)
legend(0.6,0.2,auc,title="AUC",cex=0.5)     
# conclusion 
# We have shwon how dimensionality reduction technique like 
# principal component analysis can be used to reduce a lage number 
# of highly correlated predictors to small set of linear combinations
# of those predictors.In doing so, we unveiled patterens
# in the data which led us to build a classification rule using
# linear discriminant analysis .By applaying the clssification
# rule we have constructed a diagnostic system that predicts
# malignant tumors at 99.9%

