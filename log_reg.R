#Logistic Regression functions go here.
if (!requireNamespace("nnet", quietly = TRUE)) {
  install.packages("nnet")
}
library(nnet)



mydata=read.csv('cirrhosis.csv')

log_reg_general <- function(Y,Xmatrix) {
  
  predictors=ncol(Xmatrix)
  M = multinom(Y ~ Xmatrix)
  bhat=coef(M)
  

  return(list(bhats=bhat,model=M))
  #this will be expanded
}

multinom_probs <- function(bhat,X){
  if(is.matrix(bhat)==TRUE){
    b=bhat
  }
  else{
    b=matrix(bhat,nrow=1)
  }
  
  X=cbind(1,as.matrix(X))
  xponent_term=cbind(0,X%*%t(b))
  exponentiated=exp(xponent_term)
  exponentiated/rowSums(exponentiated)
  #this just mirrors slide 4 lec12 w.r.t matrix operations
  
}



rmsefunc <- function(Y,P){
  #this entirely fucking breaks with one NA value... hahah.
  Y<- factor(Y)
  yasnum<-as.integer(Y)
  number=nlevels(Y)
  sequence=seq_len(number)
  onehot=outer(yasnum,sequence,"==")
  rmseval=sqrt(mean(onehot-P)^2)
  return(rmseval)
}


LogisticRegressionFull <- function(Y,Xmatrix){
  main=log_reg_general(Y,Xmatrix)
  bhat=main$bhat
  prob=multinom_probs(bhat,Xmatrix)
  rmses=rmsefunc(Y,prob)
  
  
  returnlist=list(model=main$model,bhatval=bhat,probabilitys=prob,rmseval=rmses)
}
  
asdf=LogisticRegressionFull(Y,X)
asdf$model$n
asdf$bhatval
asdf$rmseval


X1=mydata$N_Days
Y=mydata$Status
X2=mydata$Platelets
k=log_reg_general(Y,cbind(X1,X2))
X=cbind(X1,X2)
length(Y)
k$bhats
p=multinom_probs(k$bhats,X)
#rmsefunc(Y,p) #-- breaks w. any NA. have not done a data cleaning function or looked at how to clean data yet
#so commented this out; though it should work
