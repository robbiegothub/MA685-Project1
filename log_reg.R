#Logistic Regression functions go here.
if (!requireNamespace("nnet", quietly = TRUE)) {
  install.packages("nnet")
}
library(nnet)

#testing

mydata=read.csv('cirrhosis.csv')

X1=mydata$N_Days
Y=mydata$Status
X2=mydata$Platelets
X=cbind(X1,X2)


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
  rmseval=sqrt(mean((onehot-P)^2))
  return(rmseval)
}











null_logreg=function(Y){
  n=length(Y)
  Y=factor(Y)
  mnull=multinom(Y~1)
  bhatnull=coef(mnull)
  bhatnull=c(0,as.numeric(bhatnull))
  phatnull= exp(bhatnull)/sum(exp(bhatnull))# as.numeric(rep(1/(1+exp(-bhatnull)),n))
  lhatnull=sum(table(Y)*log(phatnull))
  retlist=list(phat=phatnull,lhat=lhatnull)
  return(retlist)
  }
null_logreg(Y)








LogisticRegressionFull <- function(Y,Xmatrix){
  main=log_reg_general(Y,Xmatrix)
  bhat=main$bhat
  prob=multinom_probs(bhat,Xmatrix)
  
  lhat=main$model$value
  rmses=rmsefunc(Y,prob)
  lhatnull=null_logreg(Y)$lhat
  rmcfad=1-(-(lhat))/(lhatnull)
  metricslist=list(AIC=main$AIC,RMSE=rmses,BIC=BIC(main$model),R_mcf=rmcfad)#can we use BIC()?
  returnlist=list(model=main$model,bhatval=bhat,probabilitys=prob,metrics=metricslist)
}



asdf=LogisticRegressionFull(Y,X)
asdf$model$value
asdf$bhatval
asdf$metrics$R_mcf
asdf$model$AIC
BIC(asdf$model)


#rmsefunc(Y,p) #-- breaks w. any NA. have not done a data cleaning function or looked at how to clean data yet
#so commented this out; though it should work
