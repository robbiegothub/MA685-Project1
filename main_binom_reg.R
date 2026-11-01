#IMPORTS
source("log_reg.R")
source('hyper_param.R')
source('ma685_functions.R')
mydata=read.csv('cirrhosis.csv')


#DATA CLEANING
cols = c("Age","Sex","Spiders","Edema","Bilirubin","Albumin","Alk_Phos","SGOT",
         "Platelets","Cholesterol","Tryglicerides","Hepatomegaly","Copper")
clean = mydata[1:312, cols]
clean = clean[complete.cases(clean), ]




#BUILDING MODEL
Y1 = as.integer(clean$Hepatomegaly == "Y")
Y2 = as.integer(clean$Copper > 100)#BINARIZE 

XSet1 = model.matrix(~ Age + Sex + Spiders + Edema, data = clean)[, -1]
XSet2 = model.matrix(~ Age + Sex + Spiders + Edema + Bilirubin + Albumin + Alk_Phos +
                       SGOT + Platelets + Cholesterol + Tryglicerides, data = clean)[, -1]

Model1Set1=LogisticRegressionFull(Y1,XSet1)
Model2Set1=LogisticRegressionFull(Y2,XSet1)

Model1Set2=LogisticRegressionFull(Y1,XSet2)
Model2Set2=LogisticRegressionFull(Y2,XSet2)


#METRICS OUTPUT
Model1Set1$metrics

Model2Set1$metrics


Model1Set2$metrics
Model2Set2$metrics


#TAU TUNING 
t11=tune_seeds(Y1, XSet1)
t11=eval_tau(Y1, XSet1, t11$tau)


t21=tune_seeds(Y2, XSet1)
t21=eval_tau(Y2, XSet1, t21$tau)

t12=tune_seeds(Y1, XSet2)
t12=eval_tau(Y1, XSet2, t12$tau)

t22=tune_seeds(Y2, XSet2)
t22= eval_tau(Y2, XSet2, t22$tau)

t11
t21
t12
t22


max(mean(Y1), 1 - mean(Y1)); max(mean(Y2), 1 - mean(Y2)) #baselines




#GINI ANALYSIS


gini_Y1 = impure(Y1)   # Hepat
gini_Y2 = impure(Y2)   #Copper piss


tau = 0.375
#this is w.r.t gini 75/25 split h0= outcome is at least as imbalanced as a 75/25 split, 
#ha=split is more balanced than 75/25


mc_Y1 = mygini_mc(Y1, tau = tau, nmc = 25000, alpha = 0.05, seed = 1234)
mc_Y2 = mygini_mc(Y2, tau = tau, nmc = 25000, alpha = 0.05, seed = 1234)
#awesome sauce making this function already. was gonna do this, andrealized you alr did it. nice!!


