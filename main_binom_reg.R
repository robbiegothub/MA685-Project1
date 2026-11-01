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
tune11 = tune_seeds(Y1, XSet1); eval11 = eval_tau(Y1, XSet1, tune11$tau)
tune21 = tune_seeds(Y2, XSet1); eval21 = eval_tau(Y2, XSet1, tune21$tau)
tune12 = tune_seeds(Y1, XSet2); eval12 = eval_tau(Y1, XSet2, tune12$tau)
tune22 = tune_seeds(Y2, XSet2); eval22 = eval_tau(Y2, XSet2, tune22$tau)



results = rbind(
  Hep_Set1 = c(tau = tune11$tau, eval11["mean", ]),
  Cu_Set1  = c(tau = tune21$tau, eval21["mean", ]),
  Hep_Set2 = c(tau = tune12$tau, eval12["mean", ]),
  Cu_Set2  = c(tau = tune22$tau, eval22["mean", ]))
round(results, 3)


max(mean(Y1), 1 - mean(Y1)); max(mean(Y2), 1 - mean(Y2)) #baselines




#GINI ANALYSIS


gini_Y1 = impure(Y1)   # Hepat
gini_Y2 = impure(Y2)   #Copper piss

gini_Y1
gini_Y2

tau = 0.375
#this is w.r.t gini 75/25 split h0= outcome is at least as imbalanced as a 75/25 split, 
#ha=split is more balanced than 75/25


mc_Y1 = mygini_mc(Y1, tau = tau, nmc = 25000, alpha = 0.05, seed = 1234)
mc_Y2 = mygini_mc(Y2, tau = tau, nmc = 25000, alpha = 0.05, seed = 1234)
#awesome sauce making this function already. was gonna do this, andrealized you alr did it. nice!!
gini_results = rbind(
  Hepatomegaly = c(Eg = mc_Y1$Eg, crit = unname(mc_Y1$gini_crit), pval = mc_Y1$pval),
  Copper       = c(Eg = mc_Y2$Eg, crit = unname(mc_Y2$gini_crit), pval = mc_Y2$pval))
round(gini_results, 4)
gini_results
