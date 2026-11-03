#IMPORTS
source("log_reg.R")
source('hyper_param.R')
source('ma685_functions.R')
mydata=read.csv('cirrhosis.csv')



#DATA CLEANING
cols = c("Age","Sex","Spiders","Edema","Bilirubin","Albumin","Alk_Phos","SGOT",
         "Platelets","Cholesterol","Tryglicerides","Hepatomegaly","Copper","Stage")
clean = mydata[1:312, cols]
clean = clean[complete.cases(clean), ]


y1 = factor(clean$Stage)

XSet1 = model.matrix(~ Age + Sex + Spiders + Edema, data = clean)[, -1]
XSet2 = model.matrix(~ Age + Sex + Spiders + Edema + Bilirubin + Albumin + Alk_Phos +
                       SGOT + Platelets + Cholesterol + Tryglicerides, data = clean)[, -1]
XSet3= model.matrix(~ Age + Sex + Spiders + Edema + Bilirubin + Albumin + Alk_Phos +
                             SGOT + Platelets + Cholesterol + Tryglicerides + Hepatomegaly + Copper, data = clean)[, -1]

Model1Set1=LogisticRegressionFull(y1,XSet1)
Model1Set2=LogisticRegressionFull(y1,XSet2)
Model1Set3=LogisticRegressionFull(y1,XSet3)


Model1Set1$metrics
Model1Set2$metrics
Model1Set3$metrics







gini_y1 = impure(y1)   # stages

gini_y1


tau = 0.5625
#this is w.r.t gini 62.5/12.5/12.5/12.5 split h0= outcome is at least as imbalanced as a 62.5/12.5/12.5/12.5 split, 
#ha=split is more balanced than 62.5/12.5/12.5/12.5


mc_Y1 = mygini_mc(y1, tau = tau, nmc = 25000, alpha = 0.05, seed = 1234)
#awesome sauce making this function already. was gonna do this, andrealized you alr did it. nice!!
gini_results = rbind(
  Stages = c(Eg = mc_Y1$Eg, crit = unname(mc_Y1$gini_crit), pval = mc_Y1$pval))
round(gini_results, 4)
gini_results

#counts and proportion. was lazy sorry.
sum(y1==4)
sum(y1==3)
sum(y1==2)
sum(y1==1)
sum(y1==4)/length(y1)
sum(y1==3)/length(y1)
sum(y1==2)/length(y1)
sum(y1==1)/length(y1)


eval_multi(y1, XSet1); eval_multi(y1, XSet2); eval_multi(y1, XSet3)



#ridge logistic regression

lams = c(0, 10^seq(-3, 2.5, 0.5))
tuneS1 = tune_decay(y1, XSet1, lams = lams)
tuneS2 = tune_decay(y1, XSet2, lams = lams)
tuneS3 = tune_decay(y1, XSet3, lams = lams)
rbind(bacc = c(tuneS1$lam, tuneS2$lam, tuneS3$lam),
      acc  = c(tuneS1$lam_acc, tuneS2$lam_acc, tuneS3$lam_acc))

c(tuneS1$lam, tuneS2$lam, tuneS3$lam)


eval_multi(y1, XSet1, decay = tuneS1$lam)
eval_multi(y1, XSet2, decay = tuneS2$lam)
eval_multi(y1, XSet3, decay = tuneS3$lam)



#ridge logistic regression with noisy bootstrap sampling

tuneA1 = tune_decay(y1, XSet1, aug = TRUE,lams = lams)            # retune lambda with augmentation
eval_multi(y1, XSet1, decay = tuneA1$lam, aug = TRUE)


tuneA2 = tune_decay(y1, XSet2, aug = TRUE,lams = lams)            # retune lambda with augmentation
eval_multi(y1, XSet2, decay = tuneA2$lam, aug = TRUE)

tuneA3 = tune_decay(y1, XSet3, aug = TRUE,lams = lams)            # retune lambda with augmentation
eval_multi(y1, XSet3, decay = tuneA3$lam, aug = TRUE)


tuneA1$lam 
tuneA2$lam 
tuneA3$lam

plot(log10(tuneA1$lams + 1e-4), tuneA1$bacc, type = 'b')
