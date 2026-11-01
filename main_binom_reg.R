source("log_reg.R")
source('hyper_param.R')
source('ma685_functions.R')
mydata=read.csv('cirrhosis.csv')

cols = c("Age","Sex","Spiders","Edema","Bilirubin","Albumin","Alk_Phos","SGOT",
         "Platelets","Cholesterol","Tryglicerides","Hepatomegaly","Copper")
clean = mydata[1:312, cols]
clean = clean[complete.cases(clean), ]

Y1 = as.integer(clean$Hepatomegaly == "Y")
Y2 = as.integer(clean$Copper > 100)

XSet1 = model.matrix(~ Age + Sex + Spiders + Edema, data = clean)[, -1]
XSet2 = model.matrix(~ Age + Sex + Spiders + Edema + Bilirubin + Albumin + Alk_Phos +
                       SGOT + Platelets + Cholesterol + Tryglicerides, data = clean)[, -1]

Model1Set1=LogisticRegressionFull(Y1,XSet1)
Model2Set1=LogisticRegressionFull(Y2,XSet1)

Model1Set2=LogisticRegressionFull(Y1,XSet2)
Model2Set2=LogisticRegressionFull(Y2,XSet2)



Model1Set1$metrics

Model2Set1$metrics


Model1Set2$metrics
Model2Set2$metrics

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


myconfusion
