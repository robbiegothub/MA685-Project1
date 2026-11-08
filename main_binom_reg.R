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


Model1Set2$bhatval

#TAU TUNING 
tune11 = tune_seeds(Y1, XSet1); eval11 = eval_tau(Y1, XSet1, tune11$tau)
tune21 = tune_seeds(Y2, XSet1); eval21 = eval_tau(Y2, XSet1, tune21$tau)
tune12 = tune_seeds(Y1, XSet2); eval12 = eval_tau(Y1, XSet2, tune12$tau)
tune22 = tune_seeds(Y2, XSet2); eval22 = eval_tau(Y2, XSet2, tune22$tau)

eval22

results = rbind(
  Hep_Set1 = c(tau = tune11$tau, eval11["mean", ]),
  Cu_Set1  = c(tau = tune21$tau, eval21["mean", ]),
  Hep_Set2 = c(tau = tune12$tau, eval12["mean", ]),
  Cu_Set2  = c(tau = tune22$tau, eval22["mean", ]))
round(results, 3)


tau_tab = rbind(Hep_Set1 = tau_spread(tune11), Cu_Set1 = tau_spread(tune21),
                Hep_Set2 = tau_spread(tune12), Cu_Set2 = tau_spread(tune22))
round(tau_tab, 3)
round(c(sd(tune11$topt), sd(tune12$topt), sd(tune21$topt), sd(tune22$topt)), 3)

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

# PCA with only Xset2
Xp = XSet2
dim(Xp)                          # 276 x 12
colnames(Xp)

stopifnot(all(apply(Xp, 2, sd) > 0))   # scale() breaks on a constant column

# pca on everything in set 2
mypca(Xp)
pca = mypca(Xp)

round(cbind(Eigenvalue = pca$values,
            Prop = pca$values / sum(pca$values),
            Flam = pca$Flam), 4)
elbow = 4

data.frame(Method = c('Kaiser', 'Elbow', 'Var Threshold'),
           PCs = c(pca$kaiser, elbow, pca$vthresh))
# use Var Thresh, keeps most features

dr = pca$vthresh
pc = mypca(Xp, dr = dr, plot = FALSE)

round(cov(pc$Xnew), 4)           # kept eigenvalues
round(pc$Lam, 3)                 # loadings

XPC = pc$Xnew
# Fit both models
Model1PC = LogisticRegressionFull(Y1, XPC)
Model2PC = LogisticRegressionFull(Y2, XPC)

cbind(Hep_Set2 = Model1Set2$metrics, Hep_PCA = Model1PC$metrics,
            Cu_Set2  = Model2Set2$metrics, Cu_PCA  = Model2PC$metrics)

tunePC1 = tune_seeds(Y1, XPC); evalPC1 = eval_tau(Y1, XPC, tunePC1$tau)
tunePC2 = tune_seeds(Y2, XPC); evalPC2 = eval_tau(Y2, XPC, tunePC2$tau)

resultsPC = rbind(
  H_PCA  = c(tau = tunePC1$tau, evalPC1['mean', ]),
  UC_PCA = c(tau = tunePC2$tau, evalPC2['mean', ]))

rbind(results[c('Hep_Set2', 'Cu_Set2'), ], resultsPC)

#VIF FOR Xset2
vif1=myvif(XSet1)
round(sort(vif1,decreasing=TRUE),2)
vif2 = myvif(XSet2)
round(sort(vif2, decreasing = TRUE), 2)

# drop the worst predictor one at a time until all VIF <= 5
vd = myvif_drop(XSet2, cutoff = 5)
vd$kept
# all kept under vif <= 5, Multicolinearity assumption satisfied?


#ALL POSSIBLE MODELS 
XS = XSet2                                 
#renaming alk_phos for cleaner look in table printout 
colnames(XS) = gsub('Alk_Phos', 'AlkPhos', colnames(XS))

#swapping edema to a factor 
# keeping all predictors except the two edema dummies
keep = setdiff(colnames(XS), c('EdemaS', 'EdemaY'))    
DF = data.frame(XS[, keep], Edema = factor(clean$Edema), Hep = Y1, Cu = Y2)

#Sanity check to see if i broke it 
levels(DF$Edema)       
# N S Y, N is the reference

vars = c(keep, 'Edema')
# 11
length(vars)                                           
2^length(vars) - 1                                     
# 2047
head(DF)
# nothing broke 

t0 = Sys.time()
tb1 = mymodeltable(DF, 'Hep', vars, type = 'binomial')
tb2 = mymodeltable(DF, 'Cu',  vars, type = 'binomial')
Sys.time() - t0

#must be 2047
nrow(tb1$table)                                       
# S and Y must both appear in the training rows
table(DF$Edema[tb1$split$itrain])                      
tb1$best['AIC']; tb1$best['BIC']
tb2$best['AIC']; tb2$best['BIC']

#table to make top 5
top5 = function(tb, by) head(tb$table[order(tb$table[[by]]), ], 5)
top5(tb1, 'AIC'); top5(tb1, 'BIC')
top5(tb2, 'AIC'); top5(tb2, 'BIC')

# models tied with the best (within 2 AIC)
sum(tb1$table$AIC - min(tb1$table$AIC) <= 2)           
sum(tb2$table$AIC - min(tb2$table$AIC) <= 2)
# 14 and 8 :( soo many ties bro

#now building out the tables for latex
pre_table = function(tb, by = 'AIC', n = 5){
  idx = order(tb$table[[by]])[1:n]
  t = tb$table[idx, ]
  data.frame(Rank = 1:n, Model = t$Model, AIC = t$AIC, BIC = t$BIC,
             RMSE = t$RMSE, McF = t[['McF R2']], row.names = NULL)
}
# first models
pre1A = pre_table(tb1, 'AIC')
pre1B = pre_table(tb1, 'BIC')
#second models
pre2A = pre_table(tb2, 'AIC')
pre2B = pre_table(tb2, 'BIC')
#printing out
pre1A; pre1B; pre2A; pre2B

#building out tables post tau tuning 
#20 total tunings, might take a bit to run

post_table = function(tb, Y, by = 'AIC', n = 5, metric = 'acc'){
  idx = order(tb$table[[by]])[1:n]
  tuned = t(sapply(idx, function(i){
    m = tb$models[[i]]
    X = XS[, names(coef(m))[-1], drop = FALSE]            # predictors of that model (Edema = both dummies)
    tn = tune_seeds(Y, X); ev = eval_tau(Y, X, tn$tau)
    c(tau = tn$tau, ev['mean', ])
  }))
  res = data.frame(Rank = 1:n, Model = tb$table$Model[idx],
                   AIC = tb$table$AIC[idx], BIC = tb$table$BIC[idx],
                   data.frame(tuned, check.names = FALSE), row.names = NULL, check.names = FALSE)
  res$NewRank = rank(-res[[metric]], ties.method = 'min')  # 1 = best tuned accuracy
  res
}

#first models post tune
post1A = post_table(tb1, Y1, 'AIC')
post1B = post_table(tb1, Y1, 'BIC')
#second models post tune
post2A = post_table(tb2, Y2, 'AIC')
post2B = post_table(tb2, Y2, 'BIC')
#print out + checks
colnames(post1A)
post1A; post1B; post2A; post2B

# this was shakey, a lot of the table and weird functions were made by claude
# it was really awkward to look at all the ties so i just decided to let claude 
# do the weird stuff. I made sure to check over it and none of the numbers look 
# like claude fudged them, then i just pasted the output to have a table made.
# feel free to change any of this section at ur discretion. 
# i just have to finish the interpretation of it.



