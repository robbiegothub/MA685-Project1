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


-Model1Set1$model$value
-Model1Set2$model$value

-Model2Set1$model$value
-Model2Set2$model$value

null_logreg(Y1)$lhat                      
null_logreg(Y2)$lhat                      


#Using Optim
myoptim_binom(XSet1,Y1)
myoptim_binom(XSet2,Y1)

myoptim_binom(XSet1,Y2)
myoptim_binom(XSet2,Y2)

myoptim_binom(XSet2[, 0], Y1)$loglik      
myoptim_binom(XSet2[, 0], Y2)$loglik  

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

#VIF for Multicolinearirty
vif2 = myvif(XSet2)
round(sort(vif2, decreasing = TRUE), 2)

# drop the worst predictor one at a time until all VIF <= 5
vd = myvif_drop(XSet2, cutoff = 5)
vd$kept
# all kept under vif <= 5, Multi-colinearity assumption satisfied?

# PCA with only lab predictors from Xset2
labs = c('Bilirubin', 'Albumin', 'Alk_Phos', 'SGOT', 'Platelets', 'Cholesterol', 'Tryglicerides')
Xp = XSet2[, labs]
dim(Xp)                          # 276 x 7
colnames(Xp)

stopifnot(all(apply(Xp, 2, sd) > 0))   # scale() breaks on a constant column

# pca on the lab values
mypca(Xp)
pca = mypca(Xp)

# figures for paper
ev = pca$values
vt = 0.85                                  # variance threshold used by mypca
pdf('pca_scree_cum.pdf', width = 10, height = 4.5)
par(mfrow = c(1, 2))
plot(ev, type = 'b', pch = 16, xlab = 'Component', ylab = 'Eigenvalue', main = 'Scree Plot')
abline(h = 1, lty = 2)
plot(cumsum(ev) / sum(ev), type = 'b', pch = 16, ylim = c(0, 1),
     xlab = 'Number of components', ylab = 'Cumulative proportion of variance',
     main = 'Cumulative Variance')
abline(h = vt, lty = 2)
dev.off()
par(mfrow = c(1, 1))
# end of figures for PCA

round(cbind(Eigenvalue = pca$values,
            Prop = pca$values / sum(pca$values),
            Flam = pca$Flam), 4)
elbow = 3                        # placeholder, reset after viewing the new scree plot

data.frame(Method = c('Kaiser', 'Elbow', 'Var Threshold'),
           PCs = c(pca$kaiser, elbow, pca$vthresh))
# use Var Thresh, keeps most of the variance

dr = pca$vthresh
pc = mypca(Xp, dr = dr, plot = FALSE)

round(cov(pc$Xnew), 4)           # kept eigenvalues
round(pc$Lam, 3)                 # loadings

# non-lab predictors (Set 1) stay as they are; only the labs are replaced by components
XPC = cbind(XSet1, pc$Xnew)
colnames(XPC)

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

###############

#ASSUMPTIONS TESTING
nums = c('Age', 'Bilirubin', 'Albumin', 'Alk_Phos', 'SGOT', 'Platelets', 'Cholesterol', 'Tryglicerides')
stopifnot(all(nums %in% names(clean)), all(nums %in% colnames(XSet2)))

#Assumption 2
#claude assisted for the squared term test along with how to make the plots
#empirical log-odds in 6 equal-count bins 
elogit = function(x, y, bins = 6){
  g = cut(x, unique(quantile(x, 0:bins / bins)), include.lowest = TRUE)
  n1 = tapply(y, g, sum); n0 = tapply(1 - y, g, sum)
  list(x = as.numeric(tapply(x, g, mean)), el = log((n1 + 0.5) / (n0 + 0.5)))
}
xval = function(v){ x = clean[[v]]; if(v == 'Age' && max(x) > 200) x / 365.25 else x }

lin_panels = function(Y){
  par(mfrow = c(2, 4), mar = c(4, 4, 2.5, 1))
  for(v in nums){
    e = elogit(xval(v), Y)
    plot(e$x, e$el, pch = 16, xlab = ifelse(v == 'Age', 'Age (years)', v),
         ylab = 'Empirical log-odds', main = paste0('r = ', round(cor(e$x, e$el), 2)))
    abline(lm(e$el ~ e$x), lty = 2)
  }
}
par(mfrow = c(1, 1))
#Assumption 2 fig 1
lin_panels(Y1)
#Assumtption 2 fig 2
lin_panels(Y2)   
#pdf('linearity_hep.pdf', width = 10, height = 5.5); lin_panels(Y1); dev.off()
#pdf('linearity_cu.pdf',  width = 10, height = 5.5); lin_panels(Y2); dev.off()

#some strong predictors but also quite a few lab tests that seem right skewed
#these are the models and which predictors no clear trends
#For Y1 Alk, Age, Chol and triglycerides
#For Y2 Alk, Age, and Platelets

# looking at the predictors to try a squared term (with all other predictors held constant)
# testing and printing pvals with a bonferroni correction at end 
lin_table = function(Y){
  d = data.frame(XSet2)
  d[nums] = scale(as.matrix(d[nums]))
  d$yy = as.numeric(Y)
  m0 = glm(yy ~ ., data = d, family = binomial)
  sapply(nums, function(v){
    e = elogit(xval(v), Y)
    m1 = update(m0, as.formula(paste('. ~ . + I(', v, '^2)')))
    c(r = cor(e$x, e$el), p_quad = anova(m0, m1, test = 'Chisq')[2, 'Pr(>Chi)'])
  })
}
# ASSUMPTION 2 TABLE 2
lin1 = lin_table(Y1); lin2 = lin_table(Y2)
print('lin_table_1')
round(t(lin1), 3) 
print('lin_table_2')
round(t(lin2), 3)

p_all = c(lin1['p_quad', ], lin2['p_quad', ])
names(p_all) = c(paste0('Hep_', nums), paste0('Cu_', nums))
round(sort(p.adjust(p_all, method = 'bonferroni'))[1:4], 3) 
#with bonferoni correction only the Bilirubin on the copper test is significant

#log transfrom the Bilirubin to see if truly requires non linear transformation
d_raw = data.frame(XSet2)
aic_cmp = function(Y, skew){
  d_log = d_raw; d_log[skew] = log(d_log[skew])
  c(raw = AIC(glm(Y ~ ., data = d_raw, family = binomial)),
    log = AIC(glm(Y ~ ., data = d_log, family = binomial)))
}

bil  = c('Bilirubin')
five = c('Bilirubin', 'Alk_Phos', 'SGOT', 'Cholesterol', 'Tryglicerides')
round(rbind(Hep_bil = aic_cmp(Y1, bil), Cu_bil = aic_cmp(Y2, bil),
            Hep_five = aic_cmp(Y1, five), Cu_five = aic_cmp(Y2, five)), 1)

XL = XSet2
XL[, c('Bilirubin')] = log(XL[, c('Bilirubin')])
tuneL1 = tune_seeds(Y1, XL); evalL1 = eval_tau(Y1, XL, tuneL1$tau)
tuneL2 = tune_seeds(Y2, XL); evalL2 = eval_tau(Y2, XL, tuneL2$tau)
# ASSUMPTION 2 TABLE 3 
round(rbind(Hep_Set2 = results['Hep_Set2', ], Hep_log = c(tau = tuneL1$tau, evalL1['mean', ]),
            Cu_Set2  = results['Cu_Set2', ],  Cu_log  = c(tau = tuneL2$tau, evalL2['mean', ])), 3)






#BOOTSTRAP CI FOR GINI IMPURITY
bs_g1 = mybootstrap(Y1, impure, nboot = 10000, seed = 123)   # hepatomegaly
bs_g2 = mybootstrap(Y2, impure, nboot = 10000, seed = 123)   # urine copper
round(rbind(Hep = c(bs_g1$est, bs_g1$ci), Cu = c(bs_g2$est, bs_g2$ci)), 4)