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

#VIF
round(sort(myvif(XSet3), decreasing=TRUE), 2)
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


#ASSUMPTION 2: linearity of log-odds, stage k vs stage 1
library(nnet)
nums3 = c('Age', 'Bilirubin', 'Albumin', 'Alk_Phos', 'SGOT', 'Platelets', 'Cholesterol', 'Tryglicerides', 'Copper')
stopifnot(all(nums3 %in% colnames(XSet3)), all(nums3 %in% names(clean)))

xval = function(v){ x = clean[[v]]; if(v == 'Age' && max(x) > 200) x / 365.25 else x }

#empirical log-odds of stage k vs stage 1 in 3 equal-count groups (only 12 patients in stage 1)
elogit_k = function(x, k, bins = 3){
  g = cut(x, unique(quantile(x, 0:bins / bins)), include.lowest = TRUE)
  nk = tapply(y1 == k, g, sum); n1 = tapply(y1 == 1, g, sum)
  list(x = as.numeric(tapply(x, g, mean)), el = log((nk + 0.5) / (n1 + 0.5)))
}

stage_panels = function(){
  par(mfrow = c(2, 5), mar = c(4, 4, 2, 1))
  for(v in nums3){
    e = lapply(2:4, function(k) elogit_k(xval(v), k))
    plot(range(unlist(lapply(e, function(a) a$x))), range(sapply(e, function(a) a$el)), type = 'n',
         xlab = ifelse(v == 'Age', 'Age (years)', v), ylab = 'Log-odds vs stage 1', main = v)
    for(i in 1:3) lines(e[[i]]$x, e[[i]]$el, type = 'b', pch = c(16, 17, 15)[i], lty = i)
  }
  plot.new()
  legend('center', legend = c('Stage 2', 'Stage 3', 'Stage 4'), pch = c(16, 17, 15), lty = 1:3, bty = 'n')
}
stage_panels()
pdf('linearity_stage.pdf', width = 10, height = 4.5); stage_panels(); dev.off()
par(mfrow = c(1, 1))

#claude assisted for the squared curvature testing, wanted a sort of quantitative method
#squared-term test: x^2 enters all three stage-vs-1 logits (3 df)
d3 = data.frame(XSet3)
d3[nums3] = scale(as.matrix(d3[nums3]))
d3$yy = y1
m0 = multinom(yy ~ ., data = d3, trace = FALSE, maxit = 500)
stage_quad = sapply(nums3, function(v){
  m1 = update(m0, as.formula(paste('. ~ . + I(', v, '^2)')))
  anova(m0, m1)[2, 'Pr(Chi)']
})
round(stage_quad, 3)
round(p.adjust(stage_quad, method = 'bonferroni'), 3)


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







#BOOTSTRAP CI FOR GINI IMPURITY
bs_g = mybootstrap(y1, impure, nboot = 10000, seed = 123)    # stage
round(c(bs_g$est, bs_g$ci), 4)
