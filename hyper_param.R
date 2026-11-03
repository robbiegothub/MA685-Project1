source('ma685_functions.R')

#claude wrote alot of the actual code for these functions, based on my initial functions, to optimize a bit.
tune_seeds = function(Y, X, nseed = 500, ptrain = 0.65){
  topt = c(); accopt = c()
  for(s in 1:nseed){
    sp = mysplit(length(Y), ptrain, seed = s)
    fit = log_reg_general(Y[sp$itrain], X[sp$itrain, ])
    phat = multinom_probs(fit$bhats, X[sp$itest, ])[, 2]   
    tu = mytunetau(phat, Y[sp$itest], plot = FALSE)
    topt = c(topt, tu$tmax); accopt = c(accopt, tu$accmax)
  }
  return(list(topt = topt, accopt = accopt, tau = mean(topt)))
}

#this isn't taught, but theres the issue of optimism, the best tau for a specific split could just be noise; the best tau over all splits could be noise, meaning we should evaluate 
#the performance of tau on new splits, not the split that generated it.
eval_tau = function(Y, X, tau, seeds = 1001:1500, ptrain = 0.65){
  res = sapply(seeds, function(s){
    sp = mysplit(length(Y), ptrain, seed = s)
    fit = log_reg_general(Y[sp$itrain], X[sp$itrain, ])
    phat = multinom_probs(fit$bhats, X[sp$itest, ])[, 2]
    yhat = myclass_binom(phat, tau); y = Y[sp$itest]
    c(acc  = mean(yhat == y),
      sens = mean(yhat[y == 1] == 1),
      spec = mean(yhat[y == 0] == 0))
  })
  return(rbind(mean = rowMeans(res), sd = apply(res, 1, sd)))
}



#claude geenrated
eval_multi = function(Y, X, seeds = 1:500, ptrain = 0.65){
  Y = factor(Y); lv = levels(Y)
  res = sapply(seeds, function(s){
    sp = mysplit(length(Y), ptrain, seed = s)
    fit = log_reg_general(Y[sp$itrain], X[sp$itrain, ])
    P = multinom_probs(fit$bhats, X[sp$itest, ])
    yhat = lv[apply(P, 1, which.max)]; y = as.character(Y[sp$itest])
    c(acc = mean(yhat == y),
      sapply(lv, function(k) mean(yhat[y == k] == k)))   # per-stage recall
  })
  rownames(res) = c('acc', paste0('recall_', lv))
  rbind(mean = rowMeans(res, na.rm = TRUE), sd = apply(res, 1, sd, na.rm = TRUE))
}

#claude generated.
tune_decay = function(Y, X, lams = c(0, 10^seq(-3, 1, 0.5)), seeds = 1:500, ptrain = 0.65){
  Y = factor(Y); lv = levels(Y)
  acc = sapply(seeds, function(s){
    sp = mysplit(length(Y), ptrain, seed = s)
    mu = colMeans(X[sp$itrain, ]); sdv = apply(X[sp$itrain, ], 2, sd)
    Xtr = scale(X[sp$itrain, ], mu, sdv); Xte = scale(X[sp$itest, ], mu, sdv)
    sapply(lams, function(l){
      fit = log_reg_general(Y[sp$itrain], Xtr, decay = l)
      P = multinom_probs(fit$bhats, Xte)
      mean(lv[apply(P, 1, which.max)] == as.character(Y[sp$itest]))
    })
  })
  m = rowMeans(acc)
  list(lam = lams[which.max(m)], acc = m, lams = lams)
}