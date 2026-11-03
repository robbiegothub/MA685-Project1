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



#claude generated. standardize X using training-split means/sds only
std_split = function(X, sp){
  mu = colMeans(X[sp$itrain, ]); sdv = apply(X[sp$itrain, ], 2, sd)
  sdv[sdv == 0] = 1   # rare dummy (e.g. EdemaY) can be constant in a split
  list(tr = scale(X[sp$itrain, ], mu, sdv), te = scale(X[sp$itest, ], mu, sdv))
}

#claude generated. smoothed (noisy) bootstrap oversampling, applied to the TRAINING split only.
# each class is resampled up to `target` rows; continuous columns get N(0, (h*sd_k)^2) noise,
# sd_k = that column's sd within class k (standardized scale). binary dummies are copied as-is.
noisy_boot = function(Xtr, Ytr, bin, h = 0.5, target = max(table(Ytr)), seed = 1){
  set.seed(seed)
  Ytr = factor(Ytr); lv = levels(Ytr); cont = !bin
  Xn = list(Xtr); Yn = as.character(Ytr)
  for(k in lv){
    idx = which(Ytr == k); m = target - length(idx)
    if(m <= 0 || length(idx) == 0) next
    Xb = Xtr[idx[sample.int(length(idx), m, replace = TRUE)], , drop = FALSE]
    sdk = apply(Xtr[idx, , drop = FALSE], 2, sd); sdk[is.na(sdk)] = 0
    Xb[, cont] = Xb[, cont] + matrix(rnorm(m*sum(cont)), m) * rep(h*sdk[cont], each = m)
    Xn[[length(Xn) + 1]] = Xb; Yn = c(Yn, rep(k, m))
  }
  list(X = do.call(rbind, Xn), Y = factor(Yn, levels = lv))
}

# fit on a (possibly augmented) training split, return test probabilities
fit_pred = function(Y, X, sp, decay, aug, h, bin, seed){
  Z = std_split(X, sp); Ytr = Y[sp$itrain]; Xtr = Z$tr
  if(aug){ B = noisy_boot(Xtr, Ytr, bin, h = h, seed = seed); Xtr = B$X; Ytr = B$Y }
  multinom_probs(log_reg_general(Ytr, Xtr, decay = decay)$bhats, Z$te)
}

#claude geenrated
eval_multi = function(Y, X, decay = 0, seeds = 1001:1500, ptrain = 0.65, aug = FALSE, h = 0.5){
  Y = factor(Y); lv = levels(Y); bin = apply(X, 2, function(c) all(c %in% c(0, 1)))
  res = sapply(seeds, function(s){
    sp = mysplit(length(Y), ptrain, seed = s)
    P = fit_pred(Y, X, sp, decay, aug, h, bin, seed = s)
    yhat = lv[apply(P, 1, which.max)]; y = as.character(Y[sp$itest])
    rec = sapply(lv, function(k) mean(yhat[y == k] == k))
    c(acc = mean(yhat == y), bacc = mean(rec, na.rm = TRUE), rec)
  })
  rownames(res) = c('acc', 'bacc', paste0('recall_', lv))
  rbind(mean = rowMeans(res, na.rm = TRUE), sd = apply(res, 1, sd, na.rm = TRUE))
}
#claude generated.
tune_decay = function(Y, X, lams = c(0, 10^seq(-3, 1, 0.5)), seeds = 1:500, ptrain = 0.65,
                      aug = FALSE, h = 0.5){
  Y = factor(Y); lv = levels(Y); bin = apply(X, 2, function(c) all(c %in% c(0, 1)))
  res = sapply(seeds, function(s){
    sp = mysplit(length(Y), ptrain, seed = s)
    y = as.character(Y[sp$itest])
    sapply(lams, function(l){
      P = fit_pred(Y, X, sp, l, aug, h, bin, seed = s)
      yhat = lv[apply(P, 1, which.max)]
      mean(sapply(lv, function(k) mean(yhat[y == k] == k)), na.rm = TRUE)  # balanced acc
    })
  })
  m = rowMeans(res)
  plot(log10(lams + 1e-4), m, type = 'b', pch = 16,
       xlab = 'log10(lambda)', ylab = 'Mean balanced accuracy')
  list(lam = lams[which.max(m)], bacc = m, lams = lams)
}