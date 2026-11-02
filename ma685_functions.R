# This is where all my functions relating to the course will be built
# Use:  source('ma685_functions.R')
# Only base R + nnet (+ CVXR for the lasso/LAD functions, as in Lecture 16/17)

library(nnet)


# ======================================================================
# PART 1: Lect 1-10
# ======================================================================

# Function 1 
# given sample s, calculate Variance matrix (Sigma_hatx), mean vector(X_bar) 
mysummary = function(X){
  X = as.matrix(X)
  xbar = colMeans(X)
  sighat_x = cov(X)           # divides by n-1
  return(list(xbar = xbar, sighat_x = sighat_x))
}



# Function 2
# given d dimensional sample, S, of size n, return Mahalanobis quantiles 
# w.r.t Sigma_hatx


# function 3 
# mymaha(X) which calculates non squared mahalanobis quantiles D_m(xk, xbar) w.r.t sighat_x
# X is d dimensional, of size n. sighat_x is the cov(X)
mymaha = function(X){
  X = as.matrix(X)
  s = mysummary(X)
  D2 = mahalanobis(X, center = s$xbar, cov = s$sighat_x)   # squared distances
  return(sqrt(D2))
}

# m2qq_ideal = (1-1/n) * qchisq((1:n)/(n+1), d, lower.tail = T)
# compare sort(function 2) to m2qq and make plot 
# Function 2a
# given d dimensional sample X of size n, compare sorted NON-squared Mahalanobis
# distances to the square roots of the ideal chi-square quantiles
mymahaqq = function(X){
  X = as.matrix(X)
  n = nrow(X); d = ncol(X)
  m = sort(mymaha(X))
  m2qq_ideal = (1 - 1/n) * qchisq((1:n)/(n+1), df = d, lower.tail = TRUE)
  plot(sqrt(m2qq_ideal), m,
       xlab = "Ideal quantiles (sqrt)", ylab = "Sorted Mahalanobis distances",
       main = "Non-squared QQ plot")
  abline(0, 1, col = "red")
  invisible(cbind(ideal = sqrt(m2qq_ideal), sample = m))
}

# Function 2b
# same comparison for the SQUARED Mahalanobis distances vs chi-square quantiles
mymaha2qq = function(X){
  X = as.matrix(X)
  n = nrow(X); d = ncol(X)
  m2 = sort(mymaha(X)^2)
  m2qq_ideal = (1 - 1/n) * qchisq((1:n)/(n+1), df = d, lower.tail = TRUE)
  plot(m2qq_ideal, m2,
       xlab = "Ideal chi-square quantiles", 
       ylab = "Sorted squared Mahalanobis distances",
       main = "Squared QQ plot")
  abline(0, 1, col = "red")
  invisible(cbind(ideal = m2qq_ideal, sample = m2))
}


# function 4 Whiten Data

whiten = function(X){
  X = as.matrix(X)
  sx = apply(X, 2, sd)
  Sig_hatx = cov(X)
  xbar = colMeans(X)
  Xc = sweep(X, 2, xbar, '-')
  Wx = Xc %*% solve(chol(Sig_hatx))
  return(Wx)
}

# function 5 Projection-Pursuit KLD
rsphere = function(m,d){
  V = matrix(rnorm(m*d), m, d)
  return(V/sqrt(rowSums(V^2)))
}

#function to take in Wx and then produce the projected data as well as 
#summing the KL divergence

myppkld = function(Wx, m = 10000){
  Wx = as.matrix(Wx)
  n = nrow(Wx); d = ncol(Wx)
  pvec = rsphere(m,d)
  pWx = Wx %*% t(pvec)
  pWx = apply(pWx, 2, sort)
  
  Fx = (1:n)/(n+1)
  qx_ideal = dnorm(qnorm(Fx))
  KL_values = c()
  
  for (k in 1:m){
    px_approx = 0*qx_ideal
    px_approx[1:(n-1)] = (1/(n+1))/(pWx[2:n, k] - pWx[1:(n-1),k])
    px_approx[n] = (1/(n+1))/(pWx[n,k] - pWx[n-1, k])
    KL_values = c(KL_values, mean(log(px_approx/qx_ideal)))
  }
  invisible(list(V = pvec, P = pWx, kl = KL_values, kl_sums = sum(KL_values)))
}



# ======================================================================
# PART 2: MULTIVARIATE LEAST SQUARES (Lecture 10)
# ======================================================================

# Function 6
# given predictors X (n x p, NO column of 1s) and responses Y (n x d), fit the
# multivariate least-squares model Y = X*beta + eps (the 1s are added for me)
# returns bhat, hat matrix H, the three sum-of-squares matrices (T = E + M),
# the estimated covariance sigh = ssige/(n-p-1), and the variance of vec(bhat)
rjv_mmlr = function(X, Y){
  Xo = as.matrix(X); Y = as.matrix(Y)
  X = cbind(1, Xo)
  n = nrow(X); p = ncol(X) - 1; d = ncol(Y)
  
  XtXi = solve(t(X) %*% X)
  bhat = XtXi %*% t(X) %*% Y
  rownames(bhat) = c('Int', if(is.null(colnames(Xo))) paste0('X', 1:p) else colnames(Xo))
  
  H = X %*% XtXi %*% t(X)
  ones = matrix(1, n, n)
  ssigt = t(Y) %*% (diag(n) - (1/n)*ones) %*% Y     # total
  ssige = t(Y) %*% (diag(n) - H) %*% Y              # error
  ssigm = t(Y) %*% (H - (1/n)*ones) %*% Y           # model
  
  sigh = ssige/(n - p - 1)
  varbhat = kronecker(sigh, XtXi)                   # V(vec(bhat)) = Sigma (x) (X'X)^-1
  
  return(list(bhat = bhat, H = H, n = n, p = p, d = d,
              ssigt = ssigt, ssige = ssige, ssigm = ssigm,
              sigh = sigh, varbhat = varbhat))
}

# Function 7
# given the output of rjv_mmlr(), compute the six test statistics for
# H0: all slopes = 0 (Roy, Hotelling-Lawley, Pillai, Wilks' chi-square,
# Wilks' F via Rao's t and nu2, and the NM statistic with constant c) 
# and the p-values for Wilks (chi-square and F)
rjv_manova = function(fit, c = 0.03){
  n = fit$n; p = fit$p; d = fit$d; g = p + 1
  E = fit$ssige; M = fit$ssigm; sst = fit$ssigt
  
  Ro = max(Re(eigen(solve(E) %*% M)$values))
  HL = sum(diag(solve(E) %*% M))
  Pi = sum(diag(solve(sst) %*% M))
  
  W = (n - 1 - (d + g)/2) * log(det(sst)/det(E))
  Lambda = det(E)/det(sst)
  den = d^2 + (g-1)^2 - 5
  tt = if(den > 0) sqrt((d^2*(g-1)^2 - 4)/den) else 1     # Rao's t
  nu2 = (n - 1 - (d + g)/2)*tt - (d*(g-1) - 2)/2
  Wm = (Lambda^(-1/tt) - 1) * nu2/(d*(g-1))
  NM = det(M)/det(M + c*E)
  
  pchi = pchisq(W, d*(g-1), lower.tail = FALSE)
  pF = pf(Wm, d*(g-1), nu2, lower.tail = FALSE)
  
  return(list(stats = c(Ro = Ro, HL = HL, Pi = Pi, W = W, Wm = Wm, NM = NM),
              tt = tt, nu2 = nu2, p_chisq = pchi, p_F = pF))
}



# ======================================================================
# PART 3: LOGISTIC / MULTINOMIAL REGRESSION (Lectures 11, 12)
# ======================================================================

# Function 8
# the logistic function (inverse logit) and the row-wise softmax
# myinvlogit takes any logits L -> probabilities
# mysoftmax takes an n x d matrix of eta's -> n x d probabilities (rows sum to 1)
#   (subtracts the row max first so exp() doesn't overflow)
myinvlogit = function(L){
  return(1/(1 + exp(-L)))
}

mysoftmax = function(eta){
  m = apply(eta, 1, max)
  ex = exp(eta - m)
  return(ex/rowSums(ex))
}

# Function 9
# given a response (any categorical vector), return the n x g matrix of 0/1
# indicators (this is what the training "Ptr" was in Lecture 12)
# levels = order of the columns (first one is the reference group)
mydummy = function(Y, levels = NULL){
  if(is.null(levels)) levels = levels(factor(Y))
  Yf = factor(Y, levels = levels)
  P = sapply(levels, function(f) ifelse(Yf == f, 1, 0))
  return(P)
}

# Function 10
# negative log-likelihood (and its gradient) for BINOMIAL logistic regression
# X here INCLUDES the column of 1s, Y is 0/1. These are what optim() minimizes.
# -l(b) = -sum[ y*eta - ln(1+e^eta) ]
mynll_binom = function(bhat, X, Y){
  eta = as.vector(X %*% bhat)
  return(-sum(Y*eta - (pmax(eta, 0) + log1p(exp(-abs(eta))))))
}

mygr_binom = function(bhat, X, Y){
  phat = myinvlogit(as.vector(X %*% bhat))
  return(-as.vector(t(X) %*% (Y - phat)))
}

# Function 11
# negative log-likelihood (cross-entropy) and gradient for MULTINOMIAL logistic regression
# bvec = vectorized (p+1) x (d-1) coefficient matrix (reference = group 1)
# X includes the 1s, P is the n x d indicator matrix from mydummy()
mynll_multi = function(bvec, X, P){
  d = ncol(P)
  B = matrix(bvec, ncol = d - 1)
  eta = cbind(0, X %*% B)
  m = apply(eta, 1, max)
  lse = m + log(rowSums(exp(eta - m)))        # log of the softmax denominator
  return(-sum(P*(eta - lse)))
}

mygr_multi = function(bvec, X, P){
  d = ncol(P)
  B = matrix(bvec, ncol = d - 1)
  phat = mysoftmax(cbind(0, X %*% B))
  return(-as.vector(t(X) %*% (P[, -1, drop = FALSE] - phat[, -1, drop = FALSE])))
}

# Function 12
# fit binomial logistic regression NUMERICALLY with optim(BFGS) and compare to glm
# X = predictors (no 1s), Y = 0/1. Returns my estimate, glm's estimate, the
# absolute difference, and the maximized log-likelihood.
# (tip: standardize X first if BFGS is slow to converge)
myoptim_binom = function(X, Y, bhat0 = NULL){
  X = cbind(1, as.matrix(X))
  if(is.null(bhat0)) bhat0 = rep(0, ncol(X))
  So = optim(bhat0, mynll_binom, mygr_binom, X = X, Y = Y, method = 'BFGS',
             control = list(reltol = 1e-12, maxit = 1000))
  bglm = glm.fit(X, Y, family = binomial())$coefficients
  return(list(bhat_optim = So$par, bhat_glm = unname(bglm),
              absdiff = sum(abs(So$par - bglm)),
              loglik = -So$value, converged = (So$convergence == 0)))
}

# Function 13
# same thing for MULTINOMIAL: optim(BFGS) vs multinom() from nnet
# Y = categorical (first level = reference); returns (p+1) x (d-1) matrices
myoptim_multi = function(X, Y, maxit = 2000){
  X = as.matrix(X)
  Yf = factor(Y)
  P = mydummy(Yf)
  Xd = cbind(1, X)
  d = ncol(P)
  b0 = rep(0, ncol(Xd)*(d - 1))
  So = optim(b0, mynll_multi, mygr_multi, X = Xd, P = P, method = 'BFGS',
             control = list(reltol = 1e-12, maxit = maxit))
  Bopt = matrix(So$par, ncol = d - 1)
  M = multinom(Yf ~ X, trace = FALSE, maxit = 500)
  Bmn = t(coef(M))
  return(list(B_optim = Bopt, B_multinom = Bmn,
              absdiff = sum(abs(Bopt - Bmn)), loglik = -So$value,
              converged = (So$convergence == 0)))
}

# Function 14
# predicted probabilities for NEW data
# binomial: bhat = coef(M) (vector), returns vector of P(Y=1)
# multinomial: bhat = coef(M) ((d-1) x (p+1) matrix); returns n x d matrix of
# probabilities with columns named by FF (default: 'Ref' + the other levels)
mypred_binom = function(bhat, Xnew){
  Xnew = cbind(1, as.matrix(Xnew))
  return(myinvlogit(as.vector(Xnew %*% bhat)))
}

mypred_multi = function(bhat, Xnew, FF = NULL){
  Xnew = cbind(1, as.matrix(Xnew))
  eta = cbind(0, Xnew %*% t(bhat))
  phat = mysoftmax(eta)
  if(is.null(FF)) FF = c('Ref', rownames(bhat))
  colnames(phat) = FF
  return(phat)
}

# Function 15
# model performance metrics from predicted probabilities
# binomial: Y is 0/1, phat a vector, k = number of estimated parameters (p+1)
# multinomial: Y categorical, phat an n x d matrix (named columns), k = (p+1)*(d-1)
# returns RMSE, log-likelihood, AIC = 2k - 2ll, BIC = k ln(n) - 2ll, and
# McFadden R^2 = 1 - ll/ll_null  (null = intercept only = overall proportions)
mymetrics_binom = function(Y, phat, k){
  n = length(Y)
  ll = sum(Y*log(phat) + (1 - Y)*log(1 - phat))
  p0 = mean(Y)
  ll0 = sum(Y*log(p0) + (1 - Y)*log(1 - p0))
  return(c(RMSE = sqrt(mean((Y - phat)^2)), logLik = ll,
           AIC = 2*k - 2*ll, BIC = k*log(n) - 2*ll, McF = 1 - ll/ll0))
}

mymetrics_multi = function(Y, phat, k){
  n = nrow(phat)
  P = mydummy(Y, levels = colnames(phat))
  ll = sum(P*log(phat))
  p0 = colMeans(P)
  ll0 = n*sum(p0*log(p0))
  return(c(RMSE = sqrt(mean((P - phat)^2)), logLik = ll,
           AIC = 2*k - 2*ll, BIC = k*log(n) - 2*ll, McF = 1 - ll/ll0))
}

# Function 16
# classification + accuracy + confusion matrix
# binary: threshold tau.  multinomial: group with the maximum predicted probability
myclass_binom = function(phat, tau = 0.5){
  return(ifelse(phat >= tau, 1, 0))
}

myclass_multi = function(phat){
  return(colnames(phat)[apply(phat, 1, which.max)])
}

myaccuracy = function(yhat, y){
  return(mean(as.character(yhat) == as.character(y)))
}

myconfusion = function(yhat, y){
  return(table(Predicted = yhat, Actual = y))
}

# Function 17
# random training/testing split (same sampling as in lecture: sample w/o replacement)
# returns the index vectors itrain, itest
mysplit = function(n, ptrain = 0.65, seed = 2026){
  set.seed(seed)
  ntrain = floor(ptrain*n)
  itrain = sample(1:n, ntrain, replace = FALSE)
  itest = setdiff(1:n, itrain)
  return(list(itrain = itrain, itest = itest, ntrain = ntrain, ntest = n - ntrain))
}



# ======================================================================
# PART 4: VARIABLE SELECTION / ALL-POSSIBLE-MODELS / TUNING (Lecture 13)
# ======================================================================

# Function 18
# builds the formula  yname ~ v1 + v2 + ...   (or  yname ~ 1  if v is empty)
myformula = function(yname, v){
  rhs = if(length(v) == 0) '1' else paste(v, collapse = ' + ')
  return(as.formula(paste(yname, '~', rhs)))
}

# Function 19
# ALL-POSSIBLE-MODELS: fit every one of the 2^p - 1 combinations of vars
# df = data frame holding the response (column yname) and the predictors (vars)
# type = 'binomial' (glm, response must be 0/1) or 'multinom' (response a factor, d >= 3)
# For each model: RMSE, logLik, AIC, BIC, McFadden (fit on ALL data)
# If itrain is given (see mysplit), the model is also refit on the training rows
# and the testing Accuracy is reported (tau used for binomial).
# returns the metrics table (one row per model) and the list of fitted models
myallsubsets = function(df, yname, vars, type = 'binomial', itrain = NULL, tau = 0.5){
  MM = list(); rows = list()
  Y = df[[yname]]
  for(k in 1:length(vars)){
    combos = combn(vars, k, simplify = FALSE)
    for(j in combos){
      nm = paste(j, collapse = '_')
      form = myformula(yname, j)
      
      if(type == 'binomial'){
        M = glm(form, data = df, family = 'binomial')
        met = mymetrics_binom(Y, fitted(M), length(coef(M)))
      } else {
        M = multinom(form, data = df, trace = FALSE, maxit = 500)
        met = mymetrics_multi(Y, fitted(M), length(coef(M)))
      }
      
      if(!is.null(itrain)){
        dtr = df[itrain, ]; dte = df[-itrain, ]
        if(type == 'binomial'){
          Mtr = glm(form, data = dtr, family = 'binomial')
          yh = ifelse(predict(Mtr, newdata = dte, type = 'response') >= tau, 1, 0)
        } else {
          Mtr = multinom(form, data = dtr, trace = FALSE, maxit = 500)
          yh = predict(Mtr, newdata = dte)
        }
        met = c(met, Accuracy = myaccuracy(yh, dte[[yname]]))
      }
      MM[[nm]] = M
      rows[[nm]] = met
    }
  }
  tab = do.call(rbind, rows)
  return(list(table = tab, models = MM))
}

# Function 20
# BACKWARD ELIMINATION by p-value (binomial glm, numeric predictors)
# repeatedly drops the predictor with the largest p-value until all p <= alpha
# returns the kept variables and the p-values at every step
mybackward_p = function(df, yname, vars, alpha = 0.05){
  steps = list(); cur = vars
  while(length(cur) > 0){
    M = glm(myformula(yname, cur), data = df, family = 'binomial')
    pv = summary(M)$coefficients[-1, 4]
    names(pv) = cur                      # (assumes numeric predictors)
    steps[[length(steps) + 1]] = list(model = cur, pvalues = round(pv, 4))
    worst = names(which.max(pv))
    if(pv[worst] > alpha){
      cur = setdiff(cur, worst)
    } else break
  }
  return(list(final = cur, steps = steps))
}

# Function 21
# BACKWARD ELIMINATION by AIC (works for binomial AND multinomial, where
# a single p-value per predictor doesn't exist): drop the predictor whose
# removal lowers AIC the most; stop when no removal lowers AIC
mybackward_aic = function(df, yname, vars, type = 'binomial'){
  fitaic = function(v){
    form = myformula(yname, v)
    M = if(type == 'binomial') glm(form, data = df, family = 'binomial') else
      multinom(form, data = df, trace = FALSE, maxit = 500)
    return(AIC(M))
  }
  cur = vars; aic_cur = fitaic(cur)
  steps = list(list(model = cur, AIC = aic_cur))
  while(length(cur) > 0){
    cands = sapply(cur, function(v) fitaic(setdiff(cur, v)))
    best = which.min(cands)
    if(cands[best] < aic_cur){
      cur = setdiff(cur, cur[best]); aic_cur = cands[best]
      steps[[length(steps) + 1]] = list(model = cur, AIC = aic_cur)
    } else break
  }
  return(list(final = cur, steps = steps))
}

# Function 22
# HYPERPARAMETER TUNING of the threshold tau for ONE train/test partition
# phat = predicted probabilities on the test set, Y = true test 0/1
# tries every tau on the grid, returns the best tau, its accuracy, and all accuracies
mytunetau = function(phat, Y, taus = seq(0.01, 0.99, 0.01), plot = TRUE){
  acc = sapply(taus, function(t) mean(ifelse(phat >= t, 1, 0) == Y))
  tmax = taus[which.max(acc)]; accmax = max(acc)
  if(plot){
    plot(taus, acc, type = 'b', pch = 16, cex = 0.5, xlab = 'Threshold Value',
         ylab = 'Accuracy', main = 'Testing Accuracy vs Threshold')
    abline(v = tmax, lty = 2); abline(h = accmax, lty = 2)
  }
  return(list(tmax = tmax, accmax = accmax, acc = acc, taus = taus))
}

# Function 23
# AGGREGATION over seeds: for seeds 1..nseed redo the split, refit the binomial
# model (form = a formula like Y ~ X2 + X3), tune tau, and record the best tau/accuracy
# returns all of them, their mean / sd / quantiles, and plots (scatter + histogram)
# -> shows how much the partition itself changes the chosen threshold
mytauseeds = function(form, df, yname, nseed = 500, ptrain = 0.65,
                      taus = seq(0.01, 0.99, 0.01), plot = TRUE){
  n = nrow(df); ntrain = floor(ptrain*n)
  mytopt = c(); myaccopt = c()
  for(s in 1:nseed){
    set.seed(s)
    itrain = sample(1:n, ntrain, replace = FALSE)
    dtr = df[itrain, ]; dte = df[-itrain, ]
    Mtr = glm(form, data = dtr, family = 'binomial')
    phat = predict(Mtr, newdata = dte, type = 'response')
    tu = mytunetau(phat, dte[[yname]], taus, plot = FALSE)
    mytopt = c(mytopt, tu$tmax); myaccopt = c(myaccopt, tu$accmax)
  }
  if(plot){
    op = par(mfrow = c(1, 2))
    plot(mytopt, myaccopt, pch = 16, cex = 0.6, xlab = 'Optimal Thresholds',
         ylab = 'Associated Optimal Accuracy',
         main = paste('Optimal Threshold and Accuracy,', nseed, 'Seeds'))
    abline(v = mean(mytopt), lty = 2); abline(h = mean(myaccopt), lty = 2)
    hist(mytopt, col = 'lightblue', border = 'darkblue', xlab = 'Optimal Thresholds',
         main = paste('Histogram of Optimal Thresholds,', nseed, 'Seeds'))
    par(op)
  }
  return(list(topt = mytopt, accopt = myaccopt,
              means = c(topt = mean(mytopt), accopt = mean(myaccopt)),
              sds = c(topt = sd(mytopt), accopt = sd(myaccopt)),
              tquant = quantile(mytopt, c(0.025, 0.5, 0.975))))
}



# ======================================================================
# PART 5: BOOTSTRAPPING, GINI IMPURITY, MONTE-CARLO (Lecture 14)
# ======================================================================

# Function 24
# generic BOOTSTRAP: resample the n observations (rows) with replacement nboot times
# and recompute stat() each time. X = vector, matrix or data frame.
# stat = function that takes the resampled data and returns a number OR a vector
#   e.g. mean:        stat = function(x) mean(x)
#   regression line:  X = data.frame(x, y), stat = function(d) coef(lm(d$y ~ d$x))
#     (resamples the PAIRS together, as in lecture)
# returns the bootstrap values, the percentile CI (alpha/2, 1-alpha/2), the
# original estimate, and a histogram with the CI for each statistic
mybootstrap = function(X, stat, nboot = 10000, alpha = 0.05, seed = 123, plot = TRUE){
  set.seed(seed)
  isvec = is.null(dim(X))
  n = if(isvec) length(X) else nrow(X)
  bvals = sapply(1:nboot, function(k){
    s = sample(1:n, n, replace = TRUE)
    stat(if(isvec) X[s] else X[s, , drop = FALSE])
  })
  if(is.null(dim(bvals))) bvals = rbind(bvals)
  ci = t(apply(bvals, 1, quantile, c(alpha/2, 1 - alpha/2)))
  if(plot){
    op = par(mfrow = c(1, nrow(bvals)))
    for(i in 1:nrow(bvals)){
      hist(bvals[i, ], breaks = 30, col = 'lightblue', border = 'darkblue',
           xlab = 'Bootstrap statistic', main = paste('Bootstrap Distribution', i))
      abline(v = ci[i, ], col = 'red', lwd = 2, lty = 2)
    }
    par(op)
  }
  est = if(isvec) stat(X) else stat(X)
  return(list(boot = bvals, ci = ci, est = est))
}

# Function 25
# Gini impurity of a categorical vector:  E_g = 1 - sum(p_k^2)
impure = function(Y){
  p = table(Y)/length(Y)
  return(1 - sum(p^2))
}

# Function 26
# Monte-Carlo test of purity:  H0: E_g(P) <= tau  vs  H1: E_g(P) > tau
# Y = observed categorical sample (size n, g groups). Builds the population whose
# impurity is exactly tau (one major class pmajor, g-1 equal minor classes),
# simulates nmc samples of size n from it, and returns the critical Gini value
# (1-alpha quantile) and the p-value  (1 + #{sim >= observed}) / (1 + nmc)
# pmajor solves the quadratic:  pM = (1 + sqrt((g-1)(g-1-g*tau)))/g
mygini_mc = function(Y, tau = 0.03, nmc = 25000, alpha = 0.05, seed = 1234, plot = TRUE){
  n = length(Y); g = length(unique(Y)); Eg = impure(Y)
  if(tau >= (g - 1)/g) stop('tau must be below the maximum impurity (g-1)/g')
  pmajor = (1 + sqrt((g - 1)*(g - 1 - g*tau)))/g
  pminor = (1 - pmajor)/(g - 1)
  
  set.seed(seed)
  groups = 1:g
  mcgini = sapply(1:nmc, function(k){
    S = sample(groups, n, replace = TRUE, prob = c(pmajor, rep(pminor, g - 1)))
    impure(S)
  })
  gini_crit = quantile(mcgini, 1 - alpha)
  pval = sum(c(Eg, mcgini) >= Eg)/(1 + nmc)
  if(plot){
    hist(mcgini, breaks = 30, col = 'lightblue', border = 'darkblue',
         xlab = 'Simulated Gini Impurity', main = 'M.C. Gini Impurities')
    abline(v = gini_crit, col = 'red', lwd = 2, lty = 2)
  }
  return(list(Eg = Eg, gini_crit = gini_crit, pval = pval,
              pmajor = pmajor, pminor = pminor, mcgini = mcgini))
}



# ======================================================================
# PART 6: DIMENSIONALITY REDUCTION (Lecture 15)
# ======================================================================

# Function 27
# PCA on the STANDARDIZED data (eigen-pairs of cov(Z) = correlation matrix of X)
# X = n x d numeric. tau = explained-variance threshold (default 0.85).
# computes the three "how many to keep" answers:
#   Kaiser (eigenvalue > 1), Elbow (read it off the scree plot!), 
#   Variance threshold (smallest k with F_k >= tau)
# dr = number of components to keep (default = the variance-threshold answer)
# returns eigenvalues, F_k, loadings (Lam), scores Xnew = Z*Lam, and the
# center/scale used (so test data can be projected with mypca_apply)
mypca = function(X, dr = NULL, tau = 0.85, plot = TRUE){
  X = as.matrix(X)
  Zx = scale(X)
  ctr = attr(Zx, 'scaled:center'); scl = attr(Zx, 'scaled:scale')
  Ez = eigen(cov(Zx))
  lam = Ez$values; d = length(lam)
  Flam = cumsum(lam)/sum(lam)
  kaiser = sum(lam > 1)
  vthresh = which(Flam >= tau)[1]
  if(is.null(dr)) dr = vthresh
  
  Lam = Ez$vectors[, 1:dr, drop = FALSE]
  rownames(Lam) = colnames(X); colnames(Lam) = paste0('PC', 1:dr)
  Xnew = Zx %*% Lam
  
  if(plot){
    plot(1:d, lam, type = 'b', pch = 16, xlab = 'Eigenvalue Rank',
         ylab = 'Eigenvalues', main = 'Eigenvalue(Z) Scree Plot')
    abline(h = 1, lty = 2)
  }
  return(list(values = lam, Flam = Flam, kaiser = kaiser, vthresh = vthresh,
              dr = dr, Lam = Lam, Xnew = Xnew, center = ctr, scale = scl))
}

# Function 28
# project NEW data (e.g. the testing set) onto the principal components found
# from the TRAINING data: standardize with the training mean/sd, then multiply by Lam
mypca_apply = function(pc, Xnew){
  Xnew = as.matrix(Xnew)
  Z = sweep(sweep(Xnew, 2, pc$center, '-'), 2, pc$scale, '/')
  return(Z %*% pc$Lam)
}

# Function 29
# LDA via the eigenpairs of S = (sSigma_E)^-1 (sSigma_M)
# X = n x d numeric, Y = class labels (g classes). 
# builds within-class means, within-class scatter (SSE), between-class scatter (SSM),
# eigen-decomposes S. At most min(d, g-1) eigenvalues are non-zero.
# dr = number of linear discriminants kept (default min(d, g-1))
# returns eigenvalues, proportion of the non-zero eigenvalues, discriminants L,
# and the scores Xlda = X*L
mylda = function(X, Y, dr = NULL, plot = TRUE){
  X = as.data.frame(X); Y = factor(Y)
  g = nlevels(Y); d = ncol(X)
  Xs = split(X, Y)
  nj = sapply(Xs, nrow)
  xbarc = sapply(Xs, colMeans)                 # d x g within-class means
  
  Sj = lapply(Xs, function(df){
    M = as.matrix(df)
    Xc = sweep(M, 2, colMeans(M), '-')
    t(Xc) %*% Xc
  })
  SSE = Reduce('+', Sj)
  
  xbar = colMeans(X)
  D = sweep(xbarc, 1, xbar, '-')
  SSM = D %*% diag(nj) %*% t(D)
  
  S = solve(SSE) %*% SSM
  ES = eigen(S)
  lam = Re(ES$values); vec = Re(ES$vectors)
  nz = min(d, g - 1)
  if(is.null(dr)) dr = nz
  
  LL = vec[, 1:dr, drop = FALSE]
  rownames(LL) = colnames(X); colnames(LL) = paste0('LD', 1:dr)
  Xlda = as.matrix(X) %*% LL
  
  if(plot){
    plot(1:d, lam, type = 'b', pch = 16, xlab = 'Eigenvalue Rank',
         ylab = 'Eigenvalues', main = 'Eigenvalue(S) Scree Plot')
  }
  return(list(values = lam, propvar = lam[1:nz]/sum(lam[1:nz]), L = LL, Xlda = Xlda,
              xbarc = xbarc, nj = nj, SSE = SSE, SSM = SSM, S = S))
}

# Function 30
# project NEW data onto the linear discriminants found from the TRAINING data
mylda_apply = function(lda, Xnew){
  return(as.matrix(Xnew) %*% lda$L)
}



# ======================================================================
# PART 7: OPTIMIZATION (Lecture 16)
# ======================================================================

# Function 31
# Newton-Raphson root finder:  x_{n+1} = x_n - f(x_n)/f'(x_n)
# f = function, fp = its derivative. To OPTIMIZE g, pass f = g' and fp = g''.
# returns the root and the whole path of iterates
mynewton = function(f, fp, x0, tol = 1e-10, maxit = 100){
  x = x0; path = x0
  for(i in 1:maxit){
    xn = x - f(x)/fp(x)
    path = c(path, xn)
    if(abs(xn - x) < tol) break
    x = xn
  }
  return(list(root = xn, path = path))
}

# Function 32
# Secant method (no derivative needed, two starting points x_-1 and x_0):
# x_{n+1} = x_n - f(x_n)*(x_n - x_{n-1})/(f(x_n) - f(x_{n-1}))
mysecant = function(f, xm1, x0, tol = 1e-10, maxit = 100){
  path = c(xm1, x0)
  for(i in 1:maxit){
    xa = path[length(path) - 1]; xb = path[length(path)]
    xn = xb - f(xb)*(xb - xa)/(f(xb) - f(xa))
    path = c(path, xn)
    if(is.nan(xn) || abs(xn - xb) < tol) break
  }
  return(list(root = xn, path = path))
}



# ======================================================================
# PART 8: MULTIVARIATE REGRESSION VIA OPTIMIZATION (Lecture 17)
# ======================================================================

# Function 33
# STANDARDIZE before penalizing: Z = scale(X), Yz = scale(Y), Xz = [1, Z]
# returns Xz, Yz, the sd's (sX, sY) needed to go back to original units, n, p, d
mystd = function(X, Y){
  X = as.matrix(X); Y = as.matrix(Y)
  Z = scale(X)[, , drop = FALSE]
  Yz = scale(Y)[, , drop = FALSE]
  Xz = cbind(1, Z)
  return(list(Xz = Xz, Yz = Yz, sX = apply(X, 2, sd), sY = apply(Y, 2, sd),
              n = nrow(Xz), p = ncol(Xz) - 1, d = ncol(Yz)))
}

# Function 34
# convert standardized slopes back to the original units:
# b_raw[i,j] = b[i,j] * sY_j / sX_i   (intercept row removed)
mybackraw = function(bhat, sX, sY){
  return(bhat[-1, , drop = FALSE] * outer(1/sX, sY))
}

# Function 35
# RIDGE regression, closed form:  bhat = (X'X + lam*D)^-1 X'Y,  D = diag(0,1,...,1)
# (the intercept is never penalized). Xz = [1, Z], Yz standardized.
myridge = function(Xz, Yz, lam){
  p = ncol(Xz) - 1
  D = diag(c(0, rep(1, p)))
  return(solve(t(Xz) %*% Xz + lam*D) %*% t(Xz) %*% Yz)
}

# Function 36
# Numerical fits with CVXR (ECOS solver) for the non-differentiable functionals
#   loss = 'ls'  & lam > 0 -> LASSO         ||Y - Xb||_F^2 + lam*sum|b|
#   loss = 'lad' & lam = 0 -> LAD           sum|Y - Xb|
#   loss = 'lad' & lam > 0 -> LAD-LASSO     sum|Y - Xb| + lam*sum|b|
# (intercept row is not penalized). Values within 1e-6 of 0 are set to exactly 0
# because the solver returns things like 1e-10 instead of 0.
myconvexfit = function(Xz, Yz, lam = 0, loss = 'ls'){
  if(!requireNamespace('CVXR', quietly = TRUE)) stop('install.packages("CVXR") first')
  library(CVXR)
  p = ncol(Xz) - 1; d = ncol(Yz)
  bhat = Variable(p + 1, d)
  pen = lam*sum(abs(bhat[2:(p+1), ]))
  obj = if(loss == 'ls') Minimize(sum_squares(Yz - Xz %*% bhat) + pen) else
    Minimize(sum(abs(Yz - Xz %*% bhat)) + pen)
  So = solve(Problem(obj), solver = 'ECOS')
  B = So$getValue(bhat)
  B[abs(B) < 1e-6] = 0
  dimnames(B) = list(colnames(Xz), colnames(Yz))
  return(list(bhat = B, value = So$value))
}

# Function 37
# value of all five functionals (LS, Ridge, Lasso, LAD, LAD-Lasso) at a coefficient
# matrix B. Use it on every fit to see that each fit wins its own functional.
myfunctionals = function(B, Yz, Xz, lam){
  R = Yz - Xz %*% B
  S2 = sum(R^2); S1 = sum(abs(R))
  L2 = sum(B[-1, ]^2); L1 = sum(abs(B[-1, ]))
  return(c(LS = S2, Ridge = S2 + lam*L2, Lasso = S2 + lam*L1,
           LAD = S1, LADLasso = S1 + lam*L1))
}

# Function 38
# tune lambda on a testing set: split once (mysplit), fit on the training rows with
# fitfun(Xz_train, Yz_train, lam) for every lam on the grid, and compute the
# testing RMSE. fitfun e.g. myridge, or function(X, Y, l) myconvexfit(X, Y, l)$bhat
# (note: Xz/Yz were standardized on all the data; fine for tuning, mention it)
mytunelam = function(Xz, Yz, fitfun, lams, ptrain = 0.65, seed = 2026, plot = TRUE){
  sp = mysplit(nrow(Xz), ptrain, seed)
  Xtr = Xz[sp$itrain, , drop = FALSE]; Ytr = Yz[sp$itrain, , drop = FALSE]
  Xte = Xz[sp$itest, , drop = FALSE];  Yte = Yz[sp$itest, , drop = FALSE]
  rmse = sapply(lams, function(l){
    B = fitfun(Xtr, Ytr, l)
    sqrt(mean((Yte - Xte %*% B)^2))
  })
  lbest = lams[which.min(rmse)]
  if(plot){
    plot(lams, rmse, type = 'b', pch = 16, cex = 0.6, xlab = 'lambda',
         ylab = 'Testing RMSE', main = 'Testing RMSE vs lambda')
    abline(v = lbest, lty = 2)
  }
  return(list(lbest = lbest, rmse = rmse, lams = lams))
}

# Function 39
# WEIGHTED least squares, a separate weight for each response:
# bhat_j = (X'W_j X)^-1 X'W_j Y_j. W is an n x d matrix (column j = weights for response j)
mywls = function(Xz, Yz, W){
  d = ncol(Yz)
  B = matrix(0, ncol(Xz), d, dimnames = list(colnames(Xz), colnames(Yz)))
  for(j in 1:d){
    Wj = diag(W[, j])
    B[, j] = solve(t(Xz) %*% Wj %*% Xz) %*% t(Xz) %*% Wj %*% Yz[, j]
  }
  return(B)
}

# Function 40
# weights 1/sigma^2 by group (e.g. job rank): E = residual matrix (n x d) from the LS fit,
# group = the label of each observation. Returns the n x d weight matrix and the
# group variances (s2)
mywts_group = function(E, group){
  s2 = apply(E, 2, function(e) tapply(e, group, var))
  return(list(W = 1/s2[as.character(group), , drop = FALSE], s2 = s2))
}

# ======================================================================
# PART 9: MULTICOLLINEARITY (VIF) FOR LOGISTIC REGRESSION
# ======================================================================

# Function 41
# Variance Inflation Factor for each predictor.  VIF_j = 1/(1 - R_j^2), where R_j^2 comes
# from regressing X_j on all the OTHER predictors. Shortcut: VIF_j is the j-th
# diagonal entry of the inverse of the correlation matrix of X.
# X = n x p numeric predictors (NO column of 1s, NO response).
# Rule of thumb: VIF > 5 is worrying, VIF > 10 is severe.
# (The response never enters, so this is the same VIF for binomial, multinomial or linear.)
myvif = function(X){
  X = as.matrix(X)
  Rx = cor(X)
  vif = diag(solve(Rx))
  names(vif) = colnames(X)
  return(vif)
}

# Function 42
# VIF straight from a fitted glm (binomial) object M: pulls the predictors out of
# the model matrix (dropping the intercept) and calls myvif()
# returns the VIFs plus the correlation matrix, so I can see WHICH pair is the problem
myvif_glm = function(M){
  X = model.matrix(M)[, -1, drop = FALSE]
  if(ncol(X) < 2) stop('need at least 2 predictors for a VIF')
  return(list(vif = myvif(X), cor = cor(X)))
}

# Function 43
# drop predictors with large VIF one at a time: remove the predictor with the largest
# VIF, recompute, and stop when every VIF <= cutoff. Returns the kept variables and
# the VIFs at every step (for the write-up)
myvif_drop = function(X, cutoff = 5){
  X = as.matrix(X)
  steps = list()
  while(ncol(X) > 1){
    v = myvif(X)
    steps[[length(steps) + 1]] = round(v, 3)
    if(max(v) <= cutoff) break
    X = X[, colnames(X) != names(which.max(v)), drop = FALSE]
  }
  return(list(kept = colnames(X), steps = steps))
}



# ======================================================================
# PART 10: ALL-POSSIBLE-MODELS SUMMARY TABLE (Lecture 13 "Result Aggregation")
# ======================================================================

# Function 44
# builds the Lecture 13 comparison table: one row per model (all 2^p - 1 combinations)
# with columns  Model | RMSE | Accuracy | AIC | BIC | McF R2
#   RMSE, AIC, BIC, McF  -> from the model fit on ALL the data
#   Accuracy             -> model fit on the TRAINING rows, classified on the TESTING rows
#                           (binomial: threshold tau; multinomial: max predicted probability)
# df = data frame with the response (column yname) and the predictors (vars)
# type = 'binomial' (Y must be 0/1) or 'multinom' (Y a factor, first level = reference)
# latex = TRUE also prints the table as a LaTeX tabular (copy it into the report)
# returns the table, the best model for each metric (lowest RMSE/AIC/BIC, highest
# Accuracy/McF), the fitted models, and the split used
# (McF always grows with more predictors, so judge it by how little it gains)
mymodeltable = function(df, yname, vars, type = 'binomial', ptrain = 0.65, seed = 2026,
                        tau = 0.5, digits = 4, latex = FALSE){
  sp = mysplit(nrow(df), ptrain, seed)
  res = myallsubsets(df, yname, vars, type, itrain = sp$itrain, tau = tau)
  
  tab = res$table[, c('RMSE', 'Accuracy', 'AIC', 'BIC', 'McF')]
  rownames(tab) = names(res$models)
  
  best = c(RMSE = rownames(tab)[which.min(tab[, 'RMSE'])],
           Accuracy = paste(rownames(tab)[tab[, 'Accuracy'] == max(tab[, 'Accuracy'])], collapse = ' / '),
           AIC = rownames(tab)[which.min(tab[, 'AIC'])],
           BIC = rownames(tab)[which.min(tab[, 'BIC'])],
           McF = rownames(tab)[which.max(tab[, 'McF'])])
  
  # rounded copy for display: RMSE/Accuracy/McF to `digits`, AIC/BIC to 1 decimal
  show = data.frame(Model = gsub('_', ', ', rownames(tab)),
                    RMSE = round(tab[, 'RMSE'], digits),
                    Accuracy = round(tab[, 'Accuracy'], digits),
                    AIC = round(tab[, 'AIC'], 1),
                    BIC = round(tab[, 'BIC'], 1),
                    McF = round(tab[, 'McF'], digits),
                    row.names = NULL)
  names(show)[6] = 'McF R2'
  
  if(latex){
    lab = sapply(rownames(tab), function(m){
      if(grepl('X[0-9]+', m)) paste0('$', gsub('X([0-9]+)', 'X_{\\1}', gsub('_', ', ', m)), '$')
      else gsub('_', '\\\\_', gsub('_', ', ', m))
    })
    cat('\\begin{tabular}{l|ccccc}\n\\hline\n')
    cat('Model & RMSE & Accuracy & AIC & BIC & $R^2_{McF}$ \\\\ \\hline\n')
    for(i in 1:nrow(show)){
      cat(lab[i], '&', sprintf(paste0('%.', digits, 'f'), show$RMSE[i]), '&',
          sprintf(paste0('%.', digits, 'f'), show$Accuracy[i]), '&',
          sprintf('%.1f', show$AIC[i]), '&', sprintf('%.1f', show$BIC[i]), '&',
          sprintf(paste0('%.', digits, 'f'), show[i, 6]), '\\\\\n')
    }
    cat('\\hline\n\\end{tabular}\n')
  }
  return(list(table = show, best = best, models = res$models, split = sp))
}

