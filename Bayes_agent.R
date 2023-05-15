weighted_bayes_f = function(previous_expectation, present_stimulus) {
  # waht would affect association?
  # RW association = inv_logit_scaled(beta_1 * logit_scaled(previous_association) + beta_2 * logit_scaled(e))
  
  present_percept = inv_logit_scaled(w_1 * logit_scaled(expectation) + w_2 * logit_scaled(present_stimulus))
  
  return(present_percept)
  
}

weighted_Bayes_f = function(parameters) {
  
  trials = c(40,40,40,40)
  bias = rep(c(0.1,0.9,0.1,0.9),trials)
  ntrials = sum(trials)
  
  cue = rbinom(sum(trials),1,0.5)
  
  stim = array(sum(trials, NA))
  for(i in 1:ntrials){
    stim[i] = ifelse(cue[i] == 1, rbinom(1,1,bias[i]), rbinom(1,1,(1-bias[i])))
    
  }
  
  u = ifelse(cue == stim, 1,0)
  dd = data.frame(stim = stim, cue = cue, u = u)
  
  dd %>% ggplot(aes(x = 1:nrow(.), y = u, col = as.factor(stim)))+geom_point()+theme_classic()
  
  percept = array(NA, ntrials)
  perceptmu = array(NA, ntrials)
  percept_bin = array(NA, ntrials)
  exp = array(NA, ntrials)
  pe = array(NA, ntrials)
  association = array(NA, ntrials)
  pred = array(NA,ntrials)
  
  association[1] = 0.5
  
  alpha = parameters$alpha
  w1 = parameters$w1
  w2 = parameters$w2
  percept_precision = parameters$percept_precision
  beta = parameters$beta
  
    
  for (i in seq(ntrials)){
    
    exp[i] = ifelse(cue[i] == 1, association[i], 1-association[i])
    
    pred[i] = rbinom(1,1,(exp[i]^beta)/((exp[i]^beta)+(1-exp[i])^(beta)))
    
    # correlated weight / Resterla Wagner paradigm
    #perceptmu[i] = ifelse(w1*stim[i]+(1-w1)*exp[i] == 0, 0.01, ifelse(w1*stim[i]+(1-w1)*exp[i] == 1, 0.99,w1*stim[i]+(1-w1)*exp[i]))
    
    # Bayesian inference paradigm: weighted Bayes ## Do we add bias here?
    stim[i] = ifelse(stim[i] == 0, 0.01, ifelse(stim[i] == 1, 0.99, stim[i]))
    exp[i] = ifelse(exp[i] == 0, 0.01, ifelse(exp[i] == 1, 0.99, exp[i]))
    perceptmu[i] = inv_logit_scaled(w1 * logit_scaled(stim[i]) + w2 * logit_scaled(exp[i]))
    # is this the rating?
    percept[i] = extraDistr::rprop(1,percept_precision,perceptmu[i])
    percept_bin[i] = rbinom(1,1,(perceptmu[i]^beta)/((perceptmu[i]^beta)+(1-perceptmu[i])^(beta)))
    
    pe[i] = ifelse(cue[i] == 1, (perceptmu[i] -exp[i]),(-(perceptmu[i] -exp[i])))
    association[i+1] = association[i]+alpha*pe[i]
      
  }
  
  df = data.frame(u = u[1:ntrials], cue = cue[1:ntrials], stim = stim[1:ntrials], exp = exp[1:ntrials], pred = pred[1:ntrials], perceptmu = perceptmu[1:ntrials], percept_bin = percept_bin[1:ntrials], percept = percept[1:ntrials], association = association[1:ntrials],
                  pe = pe[1:ntrials], alpha = alpha, w1 = w1, w2 = w2, percept_precision = percept_precision, beta = beta, x = 1:ntrials, desired = rep(bias,1))
  
  return(df)
  
}
