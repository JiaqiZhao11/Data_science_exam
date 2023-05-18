get_experiment = function(){
  trials = c(40,40,40,40)
  bias = rep(c(0.2,0.8,0.3,0.7),trials)
  
  
  cue = rbinom(sum(trials),1,0.5)
  
  stim = array(sum(trials, NA))
  for(i in 1:sum(trials)){
    stim[i] = ifelse(cue[i] == 1, rbinom(1,1,bias[i]), rbinom(1,1,(1-bias[i])))
    
  }
  
  u = ifelse(cue == stim, 1,0)
  dd = data.frame(stim = stim, cue = cue, u = u, bias = bias)
  
  return(dd)
  
  
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
  expect = array(NA, ntrials)
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
    
    expect[i] = ifelse(cue[i] == 1, association[i], 1-association[i])
    
    pred[i] = rbinom(1,1,(expect[i]^beta)/((expect[i]^beta)+(1-expect[i])^(beta)))
    
    # correlated weight / Resterla Wagner paradigm
    #perceptmu[i] = ifelse(w1*stim[i]+(1-w1)*expect[i] == 0, 0.01, ifelse(w1*stim[i]+(1-w1)*expect[i] == 1, 0.99,w1*stim[i]+(1-w1)*expect[i]))
    
    # Bayesian inference paradigm: weighted Bayes
    stim[i] = ifelse(stim[i] == 0, 0.01, ifelse(stim[i] == 1, 0.99, stim[i]))
    expect[i] = ifelse(expect[i] == 0, 0.01, ifelse(expect[i] == 1, 0.99, expect[i]))
    perceptmu[i] = inv_logit_scaled(w1 * logit_scaled(stim[i]) + w2 * logit_scaled(expect[i]))
    
    # agent's rating on the intensity of the stimulus
    percept[i] = extraDistr::rprop(1,percept_precision,perceptmu[i])
    # agent's response on burning or not burning
    percept_bin[i] = rbinom(1,1,(perceptmu[i]^beta)/((perceptmu[i]^beta)+(1-perceptmu[i])^(beta)))
    
    pe[i] = ifelse(cue[i] == 1, (perceptmu[i] -expect[i]),(-(perceptmu[i] -expect[i])))
    association[i+1] = association[i]+alpha*pe[i]
      
  }
  
  df = data.frame(u = u[1:ntrials], cue = cue[1:ntrials], stim = stim[1:ntrials], expect = expect[1:ntrials], pred = pred[1:ntrials], perceptmu = perceptmu[1:ntrials], percept_bin = percept_bin[1:ntrials], percept = percept[1:ntrials], association = association[1:ntrials],
                  pe = pe[1:ntrials], alpha = alpha, w1 = w1, w2 = w2, percept_precision = percept_precision, beta = beta, x = 1:ntrials, desired = rep(bias,1))
  
  return(df)
  
}


weighted_Bayes_f_v2 = function(alpha, w1, w2, percept_precision, beta) {
  
  dd = get_experiment()
  
  cue = dd$cue
  stim = dd$stim
  u = dd$u
  bias = dd$bias
  
  ntrials = nrow(dd)
  
  percept = array(NA, ntrials)
  perceptmu = array(NA, ntrials)
  percept_bin = array(NA, ntrials)
  expect = array(NA, ntrials)
  pe = array(NA, ntrials)
  association = array(NA, ntrials)
  pred = array(NA,ntrials)
  
  association[1] = 0.5
  
  for (i in seq(ntrials)){
    
    expect[i] = ifelse(cue[i] == 1, association[i], 1-association[i])
    
    pred[i] = rbinom(1,1,(expect[i]^beta)/((expect[i]^beta)+(1-expect[i])^(beta)))
    
    # correlated weight / Resterla Wagner paradigm
    #perceptmu[i] = ifelse(w1*stim[i]+(1-w1)*expect[i] == 0, 0.01, ifelse(w1*stim[i]+(1-w1)*expect[i] == 1, 0.99,w1*stim[i]+(1-w1)*expect[i]))
    
    # Bayesian inference paradigm: weighted Bayes 
    stim[i] = ifelse(stim[i] == 0, 0.01, ifelse(stim[i] == 1, 0.99, stim[i]))
    expect[i] = ifelse(expect[i] == 0, 0.01, ifelse(expect[i] == 1, 0.99, expect[i]))
    perceptmu[i] = inv_logit_scaled(w1 * logit_scaled(stim[i]) + w2 * logit_scaled(expect[i]))
    
    # agent's rating on the intensity of the stimulus
    percept[i] = extraDistr::rprop(1,percept_precision,perceptmu[i])
    # agent's response on burning or not burning
    percept_bin[i] = rbinom(1,1,(perceptmu[i]^beta)/((perceptmu[i]^beta)+(1-perceptmu[i])^(beta)))
    
    pe[i] = ifelse(cue[i] == 1, (perceptmu[i] -expect[i]),(-(perceptmu[i] -expect[i])))
    association[i+1] = association[i]+alpha*pe[i]
    
  }
  
  df = data.frame(u = u[1:ntrials], cue = cue[1:ntrials], stim = stim[1:ntrials], expect = expect[1:ntrials], pred = pred[1:ntrials], perceptmu = perceptmu[1:ntrials], percept_bin = percept_bin[1:ntrials], percept = percept[1:ntrials], association = association[1:ntrials],
                  pe = pe[1:ntrials], alpha = alpha, w1 = w1, w2 = w2, percept_precision = percept_precision, beta = beta, x = 1:ntrials, desired = rep(bias,1), id = rnorm(1,0,1))
  
  return(df)
  
}

hier_Bayes_f = function(parameters){
  
  nsubs = parameters$nsubs
  
  df = data.frame()
  for (s in 1:nsubs){
    
    df1 = weighted_Bayes_f_v2(w1 = extraDistr::rprop(1,parameters$kappa_w1,parameters$mu_w1),
                          w2 = extraDistr::rprop(1,parameters$kappa_w2,parameters$mu_w2),    
                          alpha = extraDistr::rprop(1,parameters$kappa_alpha,parameters$mu_alpha),
                          percept_precision = rlnorm(1,parameters$mu_percept_precision,parameters$sd_percept_precision),
                          beta = rlnorm(1,parameters$mu_beta,parameters$sd_beta)
    )
    
    df = rbind(df,df1)
    
  }
  
  
  return(list(df, parameters))
  
  
  
}
