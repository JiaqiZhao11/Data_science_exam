

rm_agent = function(bias,trials){
  u = c()
  for (i in 1:length(trials)){
    u1 = rbinom(n = trials[i], size = 1, prob = bias[i])
    u = c(u,u1)
  }
  return(u)
}

source((here("~","Advanced-cognitive-modeling","assignment2","hgf_agent.R")))

trials = c(50,50,50,50)
bias = rep(c(0.1,0.9,0.1,0.9),trials)


cue = rbinom(sum(trials),1,0.5)

stim = array(sum(trials, NA))
for(i in 1:200){
  stim[i] = ifelse(cue[i] == 1, rbinom(1,1,bias[i]), rbinom(1,1,(1-bias[i])))
  
}

u = ifelse(cue == stim, 1,0)
dd = data.frame(stim = stim, cue = cue, u = u)

dd %>% ggplot(aes(x = 1:nrow(.), y = u, col = as.factor(stim)))+geom_point()+theme_classic()

ntrials = sum(trials)

percept = array(NA, ntrials)
perceptmu = array(NA, ntrials)
exp = array(NA, ntrials)

association = array(NA, ntrials)
per_con = array(NA, ntrials)
pred = array(NA,ntrials)

association[1] = 0.5

alpha = 0.2
w1 = 2
w2 = 

for (i in seq(ntrials)){
  exp[i] = ifelse(cue[i] == 1, association[i], 1-association[i])
  
  pred[i] = rbinom(1,1,exp[i])

  # correlated weight / Resterla Wagner paradigm
  #perceptmu[i] = ifelse(w1*stim[i]+(1-w1)*exp[i] == 0, 0.01, ifelse(w1*stim[i]+(1-w1)*exp[i] == 1, 0.99,w1*stim[i]+(1-w1)*exp[i]))
  
  # Bayesian inference paradigm: weighted Bayes ## Do we add bias here?
  stim[i] = ifelse(stim[i] == 0, 0.01, ifelse(stim[i] == 1, 0.99, stim[i]))
  exp[i] = ifelse(exp[i] == 0, 0.01, ifelse(exp[i] == 1, 0.99, exp[i]))
  perceptmu[i] = inv_logit_scaled(w1 * logit_scaled(stim[i]) + w2 * logit_scaled(exp[i]))
  percept[i] = extraDistr::rprop(1,10,perceptmu[i])
  
  
  per_con[i] = ifelse(cue[i] == 1, percept[i], 1-percept[i])
  association[i+1] = association[i]+alpha*(per_con[i] -association[i])
  
}


data.frame(per_con = per_con[1:200],stim = stim[1:200], percept = percept[1:200], pred = pred[1:200],association = association[1:200]) %>% 
  ggplot(aes(x = 1:nrow(.), y = u))+geom_line(aes(y = association))+geom_point()+geom_point(aes(x = 1:200, y = per_con), col = "red")+geom_point(aes(x = 1:200, y = stim-0.1), col = "blue")


data.frame(per_con = per_con[1:200],stim = stim[1:200], percept = percept[1:200], pred = pred[1:200],association = association[1:200]) %>% 
  ggplot(aes(as.factor(stim), y = percept))+geom_boxplot()




data.frame(stim = stim[1:200], percept = percept[1:200], pred = pred[1:200],association = association[1:200], cue = cue[1:200]) %>% ggplot(aes(percept, y = association))+geom_point()+geom_smooth(method = "lm")+facet_grid(~cue)
data.frame(stim = stim[1:200], percept = percept[1:200], pred = pred[1:200],association = association[1:200], cue = cue[1:200]) %>% ggplot(aes(as.factor(pred), y = percept, col = as.factor(stim)))+geom_boxplot()+facet_grid(~cue)








