

data {
  int<lower=0> ntrials;
  vector[ntrials] cue;
  vector[ntrials] stim;
  vector[ntrials] u;
  int pred[ntrials];
  vector[ntrials] percept;
  int percept_bin[ntrials];
}

parameters {
  real<lower=0> alpha;
  real<lower=0> w1; // weight on stimulus
  real<lower=0> w2; // weight on expectation
  real<lower=0> percept_precision;
  real<lower=0> beta;
  
}


transformed parameters{
  
  vector<lower=0,upper=1>[ntrials] perceptmu; 
  vector<lower=0,upper=1>[ntrials] expect;
  vector<lower=0,upper=1>[ntrials+1] association;
  vector[ntrials] pe;
  
  association[1] = 0.5;
  
  for (t in 1:ntrials){
    
    if(cue[t] == 1){
        expect[t] = association[t];
      }else{
        expect[t] = 1-association[t];
      }
    
   
    perceptmu[t] = inv_logit(w1 * logit(stim[t]) + w2 * logit(expect[t]));
    
    
    if(cue[t] == 1){
        pe[t] = (perceptmu[t] - expect[t]);
      }else{
        pe[t] = -(perceptmu[t] - expect[t]);
      }
      
    association[t+1] = association[t] + alpha * pe[t];
  }
  
}


model {
  // priors for parameters
  target += beta_proportion_lpdf(alpha | 0.3, 1);
  target += beta_proportion_lpdf(w1 | 0.5, 1);
  target += beta_proportion_lpdf(w2 | 0.5, 1);
  target += lognormal_lpdf(percept_precision | 0, 1);
  target += lognormal_lpdf(beta | 0, 1);
  
  // generating data
  for (t in 1:ntrials){
    target += beta_proportion_lpdf(percept[t] | perceptmu[t], percept_precision);
    target += bernoulli_lpmf(percept_bin[t] | (perceptmu[t]^beta)/((perceptmu[t]^beta)+(1-perceptmu[t])^(beta)));
    target += bernoulli_lpmf(pred[t] |  (expect[t]^beta)/((expect[t]^beta)+(1-expect[t])^(beta)));
  }

}


