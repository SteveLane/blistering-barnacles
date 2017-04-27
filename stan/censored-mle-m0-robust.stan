////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Title: Censored MLE, Model 0, Group Level, Robust
// Author: Steve Lane
// Date: Wednesday, 08 March 2017
// Synopsis: Sampling statements to fit a regression with censored outcome data.
// Includes boat-level intercept, and observation level location ID.
// No boat-level predictors.
// Based off M0, but with cauchy distribution for outcome for added robustness.
// Time-stamp: <2017-04-27 10:32:10 (slane)>
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

functions{
  /* Log-cauchy pdf, returns sum of log-density */
  real logcauchy_l(real y, real mu, real sigma){
    real lpdf;
    lpdf = log(sigma) - log(pi()) - log(y) - log((log(y) - mu)^2 + sigma^2);
    return lpdf;
  }
  /* Allows vectors to come in */
  /* real logcauchy_lpdf(vector y, vector mu, real sigma){ */
  /*   vector[num_elements(y)] lpdf; */
  /*   vector[num_elements(y)] ly; */
  /*   real log_sigma; */
  /*   real log_pi; */
  /*   real sigma2; */
  /*   real lprob; */
  /*   ly = log(y); */
  /*   log_sigma = log(sigma); */
  /*   log_pi = log(pi()); */
  /*   sigma2 = sigma^2; */
  /*   for(i in 1:num_elements(y)){ */
  /*     lpdf[i] = log_sigma - log_pi - ly[i] - log((ly[i] - mu[i])^2 + sigma2); */
  /*   } */
  /*   lprob = sum(lpdf); */
  /*   return lprob; */
  /* } */
}

data{
  // Data to be supplied to sampler
  /* Number of (fully observed) observations */
  int<lower=1> N;
  /* Number of (censored) observations */
  int<lower=1> nCens;
  /* Categorical predictors, entered as matrices of indicators */
  /* Location of measurement, hull as base case */
  int<lower=1> numLoc;
  matrix[N, numLoc - 1] locID;
  matrix[nCens, numLoc - 1] locIDCens;
  /* Boat random effect */
  int<lower=1> numBoat;
  int<lower=1,upper=numBoat> boatID[N];
  int<lower=1,upper=numBoat> boatIDCens[nCens];
  /* Observed data */
  real<lower=1.5> Y[N];
  /* Truncated data (brute force, all equal */
  real<upper=min(Y)> U;
}

parameters{
  // Parameters for the model
  /* Intercept */
  real mu;
  /* Betas for categorical indicators */
  vector[numLoc - 1] betaLoc;
  /* Alphas for modelled random effect */
  vector[numBoat] alphaBoat;
  /* Errors for categorical predictors */
  real<lower=0> sigma_alphaBoat;
  /* Censored data */
  vector<lower=0,upper=U>[nCens] yCens;
  /* Error */
  real<lower=0> sigma;
}

transformed parameters{
  // Make it easier for some sampling statements (not necessary)
    /* Regression for observed data */
  vector[N] muHat;
  /* Regression for censored data */
  vector[nCens] muHatCens;
  for(i in 1:N){
    muHat[i] = mu + locID[i] * betaLoc + alphaBoat[boatID[i]];
  }
  for(j in 1:nCens){
    muHatCens[j] = mu + locIDCens[j] * betaLoc + alphaBoat[boatIDCens[j]];
  }
}

model{
  // Model sampling statements
  /* Priors for categorical indicators */
  betaLoc ~ student_t(3, 0, 1);
  /* Priors for modelled random effect */
  mu ~ normal(0, 5);
  sigma_alphaBoat ~ cauchy(0, 2.5);
  alphaBoat ~ cauchy(0, sigma_alphaBoat);
  /* Prior for observation (model) error */
  sigma ~ cauchy(0, 2.5);
  /* Observed log-likelihood */
  for(i in 1:N){
    target += logcauchy_l(Y[i], muHat[i], sigma);
  }
  /* Y ~ logcauchy(muHat, sigma); */
  /* Censored log-likelihood */
  for(j in 1:nCens){
    target += logcauchy_l(yCens[j], muHatCens[j], sigma);
  }
  /* yCens ~ logcauchy(muHatCens, sigma); */
}

generated quantities{
  /* Need the log likelihood for leave one out prediction errors */
  vector[N + nCens] log_lik;
  {
    real linPred;
    for(i in 1:N){
      linPred = mu + locID[i] * betaLoc + alphaBoat[boatID[i]];
      log_lik[i] = logcauchy_l(Y[i], linPred, sigma);
    }
    for(j in 1:nCens){
      linPred = mu + locIDCens[j] * betaLoc + alphaBoat[boatIDCens[j]];
      log_lik[N + j] = logcauchy_l(yCens[j], linPred, sigma);
    }
  }
}
