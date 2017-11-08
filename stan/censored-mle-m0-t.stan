////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Title: Censored MLE, Model 0, Group Level, Robust
// Author: Steve Lane
// Date: Wednesday, 08 March 2017
// Synopsis: Sampling statements to fit a regression with censored outcome data.
// Includes boat-level intercept, and observation level location ID.
// No boat-level predictors.
// Based off M0, but with t distribution for outcome for added robustness.
// Time-stamp: <2017-11-03 00:20:58 (overlordR)>
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

data{
  // Data to be supplied to sampler
  /* Number of (fully observed) observations */
  int<lower=1> N;
  /* Number of (censored) observations */
  int<lower=1> nCens;
  /* Categorical predictors */
  /* Location of measurement */
  int<lower=1> numLoc;
  int<lower=1,upper=numLoc> locID[N];
  int<lower=1,upper=numLoc> locIDCens[nCens];
  /* Boat random effect */
  int<lower=1> numBoat;
  int<lower=1,upper=numBoat> boatID[N];
  int<lower=1,upper=numBoat> boatIDCens[nCens];
  /* Observed data */
  real<lower=1.5> Y[N];
  /* Truncated data (brute force, all equal */
  real<upper=min(Y)> U;
}

transformed data{
  real logY[N];
  real logU;
  for(n in 1:N){
    logY[n] = log(Y[n]);
  }
  logU = log(U);
}

parameters{
  // Parameters for the model
  /* Intercept */
  real mu;
  /* Raw betas for categorical indicators */
  vector[numLoc] locRaw;
  /* Raw alphas for modelled random effect */
  vector[numBoat] alphaRaw;
  /* Errors for categorical predictors */
  real<lower=0> sigma_alphaBoat;
  real<lower=0> sigmaLoc;
  /* Error */
  real<lower=0> sigma;
  /* Degrees of freedom */
  real<lower=2> nu;
}

transformed parameters{
  /* Location intercept */
  vector[numLoc] betaLoc;
  /* Boat intercept */
  vector[numBoat] alphaBoat;
  /* Regression for observed data */
  vector[N] muHat;
  /* Regression for censored data */
  vector[nCens] muHatCens;
  betaLoc = sigmaLoc * locRaw;
  alphaBoat = sigma_alphaBoat * alphaRaw;
  for(i in 1:N){
    muHat[i] = mu + betaLoc[locID[i]] + alphaBoat[boatID[i]];
  }
  for(j in 1:nCens){
    muHatCens[j] = mu + betaLoc[locIDCens[j]] + alphaBoat[boatIDCens[j]];
  }
}

model{
  // Model sampling statements
  /* Prior for intercept */
  mu ~ normal(0, 5);
  /* Priors for categorical indicators */
  sigmaLoc ~ cauchy(0, 2.5);
  locRaw ~ student_t(3, 0, 1);
  sigma_alphaBoat ~ cauchy(0, 2.5);
  alphaRaw ~ student_t(3, 0, 1);
  /* Prior for observation (model) error */
  sigma ~ cauchy(0, 2.5);
  /* Prior for df */
  nu ~ gamma(2, 10);
  /* Observed log-likelihood */
  for(i in 1:N){
    target += student_t_lpdf(logY[i] | nu, muHat[i], sigma);
  }
  /* Censored log-likelihood */
  for(j in 1:nCens){
    target += student_t_lcdf(logU | nu, muHatCens[j], sigma);
  }
}

generated quantities{
  /* Need the log likelihood for leave one out prediction errors */
  vector[N + nCens] log_lik;
  /* Replications for posterior predictive checks */
  vector[N + nCens] y_ppc;
  for(i in 1:N){
    log_lik[i] = student_t_lpdf(logY[i] | nu, muHat[i], sigma);
    y_ppc[i] = student_t_rng(nu, muHat[i], sigma);
  }
  for(j in 1:nCens){
    log_lik[N + j] = student_t_lcdf(logU | nu, muHatCens[j], sigma);
    y_ppc[N + j] = student_t_rng(nu, muHatCens[j], sigma);
  }
}
