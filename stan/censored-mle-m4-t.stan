////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Title: Censored MLE, Model 4, Group Level
// Author: Steve Lane
// Date: Wednesday, 08 March 2017
// Synopsis: Sampling statements to fit a regression with censored outcome data.
// Includes boat-level intercept, and observation level location ID.
// Restricted model form.
// Based off M4, but with t distribution for outcome for added robustness.
// Time-stamp: <2017-10-17 22:50:49 (overlordR)>
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

data{
  // Data to be supplied to sampler
  /* Number of (fully observed) observations */
  int<lower=1> N;
  /* Number of (censored) observations */
  int<lower=1> nCens;
  /* Number of boats/vessels */
  int<lower=1> numBoat;
  /* Numeric/ordinal predictors */
  real days1[numBoat];
  real days2[numBoat];
  /* Categorical predictors */
  /* Location of measurement */
  int<lower=1> numLoc;
  int<lower=1,upper=numLoc> locID[N];
  int<lower=1,upper=numLoc> locIDCens[nCens];
  /* Boat type */
  int<lower=1> numType;
  int<lower=1,upper=numType> boatType[numBoat];
  /* Boat random effect */
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
  /* Betas for continuous */
  real betaDays1;
  real betaDays2;
  /* Betas for categorical indicators */
  real betaLoc[numLoc];
  real betaType[numType];
  /* Betas for interaction terms */
  real betaDaysType[numType];
  /* Alphas for modelled random effect */
  vector[numBoat] alphaBoat;
  /* Errors for categorical predictors */
  real<lower=0> sigma_alphaBoat;
  real<lower=0> sigmaLoc;
  real<lower=0> sigmaType;
  /* Errors for interaction terms */
  real<lower=0> sigmaDaysType;
  /* Error */
  real<lower=0> sigma;
  /* Degrees of freedom */
  real<lower=5> nu;
}

transformed parameters{
  // Make it easier for some sampling statements (not necessary)
  /* Regression for observed data */
  vector[N] muHat;
  /* Regression for censored data */
  vector[nCens] muHatCens;
  /* Regression for boat-level intercept */
  vector[numBoat] alphaHat;
  for(n in 1:numBoat){
    alphaHat[n] = alphaBoat[n] + betaDays1 * days1[n] + betaDays2 * days2[n] + betaType[boatType[n]] + betaDaysType[boatType[n]] * days1[n];
  }
  for(i in 1:N){
    muHat[i] = mu + betaLoc[locID[i]] + alphaHat[boatID[i]];
  }
  for(j in 1:nCens){
    muHatCens[j] = mu + betaLoc[locIDCens[j]] + alphaHat[boatIDCens[j]];
  }
}

model{
  // Model sampling statements
  /* Priors for intercept + continuous */
  mu ~ normal(0, 5);
  betaDays1 ~ student_t(3, 0, 1);
  betaDays2 ~ student_t(3, 0, 1);
  /* Priors for categorical indicators */
  sigmaLoc ~ cauchy(0, 2.5);
  betaLoc ~ student_t(3, 0, sigmaLoc);
  sigmaType ~ cauchy(0, 2.5);
  betaType ~ student_t(3, 0, sigmaType);
  /* Priors for interactions */
  sigmaDaysType ~ cauchy(0, 2.5);
  betaDaysType ~ student_t(3, 0, sigmaDaysType);
  /* Priors for modelled effects */
  sigma_alphaBoat ~ cauchy(0, 2.5);
  alphaBoat ~ student_t(3, 0, sigma_alphaBoat);
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
