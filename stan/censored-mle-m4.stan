////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Title: Censored MLE, Model 4, Group Level
// Author: Steve Lane
// Date: Wednesday, 08 March 2017
// Synopsis: Sampling statements to fit a regression with censored outcome data.
// Includes boat-level intercept, and observation level location ID.
// Restricted model form.
// Time-stamp: <2017-10-09 03:15:22 (overlordR)>
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
  /* Categorical predictors, entered as matrices of indicators */
  /* Location of measurement, hull as base case */
  int<lower=1> numLoc;
  matrix[N, numLoc - 1] locID;
  matrix[nCens, numLoc - 1] locIDCens;
  /* Boat type, yacht as base case */
  int<lower=1> numType;
  matrix[numBoat, numType - 1] boatType;
  /* Days1:BoatType interaction */
  matrix[numBoat, numType - 1] days1boat;
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
  vector[numLoc - 1] betaLoc;
  vector[numType - 1] betaType;
  /* Betas for interaction terms */
  vector[numType - 1] betaDaysType;
  /* Alphas for modelled random effect */
  vector[numBoat] alphaBoat;
  /* Errors for categorical predictors */
  real<lower=0> sigma_alphaBoat;
  /* Error */
  real<lower=0> sigma;
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
    alphaHat[n] = betaDays1 * days1[n] + betaDays2 * days2[n] + boatType[n] * betaType + days1boat[n] * betaDaysType;
  }
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
  betaDays1 ~ student_t(3, 0, 1);
  betaDays2 ~ student_t(3, 0, 1);
  /* Priors for categorical indicators */
  betaType ~ student_t(3, 0, 1);
  /* Priors for interactions */
  betaDaysType ~ student_t(3, 0, 1);
  /* Priors for modelled effects */
  sigma_alphaBoat ~ cauchy(0, 2.5);
  alphaBoat ~ cauchy(alphaHat, sigma_alphaBoat);
  /* Prior for observation (model) error */
  sigma ~ cauchy(0, 2.5);
  /* Observed log-likelihood */
  for(i in 1:N){
    target += normal_lpdf(logY[i] | muHat[i], sigma);
  }
  /* Censored log-likelihood */
  for(i in 1:nCens){
    target += normal_lcdf(logU | muHatCens[i], sigma);
  }
}

generated quantities{
  /* Need the log likelihood for leave one out prediction errors */
  vector[N + nCens] log_lik;
  /* Replications for posterior predictive checks */
  vector[N + nCens] y_ppc;
  for(i in 1:N){
    log_lik[i] = normal_lpdf(logY[i] | muHat[i], sigma);
    y_ppc[i] = normal_rng(muHat[i], sigma);
  }
  for(j in 1:nCens){
    log_lik[N + j] = normal_lcdf(logU | muHatCens[j], sigma);
    y_ppc[N + j] = normal_rng(muHatCens[j], sigma);
  }
}
