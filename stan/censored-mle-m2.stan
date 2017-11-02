////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Title: Censored MLE, Model 2, Group Level
// Author: Steve Lane
// Date: Wednesday, 08 March 2017
// Synopsis: Sampling statements to fit a regression with censored outcome data.
// Includes boat-level intercept, and observation level location ID.
// Removed hull surface area.
// Time-stamp: <2017-11-02 00:48:02 (overlordR)>
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
  real midTrips[numBoat];
  /* Categorical predictors */
  /* Location of measurement */
  int<lower=1> numLoc;
  int<lower=1,upper=numLoc> locID[N];
  int<lower=1,upper=numLoc> locIDCens[nCens];
  /* Paint type */
  int<lower=1> numPaint;
  int<lower=1,upper=numPaint> paintType[numBoat];
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
  real betaMidTrips;
  /* Betas for categorical indicators */
  real betaLoc[numLoc];
  real betaPaint[numPaint];
  real betaType[numType];
  /* Alphas for modelled random effect */
  vector[numBoat] alphaBoat;
  /* Errors for categorical predictors */
  real<lower=0> sigma_alphaBoat;
  real<lower=0> sigmaLoc;
  real<lower=0> sigmaPaint;
  real<lower=0> sigmaType;
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
    alphaHat[n] = alphaBoat[n] + betaDays1 * days1[n] + betaDays2 * days2[n] + betaMidTrips * midTrips[n] + betaPaint[paintType[n]] + betaType[boatType[n]];
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
  betaDays1 ~ normal(0, 1);
  betaDays2 ~ normal(0, 1);
  betaMidTrips ~ normal(0, 1);
  /* Priors for categorical indicators */
  sigmaLoc ~ cauchy(0, 2.5);
  betaLoc ~ normal(0, sigmaLoc);
  sigmaPaint ~ cauchy(0, 2.5);
  betaPaint ~ normal(0, sigmaPaint);
  sigmaType ~ cauchy(0, 2.5);
  betaType ~ normal(0, sigmaType);
  sigma_alphaBoat ~ cauchy(0, 2.5);
  alphaBoat ~ normal(0, sigma_alphaBoat);
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
