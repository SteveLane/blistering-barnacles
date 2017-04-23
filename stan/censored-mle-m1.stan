////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Title: Censored MLE, Model 1, Group Level
// Author: Steve Lane
// Date: Wednesday, 08 March 2017
// Synopsis: Sampling statements to fit a regression with censored outcome data.
// Includes boat-level intercept, and observation level location ID.
// All boat-level intercept predictors included.
// Time-stamp: <2017-04-24 09:52:50 (slane)>
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
  real hullSA[numBoat];
  /* Categorical predictors, entered as matrices of indicators */
  /* Location of measurement, hull as base case */
  int<lower=1> numLoc;
  matrix[N, numLoc - 1] locID;
  matrix[nCens, numLoc - 1] locIDCens;
  /* Paint type, ablative as base case */
  int<lower=1> numPaint;
  matrix[numBoat, numPaint - 1] paintType;
  /* Boat type, yacht as base case */
  int<lower=1> numType;
  matrix[numBoat, numType - 1] boatType;
  /* Boat random effect */
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
  /* Betas for continuous */
  real betaDays1;
  real betaDays2;
  real betaMidTrips;
  real betaHullSA;
  /* Betas for categorical indicators */
  vector[numLoc - 1] betaLoc;
  vector[numPaint - 1] betaPaint;
  vector[numType - 1] betaType;
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
  /* Regression for boat-level intercept */
  vector[numBoat] alphaHat;
  for(n in 1:numBoat){
    alphaHat[n] = mu + betaDays1 * days1[n] + betaDays2 * days2[n] + betaMidTrips * midTrips[n] + betaHullSA * hullSA[n] + paintType[n] * betaPaint + boatType[n] * betaType;
  }
  for(i in 1:N){
    muHat[i] = locID[i] * betaLoc + alphaBoat[boatID[i]];
  }
  for(j in 1:nCens){
    muHatCens[j] = locIDCens[j] * betaLoc + alphaBoat[boatIDCens[j]];
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
  betaMidTrips ~ student_t(3, 0, 1);
  betaHullSA ~ student_t(3, 0, 1);
  /* Priors for categorical indicators */
  betaPaint ~ student_t(3, 0, 1);
  betaType ~ student_t(3, 0, 1);
  sigma_alphaBoat ~ cauchy(0, 2.5);
  alphaBoat ~ cauchy(alphaHat, sigma_alphaBoat);
  /* Prior for observation (model) error */
  sigma ~ cauchy(0, 2.5);
  /* Observed log-likelihood */
  Y ~ lognormal(muHat, sigma);
  /* Censored log-likelihood */
  yCens ~ lognormal(muHatCens, sigma);
}

generated quantities{
  /* Need the log likelihood for leave one out prediction errors */
  vector[N + nCens] log_lik;
  {
    real linPred;
    for(i in 1:N){
      linPred = locID[i] * betaLoc + alphaBoat[boatID[i]];
      log_lik[i] = lognormal_lpdf(Y[i] | linPred, sigma);
    }
    for(j in 1:nCens){
      linPred = locIDCens[j] * betaLoc + alphaBoat[boatIDCens[j]];
      log_lik[N + j] = lognormal_lpdf(yCens[j] | linPred, sigma);
    }
  }
}
