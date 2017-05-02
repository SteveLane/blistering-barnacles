////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Title: Censored MLE, Model 3, Group Level
// Author: Steve Lane
// Date: Wednesday, 08 March 2017
// Synopsis: Sampling statements to fit a regression with censored outcome data.
// Includes boat-level intercept, and observation level location ID.
// Adds in some interactions terms.
// Time-stamp: <2017-05-02 11:56:00 (slane)>
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
  /* Days1:BoatType interaction */
  matrix[numBoat, numType - 1] days1boat;
  /* MidTrips:BoatType interaction */
  matrix[numBoat, numType - 1] midTripsboat;
  /* MidTrips:PaintType interaction */
  matrix[numBoat, numPaint - 1] midTripspaint;
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
  /* Betas for categorical indicators */
  vector[numLoc - 1] betaLoc;
  vector[numPaint - 1] betaPaint;
  vector[numType - 1] betaType;
  /* Betas for interaction terms */
  vector[numType - 1] betaDaysType;
  vector[numType - 1] betaTripsType;
  vector[numPaint - 1] betaTripsPaint;
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
    alphaHat[n] = betaDays1 * days1[n] + betaDays2 * days2[n] + betaMidTrips * midTrips[n] + paintType[n] * betaPaint + boatType[n] * betaType + days1boat[n] * betaDaysType + midTripsboat[n] * betaTripsType + midTripspaint[n] * betaTripsPaint;
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
  betaMidTrips ~ student_t(3, 0, 1);
  /* Priors for categorical indicators */
  betaPaint ~ student_t(3, 0, 1);
  betaType ~ student_t(3, 0, 1);
  /* Priors for interactions */
  betaDaysType ~ student_t(3, 0, 1);
  betaTripsType ~ student_t(3, 0, 1);
  betaTripsPaint ~ student_t(3, 0, 1);
  /* Priors for modelled effects */
  sigma_alphaBoat ~ cauchy(0, 2.5);
  alphaBoat ~ cauchy(alphaHat, sigma_alphaBoat);
  /* Prior for observation (model) error */
  sigma ~ cauchy(0, 2.5);
  /* Observed log-likelihood */
  Y ~ lognormal(muHat, sigma);
  /* Censored log-likelihood */
  for(i in 1:nCens){
    target += lognormal_lcdf(U | muHatCens[i], sigma);
  }
}

generated quantities{
  /* Need the log likelihood for leave one out prediction errors */
  vector[N + nCens] log_lik;
  /* Replications for posterior predictive checks */
  vector[N + nCens] y_ppc;
  for(i in 1:N){
    log_lik[i] = lognormal_lpdf(Y[i] | muHat[i], sigma);
    y_ppc[i] = lognormal_rng(muHat[i], sigma);
  }
  for(j in 1:nCens){
    log_lik[N + j] = lognormal_lcdf(U | muHatCens[j], sigma);
    y_ppc[N + j] = lognormal_rng(muHatCens[j], sigma);
  }
}
