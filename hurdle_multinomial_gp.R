library(ggpubr)
library(tidyverse)
library(rstan)
library(bayesplot)

#Red data
df <- read.csv("All_Files_For_Publication/data/Updated 2017 and 2015 Cross-sectional Data for Bayesian Analsysis_with_UTM.csv")

#Two models to choose from - thesis or publication (includes Oct-Dec rainfall and populations of sheep cattle and goats)
modelThesis <- FALSE # Model 3 in Table S1
modelPub <- FALSE # Model 4 in Table S1 
model2 <- TRUE #Model 1 in Table S1 - Selected Model
modelPres <- FALSE # Model 5 in Table S1
modelThesis_Oct <- FALSE # Model 6 in Table S1
#Note Model 2 requires the use of the hurdle-multinomial-re.stan and hurdle-multinomial-re_for_publication.r

#Organize data
df<- df%>%
  rename(Y = UTM_Lat)%>%
  rename(X = UTM_Long)

#Add column to indicate if there were any deaths or abortions
df <- df%>%
  mutate(RVF_Anydeaths = if_else(TotalRVFDeathsDomR.2010.Stand >1, 1, 0))%>%
  mutate(RVF_Anyaborts = if_else(TotalRVFAbortionsDomR.2010.Stand >1, 1, 0))


#f18.Production.15 is a categorical variable with three options, need to set them as integers so the model can index it appropriately
#relevel so that the small holders are the reference group, then turn into an integer

prod_dat <- as.factor(df$f18.Production.15)
prod_dat <- relevel(prod_dat, ref="Small holder")
df$f18.Production.15 <- prod_dat


# Scale coordinates so the maximum distance between two points is one and center
distance_scale <- max(dist(cbind(df$X, df$Y)))
locs_standardized <- scale(cbind(df$X, df$Y) / distance_scale, center = TRUE, scale = FALSE)


# Convert the model to C++ code, then compile (Not normally done, just for
# illustration. stan_model() can do this in one step)
hmgp_mod_c <- stanc("hurdle-multinomial-gp.stan", verbose = TRUE)
hmgp_mod <- stan_model(stanc_ret = hmgp_mod_c, verbose = TRUE)

# This allows us to use the functions we define in stan code, notably rsoftmax()
expose_stan_functions(hmgp_mod)
rsoftmax(c(1,1,1))
sum(rsoftmax(c(0.2, 2, 3)))



# Assemble the data.  All terms in the `data` block in the Stan code
# must be defined as part of a named list.  We use normalized values for more
# efficient sampling. This would be taken care of under the hood in `brms` models.
if(modelThesis == TRUE){
standata <- lst(
  # Matrix of the three possible outcomes for each farm
  Y = cbind(
    df$TotalPopDomR_2010 - df$TotalRVFDeathsDomR.2010 - df$TotalRVFAbortionsDomR.2010,
    df$TotalRVFDeathsDomR.2010,
    df$TotalRVFAbortionsDomR.2010
  ),
  N = nrow(Y),
  ncat = ncol(Y),

  #data for residual estimation
  outbreak = df$X3.Outbreak,
  num_animals = df$TotalPopDomR_2010,

  # Design matrix for outcomes
  X = model.matrix(~PrctSheep_Ave.2010.Stand +
                     TotalPopDomR_2010.Stand +
                     X5.RVF.Vax.before.Outbreak.2010 +
                     Num_Natural_WaterSources.Ave.Stand +
                     CumSum.2m.Rain.Stand +
                     Prct.Rain.Days.2m.Stand,
                   data = df),
  K = ncol(X),

  # Design matrix for hurdles (not inluding random effects)
  X_eta = model.matrix(~Sum.Ruminants.Gained.Ave.Stand +
                         Max.Dist.Purch.Ave.Stand +
                         f18.Production.15 +
                         f9.Pan.Num.Ave.Stand +
                         f14.Size.Farm.15.Stand +
                         f5.Wildlife.Mix.15 +
                         Ave.2m.CumSum.Stand,
                       data = df),
  K_eta = ncol(X_eta),

  # Gaussian Process
  Dgp = 2,  # Dimension
  Xgp = locs_standardized, #Values


  prior_only = FALSE # Change this to sample from the priors rather than posterior
)
}

if(modelThesis_Oct == TRUE){
  standata <- lst(
    # Matrix of the three possible outcomes for each farm
    Y = cbind(
      df$TotalPopDomR_2010 - df$TotalRVFDeathsDomR.2010 - df$TotalRVFAbortionsDomR.2010,
      df$TotalRVFDeathsDomR.2010,
      df$TotalRVFAbortionsDomR.2010
    ),
    N = nrow(Y),
    ncat = ncol(Y),

    #data for residual estimation
    outbreak = df$X3.Outbreak,
    num_animals = df$TotalPopDomR_2010,

    # Design matrix for outcomes
    X = model.matrix(~PrctSheep_Ave.2010.Stand +
                       TotalPopDomR_2010.Stand +
                       X5.RVF.Vax.before.Outbreak.2010 +
                       Num_Natural_WaterSources.Ave.Stand +
                       CumSum.2m.Rain.Stand +
                       Prct.Rain.Days.2m.Stand,
                     data = df),
    K = ncol(X),

    # Design matrix for hurdles (not inluding random effects)
    X_eta = model.matrix(~Sum.Ruminants.Gained.Ave.Stand +
                           Max.Dist.Purch.Ave.Stand +
                           f18.Production.15 +
                           f9.Pan.Num.Ave.Stand +
                           f14.Size.Farm.15.Stand +
                           f5.Wildlife.Mix.15 +
                           CumSum.2m.Rain_Oct.Stand,
                         data = df),
    K_eta = ncol(X_eta),

    # Gaussian Process
    Dgp = 2,  # Dimension
    Xgp = locs_standardized, #Values


    prior_only = FALSE # Change this to sample from the priors rather than posterior
  )
}

if(modelPub == TRUE){
  standata <- lst(
    # Matrix of the three possible outcomes for each farm
    Y = cbind(
      df$TotalPopDomR_2010 - df$TotalRVFDeathsDomR.2010 - df$TotalRVFAbortionsDomR.2010,
      df$TotalRVFDeathsDomR.2010,
      df$TotalRVFAbortionsDomR.2010
    ),
    N = nrow(Y),
    ncat = ncol(Y),

    #data for residual estimation
    outbreak = df$X3.Outbreak,
    num_animals = df$TotalPopDomR_2010,

    # Design matrix for outcomes
    X = model.matrix(~ X2a.Total.Cattle.2010.Stand +
                       X2b.Total.Sheep.2010.Stand +
                       X2c.Total.Goats.2010.Stand +
                       X5.RVF.Vax.before.Outbreak.2010 +
                       Num_Natural_WaterSources.Ave.Stand +
                       CumSum.2m.Rain.Stand +
                       Prct.Rain.Days.2m.Stand,
                     data = df),
    K = ncol(X),

    # Design matrix for hurdles (not inluding random effects)
    X_eta = model.matrix(~Sum.Ruminants.Gained.Ave.Stand +
                           Max.Dist.Purch.Ave.Stand +
                           f18.Production.15 +
                           f9.Pan.Num.Ave.Stand +
                           f14.Size.Farm.15.Stand +
                           f5.Wildlife.Mix.15 +
                           CumSum.2m.Rain_Oct.Stand,
                         data = df),
    K_eta = ncol(X_eta),

    # Gaussian Process
    Dgp = 2,  # Dimension
    Xgp = locs_standardized, #Values


    prior_only = FALSE # Change this to sample from the priors rather than posterior
  )
}


if(model2 == TRUE){#Can't do this - %Cattle and %sheep are collinear
  standata <- lst(
    # Matrix of the three possible outcomes for each farm
    Y = cbind(
      df$TotalPopDomR_2010 - df$TotalRVFDeathsDomR.2010 - df$TotalRVFAbortionsDomR.2010,
      df$TotalRVFDeathsDomR.2010,
      df$TotalRVFAbortionsDomR.2010
    ),
    N = nrow(Y),
    ncat = ncol(Y),

    #data for residual estimation
    outbreak = df$X3.Outbreak,
    num_animals = df$TotalPopDomR_2010,

    # Design matrix for outcomes
    X = model.matrix(~ PrctSheep_Ave.2010.Stand +
                       PrctGoats_Ave.2010.Stand +
                       TotalPopDomR_2010.Stand+
                       X5.RVF.Vax.before.Outbreak.2010 +
                       Num_Natural_WaterSources.Ave.Stand +
                       CumSum.2m.Rain.Stand +
                       Prct.Rain.Days.2m.Stand,
                     data = df),
    K = ncol(X),

    # Design matrix for hurdles (not inluding random effects)
    X_eta = model.matrix(~Sum.Ruminants.Gained.Ave.Stand +
                           Max.Dist.Purch.Ave.Stand +
                           f18.Production.15 +
                           f9.Pan.Num.Ave.Stand +
                           f14.Size.Farm.15.Stand +
                           f5.Wildlife.Mix.15 +
                           CumSum.2m.Rain_Oct.Stand,
                         data = df),
    K_eta = ncol(X_eta),

    # Gaussian Process
    Dgp = 2,  # Dimension
    Xgp = locs_standardized, #Values


    prior_only = FALSE # Change this to sample from the priors rather than posterior
  )
}

if(modelPres == TRUE){
  standata <- lst(
    # Matrix of the three possible outcomes for each farm
    Y = cbind(
      df$TotalPopDomR_2010 - df$TotalRVFDeathsDomR.2010 - df$TotalRVFAbortionsDomR.2010,
      df$TotalRVFDeathsDomR.2010,
      df$TotalRVFAbortionsDomR.2010
    ),
    N = nrow(Y),
    ncat = ncol(Y),

    #data for residual estimation
    outbreak = df$X3.Outbreak,
    num_animals = df$TotalPopDomR_2010,

    # Design matrix for outcomes
    X = model.matrix(~ Sheep.Present.2010 +
                       Cattle.Present.2010 +
                       Goats.Present.2010 +
                       TotalPopDomR_2010.Stand+
                       X5.RVF.Vax.before.Outbreak.2010 +
                       Num_Natural_WaterSources.Ave.Stand +
                       CumSum.2m.Rain.Stand +
                       Prct.Rain.Days.2m.Stand,
                     data = df),
    K = ncol(X),

    # Design matrix for hurdles (not inluding random effects)
    X_eta = model.matrix(~Sum.Ruminants.Gained.Ave.Stand +
                           Max.Dist.Purch.Ave.Stand +
                           f18.Production.15 +
                           f9.Pan.Num.Ave.Stand +
                           f14.Size.Farm.15.Stand +
                           f5.Wildlife.Mix.15 +
                           CumSum.2m.Rain_Oct.Stand,
                         data = df),
    K_eta = ncol(X_eta),

    # Gaussian Process
    Dgp = 2,  # Dimension
    Xgp = locs_standardized, #Values


    prior_only = FALSE # Change this to sample from the priors rather than posterior
  )
}
# Sample from the model.  Use cores as you have available.
hmgp_fit <- sampling(
  hmgp_mod,
  data = standata,
  chains = 4,
  iter = 5000,
  cores = 4)

# Examine the model parameters and diagnostic values.
# Note there are no variable names, just matrix coordinates.
# b_1 and b_2 are coefficients for the multinomial outcome, b_eta for the
# hurdle.
# lp__ is model likelihood, but only relative within samples.

#Now there are a larger number of parameters for theta and eta
#hmgp_fit

#Convert to mcmc object and save
if(standata$prior_only == TRUE){
  hm_mod.mcmc <- As.mcmc.list(hmgp_fit)
  save(hm_mod.mcmc, file = "All_Files_For_Publication/NR_hurdle_GP_mcmc_newParam2prct_PRIORS_ONLY.rdata")
}else{
if(modelThesis == TRUE){
  hm_mod.mcmc <- As.mcmc.list(hmgp_fit)
  save(hm_mod.mcmc, file = "All_Files_For_Publication/NR_hurdle_GP_mcmc_thesis.rdata")
}
if(modelPub == TRUE){
  hm_mod.mcmc <- As.mcmc.list(hmgp_fit)
  save(hm_mod.mcmc, file = "All_Files_For_Publication/NR_hurdle_GP_mcmc_newParam.rdata")
}
if(model2 == TRUE){
  hm_mod.mcmc <- As.mcmc.list(hmgp_fit)
  save(hm_mod.mcmc, file = "All_Files_For_Publication/NR_hurdle_GP_mcmc_newParam2prct.rdata")
}
if(modelPres == TRUE){
  hm_mod.mcmc <- As.mcmc.list(hmgp_fit)
  save(hm_mod.mcmc, file = "All_Files_For_Publication/NR_hurdle_GP_mcmc_newParam_present.rdata")
}
if(modelThesis_Oct == TRUE){
  hm_mod.mcmc <- As.mcmc.list(hmgp_fit)
  save(hm_mod.mcmc, file = "All_Files_For_Publication/NR_hurdle_GP_mcmc_Thesis_Oct.rdata")
}
}

# Dot plots of variables
stan_plot(hmgp_fit, "b_1")
stan_plot(hmgp_fit, "b_2")
stan_plot(hmgp_fit, "b_eta")

# GP parameters
stan_plot(hmgp_fit, c("sdgp", "lscale"))
# GP latent variables (correlated random effects of all points)
stan_plot(hmgp_fit, c("gp_eff"))


# Look at some chain diagnostics. Are they well mixed?
#b_1
bayesplot::mcmc_trace(hmgp_fit, vars(matches("b_1")))
bayesplot::mcmc_dens_overlay(hmgp_fit, vars(matches("b_1")))
# These rank-histogram should be approximately uniform for well-mixed chains.
# See https://arxiv.org/abs/1903.08008
bayesplot::mcmc_rank_hist(hmgp_fit, vars(matches("b_1")))
#b_2
bayesplot::mcmc_trace(hmgp_fit, vars(matches("b_2")))
bayesplot::mcmc_dens_overlay(hmgp_fit, vars(matches("b_2")))
bayesplot::mcmc_rank_hist(hmgp_fit, vars(matches("b_2")))
#b_eta
bayesplot::mcmc_trace(hmgp_fit, vars(matches("b_eta")))
bayesplot::mcmc_dens_overlay(hmgp_fit, vars(matches("b_eta")))
bayesplot::mcmc_rank_hist(hmgp_fit, vars(matches("b_eta")))
#GP
bayesplot::mcmc_rank_hist(hmgp_fit, regex_pars = "(sdgp|lscale)")
# Look at 10 random latent variables
bayesplot::mcmc_rank_hist(hmgp_fit, pars = vars(param_range("gp_eff", sample(nrow(df), 10))))

# We can extract the theta and eta linear predictor matrices.
etas <- extract(hmgp_fit, "eta")[[1]]
thetas <- extract(hmgp_fit, "theta")[[1]]

# And convert to the outcome scale
rift <- plogis(etas)
outcome <- apply(thetas, 1:2, rsoftmax)#1:2 break it up by both farm and sample

# The expected number of deaths and mortalities are the probability of rift
# times the outcome probabilities of each times the number of animals
expected_deaths <- t(rift * outcome[2,,]) * df$TotalPopDomR_2010
expected_abortions <- t(rift * outcome[3,,]) * df$TotalPopDomR_2010

predicted_rift <- (apply(rift, 1, function(.) as.integer(rbernoulli(length(.), .))))
predicted_deaths <-
  apply(outcome[2,,], 1, function(.) rbinom(length(.), df$TotalPopDomR_2010, .)) *
  predicted_rift

predicted_abortions <-
  apply(outcome[3,,], 1, function(.) rbinom(length(.), df$TotalPopDomR_2010, .)) *
          predicted_rift

predicted_nothing <-
  apply(outcome[1,,], 1, function(.) rbinom(length(.), df$TotalPopDomR_2010, .)) *
  predicted_rift

outbreak_order <- order(df$X3.Outbreak)
deaths_order <- order(df$TotalRVFDeathsDomR.2010)
abortion_order <- order(df$TotalRVFAbortionsDomR.2010)

ppc_bars(df$X3.Outbreak[outbreak_order],
              t(predicted_rift[outbreak_order,]))

ppc_intervals(df$TotalRVFDeathsDomR.2010[deaths_order],
              t(predicted_deaths[deaths_order, ]))

ppc_intervals(df$TotalRVFAbortionsDomR.2010[abortion_order],
              t(predicted_abortions[abortion_order, ]))

#Save data
if(standata$prior_only == TRUE){
  hmgp_fit@stanmodel@dso <- new("cxxdso")
  saveRDS(hmgp_fit,  "All_Files_For_Publication/NR_hurdle_GP_newParam2prct_PRIORS_ONLY.rds")
}else{
if(modelThesis == TRUE){
  hmgp_fit@stanmodel@dso <- new("cxxdso")
  saveRDS(hmgp_fit, "All_Files_For_Publication/NR_hurdle_GP_Thesis.rds")}
if(modelPub == TRUE){
  hmgp_fit@stanmodel@dso <- new("cxxdso")
  saveRDS(hmgp_fit,  "All_Files_For_Publication/NR_hurdle_GP_newParam.rds")
}
if(model2 == TRUE){
  hmgp_fit@stanmodel@dso <- new("cxxdso")
  saveRDS(hmgp_fit,  "All_Files_For_Publication/NR_hurdle_GP_newParam2prct.rds")
}
if(modelPres == TRUE){
  hmgp_fit@stanmodel@dso <- new("cxxdso")
  saveRDS(hmgp_fit,  "All_Files_For_Publication/NR_hurdle_GP_newParam_present.rds")
}

if(modelThesis_Oct == TRUE){
  hmgp_fit@stanmodel@dso <- new("cxxdso")
  saveRDS(hmgp_fit,  "All_Files_For_Publication/NR_hurdle_GP_Thesis_Oct.rds")
}
}
