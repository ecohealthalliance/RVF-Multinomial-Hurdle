#'Bayesian Model checking, assessment and plotting
#'author: Mindy Rostal
#'date: November 2, 2020
#'

#Load libraries
library(ggpubr)#note this masks extract from rstan so need to open this library first, then rstan
library(tidyverse)
library(ggmcmc)
library(coda)
library(lattice)
library(stringr)
library(tidybayes)
library(ggplot2)
library(ggstance)
library(dplyr)
library(pROC)
library(egg)
library(bayesplot)
library(rstan)
library(RColorBrewer)


#Select if you are analyzing the random effects model or the model with the Gaussian Process
modelThesis <- FALSE # Model 3 in Table S1
modelPub <- FALSE # Model 4 in Table S1
model2 <- TRUE #Model 1 in Table S1 - Selected Model
modelPres <- FALSE # Model 5 in Table S1
modelThesisOct <- FALSE # Model 6 in Table S1
#Note Model 2 requires the use of the hurdle-multinomial-re.stan and hurdle-multinomial-re_for_publication.r
All.One.2 <- FALSE #Model 7 All variables from model2 in both the hurdle and multinomial linear predictors
ModelRE <- FALSE
ModelGP <- TRUE
PRIORS_ONLY <- FALSE
Intercept_only <- FALSE

#Source analysis
#source("RVF_model_script_MV.R")
source("All_Files_For_Publication/Bayesian Assessment Functions.R")

# This allows us to use the functions we define in stan code, notably rsoftmax()
# Convert the model to C++ code, then compile (Not normally done, just for
# illustration. stan_model() can do this in one step)
if(ModelGP == TRUE){
  hmgp_mod_c <- stanc("All_Files_For_Publication/hurdle-multinomial-gp.stan", verbose = TRUE)
  hmgp_mod <- stan_model(stanc_ret = hmgp_mod_c, verbose = TRUE)
}
if(ModelRE == TRUE){
  hmgp_mod_c <- stanc("All_Files_For_Publication/hurdle-multinomial-re.stan", verbose = TRUE)
  hmgp_mod <- stan_model(stanc_ret = hmgp_mod_c, verbose = TRUE)
}
if(Intercept_only == TRUE){
  hmgp_mod_c <- stanc("All_Files_For_Publication/hurdle-multinomial-intercept-only.stan", verbose = TRUE)
  hmgp_mod <- stan_model(stanc_ret = hmgp_mod_c, verbose = TRUE)
}
if(ModelGP == FALSE & Intercept_only == FALSE & ModelRE == FALSE){
  hmgp_mod_c <- stanc("All_Files_For_Publication/hurdle-multinomial.stan", verbose = TRUE)
  hmgp_mod <- stan_model(stanc_ret = hmgp_mod_c, verbose = TRUE)
}
if(ModelRE == FALSE & Intercept_only == FALSE & ModelGP == FALSE){
  hmgp_mod_c <- stanc("All_Files_For_Publication/hurdle-multinomial.stan", verbose = TRUE)
  hmgp_mod <- stan_model(stanc_ret = hmgp_mod_c, verbose = TRUE)
}
expose_stan_functions(hmgp_mod)
rsoftmax(c(1,1,1))
sum(rsoftmax(c(0.2, 2, 3)))

#Modify data
mydata<- read.csv("All_Files_For_Publication/data/Updated 2017 and 2015 Cross-sectional Data for Bayesian Analsysis_with_UTM.csv", header=T)
mydata<- mydata%>%
  rename(Y = UTM_Lat)%>%
  rename(X = UTM_Long)
#Add column to indicate if there were any deaths or abortions
mydata <- mydata%>%
  mutate(RVF_Anydeaths = if_else(TotalRVFDeathsDomR.2010.Stand >1, 1, 0))%>%
  mutate(RVF_Anyaborts = if_else(TotalRVFAbortionsDomR.2010.Stand >1, 1, 0))%>%
  mutate(ID = seq(1,120,1))%>%
  mutate(ID = str_pad(ID, 3, "left", pad = "0"))%>%
  mutate_at(vars(ID), list(~as.character(.)))
#relevel so that the small holders are the reference group, then turn into an integer
prod_dat <- as.factor(mydata$f18.Production.15)
prod_dat <- relevel(prod_dat, ref="Small holder")
mydata$f18.Production.15 <- prod_dat


WAIC_comp <- read.csv("./Publication_Figures/WAIC comparison between RE and GP models_for_publication.csv", colClasses=c("Date"="Date"))
Dev_Comp <- read.csv("Publication_Figures/Comparison of Deviance between Models_for_publication.csv", stringsAsFactors=F)

if(Intercept_only == TRUE){
  #Only thing we will calculate for the intercept only is deviance, skip the rest of the anlayses
  #load data
  load("All_Files_For_Publication/NR_hurdle_intercept_only_mcmc.rdata")
  hm_fit <- readRDS("All_Files_For_Publication/NR_hurdle_Intercept_only.rds")
  #Calculate the deviance from the LL
  LL.ggs <- ggs(hm_mod.mcmc_int, family = "log_lik")
  LL.est <- mean_qi(LL.ggs$value, .width = .95)
  LL <- LL.est[[1]]
  LL.lo <- LL.est[[2]]
  LL.hi <- LL.est[[3]]

  Deviance.int <- -2*LL
  Deviance.int.hi <- -2*LL.lo
  Deviance.int.lo <- -2*LL.hi

  Dev_Comp[4,2] <- round(Deviance.int,4)
  Dev_Comp[4,3] <- round(Deviance.int.lo,4)
  Dev_Comp[4,4] <- round(Deviance.int.hi,4)
  write.csv(Dev_Comp, "Publication_Figures/Comparison of Deviance between Models.csv", row.names = FALSE)
}else{

  if(ModelRE == FALSE & Intercept_only == FALSE & ModelGP == FALSE){
    #load data
    load("All_Files_For_Publication/NR_hurdle_mcmc_no_spatial_effect.rdata")
    hm_fit <- readRDS("All_Files_For_Publication/NR_hurdle_no_spatial_effect.rds")
    #File name to save plots and data
    #Make model name
    model.file.name <- "All_Files_For_Publication/NR_hurdle_no_spatial_effect_2prct_for_publication"
    #Calculate the deviance from the LL
    LL.ggs <- ggs(hm_mod.mcmc, family = "log_lik")
    LL.est <- mean_qi(LL.ggs$value, .width = .95)
    LL <- LL.est[[1]]
    LL.lo <- LL.est[[2]]
    LL.hi <- LL.est[[3]]

    Deviance.mod <- -2*LL
    Deviance.mod.hi <- -2*LL.lo
    Deviance.mod.lo <- -2*LL.hi

    Dev_Comp[1,2] <- round(Deviance.mod,4)
    Dev_Comp[1,3] <- round(Deviance.mod.lo,4)
    Dev_Comp[1,4] <- round(Deviance.mod.hi,4)
    write.csv(Dev_Comp, "Publication_Figures/Comparison of Deviance between Models_for_publication.csv", row.names = FALSE)
  }else{

    #Name files
    if(ModelRE == TRUE){
      load("All_Files_For_Publication/NR_hurdle_RE_mcmc_for_publication.rdata")
      hm_fit <- readRDS("All_Files_For_Publication/NR_hurdle_RE_for_publication.rds")
      #File name to save plots and data
      #Make model name
      model.file.name <- "NR_hurdle_RE_for_publication_Pubmodel_2prct"
    }else{
      if(ModelGP == TRUE){
        if(model2 == TRUE){
          load("All_Files_For_Publication/NR_hurdle_GP_mcmc_newParam2prct.rdata")
          hm_fit_priors.only <- readRDS("All_Files_For_Publication/NR_hurdle_GP_newParam2prct_PRIORS_ONLY.rds")
          hm_fit <- readRDS("All_Files_For_Publication/NR_hurdle_GP_newParam2prct.rds")
          #Make model name to save plots and data
          model.file.name <- "NR_hurdle_GP_for_publication_Pubmodel_2prct"
        }
        if(All.One.2 == TRUE){
          load("All_Files_For_Publication/NR_hurdle_GP_mcmc_AllOne2prct.rdata")
          hm_fit_priors.only <- readRDS("All_Files_For_Publication/NR_hurdle_GP_mcmc_AllOne2prct_PRIORS_ONLY.rds")
          hm_fit <- readRDS("All_Files_For_Publication/NR_hurdle_GP_mcmc_AllOne2prct.rds")
          #Make model name to save plots and data
          model.file.name <- "NR_hurdle_GP_for_publication_AllOne_2prct"
        }
        if(PRIORS_ONLY == TRUE){
          load("All_Files_For_Publication/NR_hurdle_GP_mcmc_newParam2prct_PRIORS_ONLY.rdata")
          model.file.name <- "NR_hurdle_GP_for_publication_PRIORS_ONLY_2prct"
        }
        if(modelPub == TRUE){
          load("All_Files_For_Publication/NR_hurdle_GP_mcmc_newParam.rdata")
          hm_fit_priors.only <- readRDS("All_Files_For_Publication/NR_hurdle_GP_priors_only_newParam.rds")
          hm_fit <- readRDS("All_Files_For_Publication/NR_hurdle_GP_newParam.rds")
          #File name to save plots and data
          #Make model name
          model.file.name <- "NR_hurdle_GP_for_publication_Pubmodel"
        }

        if(modelThesis == TRUE){
          load("All_Files_For_Publication/NR_hurdle_GP_mcmc_Thesis.rdata")
          hm_fit_priors.only <- readRDS("All_Files_For_Publication/NR_hurdle_GP_priors_only.rds")
          hm_fit <- readRDS("All_Files_For_Publication/NR_hurdle_GP_Thesis.rds")
          #File name to save plots and data
          #Make model name
          model.file.name <- "NR_hurdle_GP_for_publication_Thesismodel"
        }
        if(modelPres == TRUE){
          load("All_Files_For_Publication/NR_hurdle_GP_mcmc_newParam_present.rdata")
          hm_fit_priors.only <- readRDS("All_Files_For_Publication/NR_hurdle_GP_priors_only.rds")
          hm_fit <- readRDS("All_Files_For_Publication/NR_hurdle_GP_newParam_present.rds")
          #File name to save plots and data
          #Make model name
          model.file.name <- "NR_hurdle_GP_for_publication_NewParamPresence"
        }
        if(modelThesisOct == TRUE){
          load("All_Files_For_Publication/NR_hurdle_GP_mcmc_Thesis_Oct.rdata")
          hm_fit_priors.only <- readRDS("All_Files_For_Publication/NR_hurdle_GP_priors_only.rds")
          hm_fit <- readRDS("All_Files_For_Publication/NR_hurdle_GP_Thesis_Oct.rds")
          #File name to save plots and data
          #Make model name
          model.file.name <- "NR_hurdle_GP_for_publication_Thesis_Oct"
        }
      }
    }

    #Date for file names
    datestamp <- Sys.Date()

    #Set ggplot theme
    pub_theme <-   theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"), text = element_text(size=20, color = "black"),
            axis.text.x = element_text(color = "black"), axis.text.y = element_text(color = "black"),
            legend.position="none")
    #############################################
    #Prepare dataframes
    hm_mod.ggs <- ggs(hm_mod.mcmc)

    #all params
    if(model2 == TRUE | ModelRE == TRUE | PRIORS_ONLY == TRUE){
      params.ggs <- hm_mod.ggs%>%
        filter(  Parameter == "b_1[1]"   | Parameter == "b_1[2]"   | Parameter == "b_1[3]"   | Parameter == "b_1[4]"   | Parameter == "b_1[5]"   | Parameter == "b_1[6]"   | Parameter == "b_1[7]"   | Parameter == "b_1[8]"   |
                 Parameter == "b_2[1]"   | Parameter == "b_2[2]"   | Parameter == "b_2[3]"   | Parameter == "b_2[4]"   | Parameter == "b_2[5]"   | Parameter == "b_2[6]"   | Parameter == "b_2[7]"   | Parameter == "b_2[8]"   |
                 Parameter == "b_eta[1]" | Parameter == "b_eta[2]" | Parameter == "b_eta[3]" | Parameter == "b_eta[4]" | Parameter == "b_eta[5]" | Parameter == "b_eta[6]" | Parameter == "b_eta[7]" | Parameter == "b_eta[8]" | Parameter == "b_eta[9]")
    }

    #all params
    if(All.One.2 == TRUE){
      params.ggs <- hm_mod.ggs%>%
        filter(  Parameter == "b_1[1]"   | Parameter == "b_1[2]"   | Parameter == "b_1[3]"   | Parameter == "b_1[4]"   | Parameter == "b_1[5]"   | Parameter == "b_1[6]"    | Parameter == "b_1[7]"   | Parameter == "b_1[8]"   | Parameter == "b_1[9]"   | Parameter == "b_1[10]"   | Parameter == "b_1[11]"  | Parameter == "b_1[12]"  | Parameter == "b_1[13]"  | Parameter == "b_1[14]"   | Parameter == "b_1[15]"   | Parameter == "b_1[16]" |
                 Parameter == "b_2[1]"   | Parameter == "b_2[2]"   | Parameter == "b_2[3]"   | Parameter == "b_2[4]"   | Parameter == "b_2[5]"   | Parameter == "b_2[6]"    | Parameter == "b_2[7]"   | Parameter == "b_2[8]"   | Parameter == "b_1[9]"   | Parameter == "b_2[10]"   | Parameter == "b_2[11]"  | Parameter == "b_2[12]"  | Parameter == "b_2[13]"  | Parameter == "b_2[14]"   | Parameter == "b_2[15]"   | Parameter == "b_2[16]" |
                 Parameter == "b_eta[1]" | Parameter == "b_eta[2]" | Parameter == "b_eta[3]" | Parameter == "b_eta[4]" | Parameter == "b_eta[5]" | Parameter == "b_eta[6]"  | Parameter == "b_eta[7]" | Parameter == "b_eta[8]" | Parameter == "b_eta[9]" | Parameter == "b_eta[10]" | Parameter == "b_eta[11]"| Parameter == "b_eta[12]"| Parameter == "b_eta[13]"| Parameter == "b_eta[14]" | Parameter == "b_eta[15]" | Parameter == "b_eta[16]")
    }

    #Eta
    eta.ggs <- ggs(hm_mod.mcmc, family = "b_eta")

    eta.ggs <- eta.ggs%>%
      mutate(Param = str_sub(Parameter, -2, -2))%>%
      mutate(Param = paste0("RVF_", Param))%>%
      select(-Parameter)

    #death
    death.ggs <- ggs(hm_mod.mcmc, family = "b_1")

    death.ggs <- death.ggs%>%
      mutate(Param = str_sub(Parameter, -2, -2))%>%
      mutate(Param = paste0("d_", Param))%>%
      select(-Parameter)

    #abortion
    abort.ggs <- ggs(hm_mod.mcmc, family = "b_2")

    abort.ggs <- abort.ggs%>%
      mutate(Param = str_sub(Parameter, -2, -2))%>%
      mutate(Param = paste0("a_", Param))%>%
      select(-Parameter)

    if(All.One.2 == TRUE){
      #Eta
      eta.ggs <- ggs(hm_mod.mcmc, family = "b_eta")

      eta.ggs <- eta.ggs%>%
        mutate(Param = str_sub(Parameter, -3, -2))%>%
        mutate(Param = if_else(str_detect(Param, "\\["), str_sub(Param, -1), Param))%>%
        mutate(Param = paste0("RVF_", Param))%>%
        select(-Parameter)

      #death
      death.ggs <- ggs(hm_mod.mcmc, family = "b_1")

      death.ggs <- death.ggs%>%
        mutate(Param = str_sub(Parameter, -3, -2))%>%
        mutate(Param = if_else(str_detect(Param, "\\["), str_sub(Param, -1), Param))%>%
        mutate(Param = paste0("d_", Param))%>%
        select(-Parameter)

      #abortion
      abort.ggs <- ggs(hm_mod.mcmc, family = "b_2")

      abort.ggs <- abort.ggs%>%
        mutate(Param = str_sub(Parameter, -3, -2))%>%
        mutate(Param = if_else(str_detect(Param, "\\["), str_sub(Param, -1), Param))%>%
        mutate(Param = paste0("a_", Param))%>%
        select(-Parameter)
    }


    Diag.file.name <- paste("Publication_Figures/", model.file.name, "Density_running_and_Traceplots_of_Parameters", datestamp, ".pdf", sep ="_")

    #Convergence and diagnostics
    #Save PDF of Density, running and traceplots
    if(PRIORS_ONLY == FALSE){
      ggmcmc(params.ggs, plot = c("traceplot", "density", "running", "autocorrelation", "crosscorrelation"), file = Diag.file.name)

      #Potential Scale Reduction Factor (R̂ , “Rhat”). It is a weighted average of the Between-chain and Within-chain variances.The most prevalent version, though is the one in the second version of the second edition of Bayesian Data Analysis (let’s call this the “BDA2” version), this is the default for ggmcmc (http://blog.xavier-fim.net/2019/03/comparison-of-rhat-versions-clarifying-formulas-for-the-potential-scale-reduction-factor-and-its-implications/)
      ggs_Rhat(params.ggs)

      ggs_geweke(params.ggs)

      #raftery.diag(hm_mod.mcmc, q=0.025, r=0.005, s=0.95, converge.eps=0.001)

      head(effectiveSize(hm_mod.mcmc),25)

      ############Residuals
      df_resid_theta <- data.frame(yd = mydata$TotalRVFDeathsDomR.2010,
                                   ya = mydata$TotalRVFAbortionsDomR.2010,
                                   y_hat = summary(hm_fit, pars = 'theta_hat')$summary[, 'mean'],
                                   residual_d = summary(hm_fit, pars = 'resid_theta_2')$summary[, 'mean'],
                                   residual_a = summary(hm_fit, pars = 'resid_theta_3')$summary[, 'mean'])


      df_resid_theta$names<- rownames(df_resid_theta)

      df_resid_theta <- df_resid_theta%>%
        mutate(Theta_hat =str_sub(names, -4,-2))%>%
        mutate(Farm = str_sub(Theta_hat, 1,1))%>%
        mutate(Outcome = str_sub(Theta_hat, -1, -1))%>%
        select(-names, -Theta_hat)

      death_resid <- filter(df_resid_theta, Outcome == "2")
      abort_resid <- filter(df_resid_theta, Outcome == "3")

      UTM <- mydata%>%
        select(X, Y)%>%
        mutate(Farm = seq(1, nrow(mydata), by = 1))%>%
        mutate_at(vars(Farm), list(~as.character(.)))

      #Residual plots
      #Death
      resid_plot_d <- ggplot(death_resid , aes(x = y_hat, y = residual_d, color = residual_d)) +
        geom_point(alpha = .75) +
        scale_color_continuous(low = 'darkred', high = 'lightsalmon') +
        geom_hline(yintercept = 0) +
        labs(title = "Plot of Death Residuals", x = "Model Estimate", y = "Residuals")+
        pub_theme+
        theme(legend.position = 'right',
              plot.background = element_blank(),
              panel.border = element_blank())

      #resid_plot_d #We expect clustering around zero due to the hurdle

      pred_plot <- ggplot(death_resid, aes(x = y_hat, y = yd, color = residual_d)) +
        geom_point(alpha = .75) +
        scale_color_continuous(low = 'darkred', high = 'lightsalmon') +
        geom_abline(slope = 1, intercept = 0) +
        labs(title = "Plot of Death Residuals", x = "Model Estimate", y = "Actual Number of Deaths")+
        pub_theme+
        theme(legend.position = 'right',
              plot.background = element_blank(),
              panel.border = element_blank())

      #pred_plot
      #Map residuals
      df_resid_death <- full_join(death_resid, UTM, by = "Farm")

      resid_plot_death_map <- ggplot(df_resid_death, aes(x = X, y = Y , color = residual_d))+
        geom_point()+
        labs( x = "UTM Longitude", y = "UTM Latitude",  color='Residuals of Number \nof Deaths')+#title = "Map of Eta Residuals",
        pub_theme+
        theme(legend.position = 'right',
              plot.background = element_blank(),
              panel.border = element_blank(),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

      #resid_plot_death_map

      D_eta_resid_h <- ggplot(death_resid, aes(x = residual_d))+
        geom_histogram()+
        labs(title = "Histogram of Death Residuals", x = "Residuals", y = "Count")

      #D_eta_resid_h

      #Abortion - theta
      resid_plot_a <- ggplot(abort_resid , aes(x = y_hat, y = residual_a, color = residual_a)) +
        geom_point(alpha = .75) +
        scale_color_continuous(low = 'darkred', high = 'lightsalmon') +
        geom_hline(yintercept = 0) +
        labs(title = "Plot of Abortion Residuals", x = "Model Estimate", y = "Residuals")+
        pub_theme+
        theme(legend.position = 'right',
              plot.background = element_blank(),
              panel.border = element_blank())

      #resid_plot_a
      A_eta_resid_h <- ggplot(abort_resid, aes(x = residual_a))+
        geom_histogram()+
        labs(title = "Histogram of Abortion Residuals", x = "Residuals", y = "Count")

      #A_eta_resid_h
      #Map residuals
      df_abort_resid <- full_join(abort_resid, UTM, by = "Farm")

      resid_plot_abort_map <- ggplot(df_abort_resid, aes(x = X, y = Y , color = residual_a))+
        geom_point()+
        labs( x = "UTM Longitude", y = "UTM Latitude", color='Residuals of Number \nof Abortions')+#title = "Map of Eta Residuals",
        pub_theme+
        theme(legend.position = 'right',
              plot.background = element_blank(),
              panel.border = element_blank(),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

      #resid_plot_abort_map

      #Save plots to
      if(model2 == TRUE){
        Fig_S6_resid_map<- ggarrange(resid_plot_death_map, resid_plot_abort_map, nrow = 2, ncol = 1, labels = c("(a)", "(b)"))
        FigAbort.title <- paste("Publication_Figures/Fig S7 Spatial Struture of Residuals ", model.file.name, ".png", sep = "")
        ggexport(Fig_S6_resid_map, filename = FigAbort.title, width = 812, height = 900)
      }

      #Eta
      df_resid_eta <- data.frame(y = mydata$X3.Outbreak,
                                 y_hat = summary(hm_fit, pars = 'eta_hat')$summary[, 'mean'],
                                 residual = summary(hm_fit, pars = 'resid_eta')$summary[, 'mean'])

      df_resid_eta$names<- rownames(df_resid_eta)

      df_resid_eta <- df_resid_eta%>%
        mutate(Farm = str_sub(names, -3,-2))%>%
        mutate(Farm = str_replace(Farm, "\\[", ""))

      df_resid_eta <- full_join(df_resid_eta, UTM, by = "Farm")

      resid_plot_eta <- ggplot(df_resid_eta , aes(x = y_hat, y = residual, color = residual)) +
        geom_point(alpha = .75) +
        scale_color_continuous(low = 'darkred', high = 'lightsalmon') +
        labs(title = "Plot of Eta Residuals", x = "Model Estimate", y = "Residuals")+
        geom_hline(yintercept = 0) +
        pub_theme+
        theme(legend.position = 'right',
              plot.background = element_blank(),
              panel.border = element_blank())

      #resid_plot_eta

      h_eta_resid_h <- ggplot(df_resid_eta, aes(x = residual))+
        geom_histogram()+
        labs(title = "Histogram of Eta Residuals", x = "Residuals", y = "Count")

      #h_eta_resid_h



      #######################################################################################
      ##############################################################################
      ##################################################################
      #Posterior predictive checks

      # Extract the theta and eta linear predictor matrices.
      etas <- rstan::extract(hm_fit, "eta")[[1]]
      thetas <- rstan::extract(hm_fit, "theta")[[1]]

      #Set order of farms/outbreak responses
      outbreak_order <- order(mydata$X3.Outbreak)

      # And convert to the outcome scale
      rift <- plogis(etas)
      outcome <- apply(thetas, 1:2, rsoftmax)
      if(ModelGP == TRUE & model2 == TRUE){
        etas_noGP <- rstan::extract(hm_fit, "eta_noGP")[[1]]
        rift_noGP <- plogis(etas_noGP)
        predicted_rift_noGP <- (apply(rift_noGP, 1, function(.) as.integer(rbernoulli(length(.), .))))
        prob_rifta_noGP <- apply(rift_noGP, 2, function(.) mean(.))
        RVFprob_vec_noGp <- as.data.frame(prob_rifta_noGP[outbreak_order])
      }

      ############################
      ##Hurdle part of model - eta
      ############################
      #RVF occurred

      #####################################################################
      ####################################################################
      #Assess prediction capability with GP
      predicted_rift <- (apply(rift, 1, function(.) as.integer(rbernoulli(length(.), .))))


      ################################
      ################################
      #Abortions and Deaths - including hurdle model
      # The expected number of deaths and mortalities are the probability of rift
      # times the outcome probabilities of each times the number of animals
      expected_deaths <- t(rift * outcome[2,,]) * mydata$TotalPopDomR_2010
      expected_abortions <- t(rift * outcome[3,,]) * mydata$TotalPopDomR_2010

      #Predicted number of deaths without the hurdle
      predicted_deaths <-
        apply(outcome[2,,], 1, function(.) rbinom(length(.), mydata$TotalPopDomR_2010, .))

      #Predicted number of deaths including the hurdle
      predicted_deaths_HM <-
        apply(outcome[2,,], 1, function(.) rbinom(length(.), mydata$TotalPopDomR_2010, .)) *
        predicted_rift

      #Predicted number of abortions without the hurdle
      predicted_abortions<-
        apply(outcome[3,,], 1, function(.) rbinom(length(.), mydata$TotalPopDomR_2010, .))

      #Predicted number of abortions including the hurdle
      predicted_abortions_HM <-
        apply(outcome[3,,], 1, function(.) rbinom(length(.), mydata$TotalPopDomR_2010, .)) *
        predicted_rift

      #Set order if needed
      deaths_order <- order(mydata$TotalRVFDeathsDomR.2010)
      abortion_order <- order(mydata$TotalRVFAbortionsDomR.2010)
      Farm_ID_order <- order(mydata$Farm_ID)

      #PPC for death
      #Change the predicted deaths to a dataframe
      n.death <- as.data.frame(as.data.frame(t(predicted_deaths)))
      nums <- seq(1,120, 1)#Number farms 1-120
      nums <- str_pad(nums, 2, "left", pad = "0")
      ids <- paste0("Farm_", nums)#Farm Numbers

      names(n.death) <- ids # label columns with Farm Numbers

      #transform the dataset to long
      Num.Death <- gather(n.death, "Farm_Num", "Estimate", Farm_01:Farm_120)

      #Histogram of the number of deaths
      hist_Deaths <- ggplot(Num.Death, aes(x = Estimate ))+
        geom_histogram()

      #hist_Deaths

      #Deaths
      ###No hurdle
      #Create summary data get median and HDI
      n_death.sum <- Num.Death%>%
        group_by(Farm_Num)%>%
        summarise(med_num_d = median(Estimate),#mean estimate of number of deaths
                  pdeath_hdi.low.1 = unlist(median_hdi(Estimate, 2.5)[["ymin"]][1]),
                  pdeath_hdi.low.2.ifbimodal = unlist(median_hdi(Estimate, 2.5)[["ymin"]][2]),
                  pdeath_hdi.hi.1 = unlist(median_hdi(Estimate, 97.5)[["ymax"]][1]),
                  pdeath_hdi.hi.2.ifbimodal = unlist(median_hdi(Estimate, 97.5)[["ymax"]][2]))

      #remove empty rows if the hdi isn't bimodal
      if(all(is.na(n_death.sum$pdeath_hdi.hi.2.ifbimodal))){
        n_death.sum <- select(n_death.sum, -pdeath_hdi.hi.2.ifbimodal)
      }
      if(all(is.na(n_death.sum$pdeath_hdi.low.2.ifbimodal))){
        n_death.sum <- select(n_death.sum, -pdeath_hdi.low.2.ifbimodal)
      }

      #Extract actual data on outbreak, number and percent of deaths
      data_deaths <- mydata%>%
        select(X3.Outbreak, TotalRVFDeathsDomR.2010, PrctRVFDeath.2010, TotalPopDomR_2010, Farm_ID)%>%
        mutate(Farm_Num = seq(1,120, 1))%>%#Make the IDs the same
        mutate(Farm_Num = str_pad(Farm_Num, 2, pad = "0"))%>%
        mutate(Farm_Num = paste("Farm", Farm_Num, sep = "_"))#Label with same IDs/Numbers

      #Join the dataframes of the predicted and actual deaths
      n_death.dat <- full_join(n_death.sum, data_deaths)

      #select the columns of the medians we want to add to the dataframe with the iterations
      ndmeans <- select(n_death.dat, Farm_Num, med_num_d, TotalRVFDeathsDomR.2010, TotalPopDomR_2010, Farm_ID)

      #Join the two databases, the medians with the raw data by farm Number
      nded <-full_join(Num.Death, ndmeans, by = "Farm_Num")

      #Arrange based on the median number of deaths predicted
      nded <- nded %>%
        arrange(med_num_d)

      #Set the levels so that ggplot doesn't change them
      nded$Farm_Num <- factor(nded$Farm_Num, levels=unique(nded$Farm_Num))

      #Boxplots
      #No hurdle model
      FigDead1 <- ggplot(nded, aes(x = Farm_Num, y = Estimate))+
        geom_boxplot(fill = "#1B9E77", alpha = .4, outlier.color = "#666666", outlier.alpha = .3, outlier.size = .5)+
        stat_summary(fun="median", geom="point", size=.3, color="#1B9E77")+
        geom_hline(yintercept = 0, color = "black") +
        geom_point(aes(y = TotalRVFDeathsDomR.2010), color  = "black", size = 1.2)+
        labs(x = "Individual Farms", y = "Number of 2010 Deaths \n(Unconditional on RVF \nPresence on Farm)")+
        pub_theme+
        theme(#text = element_text(size=16),
          axis.text.x = element_text(size = 4, colour = "black", angle = 90)
        )

      #FigDead1

      #Model With hurdle
      #Comments are same as above, except starting from the predicted deaths multiplied by the result of the hurdle model
      predicted_deaths_HM <- as.data.frame(t(as.data.frame(predicted_deaths_HM)))

      names(predicted_deaths_HM) <- ids

      Pred.Dead.HM <- gather(predicted_deaths_HM, "Farm_Num", "Estimate", Farm_01:Farm_120)

      p_death.sum.hm <- Pred.Dead.HM%>%
        group_by(Farm_Num)%>%
        summarise(med_num_d = median(Estimate),
                  ndeath_hdi.low.1 = unlist(median_hdi(Estimate, 2.5)[["ymin"]][1]),
                  ndeath_hdi.low.2.ifbimodal = unlist(median_hdi(Estimate, 2.5)[["ymin"]][2]),
                  ndeath_hdi.hi.1 = unlist(median_hdi(Estimate, 97.5)[["ymax"]][1]),
                  ndeath_hdi.hi.2.ifbimodal = unlist(median_hdi(Estimate, 97.5)[["ymax"]][2]))

      if(all(is.na(p_death.sum.hm$ndeath_hdi.hi.2.ifbimodal))){
        p_death.sum.hm <- select(p_death.sum.hm, -ndeath_hdi.hi.2.ifbimodal)
      }
      if(all(is.na(p_death.sum.hm$ndeath_hdi.low.2.ifbimodal))){
        p_death.sum.hm <- select(p_death.sum.hm, -ndeath_hdi.low.2.ifbimodal)
      }

      n_death.dat.hm <- full_join(p_death.sum.hm, data_deaths)

      ndmeans <- select(n_death.dat.hm, Farm_Num, med_num_d, TotalRVFDeathsDomR.2010, TotalPopDomR_2010, "Farm_ID")#, RVFYN)

      nded.hm <-full_join(Pred.Dead.HM, ndmeans, by = "Farm_Num")

      #Want to set the outliers for the estimates that the hurdle says are zero to zero so that they don't show up on the plot
      nded.hm.outs <- nded.hm %>%
        group_by(Farm_Num)%>%
        summarise(med = median(Estimate))

      nded.hm <- full_join(nded.hm, nded.hm.outs)

      nded.hm <- nded.hm %>%
        mutate(Estimate = if_else(med == 0, as.integer(0), Estimate))%>%
        arrange(med_num_d)

      nded.hm$Farm_Num <- factor(nded.hm$Farm_Num, levels=unique(nded.hm$Farm_Num))

      FigDead2 <- ggplot(nded.hm, aes(x = Farm_Num, y = Estimate))+
        geom_boxplot(fill = "#1B9E77", alpha = 0.4, outlier.color = "#666666", outlier.alpha = .3, outlier.size = .5)+
        stat_summary(fun="median", geom="point", size=.3, color="#1B9E77")+
        geom_hline(yintercept = 0, color = "black") +
        geom_point(aes(y = TotalRVFDeathsDomR.2010), color  = "black", size = 1.2)+
        labs(x = "Individual Farms", y = "Number of 2010 Deaths \n(Conditional on RVF \nPresence on Farm)")+
        pub_theme+
        theme(#text = element_text(size=16),
          axis.text.x = element_text(size = 4, angle = 90))

      #FigDead2
      #Plot the accuracy - with hurdle model
      FigDead3 <- ggplot(n_death.dat.hm, aes(x = med_num_d, y =  TotalRVFDeathsDomR.2010))+
        geom_point()+
        labs(x = "Predicted Number of Deaths", y = "Reported Number of \nDeaths During 2010")+
        pub_theme+
        theme(#text = element_text(size=16),
          axis.text.x = element_text(size = 4, angle = 90))

      #FigDead3

      #Save plots to
      if(model2 == TRUE){
        Fig_pDead <- ggarrange(FigDead1, FigDead2, FigDead3, nrow = 3, ncol = 1, labels = c("(a)", "(b)", "(c)"))

        FigDead.title <- paste("Publication_Figures/Fig S4 PPC of Death - A excludes hurdle ", model.file.name, ".png", sep = "")
        ggexport(Fig_pDead, filename = FigDead.title, width = 812, height = 900)
      }

      ##############################################
      #Abortions
      #The comments are the same as given for the death code
      n.abort <- as.data.frame(as.data.frame(t(predicted_abortions)))
      nums <- seq(1,120, 1)
      nums <- str_pad(nums, 2, "left", pad = "0")
      ids <- paste0("Farm_", nums)

      names(n.abort) <- ids

      Num.Abort <- gather(n.abort, "Farm_Num", "Estimate", Farm_01:Farm_120)

      #Histogram of the probability of abortions
      hist_Aborts <- ggplot(Num.Abort, aes(x = Estimate ))+
        geom_histogram()

      #hist_Aborts

      #Abortions
      ###No hurdle
      n_abort.sum <- Num.Abort%>%
        group_by(Farm_Num)%>%
        summarise(med_num_a = median(Estimate),
                  pdeath_hdi.low.1 = unlist(mean_hdi(Estimate, 2.5)[["ymin"]][1]),
                  pdeath_hdi.low.2.ifbimodal = unlist(mean_hdi(Estimate, 2.5)[["ymin"]][2]),
                  pdeath_hdi.hi.1 = unlist(mean_hdi(Estimate, 97.5)[["ymax"]][1]),
                  pdeath_hdi.hi.2.ifbimodal = unlist(mean_hdi(Estimate, 97.5)[["ymax"]][2]))

      if(all(is.na(n_abort.sum$pdeath_hdi.hi.2.ifbimodal))){
        n_abort.sum <- select(n_abort.sum, -pdeath_hdi.hi.2.ifbimodal)
      }
      if(all(is.na(n_abort.sum$pdeath_hdi.low.2.ifbimodal))){
        n_abort.sum <- select(n_abort.sum, -pdeath_hdi.low.2.ifbimodal)
      }

      data_abort <- mydata%>%
        select(X3.Outbreak, TotalRVFAbortionsDomR.2010, PrctRVFAbort.2010, TotalPopDomR_2010, Farm_ID)%>%
        mutate(Farm_Num = seq(1,120, 1))%>%
        mutate(Farm_Num = str_pad(Farm_Num, 2, pad = "0"))%>%
        mutate(Farm_Num = paste("Farm", Farm_Num, sep = "_"))

      n_abort.dat <- full_join(n_abort.sum, data_abort)

      n_abort.dat <- n_abort.dat%>%
        mutate(above0 = if_else(med_num_a >=1 & TotalRVFAbortionsDomR.2010 >0, "TP", if_else(med_num_a >=1 & TotalRVFAbortionsDomR.2010 == 0, "FP", if_else(med_num_a <1 & TotalRVFAbortionsDomR.2010 >0, "FN", "TN" ))))%>%
        mutate(AnyAbortReported = if_else(TotalRVFAbortionsDomR.2010 >0, 1, 0))%>%
        mutate(AnyAbortPredicted = if_else(med_num_a >=1, 1, 0))

      nameans <- select(n_abort.dat, Farm_Num, med_num_a, TotalRVFAbortionsDomR.2010, TotalPopDomR_2010, Farm_ID)#, RVFYN)

      nabort <-full_join(Num.Abort, nameans, by = "Farm_Num")

      nabort <- nabort %>%
        arrange(med_num_a)

      nabort$Farm_Num <- factor(nabort$Farm_Num, levels=unique(nabort$Farm_Num))

      #Boxplots
      #Plot without hurdle model
      FigAbort1 <- ggplot(nabort, aes(x = Farm_Num, y = Estimate))+
        geom_boxplot(fill = "#7570B3", alpha = 0.4, outlier.color = "#666666", outlier.alpha = .3, outlier.size = .5)+
        stat_summary(fun="median", geom="point", size=.3, color="#7570B3")+
        geom_hline(yintercept = 0, color = "black") +
        geom_point(aes(y = TotalRVFAbortionsDomR.2010), color  = "black", size = 1.2)+
        labs(x = "Individual Farms", y = "Number of 2010 Abortions \n(Unconditional on RVF \nPresence on Farm)")+
        pub_theme+
        theme(#text = element_text(size=16, colour = "black"),
          axis.text.x = element_text(size = 4, angle = 90))

      #FigAbort1
      #Model With hurdle
      predicted_abortions_HM  <- as.data.frame(t(as.data.frame(predicted_abortions_HM )))

      names(predicted_abortions_HM ) <- ids

      Pred.Abort.HM <- gather(predicted_abortions_HM , "Farm_Num", "Estimate", Farm_01:Farm_120)

      n_Abort.sum.hm <- Pred.Abort.HM%>%
        group_by(Farm_Num)%>%
        summarise(med_num_a = median(Estimate),
                  ndeath_hdi.low.1 = unlist(mean_hdi(Estimate, 2.5)[["ymin"]][1]),
                  ndeath_hdi.low.2.ifbimodal = unlist(mean_hdi(Estimate, 2.5)[["ymin"]][2]),
                  ndeath_hdi.hi.1 = unlist(mean_hdi(Estimate, 97.5)[["ymax"]][1]),
                  ndeath_hdi.hi.2.ifbimodal = unlist(mean_hdi(Estimate, 97.5)[["ymax"]][2]))

      if(all(is.na(n_Abort.sum.hm$ndeath_hdi.hi.2.ifbimodal))){
        n_Abort.sum.hm <- select(n_Abort.sum.hm, -ndeath_hdi.hi.2.ifbimodal)
      }
      if(all(is.na(p_death.sum.hm$n_Abort.sum.hm))){
        n_Abort.sum.hm <- select(n_Abort.sum.hm, -ndeath_hdi.low.2.ifbimodal)
      }

      n_abort.dat.hm <- full_join(n_Abort.sum.hm, data_abort)

      n_abort.dat.hm <- n_abort.dat.hm%>%
        mutate(above0 = if_else(med_num_a >=1 & TotalRVFAbortionsDomR.2010 >0, "TP", if_else(med_num_a >=1 & TotalRVFAbortionsDomR.2010 == 0, "FP", if_else(med_num_a <1 & TotalRVFAbortionsDomR.2010 >0, "FN", "TN" ))))%>%
        mutate(AnyAbortReported = if_else(TotalRVFAbortionsDomR.2010 >0, 1, 0))%>%
        mutate(AnyAbortPredicted = if_else(med_num_a >=1, 1, 0))

      nameans <- select(n_abort.dat.hm, Farm_Num, med_num_a, TotalRVFAbortionsDomR.2010, TotalPopDomR_2010, Farm_ID)#, RVFYN)
      nabort.hm <-full_join(Pred.Abort.HM, nameans, by = "Farm_Num")

      #Want to set the outliers for the estimates that the hurdle says are zero to zero so that they don't show up on the plot
      nabort.hm.outs <- nabort.hm %>%
        group_by(Farm_Num)%>%
        summarise(med = median(Estimate))

      nabort.hm <- full_join(nabort.hm, nabort.hm.outs)

      nabort.hm <- nabort.hm %>%
        mutate(Estimate = if_else(med == 0, as.integer(0), Estimate))%>%
        arrange(med_num_a)

      nabort.hm$Farm_Num <- factor(nabort.hm$Farm_Num, levels=unique(nabort.hm$Farm_Num))

      FigAbort2 <- ggplot(nabort.hm, aes(x = Farm_Num, y = Estimate))+
        geom_boxplot(fill = "#7570B3", alpha = 0.4,  outlier.color =  "#666666", outlier.alpha = .3, outlier.size = .5)+
        stat_summary(fun="median", geom="point", size=.3, color="#7570B3")+
        geom_hline(yintercept = 0, color = "black") +
        geom_point(aes(y = TotalRVFAbortionsDomR.2010), colour = "black", size = 1.2)+
        labs(x = "Individual Farms", y = "Number of 2010 Abortions \n(Conditional on RVF \nPresence on Farm)")+
        pub_theme+
        theme(#text = element_text(size=16, colour = "black"),
          axis.text.x = element_text(size = 4, angle = 90))

      #FigAbort2

      #Plot of model accuracy with hurdle model
      FigAbort3 <- ggplot(n_abort.dat.hm, aes(x = med_num_a, y =  TotalRVFAbortionsDomR.2010))+#, colour = above0
        geom_point()+
        labs(x = "Predicted Number of Abortions", y = "Reported Number of \nAbortions During \n 2010 Outbreak")+
        pub_theme+
        theme(#text = element_text(size=16, colour = "black"),
          axis.text.x = element_text(size = 4, angle = 90))

      #FigAbort3

      #Save plots to
      if(model2 == TRUE){
        Fig_nAbort<- ggarrange(FigAbort1, FigAbort2, FigAbort3, nrow = 3, ncol = 1, labels = c("(a)", "(b)", "(c)"))

        FigAbort.title <- paste("Publication_Figures/Fig S5 PPC of Abortion - A excludes hurdle ", model.file.name, ".png", sep = "")
        ggexport(Fig_nAbort, filename = FigAbort.title, width = 812, height = 900)
      }

      #Percent error
      n_abort.dat.hm <- n_abort.dat.hm%>%
        mutate(prct.error = if_else(med_num_a== 0 & TotalRVFAbortionsDomR.2010 == 0, 0, med_num_a/TotalRVFAbortionsDomR.2010))%>%# if_else( med_num_a > 0 & TotalRVFAbortionsDomR.2010 == 0, NA,)
        mutate( prct.error = replace(prct.error, prct.error=="Inf", NA))

            ################################
      #Save predicted results and original data together to make plot 2 in the descriptive stats file
      #Estimate the mean of each posterior predictor distribution
      posteriors <- c("nded.hm", "nabort.hm")#sequence with df names

      nded.hm$Farm_ID <- as.character(nded.hm$Farm_ID)
      nabort.hm$Farm_ID <- as.character(nabort.hm$Farm_ID)

      #Create empty dfs
      df.predicted.d <- data.frame(matrix(ncol = 5, nrow = 120))
      df.predicted.p <- data.frame(matrix(ncol = 5, nrow = 120))

      #set counter
      j <- 1
      for(post in posteriors){
        temp <- get(post)#assign df of interest to temporary name
        col.dat <- unique(colnames(temp)[4])#Get the column (is it deaths or abortions?)
        j <- 1#Reset j to 1 for each df
        for(fid in unique(temp$Farm_Num)){#For each unique ID
          filt <-  temp[temp$Farm_Num == fid,]#Filter based on the ID
          med.pred <- median_qi(filt$Estimate, .width = .95)#Take median of predicted estimate
          dat.num <- unique(filt[col.dat])#Get the number of animals on the farm that died/aborted
          dat.tot <- unique(filt["TotalPopDomR_2010"])#Get the number of animals on the farm in 2010
          farm.id <- unique(filt["Farm_ID"])#Get the
          if(col.dat == "TotalRVFDeathsDomR.2010"){#For deaths add the results to the deaths df
            df.predicted.d[j,1] <- farm.id[[1]] #add id
            df.predicted.d[j,2] <- med.pred[1] #add median estimate
            df.predicted.d[j,3] <- dat.num[[1]] #add number of deaths
            df.predicted.d[j,4] <- dat.tot[[1]] #add number of animals
            df.predicted.d[j,5] <- fid #add id
          }
          if(col.dat == "TotalRVFAbortionsDomR.2010"){#repeat for abortions
            df.predicted.p[j,1] <- farm.id[[1]]
            df.predicted.p[j,2] <- med.pred[1]
            df.predicted.p[j,3] <- dat.num[[1]]
            df.predicted.p[j,4] <- dat.tot[[1]] #add number of animals
            df.predicted.p[j,5] <- fid
          }
          j <- j+1 #Add one to counter
        }
      }

      #Name columns
      names(df.predicted.d) <- c("Farm_ID",  "Median.Predicted.Deaths", "Actual.Data.Deaths", "TotalPopDomR_2010", "Farm_Num")
      names(df.predicted.p) <- c("Farm_ID", "Median.Predicted.Abortions", "Actual.Data.Abortions", "TotalPopDomR_2010", "Farm_Num")

      #Combine data frames
      df.predicted <- full_join(df.predicted.d, df.predicted.p)

      #write csv to be manipulated in Descriptive Stats
      if(ModelGP == TRUE & model2 == TRUE){#Only want to save if it's for the GP model
        filname <- paste0("All_Files_For_Publication/data/Predicted number of deaths and abortions and original data_for_publication_", model.file.name, ".csv")
        write.csv(df.predicted, filname, row.names = FALSE)
      }
      ##################################################################
      ##################################################################

      #######################################################################
      #######################################################################
      #Coefficients
      #Death
      if(model2 == TRUE){
        FigCoef_death<- ggplot(death.ggs, aes(y = factor(Param), x = value)) +
          geom_vline(xintercept = 0, color = "red")+
          stat_pointinterval(color = "#1B9E77",.width = c(.95, .66), orientation = "horizontal")+
          geom_violinh(fill = "#1B9E77", alpha = 0.4, scale = "width") +
          scale_y_discrete(breaks = c('d_1', 'd_2', 'd_3', 'd_4', 'd_5', 'd_6', 'd_7', 'd_8'),
                           labels = c("Intercept", "Proportion of sheep", "Proportion of goats", "Total number of ruminants", "Vaccinated for RVFV before outbreak",
                                      "Number of water sources", "2-Month cumulative rainfall (mm) through mid-Mar 2010", "Proportion of days with rain through mid-Mar 2010"))+
          xlim(c(-3.5,2))+
          labs(x = "Parameter", y = "Predictors of RVF Deaths\n ")+
          pub_theme
        #FigCoef_death

        #Abortion
        FigCoef_abort<- ggplot(abort.ggs, aes(y = factor(Param), x = value)) +
          geom_vline(xintercept = 0, color = "red")+
          geom_violinh(fill = "#7570B3", alpha = 0.4, scale = "width") +
          stat_pointinterval(color = "#7570B3",.width = c(.95, .66), orientation = "horizontal")+
          scale_y_discrete(breaks = c('a_1', 'a_2', 'a_3', 'a_4', 'a_5', 'a_6', 'a_7', 'a_8'),
                           labels = c("Intercept", "Proportion of sheep", "Proportion of goats", "Total number of ruminants", "Vaccinated for RVFV before outbreak",
                                      "Number of water sources", "2-Month cumulative rainfall (mm) through mid-Mar 2010", "Proportion of days with rain through mid-Mar 2010"))+
          xlim(c(-3.5,2))+
          labs(x = "Parameter", y = "Predictors of RVF \nAbortions")+
          pub_theme
        #FigCoef_abort

        #RVF_eta
        FigCoef_RVF<- ggplot(eta.ggs, aes(y = factor(Param), x = value)) +
          geom_vline(xintercept = 0, color = "red")+
          geom_violinh(fill = "#D95F02", alpha = 0.4, scale = "width") +
          stat_pointinterval(color = "#D95F02",.width = c(.95, .66), orientation = "horizontal")+
          scale_y_discrete(breaks = c('RVF_1', 'RVF_2', 'RVF_3', 'RVF_4', 'RVF_5', 'RVF_6', 'RVF_7', 'RVF_8', 'RVF_9'),
                           labels = c("Intercept", "Number of ruminants purchased", "Distance ruminants purchased from", "Commericial farm", "Semi-commercial farm",
                                      "Number of pans", "Farm size", "Mix with wildlife", "2-Month cumulative rainfall (mm) through mid-Dec 2009"))+
          labs(x = "Parameter", y = "Predictors of RVF\n ")+
          pub_theme+
          theme(axis.text.y = element_text(margin = margin(t = 10, r = 0, b = 10, l = 0)))
        #FigCoef_RVF

        ###Percent of distribution above and below zero for eta factors
        Eta.prcts <- eta.ggs%>%
          group_by(Param)%>%
          summarise(Prct.Above = sum(value > 0)/length(value),
                    Prct.Below = sum(value < 0)/length(value),)%>%
          mutate(Param = if_else(Param == "RVF_1", "Intercept", if_else(Param == 'RVF_2', "Number of Ruminants Purchased", if_else(Param == 'RVF_3', "Distance Ruminants Purchased From",
                                                                                                                                   if_else(Param == 'RVF_4', "Commericial Farm", if_else(Param == 'RVF_5', "Semi-Commercial Farm", if_else(Param == 'RVF_6', "Number of Pans",
                                                                                                                                                                                                                                           if_else(Param == 'RVF_7', "Farm Size", if_else(Param == 'RVF_8', "Mix with wildlife", if_else(Param == 'RVF_9', "2-Month cumulative rainfall (mm) through mid-December 2009", "NA"))))))))))

        #Save plots to

        Fig_Coef <- ggarrange(FigCoef_death, FigCoef_abort, FigCoef_RVF, nrow = 3, ncol = 1, labels = c("(a)", "(b)", "(c)"))

        FigCoef.title <- paste("Publication_Figures/Fig 2 Coefficient Estimates ", model.file.name, ".png", sep = "")
        ggexport(Fig_Coef, filename = FigCoef.title, width = 812, height = 900)

        write.csv(Eta.prcts, "Publication_Figures/CSV of the percent of eta coefficients above or below zero.csv", row.names = FALSE)
      }
    }#End of if PRIORS_ONLY == FALSE


    #Death
    if(All.One.2 == TRUE){

      death.ggs$Param <- factor(death.ggs$Param, levels=unique(death.ggs$Param))

      FigCoef_death<- ggplot(death.ggs, aes(y = factor(Param), x = value)) +
        geom_vline(xintercept = 0, color = "red")+
        stat_pointinterval(color = "#1B9E77",.width = c(.95, .66), orientation = "horizontal")+
        geom_violinh(fill = "#1B9E77", alpha = 0.4, scale = "width") +
        scale_y_discrete(breaks = c('d_1', 'd_2', 'd_3', 'd_4', 'd_5', 'd_6', 'd_7', 'd_8', 'd_9', 'd_10', 'd_11', 'd_12', 'd_13', 'd_14', 'd_15', 'd_16'),
                         labels = c("Intercept", "Proportion of sheep", "Proportion of goats", "Total number of ruminants", "Vaccinated for RVFV before outbreak",
                                    "Number of water sources", "2-Month cumulative rainfall (mm) through mid-Mar 2010", "Proportion of days with rain through mid-Mar 2010",
                                    "Number of ruminants purchased", "Distance ruminants purchased from", "Commericial farm", "Semi-commercial farm",
                                    "Number of pans", "Farm size", "Mix with wildlife", "2-Month cumulative rainfall (mm) through mid-Dec 2009"))+
        xlim(c(-5,15))+
        labs(x = "Parameter", y = "Predictors of RVF Deaths\n ")+
        pub_theme
      #FigCoef_death

      #Abortion
      abort.ggs$Param <- factor(abort.ggs$Param, levels=unique(abort.ggs$Param))
      FigCoef_abort<- ggplot(abort.ggs, aes(y = factor(Param), x = value)) +
        geom_vline(xintercept = 0, color = "red")+
        geom_violinh(fill = "#7570B3", alpha = 0.4, scale = "width") +
        stat_pointinterval(color = "#7570B3",.width = c(.95, .66), orientation = "horizontal")+
        scale_y_discrete(breaks = c('a_1', 'a_2', 'a_3', 'a_4', 'a_5', 'a_6', 'a_7', 'a_8', 'a_9', 'a_10', 'a_11', 'a_12', 'a_13', 'a_14', 'a_15', 'a_16'),
                         labels = c("Intercept", "Proportion of sheep", "Proportion of goats", "Total number of ruminants", "Vaccinated for RVFV before outbreak",
                                    "Number of water sources", "2-Month cumulative rainfall (mm) through mid-Mar 2010", "Proportion of days with rain through mid-Mar 2010",
                                    "Number of ruminants purchased", "Distance ruminants purchased from", "Commericial farm", "Semi-commercial farm",
                                    "Number of pans", "Farm size", "Mix with wildlife", "2-Month cumulative rainfall (mm) through mid-Dec 2009"))+
        #xlim(c(-3.5,2))+
        labs(x = "Parameter", y = "Predictors of RVF \nAbortions")+
        pub_theme
      #FigCoef_abort

      #RVF_eta
      eta.ggs$Param <- factor(eta.ggs$Param, levels=unique(eta.ggs$Param))
      FigCoef_RVF<- ggplot(eta.ggs, aes(y = factor(Param), x = value)) +
        geom_vline(xintercept = 0, color = "red")+
        geom_violinh(fill = "#D95F02", alpha = 0.4, scale = "width") +
        stat_pointinterval(color = "#D95F02",.width = c(.95, .66), orientation = "horizontal")+
        scale_y_discrete(breaks = c('RVF_1', 'RVF_2', 'RVF_3', 'RVF_4', 'RVF_5', 'RVF_6', 'RVF_7', 'RVF_8', 'RVF_9', 'RVF_10', 'RVF_11', 'RVF_12', 'RVF_13', 'RVF_14', 'RVF_15', 'RVF_16'),
                         labels = c("Intercept", "Proportion of sheep", "Proportion of goats", "Total number of ruminants", "Vaccinated for RVFV before outbreak",
                                    "Number of water sources", "2-Month cumulative rainfall (mm) through mid-Mar 2010", "Proportion of days with rain through mid-Mar 2010",
                                    "Number of ruminants purchased", "Distance ruminants purchased from", "Commericial farm", "Semi-commercial farm",
                                    "Number of pans", "Farm size", "Mix with wildlife", "2-Month cumulative rainfall (mm) through mid-Dec 2009"))+
        labs(x = "Parameter", y = "Predictors of RVF\n ")+
        pub_theme+
        theme(axis.text.y = element_text(margin = margin(t = 10, r = 0, b = 10, l = 0)))
      #FigCoef_RVF

      ###Percent of distribution above and below zero for eta factors
      Eta.prcts <- eta.ggs%>%
        group_by(Param)%>%
        summarise(Prct.Above = sum(value > 0)/length(value),
                  Prct.Below = sum(value < 0)/length(value),)%>%
        mutate(Param = if_else(Param == "RVF_1", "Intercept", if_else(Param == 'RVF_2', "Number of Ruminants Purchased", if_else(Param == 'RVF_3', "Distance Ruminants Purchased From",
                                                                                                                                 if_else(Param == 'RVF_4', "Commericial Farm", if_else(Param == 'RVF_5', "Semi-Commercial Farm", if_else(Param == 'RVF_6', "Number of Pans",
                                                                                                                                                                                                                                         if_else(Param == 'RVF_7', "Farm Size", if_else(Param == 'RVF_8', "Mix with wildlife", if_else(Param == 'RVF_9', "2-Month cumulative rainfall (mm) through mid-December 2009", "NA"))))))))))

      #Save plots to

      Fig_Coef <- ggarrange(FigCoef_death, FigCoef_abort, FigCoef_RVF, nrow = 3, ncol = 1, labels = c("(a)", "(b)", "(c)"))

      FigCoef.title <- paste("Publication_Figures/", model.file.name, "Fig 2 Coefficient Estimates ", ".png", sep = "")
      ggexport(Fig_Coef, filename = FigCoef.title, width = 812, height = 1200)


  }#End of if PRIORS_ONLY == FALSE
    ################################################################################
    #Histograms
    #varnames(hm_mod.mcmc)
    if(ModelRE == TRUE){
      #Get datasets for Random effects parameters
      z.ggs <- ggs(hm_mod.mcmc, family = "z_1")

      RE_z <- z.ggs%>%
        mutate(Param =  "Z")%>%
        select(-Parameter)

      RE_SD.ggs <- ggs(hm_mod.mcmc, family = "sd_1")

      RE_SD <- RE_SD.ggs%>%
        mutate(Param =  "SD")%>%
        select(-Parameter)

      Farm_RE.ggs <- ggs(hm_mod.mcmc, family = "r_1_eta_1")

      Farm_RE <- Farm_RE.ggs%>%
        mutate(Param =  "Farm_RE")%>%
        select(-Parameter)
    }

    if(ModelGP == TRUE){
      #Get datasets for Gaussian process parameters
      gp.ggs <- ggs(hm_mod.mcmc, family = "gp_eff")

      GP <- gp.ggs%>%
        mutate(Param =  "GP_eff")%>%
        select(-Parameter)

      lscale.ggs <- ggs(hm_mod.mcmc, family = "lscale")

      LS <- lscale.ggs%>%
        mutate(Param =  "LS")%>%
        select(-Parameter)

      zgp.ggs <- ggs(hm_mod.mcmc, family = "zgp")

      ZGP <- zgp.ggs%>%
        mutate(Param =  "Z_GP")%>%
        select(-Parameter)

      sdgp.ggs <- ggs(hm_mod.mcmc, family = "sdgp")

      sd_GP <- sdgp.ggs%>%
        mutate(Param =  "SD_GP")%>%
        select(-Parameter)

    }
    #Join datasets
    params <- full_join(abort.ggs, death.ggs)
    params <- full_join(params, eta.ggs)

    if(ModelRE){
      params <- full_join(params, RE_z)
      params <- full_join(params, RE_SD)
      params <- full_join(params, Farm_RE)
    }else{
      if(ModelGP == TRUE){
        params <- full_join(params, GP)
        params <- full_join(params, ZGP)
        params <- full_join(params, sd_GP)
        params <- full_join(params, LS)
      }
    }

    #Summarise all
    all_params_medians <- params%>%
      group_by(Param)%>%
      summarise(median = median(value),
                mean = mean(value))

    #Make the variable a factor to use in facet
    all_params_medians$Param <- as.factor(all_params_medians$Param)


    Plot_labeller <- function(variable,value){
      return(facet.lab[value])
    }

    if(ModelRE == TRUE){
      facet.lab <- c('a_1' = "Intercept (A)", 'a_2' = "Proportion Sheep (A)", 'a_3' ="Total Ruminants (A)", 'a_4' = "RVFV Vaccinated (A)", 'a_5' = "Number water sources (A)", 'a_6' = "2m Cumulative rainfall (A)", 'a_7' = "Percent Rainy Days (A)",
                     'd_1' = "Intercept (D)", 'd_2' = "Proportion Sheep (D)", 'd_3' ="Total Ruminants (D)", 'd_4' = "RVFV Vaccinated (D)", 'd_5' = "Number water sources (D)", 'd_6' = "2m Cumulative rainfall (D)", 'd_7' = "Percent Rainy Days (D)",
                     'RVF_1' = "Intercept (RVF)", 'RVF_2' = "Ruminants Purchased (R)", 'RVF_3' = "Purchase Distance (R)", 'RVF_4' = "Commercial Farm", 'RVF_5' = "Semi-Commercial Farm", 'RVF_5' = "Accessible Pans (R)", 'RVF_6' = "Farm Size (R)", 'RVF_7' = "Wildlife (R)", 'RVF_8' = "Average Rainfall (R)",
                     "Z" = "Partial Individual RE", "SD" = "Random Effect SD", "Farm_RE" = "Farm Random Effects")
    }else{
      if(ModelGP == TRUE){
        if(model2 == TRUE){
          facet.lab <- c('a_1' = "Intercept (A)", 'a_2' = "Proportion sheep (A)", 'a_3' ="Proportion goats (A)", 'a_4' = "Total ruminants (A)", 'a_5' = "RVFV vaccinated (A)", 'a_6' = "Number water sources (A)", 'a_7' = "Cumulative rainfall Mar (A)", 'a_8' = "Percent rainy days Mar (A)",
                         'd_1' = "Intercept (D)", 'd_2' = "Proportion sheep (D)", 'd_3' ="Proportion goats (D)", 'd_4' = "Total ruminants (D)", 'd_5' = "RVFV vaccinated (D)", 'd_6' = "Number water sources (D)", 'd_7' = "Cumulative rainfall Mar (D)", 'd_8' = "Percent rainy days Mar (D)",
                         'GP_eff' = "Gaussian Process Effect", 'LS' = "Length-Scale",
                         'RVF_1' = "Intercept (R)", 'RVF_2' = "Ruminants Purchased (R)", 'RVF_3' = "Purchase Distance (R)", 'RVF_4' = "Commercial Farm", 'RVF_5' = "Semi-Commercial Farm" , 'RVF_6' = "Accessible Pans (R)" , 'RVF_7' = "Farm Size (R)" , 'RVF_8' = "Wildlife (R)", 'RVF_9' = "Cumulative rainfall Dec (R)",
                         'SD_GP' = "Gaussian Process SD", 'Z_GP' = "Latent Gaussian Process")
        }
      }
    }



      if(exists("facet.lab")){
        #Plot all histograms
        FigAllHists <- ggplot(params, aes(x = value))+
          geom_histogram()+
          facet_wrap(~Param, scales = "free", labeller = labeller(Param = facet.lab))+#
          geom_vline(data = all_params_medians, aes(xintercept = median),col='red')+
          geom_vline(data = all_params_medians, aes(xintercept = mean),col='blue')+
          labs(x = "Value", y = "Count")+#, title = "Histogram of parameters with median in red and mean in blue")+
          pub_theme

        #FigAllHists

        #Save plots to
        Fig_Hist <- ggarrange(FigAllHists, nrow = 1, ncol = 1)

        FigHist.title <- paste("Publication_Figures/SX Histograms of all parameters median in red and mean in blue ", model.file.name, ".png", sep = "")
        ggexport(Fig_Hist, filename = FigHist.title, width = 1500, height = 900)

        #Traceplots for supplement

        facet.lab.d <- c('b_1[1]' = "Intercept (D)", 'b_1[2]' = "Proportion sheep (D)", 'b_1[3]' ="Proportion goats (D)", 'b_1[4]' = "Total ruminants (D)", 'b_1[5]' = "RVFV vaccinated (D)", 'b_1[6]' = "Number water sources (D)", 'b_1[7]' = "Cumulative rainfall Mar (D)", 'b_1[8]' = "Proportion rainy days Mar (D)")
        facet.lab.a <- c('b_2[1]' = "Intercept (A)", 'b_2[2]' = "Proportion sheep (A)", 'b_2[3]' ="Proportion goats (A)", 'b_2[4]' = "Total ruminants (A)", 'b_2[5]' = "RVFV vaccinated (A)", 'b_2[6]' = "Number water sources (A)", 'b_2[7]' = "Cumulative rainfall Mar (A)", 'b_2[8]' = "Proportion rainy days Mar (A)")
        facet.lab.rvf <- c('b_eta[1]' = "Intercept (R)", 'b_eta[2]' = "Ruminants Purchased (R)", 'b_eta[3]' = "Purchase Distance (R)", 'b_eta[4]' = "Commercial Farm", 'b_eta[5]' = "Semi-Commercial Farm" , 'b_eta[6]' = "Accessible Pans (R)" , 'b_eta[7]' = "Farm Size (R)" , 'b_eta[8]' = "Wildlife (R)", 'b_eta[9]' = "Cumulative rainfall Dec (R)")

        #Death traceplots
        death.trace <- ggs_traceplot(params.ggs, family = "b_1")+
          facet_wrap(~Parameter, scales = "free", labeller = labeller(Parameter = facet.lab.d))
        abort.trace <- ggs_traceplot(params.ggs, family = "b_2")+
          facet_wrap(~Parameter, scales = "free", labeller = labeller(Parameter = facet.lab.a))
        rvf.trace <- ggs_traceplot(params.ggs, family = "b_eta")+
          facet_wrap(~Parameter, scales = "free", labeller = labeller(Parameter = facet.lab.rvf))

        rvf.trace <- ggarrange(rvf.trace, nrow = 1, ncol = 1)

        Figtrace.title <- paste("Publication_Figures/S3 Traceplots of rvf hurdle params", model.file.name, ".png", sep = "")
        ggexport(rvf.trace, filename = Figtrace.title, width = 1500, height = 900)
      }
      ##########################################
      #Table 2 with unstandardized parameters - including lp and RE
      if(model2 ==TRUE){
        #Unstandardize the standardized variables
        at <- c("a_2", #PrctSheep_Ave.2010.Stand",
                "a_3",#"PrctGoats_Ave.2010.Stand",
                "a_4", #TotalPopDomR_2010.Stand
                "a_6", #"Num_Natural_WaterSources.Ave.Stand",
                "a_7",#"CumSum.2m.Rain.Stand",
                "a_8", #"Prct.Rain.Days.2m.Stand",
                "d_2", #PrctSheep_Ave.2010.Stand",
                "d_3",#"PrctGoats_Ave.2010.Stand",
                "d_4", #TotalPopDomR_2010.Stand
                "d_6", #"Num_Natural_WaterSources.Ave.Stand",
                "d_7",#"CumSum.2m.Rain.Stand",
                "d_8", #"Prct.Rain.Days.2m.Stand",
                "RVF_2", #"Sum.Ruminants.Gained.Ave.Stand",
                "RVF_3", #"Max.Dist.Purch.Ave.Stand",
                "RVF_6",#"f9.Pan.Num.Ave.Stand")
                "RVF_7",#f14.Size.Farm.15.Stand",
                "RVF_9")#"Ave.2m.CumSum_Oct.Stand",

        key <- cbind("a_2" = "PrctSheep_Ave.2010", #as.data.frame(
                     "a_3" = "PrctGoats_Ave.2010",
                     "a_4" = "TotalPopDomR_2010",
                     "a_6" = "Num_Natural_WaterSources.Ave",
                     "a_7" = "CumSum.2m.Rain",
                     "a_8" = "Prct.Rain.Days.2m",
                     "d_2" = "PrctSheep_Ave.2010",
                     "d_3" = "PrctGoats_Ave.2010",
                     "d_4" = "TotalPopDomR_2010",
                     "d_6" = "Num_Natural_WaterSources.Ave",
                     "d_7" = "CumSum.2m.Rain",
                     "d_8" = "Prct.Rain.Days.2m",
                     "RVF_2" = "Sum.Ruminants.Gained.Ave",
                     "RVF_3" = "Max.Dist.Purch.Ave",
                     "RVF_6" = "f9.Pan.Num.Ave",
                     "RVF_7" = "f14.Size.Farm.15",
                     "RVF_9" = "Ave.2m.CumSum_Oct")

        dfs <- c("death.ggs", "abort.ggs", "eta.ggs")
        Est <- as.data.frame(matrix(ncol = 3, nrow = 17))
        j <- 1
        for(ggs in dfs){#For each dataset
          for(vr in at){#and each variable
            if(ggs == "death.ggs"){
              filt <- filter(death.ggs, Param == vr)#filter by the variable
            }else{
              if(ggs == "abort.ggs"){
                filt <- filter(abort.ggs, Param == vr)
              }else{
                if(ggs == "eta.ggs"){
                  filt <- filter(eta.ggs, Param == vr)
                }
              }
            }
            if(length(filt$Param)==0){#if the variable is 0 then skip
              print(paste("skip", vr))
            }else{
              cof <- median_qi(filt$value, .width = .95)#Otherwise calculate the median
              col.nme <- key[, vr]#name the columns
              dat_col<-mydata[col.nme]#Get the data column with the original data's mean and sd
              Est[j, 1] <- col.nme
              dat_col <- unlist(dat_col)
              unstand.median <- cof[1]*sd(dat_col) + mean(dat_col)#Unstandardize the median
              unstand.median.ci.hi <- cof[3]*sd(dat_col) + mean(dat_col)#Unstandardize the median confidence interval
              unstand.median.ci.lo <- cof[2]*sd(dat_col) + mean(dat_col)
              Est[j,2] <- round(unstand.median,2)
              Est[j,3]<- paste0("[", round(unstand.median.ci.lo,2),",", round(unstand.median.ci.hi,2), "]")#Write the results for table
              j <- j+1
            }
          }
        }

        #Correct Variable names
        var.names <- c("Percent of sheep (D)", "Percent of goats (D)", "Number of domestic Ruminants (D)", "Number water sources (D)", "2-Month cumulative rainfall through mid-March 2010 (D)", "Percent of day with rain through mid-March 2010 (D)",
                       "Percent of sheep (A)", "Percent of goats (A)", "Number of domestic Ruminants (A)", "Number water sources (A)", "2-Month cumulative rainfall through mid-March 2010 (A)", "Percent of day with rain through mid-March 2010 (A)",
                       "Ruminants purchased (R)", "Purchase distance (R)", "Accessible pans (R)", "Farm size (R)", "2-Month cumulative rainfall through mid-Dec 2009 (R)")

        Est[,1] <- var.names

        #Categorical variables
        at.cat <- c("a_1", #Intercept",
                    "a_5",#"X5.RVF.Vax.before.Outbreak.2010",
                    "d_1", #Intercept",
                    "d_5",#"X5.RVF.Vax.before.Outbreak.2010",
                    "RVF_1", #"Intercept",
                    "RVF_4", #"f18.Production.15Com",
                    "RVF_5",#"f18.Production.15SC")
                    "RVF_8")#"f5.Wildlife.Mix.15",

        key.cat <- cbind("a_1" = "Intercept (A)", #as.data.frame(
                         "a_5" = "RVFV Vaccinated (A)",
                         "d_1" = "Intercept (D)",
                         "d_5" = "RVFV Vaccinated (D)",
                         "RVF_1" = "Intercept (R)",
                         "RVF_4" = "Commercial Farm (R)",
                         "RVF_5" = "Semi-Commercial Farm (R)",
                         "RVF_8" = "Wildlife (R)")

        #Categorical variables - similar code to above, see comments there.
        Est.cat <- as.data.frame(matrix(ncol = 3, nrow = 8))
        j <- 1
        for(ggs in dfs){
          for(vr in at.cat){
            if(ggs == "death.ggs"){
              filt <- filter(death.ggs, Param == vr)
            }else{
              if(ggs == "abort.ggs"){
                filt <- filter(abort.ggs, Param == vr)
              }else{
                if(ggs == "eta.ggs"){
                  filt <- filter(eta.ggs, Param == vr)
                }
              }
            }
            if(length(filt$Param)==0){
              print(paste("skip", vr))
            }else{
              cof.cat <- median_qi(filt$value, .width = .95)
              col.nme <- key.cat[, vr]
              Est.cat[j, 1] <- col.nme
              Est.median <- cof.cat[1]
              median.ci.hi <- cof.cat[3]
              median.ci.lo <- cof.cat[2]
              Est.cat[j,2] <- round(Est.median,2)
              Est.cat[j,3]<- paste0("[", round(median.ci.lo,2),",", round(median.ci.hi,2), "]")
              j <- j+1
            }
          }
        }
      }
      #Add GP and RE estimates
      if(ModelGP == TRUE & model2 == TRUE){
        Est.gp <- data.frame()
        #Add length-scale
        filt <- filter(LS, Param == "LS")
        cof.gp <- median_qi(filt$value, .width = .95)
        col.nme <- "GP Length Scale (R)"
        Est.gp[1, 1] <- col.nme
        Est.median <- cof.gp[1]
        median.ci.hi <- cof.gp[3]
        median.ci.lo <- cof.gp[2]
        Est.gp[1,2] <- round(Est.median,2)
        Est.gp[1,3]<- paste0("[", round(median.ci.lo,2),",", round(median.ci.hi,2), "]")
        #Add variance (SD^2)
        filt <- filter(sd_GP, Param == "SD_GP")
        cof.gp_SD <- median_qi(filt$value, .width = .95)
        col.nme <- "GP Variance (R)"
        Est.gp[2, 1] <- col.nme
        Est.median <- (cof.gp_SD[1])^2
        median.ci.hi <- (cof.gp_SD[3])^2
        median.ci.lo <- (cof.gp_SD[2])^2
        Est.gp[2,2] <- round(Est.median,2)
        Est.gp[2,3]<- paste0("[", round(median.ci.lo,2),",", round(median.ci.hi,2), "]")
        names(Est.gp) <- c("Variable", "Unstandardized Median Estimate", "Unstandardized 95% CI")
        est.tab.name <- paste0("Publication_Figures/Temp GP Param Median and CI for Tab 7_", model.file.name, ".csv")
        write.csv(Est.gp, est.tab.name, row.names = FALSE)

      }

      #Put table in order
      if(model2 == TRUE){
        Est.Tab <- data.frame(matrix(ncol = 3, nrow = 29))
        Est.Tab[1,] <- c("RVF on the Farm Predictors (Hurdle)", "","")
        Est.Tab[2,] <- Est.cat[5,]
        Est.Tab[3,] <- Est[13,]
        Est.Tab[4,] <- Est[14,]
        Est.Tab[5,] <- Est.cat[6,]
        Est.Tab[6,] <- Est.cat[7,]
        Est.Tab[7,] <- Est[15,]
        Est.Tab[8,] <- Est[16,]
        Est.Tab[9,] <- Est.cat[8,]
        Est.Tab[10,] <- Est[17,]
        Est.Tab[11,] <- Est.gp[1,]
        Est.Tab[12,] <- Est.gp[2,]
        #Est.Tab[13,] <- Est.GP.RE[3,]
        Est.Tab[13,] <- c("Abortion Predictors", "","")
        Est.Tab[14,] <- Est.cat[3,]
        Est.Tab[15,] <- Est[7,]
        Est.Tab[16,] <- Est[8,]
        Est.Tab[17,] <- Est[9,]
        Est.Tab[18,] <- Est.cat[4,]
        Est.Tab[19,] <- Est[10,]
        Est.Tab[20,] <- Est[11,]
        Est.Tab[21,] <- Est[12,]
        Est.Tab[22,] <- c("Death Predictors", "","")
        Est.Tab[23,] <- Est.cat[1,]
        Est.Tab[24,] <- Est[1,]
        Est.Tab[25,] <- Est[2,]
        Est.Tab[26,] <- Est[3,]
        Est.Tab[27,] <- Est.cat[2,]
        Est.Tab[28,] <- Est[4,]
        Est.Tab[29,] <- Est[5,]
        Est.Tab[30,] <- Est[6,]


        names(Est.Tab) <- c("Variable", "Unstandardized Median Estimate", "Unstandardized 95% CI")
        fil.name.coef.tab <- paste0("Publication_Figures/Table 7 Unstandardized Coefficient Estimates and Credible Intervals ", model.file.name, ".csv")
        write.csv(Est.Tab, fil.name.coef.tab, row.names = FALSE)
      }
      ######################################
      ##Variance accounted for by random effect
      if(ModelRE == TRUE){
        #Get RE data
        SD.ggs <- ggs(hm_mod.mcmc, family = "sd_1")
        RE_SD <- mean(SD.ggs$value)
        RE_Var <- as.data.frame(RE_SD^2)#Calculate variance
        names(RE_Var) <- "Variance_of_RE"

        RE.file.name <- paste("Publication_Figures/RE Variance Calculation of overall model ", model.file.name, ".csv", sep = "")
        write.csv(RE_Var, RE.file.name, row.names = FALSE)

        #Individual district random effects
        RE_district <- unique(mydata$District)
        RE_district <- as.character(RE_district)
        z.ggs$District <- "NA"

        #Assign correct district per Z_1 number
        for(i in 1:length(RE_district)){
          for(j in 1:nrow(z.ggs)){
            re <- paste0("z_1[1,", i, "]")
            if(z.ggs$Parameter[j] == re){
              z.ggs$District[j] <- RE_district[i]
            }
          }
        }
        #Take the median of each district to get estimate and CIs
        RE.dist_Tab <- as.data.frame(matrix(ncol = 3, nrow = 20))
        k <- 1
        for(re in RE_district){
          #for(o in 1:nrow(z.ggs)){
          #if(z.ggs$District[o] == RE_district[k]){
          filt <- filter(z.ggs, District == RE_district[k])
          #Make table
          cof.re <- median_qi(filt$value, .width = .95)
          col.nme_re <- RE_district[k]
          RE.dist_Tab[k, 1] <- col.nme_re
          Est.median.re <- cof.re[1]
          median.re.ci.hi <- cof.re[3]
          median.re.ci.lo <- cof.re[2]
          RE.dist_Tab[k,2] <- round(Est.median.re,2)
          RE.dist_Tab[k,3]<- paste0("[", round(median.re.ci.lo,2),",", round(median.re.ci.hi,2), "]")
          RE.dist_Tab[k,4] <- round(median.re.ci.lo,2)
          RE.dist_Tab[k,5] <- round(median.re.ci.hi,2)

          k <- k+1
        }

        names(RE.dist_Tab) <- c("District", "Median_Variance", "95 CI", "Low 95CI", "High 95CI")



        write.csv(RE.dist_Tab, "Publication_Figures/Random Effect Variance by Farm for publication.csv", row.names = FALSE)
      }
      ###############################
      #GP plot
      if(ModelGP == TRUE){
        #Exponential curve
        ls <- median_qi(LS$value, .width = .95)
        lscale <- ls[[1]]
        lscale.low <- ls[[2]]
        lscale.hi <- ls[[3]]

        sd1 <- median_qi(sd_GP$value, .width = .95)
        SDGP <- sd1[[1]]
        SDGP.lo <- sd1[[2]]
        SDGP.hi <- sd1[[3]]

        # Scale coordinates so the maximum distance between two points is one and center
        distance_scale <- max(dist(cbind(mydata$X, mydata$Y)))

        #Change to km instead of m
        distance_scale <- distance_scale/1000

        #Bind data
        locs <- cbind(mydata$X, mydata$Y)
        #Get data in km
        locs[,1] <- locs[,1]/1000
        locs[,2] <- locs[,2]/1000

        #Calculate all parawise distances
        dists <- dist(locs)
        #Make into matrix
        dists2 <- as.matrix(dists)

        dists.v <- unique(dists)
        #Scale the distances using the distance scale
        s.dists.v <- dists.v/distance_scale

        #Distance decay function
        dist_cov <- exp(-s.dists.v^2/(2*lscale^2))#SDGP^2*
        dist_cov.lo <- exp(-s.dists.v^2/(2*lscale.low^2))#SDGP.lo^2*
        dist_cov.hi <- exp(-s.dists.v^2/(2*lscale.hi^2))#SDGP.hi^2*

        #Plot the decay rate
        dist_plot <- cbind(s.dists.v, dist_cov, dist_cov.lo, dist_cov.hi)

        dist_plot <- as.data.frame(dist_plot)

        names(dist_plot) <- c("Scaled_Dist", "Covariance", "Covariance.lo", "Covariance.hi")

        #Unscale the distances
        dist_plot$Dist.us <- dist_plot$Scaled_Dist* distance_scale

        #Distances which explain 99% of the covariance
        dist_check <- filter(dist_plot, Covariance > 0.01)

        ind <- which.min(dist_check$Covariance)

        #Distance with lowest amount of covariance (e.g. 99% the covariance is explained at smaller of this value)
        Dist.Lo.Cov <- dist_check$Dist.us[ind]

        #Pecenct of the longest distance that explains covariance
        prct.dist <- Dist.Lo.Cov/max(dist_plot$Dist.us)

        #Mean and median unscaled distances between farms
        mean.dist <- mean(dist_plot$Dist.us)
        median.dist <- median(dist_plot$Dist.us)
        max.dist <- max(dist_plot$Dist.us)

        #number of distances that explain 99% of the covariance
        n.less.Dist.Lo.Cov <- nrow(filter(dist_plot, Dist.us <= Dist.Lo.Cov))

        #Percent of total number of points that explain 99% of the covariance
        prct.less.Dist.Lo.Cov <- n.less.Dist.Lo.Cov/nrow(dist_plot)


        #Set vector for tick marks on x axis 290
        xtick <- seq(0,290, 20)

        #Figure S1
        GP_assessment <- ggplot(dist_plot, aes(x = Dist.us, y = Covariance))+
          geom_line()+
          geom_point(aes(y=0), shape = "|", color = "purple")+
          geom_line(aes(y = Covariance.hi), color = "green")+
          geom_line(aes(y = Covariance.lo), color = "green")+
          geom_vline(xintercept = Dist.Lo.Cov, color = "pink")+
          labs(x = "Distance (Km)", y = "Covariance")+
          scale_x_continuous(breaks = xtick, labels = xtick)+
          pub_theme

        #Save plots to
        if(model2 == TRUE){
          GP_assessment <- ggarrange(GP_assessment, nrow = 1, ncol = 1)

          FigGP_assessment.title <- paste("Publication_Figures/S2 Gaussian Process Covariance by Distance in KM ", model.file.name, ".png", sep = "")
          ggexport(GP_assessment, filename = FigGP_assessment.title, width = 615, height = 377)
        }
        #############
        #Prior only model
        if(PRIORS_ONLY == TRUE){
          GP_assessment <- ggarrange(GP_assessment, nrow = 1, ncol = 1)

          FigGP_assessment.title.priors.only <- paste("Publication_Figures/SX Gaussian Process Covariance by Distance in KM ", model.file.name, "_PRIORS_ONLY.png", sep = "")
          ggexport(GP_assessment, filename = FigGP_assessment.title.priors.only, width = 615, height = 377)

        }

      }


      ##############################
      #Estimate WAIC
      if(ModelRE == TRUE){
        WAIC.RE <- waic(hm_fit)
        #returns the waic and SE
        #pwaic estimated number of parameters and SE
        #lppd is the expected log pointwise predictive density (or -waic/2)
        WAIC_comp[1,2] <- WAIC.RE$total["waic"]
        WAIC_comp[1,3] <- WAIC.RE$se["waic"]
        WAIC_comp[1,4] <- as.Date(datestamp)

        #Weight the WAIC
        #sum of exponents
        #sum_exp_waic <- exp(WAIC_comp$WAIC[1])+exp(WAIC_comp$WAIC[2])

        write.csv(WAIC_comp, "Publication_Figures/WAIC comparison between RE and GP models_for_publication.csv", row.names = FALSE)
      }else{
        if(ModelGP == TRUE & PRIORS_ONLY == FALSE){
          WAIC <- waic(hm_fit)
          #returns the waic and SE
          #pwaic estimated number of parameters and SE
          #lppd is the expected log pointwise predictive density (or -waic/2)
          if(modelThesis == TRUE){
            WAIC_comp[2,2] <- WAIC$total["waic"]
            WAIC_comp[2,3] <- WAIC$se["waic"]
            WAIC_comp[2,4] <- as.Date(datestamp)}
          if(modelPub == TRUE){
            WAIC_comp[3,2] <- WAIC$total["waic"]
            WAIC_comp[3,3] <- WAIC$se["waic"]
            WAIC_comp[3,4] <- as.Date(datestamp)
          }
          if(model2 == TRUE){
            WAIC_comp[6,2] <- WAIC$total["waic"]
            WAIC_comp[6,3] <- WAIC$se["waic"]
            WAIC_comp[6,4] <- as.Date(datestamp)
          }
          if(modelPres == TRUE){
            WAIC_comp[4,2] <- WAIC$total["waic"]
            WAIC_comp[4,3] <- WAIC$se["waic"]
            WAIC_comp[4,4] <- as.Date(datestamp)
          }

          if(modelThesisOct == TRUE){
            WAIC_comp[5,2] <- WAIC$total["waic"]
            WAIC_comp[5,3] <- WAIC$se["waic"]
            WAIC_comp[5,4] <- as.Date(datestamp)
          }
          if(All.One.2 == TRUE){
            WAIC_comp[7,2] <- WAIC$total["waic"]
            WAIC_comp[7,3] <- WAIC$se["waic"]
            WAIC_comp[7,4] <- as.Date(datestamp)
          }


          write.csv(WAIC_comp, "Publication_Figures/WAIC comparison between RE and GP models_for_publication.csv", row.names = FALSE)

          WAIC_comp <- WAIC_comp%>%
            mutate(SE_hi = WAIC + SE)%>%
            mutate(SE_lo = WAIC - SE)

          WAIC_comp[,1] <- c("Mod2", "Mod3", "Mod4", "Mod5", "Mod6", "Mod1", "Mod7")

          WAIC_comp$Model <- factor(WAIC_comp$Model, levels=sort(unique(WAIC_comp$Model)))

          Mod_labs <- c("Model 1", "Model 2", "Model 3", "Model 4", "Model 5", "Model 6", "Model 7")

          #Figure S4
          WAIC.plot <- ggplot(WAIC_comp, aes(x = Model, y = WAIC))+
            geom_point()+
            geom_errorbar(mapping=aes(x = Model, ymin = SE_lo, ymax = SE_hi), width=0.2)+
            scale_x_discrete(labels= Mod_labs)+
            pub_theme+
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

          if(model2 == TRUE){
            Fig_S7_WAIC<- ggarrange(WAIC.plot)
            FigWAIC.title <- paste("Publication_Figures/Fig S7 WAIC comparison", ".png", sep = "")
            ggexport(Fig_S7_WAIC, filename = FigWAIC.title, width = 623, height = 434)
          }

        }
      }


    ####################################Deviance assessment
    if(ModelGP == TRUE & PRIORS_ONLY == FALSE){
      #Calculate the deviance from the Log Likelihood
      LL.ggs <- ggs(hm_mod.mcmc, family = "log_lik")
      LL.est <- mean_qi(LL.ggs$value, .width = .95)
      LL <- LL.est[[1]]
      LL.lo <- LL.est[[2]]
      LL.hi <- LL.est[[3]]

      Deviance.gp <- -2*LL
      Deviance.gp.hi <- -2*LL.lo
      Deviance.gp.lo <- -2*LL.hi

      if(modelThesis == TRUE){
        Dev_Comp[2,2] <- round(Deviance.gp,4)
        Dev_Comp[2,3] <- round(Deviance.gp.lo,4)
        Dev_Comp[2,4] <- round(Deviance.gp.hi,4)
      }
      if(modelPub == TRUE){
        Dev_Comp[5,2] <- round(Deviance.gp,4)
        Dev_Comp[5,3] <- round(Deviance.gp.lo,4)
        Dev_Comp[5,4] <- round(Deviance.gp.hi,4)
      }
      if(model2 == TRUE){
        Dev_Comp[6,2] <- round(Deviance.gp,4)
        Dev_Comp[6,3] <- round(Deviance.gp.lo,4)
        Dev_Comp[6,4] <- round(Deviance.gp.hi,4)
      }
      write.csv(Dev_Comp, "Publication_Figures/Comparison of Deviance between Models_for_publication.csv", row.names = FALSE)
    }

    if(ModelRE == TRUE){
      #Calculate the deviance from the LL
      LL.ggs <- ggs(hm_mod.mcmc, family = "log_lik")
      LL.est <- mean_qi(LL.ggs$value, .width = .95)
      LL <- LL.est[[1]]
      LL.lo <- LL.est[[2]]
      LL.hi <- LL.est[[3]]

      Deviance.re <- -2*LL
      Deviance.re.hi <- -2*LL.lo
      Deviance.re.lo <- -2*LL.hi

      Dev_Comp[3,2] <- round(Deviance.re,4)
      Dev_Comp[3,3] <- round(Deviance.re.lo,4)
      Dev_Comp[3,4] <- round(Deviance.re.hi,4)
      write.csv(Dev_Comp, "Publication_Figures/Comparison of Deviance between Models_for_publication.csv", row.names = FALSE)
    }

    #This is the end of what to skip if this is an intercept only model or model with no spatial effect
  }
}
Dev_Comp$Deviance <- as.numeric(Dev_Comp$Deviance)
Int.Dev <- Dev_Comp$Deviance[4]


#Calculate deviance explained
Dev_Comp$Dev.Explained <- 1-Dev_Comp$Deviance/Int.Dev


write.csv(Dev_Comp, "Publication_Figures/Comparison of Deviance between Models_for_publication.csv", row.names = FALSE)
