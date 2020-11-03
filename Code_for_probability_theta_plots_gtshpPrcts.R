library(egg)
library(ggpubr)
library(tidyverse)
library(rstan)
library(bayesplot)
library(tidybayes)
library(scales)

# Convert the model to C++ code, then compile
hmgp_mod_c <- stanc("All_Files_For_Publication/hurdle-multinomial-gp.stan", verbose = TRUE)
hmgp_mod <- stan_model(stanc_ret = hmgp_mod_c, verbose = TRUE)

#Load data
df <- read.csv("All_Files_For_Publication/data/Updated 2017 and 2015 Cross-sectional Data for Bayesian Analsysis_with_UTM.csv")
un.stand <- read.csv("All_Files_For_Publication/data/Unstandardization_Key_for_publication.csv")

#Set up data for
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


# This allows us to use the functions we define in stan code, notably rsoftmax()
expose_stan_functions(hmgp_mod)
rsoftmax(c(1,1,1))
sum(rsoftmax(c(0.2, 2, 3)))

#Set variables that we want to save data from
theta.cont.vars <- c("PrctSheep_Ave.2010.Stand", "PrctGoats_Ave.2010.Stand", "TotalPopDomR_2010.Stand",  "Num_Natural_WaterSources.Ave.Stand", "CumSum.2m.Rain.Stand","Prct.Rain.Days.2m.Stand")
                    #Don't include eta variables because they are binomial and not multinomial. "Sum.Ruminants.Gained.Ave.Stand", "Max.Dist.Purch.Ave.Stand", "f9.Pan.Num.Ave.Stand", "f14.Size.Farm.15.Stand", "CumSum.2m.Rain_Oct.Stand")


#Create blank dataframe and list
output <- data.frame(matrix(ncol = 4, nrow = 5))
all.data <- list()
j <- 1

#Set loop based on variables we are intereted in changing
for(v in theta.cont.vars){
  #Set each column to the mean or the reference value
  dat <- df%>%
    select(TotalPopDomR_2010, TotalRVFDeathsDomR.2010, TotalRVFAbortionsDomR.2010, X3.Outbreak,
           PrctSheep_Ave.2010.Stand, PrctGoats_Ave.2010.Stand, TotalPopDomR_2010.Stand,
           X5.RVF.Vax.before.Outbreak.2010, Num_Natural_WaterSources.Ave.Stand, CumSum.2m.Rain.Stand, Prct.Rain.Days.2m.Stand,
           Sum.Ruminants.Gained.Ave.Stand, Max.Dist.Purch.Ave.Stand, f18.Production.15, f9.Pan.Num.Ave.Stand, f14.Size.Farm.15.Stand,  f5.Wildlife.Mix.15, CumSum.2m.Rain_Oct.Stand)%>%
    mutate(PrctSheep_Ave.2010.Stand = mean(PrctSheep_Ave.2010.Stand))%>%
    mutate(PrctGoats_Ave.2010.Stand = mean(PrctGoats_Ave.2010.Stand))%>%
    mutate(TotalPopDomR_2010.Stand = mean(TotalPopDomR_2010.Stand))%>%
    mutate(X5.RVF.Vax.before.Outbreak.2010 = 0)%>%
    mutate(X3.Outbreak = 0)%>%
    mutate(Num_Natural_WaterSources.Ave.Stand = mean(Num_Natural_WaterSources.Ave.Stand))%>%
    mutate(CumSum.2m.Rain.Stand = mean(CumSum.2m.Rain.Stand))%>%
    mutate(Prct.Rain.Days.2m.Stand = mean(Prct.Rain.Days.2m.Stand))%>%
    mutate(Sum.Ruminants.Gained.Ave.Stand = mean(Sum.Ruminants.Gained.Ave.Stand))%>%
    mutate(Max.Dist.Purch.Ave.Stand = mean(Max.Dist.Purch.Ave.Stand))%>%
    mutate(f18.Production.15 = 0)%>%#Is this ok instead of Small Holder? If I put Small Holder I get an error that contrasts need at least two factors...
    mutate(f9.Pan.Num.Ave.Stand = mean(f9.Pan.Num.Ave.Stand))%>%
    mutate(f14.Size.Farm.15.Stand = mean(f14.Size.Farm.15.Stand))%>%
    mutate(f5.Wildlife.Mix.15 = 0)%>%
    mutate(CumSum.2m.Rain_Oct.Stand = mean(CumSum.2m.Rain_Oct.Stand))

  #Change our one variable back to the normal data
  dat[v] <- df[v]

# Assemble the data. Use dat instead of df
  standata <- lst(
    # Matrix of the three possible outcomes for each farm
    Y = cbind(
      dat$TotalPopDomR_2010 - dat$TotalRVFDeathsDomR.2010 - dat$TotalRVFAbortionsDomR.2010,
      dat$TotalRVFDeathsDomR.2010,
      dat$TotalRVFAbortionsDomR.2010
    ),
    N = nrow(Y),
    ncat = ncol(Y),

    #data for residual estimation
    outbreak = dat$X3.Outbreak,
    num_animals = dat$TotalPopDomR_2010,

    # Design matrix for outcomes
    X = model.matrix(~ PrctSheep_Ave.2010.Stand +
                       PrctGoats_Ave.2010.Stand +
                       TotalPopDomR_2010.Stand +
                       X5.RVF.Vax.before.Outbreak.2010 +
                       Num_Natural_WaterSources.Ave.Stand +
                       CumSum.2m.Rain.Stand +
                       Prct.Rain.Days.2m.Stand,
                     data = dat),
    K = ncol(X),

    # Design matrix for hurdles
    X_eta = model.matrix(~Sum.Ruminants.Gained.Ave.Stand +
                           Max.Dist.Purch.Ave.Stand +
                           f18.Production.15 +
                           f9.Pan.Num.Ave.Stand +
                           f14.Size.Farm.15.Stand +
                           f5.Wildlife.Mix.15 +
                           CumSum.2m.Rain_Oct.Stand,
                         data = dat),
    K_eta = ncol(X_eta),


    # Gaussian Process
    Dgp = 2,  # Dimension
    Xgp = locs_standardized, #Values


    prior_only = FALSE # Change this to sample from the priors rather than posterior
  )

# Sample from the model.  Use cores as you have available.
hmgp_fit <- sampling(
  hmgp_mod,
  data = standata,
  chains = 4,
  iter = 5000,
  thin = 1, #default is treedepth of 10
  cores = 4)

#Extract thetas
thetas <- extract(hmgp_fit, "theta")[[1]]

#Convert to probability scale
outcome <- apply(thetas, 1:2, rsoftmax)

#Get the probabilities for each theta outcome
prob_theta_1 <- as.data.frame(outcome[1,,])
prob_theta_2 <- as.data.frame(outcome[2,,])
prob_theta_3 <- as.data.frame(outcome[3,,])

#Set up blank dataframes
df.result <- data.frame(matrix(ncol = 7))
names(df.result) <- c( "Temp", ".lower", ".upper", ".width", ".point", ".interval",  "Farm")
df.result <- df.result[-1,]

df.result.all.thetas <- data.frame(matrix(ncol = 10))
names(df.result.all.thetas) <- c( "Median_Prob_Estimate", ".lower", ".upper", ".width", ".point", ".interval",  "Farm", "Variable",  "Theta", "data")
df.result.all.thetas <- df.result.all.thetas[-1,]


#for Theta1 - calculate median
for(i in colnames(prob_theta_1)){
  df.result.temp <- median_qi(prob_theta_1[i])

  df.result.temp$Farm <- colnames(df.result.temp[1])
  names(df.result.temp) <- c( "Median_Prob_Estimate", ".lower", ".upper", ".width", ".point", ".interval",  "Farm")
  df.result <- rbind(df.result, df.result.temp)
}

df.result <- df.result%>%
  mutate(Farm = str_replace(Farm, "V", ""))%>%
  mutate(Farm = str_pad(Farm, 3, "left", pad = "0"))%>%
  mutate(Variable = v)%>%
  mutate(Theta = "No_deaths")

df.result["data"] <- dat[v]

df.result.all.thetas <- rbind(df.result.all.thetas,df.result)

#theta2 - calculate median probability

df.result <- data.frame(matrix(ncol = 7))
names(df.result) <- c( "Temp", ".lower", ".upper", ".width", ".point", ".interval",  "Farm")
df.result <- df.result[-1,]

for(i in colnames(prob_theta_2)){
  df.result.temp <- median_qi(prob_theta_2[i])

  df.result.temp$Farm <- colnames(df.result.temp[1])
  names(df.result.temp) <- c( "Median_Prob_Estimate", ".lower", ".upper", ".width", ".point", ".interval",  "Farm")
  df.result <- rbind(df.result, df.result.temp)
}

#Get farm names correct
df.result <- df.result%>%
  mutate(Farm = str_replace(Farm, "V", ""))%>%
  mutate(Farm = str_pad(Farm, 3, "left", pad = "0"))%>%
  mutate(Variable = v)%>%
  mutate(Theta = "Deaths")

df.result["data"] <- dat[v]

df.result.all.thetas <- rbind(df.result.all.thetas,df.result)

#theta3 - median probabilities

df.result <- data.frame(matrix(ncol = 7))
names(df.result) <- c( "Temp", ".lower", ".upper", ".width", ".point", ".interval",  "Farm")
df.result <- df.result[-1,]

for(i in colnames(prob_theta_3)){#Colnames are the different farms
  df.result.temp <- median_qi(prob_theta_3[i])

  df.result.temp$Farm <- colnames(df.result.temp[1])
  names(df.result.temp) <- c( "Median_Prob_Estimate", ".lower", ".upper", ".width", ".point", ".interval",  "Farm")
  df.result <- rbind(df.result, df.result.temp)
}

df.result <- df.result%>%
  mutate(Farm = str_replace(Farm, "V", ""))%>%
  mutate(Farm = str_pad(Farm, 3, "left", pad = "0"))%>%
  mutate(Variable = v)%>%
  mutate(Theta = "Abortions")

df.result["data"] <- dat[v]

#Combine all thetas
df.result.all.thetas <- rbind(df.result.all.thetas,df.result)

#Make a column that is unstandardized for plotting
df.result.all.thetas <- left_join(df.result.all.thetas, un.stand)#Join by Variable used in the analysis

#Unstandardize the data
df.result.all.thetas <- df.result.all.thetas%>%
  mutate(Un.Stand.Data = (data * SD_Used_To_Standardize) + Mean_Used_To_Standardize)

#Filter datasets by outcome
filt1 <- filter(df.result.all.thetas, Theta == "No_deaths")
filt2 <- filter(df.result.all.thetas, Theta == "Deaths")
filt3 <- filter(df.result.all.thetas, Theta == "Abortions")

filt1 <- arrange(filt1, Un.Stand.Data, Farm)
filt2 <- arrange(filt2, Un.Stand.Data, Farm)
filt3 <- arrange(filt3, Un.Stand.Data, Farm)

#Figures 3 and S4
if(j==1){

  #percents 1e-15 are messing up the scale conversion so change them to zero if <1%
  filt1 <- filt1%>%
    mutate(Un.Stand.Data = if_else(Un.Stand.Data < 0.01 & Un.Stand.Data > -0.01, 0, Un.Stand.Data))
  filt2 <- filt2%>%
    mutate(Un.Stand.Data = if_else(Un.Stand.Data < 0.01 & Un.Stand.Data > -0.01, 0, Un.Stand.Data))
  filt3 <- filt3%>%
    mutate(Un.Stand.Data = if_else(Un.Stand.Data < 0.01 & Un.Stand.Data > -0.01, 0, Un.Stand.Data))

  #Calculate the probability of death with <10% sheep and 100% sheep
  sheep.prob.less.5 <- filter(filt2, Un.Stand.Data <=0.5 & Theta == "Deaths")
  sheep.prob.100 <- filter(filt2, Un.Stand.Data >0.995 & Theta == "Deaths")

plot1 <- ggplot(df.result.all.thetas, aes(x = Un.Stand.Data))+
  geom_line(data = filt1, aes(y = Median_Prob_Estimate, colour = "Prob_Nothing"))+
  geom_ribbon(data = filt1, aes(ymin = .lower, ymax = .upper), fill = "#666666", alpha = .4)+
  geom_line(data = filt2, aes(y = Median_Prob_Estimate, colour = "Prob_Death"))+
  geom_ribbon(data = filt2, aes(ymin = .lower, ymax = .upper), fill = "#1B9E77", alpha = .5)+
  geom_line(data = filt3, aes(y = Median_Prob_Estimate, colour = "Prob_Abort"))+
  geom_ribbon(data = filt3, aes(ymin = .lower, ymax = .upper), fill = "#7570B3", alpha = .4)+
  ylim(c(0, 1))+
  scale_color_manual(name = "RVF Clinical Signs", values = c("Prob_Death" = "#1B9E77", "Prob_Abort" = "#7570B3", "Prob_Nothing" = "#666666"), labels = c("Abortions", "Deaths", "No Signs"))+
  labs( x = "Proportion Sheep", y = "Probability", color = "Legend Title\n") +
  theme_classic()+
  theme(axis.text = element_text(size = 10, colour = "#666666"),
        axis.title = element_text(size = rel(1.5)),
        plot.margin = margin(t=.3, r = 1, b = .3, l = .3, unit = "cm"),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size = 10),
        legend.position = "none")


filname <- paste0("Publication_Figures/Probability of ", v, " for all thetas_2prcts.png")
ggsave(filename = filname)

#And zoom in with only clinical signs
plot1.zoom <- ggplot(df.result.all.thetas, aes(x = Un.Stand.Data))+
  geom_line(data = filt2, aes(y = Median_Prob_Estimate, colour = "Prob_Death"))+
  geom_ribbon(data = filt2, aes(ymin = .lower, ymax = .upper), fill = "#1B9E77", alpha = .5)+
  geom_line(data = filt3, aes(y = Median_Prob_Estimate, colour = "Prob_Abort"))+
  geom_ribbon(data = filt3, aes(ymin = .lower, ymax = .upper), fill = "#7570B3", alpha = .4)+
  ylim(c(0, 0.2))+
  scale_color_manual(name = "RVF Clinical Signs", values = c("Prob_Death" = "#1B9E77", "Prob_Abort" = "#7570B3"), labels = c("Abortions", "Deaths"))+
  labs( x = "Proportion Sheep", y = "Probability", color = "Legend Title\n") +
  theme_classic()+
  theme(axis.text = element_text(size = 10, colour = "#666666"),
        axis.title = element_text(size = rel(1.5)),
        plot.margin = margin(t=.3, r = 1, b = .3, l = .3, unit = "cm"),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size = 10),
        legend.position = "none")


filname.zoom <- paste0("Publication_Figures/Probability of ", v, " for deaths and abortions only_2prcts.png")
ggsave(filename = filname.zoom)

}else{
  if(j==2){
    plot2 <- ggplot(df.result.all.thetas, aes(x = Un.Stand.Data))+
      geom_line(data = filt1, aes(y = Median_Prob_Estimate, colour = "Prob_Nothing"))+
      geom_ribbon(data = filt1, aes(ymin = .lower, ymax = .upper), fill = "#666666", alpha = .4)+
      geom_line(data = filt2, aes(y = Median_Prob_Estimate, colour = "Prob_Death"))+
      geom_ribbon(data = filt2, aes(ymin = .lower, ymax = .upper), fill = "#1B9E77", alpha = .5)+
      geom_line(data = filt3, aes(y = Median_Prob_Estimate, colour = "Prob_Abort"))+
      geom_ribbon(data = filt3, aes(ymin = .lower, ymax = .upper), fill = "#7570B3", alpha = .4)+
      scale_color_manual(name = "RVF Clinical Signs", values = c("Prob_Death" = "#1B9E77", "Prob_Abort" = "#7570B3", "Prob_Nothing" = "#666666"), labels = c("Abortions", "Deaths", "No Signs"))+
      ylim(c(0, 1))+
      labs( x = "Proportion Goats", y = "Probability", color = "Legend Title\n") +
      theme_classic()+
      theme(axis.text = element_text(size = 10, colour = "#666666"),
            axis.title = element_text(size = rel(1.5)),
            plot.margin = margin(t=.3, r = 1, b = .3, l = .3, unit = "cm"),
            legend.key.size = unit(1, "cm"),
            legend.text = element_text(size = 10),
            legend.position = "none")


    filname <- paste0("Publication_Figures/Probability of ", v, " for all thetas_2prcts.png")
    ggsave(filename = filname)

    #And zoom in with only clinical signs
    plot2.zoom <- ggplot(df.result.all.thetas, aes(x = Un.Stand.Data))+
      geom_line(data = filt2, aes(y = Median_Prob_Estimate, colour = "Prob_Death"))+
      geom_ribbon(data = filt2, aes(ymin = .lower, ymax = .upper), fill = "#1B9E77", alpha = .5)+
      geom_line(data = filt3, aes(y = Median_Prob_Estimate, colour = "Prob_Abort"))+
      geom_ribbon(data = filt3, aes(ymin = .lower, ymax = .upper), fill = "#7570B3", alpha = .4)+
      ylim(c(0, 0.75))+
      scale_color_manual(name = "RVF Clinical Signs", values = c("Prob_Death" = "#1B9E77", "Prob_Abort" = "#7570B3"), labels = c("Abortions", "Deaths"))+
      labs( x = "Proportion Goats", y = "Probability", color = "Legend Title\n") +
      theme_classic()+
      theme(axis.text = element_text(size = 10, colour = "#666666"),
            axis.title = element_text(size = rel(1.5)),
            plot.margin = margin(t=.3, r = 1, b = .3, l = .3, unit = "cm"),
            legend.key.size = unit(1, "cm"),
            legend.text = element_text(size = 10),
            legend.position = "none")


    filname.zoom <- paste0("Publication_Figures/Probability of ", v, " for deaths and abortions only_2prcts.png")
    ggsave(filename = filname.zoom)
  }else{
    if(j==3){
      plot3 <- ggplot(df.result.all.thetas, aes(x = Un.Stand.Data))+
        geom_line(data = filt1, aes(y = Median_Prob_Estimate, colour = "Prob_Nothing"))+
        geom_ribbon(data = filt1, aes(ymin = .lower, ymax = .upper), fill = "#666666", alpha = .4)+
        geom_line(data = filt2, aes(y = Median_Prob_Estimate, colour = "Prob_Death"))+
        geom_ribbon(data = filt2, aes(ymin = .lower, ymax = .upper), fill = "#1B9E77", alpha = .5)+
        geom_line(data = filt3, aes(y = Median_Prob_Estimate, colour = "Prob_Abort"))+
        geom_ribbon(data = filt3, aes(ymin = .lower, ymax = .upper), fill = "#7570B3", alpha = .4)+
        ylim(c(0, 1))+
        scale_color_manual(name = "RVF Clinical Signs", values = c("Prob_Death" = "#1B9E77", "Prob_Abort" = "#7570B3", "Prob_Nothing" = "darkcyan"), labels = c("Abortions", "Deaths", "No Signs"))+
        labs( x = "Number of Domestic Ruminants", y = "Probability", color = "Legend Title\n") +
        theme_classic()+
        theme(axis.text = element_text(size = 10, colour = "#666666"),
              axis.title = element_text(size = rel(1.5)),
              plot.margin = margin(t=.3, r = 1, b = .3, l = .3, unit = "cm"),
              legend.key.size = unit(1, "cm"),
              legend.text = element_text(size = 10),
              legend.position = "none")


      filname <- paste0("Publication_Figures/Probability of ", v, " for all thetas_2prcts.png")
      ggsave(filename = filname)

      #And zoom in with only clinical signs
      plot3.zoom <- ggplot(df.result.all.thetas, aes(x = Un.Stand.Data))+
        geom_line(data = filt2, aes(y = Median_Prob_Estimate, colour = "Prob_Death"))+
        geom_ribbon(data = filt2, aes(ymin = .lower, ymax = .upper), fill = "#1B9E77", alpha = .5)+
        geom_line(data = filt3, aes(y = Median_Prob_Estimate, colour = "Prob_Abort"))+
        geom_ribbon(data = filt3, aes(ymin = .lower, ymax = .upper), fill = "#7570B3", alpha = .4)+
        #ylim(c(0, 1))+
        scale_color_manual(name = "RVF Clinical Signs", values = c("Prob_Death" = "#1B9E77", "Prob_Abort" = "#7570B3"), labels = c("Abortions", "Deaths"))+
        labs( x = "Number of Domestic Ruminants", y = "Probability", color = "Legend Title\n") +
        theme_classic()+
        theme(axis.text = element_text(size = 10, colour = "#666666"),
              axis.title = element_text(size = rel(1.5)),
              plot.margin = margin(t=.3, r = 1, b = .3, l = .3, unit = "cm"),
              legend.key.size = unit(1, "cm"),
              legend.text = element_text(size = 10),
              legend.position = "none")


      filname.zoom <- paste0("Publication_Figures/Probability of ", v, " for deaths and abortions only_2prcts.png")
      ggsave(filename = filname.zoom)
    }else{
      if(j==4){
        plot4 <- ggplot(df.result.all.thetas, aes(x = Un.Stand.Data))+
          geom_line(data = filt1, aes(y = Median_Prob_Estimate, colour = "Prob_Nothing"))+
          geom_ribbon(data = filt1, aes(ymin = .lower, ymax = .upper), fill = "#666666", alpha = .4)+
          geom_line(data = filt2, aes(y = Median_Prob_Estimate, colour = "Prob_Death"))+
          geom_ribbon(data = filt2, aes(ymin = .lower, ymax = .upper), fill = "#1B9E77", alpha = .5)+
          geom_line(data = filt3, aes(y = Median_Prob_Estimate, colour = "Prob_Abort"))+
          geom_ribbon(data = filt3, aes(ymin = .lower, ymax = .upper), fill = "#7570B3", alpha = .4)+
          scale_color_manual(name = "RVF Clinical Signs", values = c("Prob_Death" = "#1B9E77", "Prob_Abort" = "#7570B3", "Prob_Nothing" = "#666666"), labels = c("Abortions", "Deaths", "No Signs"))+
          labs( x = "Number of Water Sources", y = "Probability", color = "Legend Title\n") +
          theme_classic()+
          theme(axis.text = element_text(size = 10, colour = "#666666"),
                axis.title = element_text(size = rel(1.5)),
                plot.margin = margin(t=.3, r = 1, b = .3, l = .3, unit = "cm"),
                legend.key.size = unit(1, "cm"),
                legend.text = element_text(size = 10),
                legend.position = "none")


        filname <- paste0("Publication_Figures/Probability of ", v, " for all thetas_2prcts.png")
        ggsave(filename = filname)

        #And zoom in with only clinical signs
        plot4.zoom <- ggplot(df.result.all.thetas, aes(x = Un.Stand.Data))+
          geom_line(data = filt2, aes(y = Median_Prob_Estimate, colour = "Prob_Death"))+
          geom_ribbon(data = filt2, aes(ymin = .lower, ymax = .upper), fill = "#1B9E77", alpha = .5)+
          geom_line(data = filt3, aes(y = Median_Prob_Estimate, colour = "Prob_Abort"))+
          geom_ribbon(data = filt3, aes(ymin = .lower, ymax = .upper), fill = "#7570B3", alpha = .4)+
          ylim(c(0, 0.2))+
          scale_color_manual(name = "RVF Clinical Signs", values = c("Prob_Death" = "#1B9E77", "Prob_Abort" = "#7570B3"), labels = c("Abortions", "Deaths"))+
          labs( x = "Number of Water Sources", y = "Probability", color = "Legend Title\n") +
          theme_classic()+
          theme(axis.text = element_text(size = 10, colour = "#666666"),
                axis.title = element_text(size = rel(1.5)),
                plot.margin = margin(t=.3, r = 1, b = .3, l = .3, unit = "cm"),
                legend.key.size = unit(1, "cm"),
                legend.text = element_text(size = 10),
                legend.position = "none")


        filname.zoom <- paste0("Publication_Figures/Probability of ", v, " for deaths and abortions only_2prcts.png")
        ggsave(filename = filname.zoom)
      }else{
        if(j==5){
        plot5 <- ggplot(df.result.all.thetas, aes(x = Un.Stand.Data))+
          geom_line(data = filt1, aes(y = Median_Prob_Estimate, colour = "Prob_Nothing"))+
          geom_ribbon(data = filt1, aes(ymin = .lower, ymax = .upper), fill = "#666666", alpha = .4)+
          geom_line(data = filt2, aes(y = Median_Prob_Estimate, colour = "Prob_Death"))+
          geom_ribbon(data = filt2, aes(ymin = .lower, ymax = .upper), fill = "#1B9E77", alpha = .5)+
          geom_line(data = filt3, aes(y = Median_Prob_Estimate, colour = "Prob_Abort"))+
          geom_ribbon(data = filt3, aes(ymin = .lower, ymax = .upper), fill = "#7570B3", alpha = .4)+
          scale_color_manual(name = "RVF Clinical Signs", values = c("Prob_Death" = "#1B9E77", "Prob_Abort" = "#7570B3", "Prob_Nothing" = "#666666"), labels = c("Abortions", "Deaths", "No Signs"))+
          labs( x = "2-Month Cumulative Rainfall (mm)\nMarch 2010", y = "Probability", color = "Legend Title\n") +
          theme_classic()+
          theme(axis.text = element_text(size = 10, colour = "#666666"),
                axis.title = element_text(size = rel(1.5)),
                legend.text = element_text(size = 10),
                legend.key.size = unit(0.5, "cm"),
                plot.margin = margin(t=.3, r = 1, b = .3, l = .3, unit = "cm"),
                legend.position = "none")


        filname <- paste0("Publication_Figures/Probability of ", v, " for all thetas_2prcts.png")
        ggsave(filename = filname)

        #And zoom in with only clinical signs
        plot5.zoom <- ggplot(df.result.all.thetas, aes(x = Un.Stand.Data))+
          geom_line(data = filt2, aes(y = Median_Prob_Estimate, colour = "Prob_Death"))+
          geom_ribbon(data = filt2, aes(ymin = .lower, ymax = .upper), fill = "#1B9E77", alpha = .5)+
          geom_line(data = filt3, aes(y = Median_Prob_Estimate, colour = "Prob_Abort"))+
          geom_ribbon(data = filt3, aes(ymin = .lower, ymax = .upper), fill = "#7570B3", alpha = .4)+
          ylim(c(0, 0.2))+
          scale_color_manual(name = "RVF Clinical Signs", values = c("Prob_Death" = "#1B9E77", "Prob_Abort" = "#7570B3"), labels = c("Abortions", "Deaths"))+
          labs( x = "2-Month Cumulative Rainfall (mm)\nMarch 2010", y = "Probability", color = "Legend Title\n") +
          theme_classic()+
          theme(axis.text = element_text(size = 10, colour = "#666666"),
                axis.title = element_text(size = rel(1.5)),
                legend.key.size = unit(1, "cm"),
                legend.text = element_text(size = 10),
                plot.margin = margin(t=.3, r = 1, b = .3, l = .3, unit = "cm"),
                legend.position = "none")


        filname.zoom <- paste0("Publication_Figures/Probability of ", v, " for deaths and abortions only_2prcts.png")
        ggsave(filename = filname)
      }else{
        if(j==6){

          plot6 <- ggplot(df.result.all.thetas, aes(x = Un.Stand.Data))+
            geom_line(data = filt1, aes(y = Median_Prob_Estimate, colour = "Prob_Nothing"))+
            geom_ribbon(data = filt1, aes(ymin = .lower, ymax = .upper), fill = "#666666", alpha = .4)+
            geom_line(data = filt2, aes(y = Median_Prob_Estimate, colour = "Prob_Death"))+
            geom_ribbon(data = filt2, aes(ymin = .lower, ymax = .upper), fill = "#1B9E77", alpha = .5)+
            geom_line(data = filt3, aes(y = Median_Prob_Estimate, colour = "Prob_Abort"))+
            geom_ribbon(data = filt3, aes(ymin = .lower, ymax = .upper), fill = "#7570B3", alpha = .4)+
            scale_color_manual(name = "RVF Clinical Signs", values = c("Prob_Death" = "#1B9E77", "Prob_Abort" = "#7570B3", "Prob_Nothing" = "#666666"), labels = c("Abortions", "Deaths", "No Signs"))+
            labs( x = "Proportion of Days with Rain \nMarch 2010", y = "Probability", color = "Legend Title\n") +
            theme_classic()+
            theme(axis.text = element_text(size = 10, colour = "#666666"),
                  axis.title = element_text(size = rel(1.5)),
                  legend.key.size = unit(1, "cm"),
                  legend.text = element_text(size = 10),
                  plot.margin = margin(t=.3, r = 1, b = .3, l = .3, unit = "cm"),
                  legend.position = "none")


          filname <- paste0("Publication_Figures/Probability of ", v, " for all thetas_2prcts.png")
          ggsave(filename = filname)

          #And zoom in with only clinical signs
          plot6.zoom<- ggplot(df.result.all.thetas, aes(x = Un.Stand.Data))+
            geom_line(data = filt2, aes(y = Median_Prob_Estimate, colour = "Prob_Death"))+
            geom_ribbon(data = filt2, aes(ymin = .lower, ymax = .upper), fill = "#1B9E77", alpha = .5)+
            geom_line(data = filt3, aes(y = Median_Prob_Estimate, colour = "Prob_Abort"))+
            geom_ribbon(data = filt3, aes(ymin = .lower, ymax = .upper), fill = "#7570B3", alpha = .4)+
            ylim(c(0, 0.2))+
            scale_color_manual(name = "RVF Clinical Signs", values = c("Prob_Death" = "#1B9E77", "Prob_Abort" = "#7570B3"), labels = c("Abortions", "Deaths"))+
            labs( x = "Proportion of Days with Rain \nMarch 2010", y = "Probability", color = "Legend Title\n") +
            theme_classic()+
            theme(axis.text = element_text(size = 10, colour = "#666666"),
                  axis.title = element_text(size = rel(1.5)),
                  legend.key.size = unit(1, "cm"),
                  legend.text = element_text(size = 10),
                  plot.margin = margin(t=.3, r = 1, b = .3, l = .3, unit = "cm"),
                  legend.position = "none")


          filname <- paste0("Publication_Figures/Probability of ", v, " for deaths and abortions only_2prcts.png")
          ggsave(filename = filname)
        }}}}}}
all.data[[j]] <- df.result.all.thetas

j <- j+1
}


#Combine plots
FigProbs <- ggarrange(plot1, plot2, plot3, plot4, plot5, plot6, common.legend = TRUE, ncol = 2, nrow = 3, labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"), widths = c(350,350))

ggexport(FigProbs, filename = "Publication_Figures/Figure S1 How probability changes with continuous parameters_2prct.png", ncol = 2, nrow = 3,
         width = 800, height = 1130)

FigProbsZoom <- ggarrange(plot1.zoom, plot2.zoom, plot3.zoom, plot4.zoom, plot5.zoom, plot6.zoom, common.legend = TRUE, ncol = 2, nrow = 3, labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"), widths = c(350,350))

ggexport(FigProbsZoom, filename = "Publication_Figures/Figure 3 How probability for deaths and abortions changes with continuous parameters_2prct.png", ncol = 2, nrow = 3,
         width = 800, height = 1130)

saveRDS(all.data, "Publication_Figures/data to make probabilty change figure for publication_2prct.Rds")

print("completed and saved data and figures")

print(paste0("The mean probability of deaths for percent of sheep <5% is ", round(mean(sheep.prob.less.5$Median_Prob_Estimate),3)))
print(paste0("The mean probability of deaths for percent of sheep 100% is ", round(mean(sheep.prob.100$Median_Prob_Estimate),3)))
