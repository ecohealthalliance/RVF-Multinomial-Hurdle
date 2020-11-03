#'Questionnaire summary
#'author: Mindy Rostal
#'date: November 2, 2020
#'

#Load libraries
library(ggpubr)#note this masks extract from rstan so need to open this library first, then rstan
library(stringr)
library(ggplot2)
library(dplyr)
library(egg)
library(survey)
library(binom)
library(tidyr)
library(DescTools)


#Load data
mydata<- read.csv("All_Files_For_Publication/data/Updated 2017 and 2015 Cross-sectional Data for Bayesian Analsysis_with_UTM.csv", header=T)
predicted.dat <- read.csv("All_Files_For_Publication/data/Predicted number of deaths and abortions and original data_for_publication_NR_hurdle_GP_for_publication_Pubmodel_2prct.csv") #Must run the Bayesian Model Assessment code first to make this file

#Must have a folder called: Publication_Figures to write to

#Make new categories/regname columns
mydata<- mydata%>%
  rename(Y = UTM_Lat)%>%
  rename(X = UTM_Long)%>%
  mutate(RVF_Anydeaths = if_else(TotalRVFDeathsDomR.2010.Stand >1, 1, 0))%>%
  mutate(RVF_Anyaborts = if_else(TotalRVFAbortionsDomR.2010.Stand >1, 1, 0))%>%
  mutate(X4b.All.Cattle.Die = X4b.AdultCattle.Die.2010 + X4b.Calves.Die.2010)%>%
  mutate(X4b.All.Sheep.Die = X4b.AdultSheep.Die.2010 + X4b.Lambs.Die.2010)%>%
  mutate(X4b.All.Goats.Die = X4b.AdultGoats.Die.2010 + X4b.GoatKids.Die.2010)%>%
  mutate(PrctCattle.Abort.2010 = X4a.Cows.Abort.2010/X2a.Total.Cattle.2010)%>%
  mutate(PrctGoats.Abort.2010 = X4a.Does.Abort.2010/X2c.Total.Goats.2010)%>%
  mutate(PrctSheep.Abort.2010 = X4a.Ewes.Abort.2010/X2b.Total.Sheep.2010)%>%
  mutate(PrctCattle.Die.2010 = X4b.All.Cattle.Die/X2a.Total.Cattle.2010)%>%
  mutate(PrctGoats.Die.2010 = X4b.All.Goats.Die/X2c.Total.Goats.2010)%>%
  mutate(PrctSheep.Die.2010 = X4b.All.Sheep.Die/X2b.Total.Sheep.2010)

# of surveys
n.surveys <- nrow(mydata)

#Median and range of various continuous variables
median.prct.rain <- median(mydata$Prct.Rain.Days.2m)
range.prct.rain <- range(mydata$Prct.Rain.Days.2m)

median.cum.rain.2010 <- median(mydata$CumSum.2m.Rain)
range.cum.rain.2010 <- range(mydata$CumSum.2m.Rain)

median.Ave.cum.rain <- median(mydata$Ave.2m.CumSum)
range.Ave.cum.rain <- range(mydata$Ave.2m.CumSum)

median.cum.rain.2009.Oct <- median(mydata$CumSum.2m.Rain_Oct)
range.cum.rain.2009.Oct <- range(mydata$CumSum.2m.Rain_Oct)

#Environment
pans <- filter(mydata, f9.Pan.Num.Ave >0)
median.num.pans <- median(pans$f9.Pan.Num.Ave)
range.num.pans <- range(pans$f9.Pan.Num.Ave)

nat <- filter(mydata, Num_Natural_WaterSources.Ave >0 )
median.water.source <- median(nat$Num_Natural_WaterSources.Ave)
range.water.source <- range(nat$Num_Natural_WaterSources.Ave)

median.farm.size <- median(mydata$f14.Size.Farm.15)
range.farm.size <- range(mydata$f14.Size.Farm.15)

#Animals
median.Tot.2010 <- median(mydata$TotalPopDomR_2010)
range.Tot.2010 <- range(mydata$TotalPopDomR_2010)

catt <- filter(mydata, X2a.Total.Cattle.2010>0)
median.Catt.2010 <- median(catt$X2a.Total.Cattle.2010)
range.Catt.2010 <- range(catt$X2a.Total.Cattle.2010)

sheep <- filter(mydata, X2b.Total.Sheep.2010>0)
median.Sheep.2010 <- median(sheep$X2b.Total.Sheep.2010)
range.Sheep.2010 <- range(sheep$X2b.Total.Sheep.2010)

median.Prct.Sheep <- median(sheep$PrctSheep_Ave.2010)
range.Prct.Sheep <- range(sheep$PrctSheep_Ave.2010)

goats <- filter(mydata, X2c.Total.Goats.2010>0)
median.Goat.2010 <- median(goats$X2c.Total.Goats.2010)
range.Goat.2010 <- range(goats$X2c.Total.Goats.2010)

median.Prct.Goats <- median(goats$PrctGoats_Ave.2010)
range.Prct.Goats <- range(goats$PrctGoats_Ave.2010)

purch <- filter(mydata, Sum.Ruminants.Gained.Ave>0)
median.Rum.Purch <- median(purch$Sum.Ruminants.Gained.Ave)
range.Rum.Purch <- range(purch$Sum.Ruminants.Gained.Ave)

dist.purch <- filter(mydata, Max.Dist.Purch.Ave>0)
median.Dist.Purch <- median(dist.purch$Max.Dist.Purch.Ave)
range.Dist.Purch <- range(dist.purch$Max.Dist.Purch.Ave)

affected <- filter(mydata, X3.Outbreak == 1)
median.Prct.Abort <- median(affected$PrctRVFAbort.2010)
range.Prct.Abort <- range(affected$PrctRVFAbort.2010)
median.Prct.Death <- median(affected$PrctRVFDeath.2010)
range.Prct.Death <- range(affected$PrctRVFDeath.2010)

sheep.affected <- mydata%>%
  filter( X3.Outbreak == 1)%>%
  filter(X4a.Ewes.Abort.2010 >1 | X4b.All.Sheep.Die > 1)

median.Num.Abort.Sheep <- median(sheep.affected$X4a.Ewes.Abort.2010)
range.Num.Abort.Sheep <- range(sheep.affected$X4a.Ewes.Abort.2010)
median.Prct.Abort.Sheep <- median(sheep.affected$PrctSheep.Abort.2010)
median.Num.Death.All.Sheep <- median(sheep.affected$X4b.All.Sheep.Die)
range.Num.Death.All.Sheep <- range(sheep.affected$X4b.All.Sheep.Die)
median.Prct.Death.All.Sheep <- median(sheep.affected$PrctSheep.Die.2010)
median.Num.Death.Adult.Sheep <- median(sheep.affected$X4b.AdultSheep.Die.2010)
range.Num.Death.Adult.Sheep <- range(sheep.affected$X4b.AdultSheep.Die.2010)
median.Num.Death.Lambs <- median(sheep.affected$X4b.Lambs.Die.2010)
range.Num.Death.Lambs <- range(sheep.affected$X4b.Lambs.Die.2010)

#####Percents
svy_design <- svydesign(id = ~1, data = mydata)

#Get Percents and confidence intervals
Wildlife.prct.ci<- svyciprop(~f5.Wildlife.Mix.15, design= svy_design, method = c("likelihood"), level = 0.95)

Vax.Before.prct.ci<- svyciprop(~X5.RVF.Vax.before.Outbreak.2010, design= svy_design, method = c("likelihood"), level = 0.95)

Vax.During.prct.ci<- svyciprop(~X6.RVF.Vax.During.Outbreak.2010, design= svy_design, method = c("likelihood"), level = 0.95)

RVFOutbreak.ci<- svyciprop(~X3.Outbreak, design= svy_design, method = c("likelihood"), level = 0.95)

n.Report.RVF <- sum(mydata$X3.Outbreak == 1)

AnyDeath.ci<- svyciprop(~RVF_Anydeaths, design= svy_design, method = c("likelihood"), level = 0.95)

AnyAbort.ci<- svyciprop(~RVF_Anyaborts, design= svy_design, method = c("likelihood"), level = 0.95)

farm.type <- c(sum(mydata$f18.Production.15 == "Small holder"), sum(mydata$f18.Production.15 == "Commercial"),sum(mydata$f18.Production.15 == "Semi-Commercial"))

farm.type.ci <- MultinomCI(farm.type,
                           conf.level=0.95,
                           method="wald")

##Table 1 Means and ranges
Tab2.Medians <- data.frame(matrix(ncol = 4, nrow = 13))

Tab2.Medians[1,] <- c("Number of cattle in 2010", median.Catt.2010, nrow(catt), paste("[", range.Catt.2010[1], ", ", range.Catt.2010[2], "]", sep = ""))
Tab2.Medians[2,] <- c("Number of sheep in 2010", median.Sheep.2010, nrow(sheep), paste("[", range.Sheep.2010[1], ", ", range.Sheep.2010[2], "]", sep = ""))
Tab2.Medians[3,] <- c("Number of goats in 2010", median.Goat.2010, nrow(goats), paste("[", range.Goat.2010[1], ", ", range.Goat.2010[2], "]", sep = ""))
Tab2.Medians[4,] <- c("Total Number of ruminants in 2010", median.Tot.2010, nrow(mydata), paste("[", range.Tot.2010[1], ", ", range.Tot.2010[2], "]", sep = ""))
Tab2.Medians[5,] <- c("Percent of sheep", round(median.Prct.Sheep*100,1), nrow(sheep), paste("[", round(range.Prct.Sheep[1]*100,1), ", ", round(range.Prct.Sheep[2]*100,1), "]", sep = ""))
Tab2.Medians[6,] <- c("Percent of goats", round(median.Prct.Goats*100,1), nrow(goats), paste("[", round(range.Prct.Goats[1]*100,1), ", ", round(range.Prct.Goats[2]*100,1), "]", sep = ""))
Tab2.Medians[7,] <- c("Number of ruminants purchased", median.Rum.Purch, nrow(purch), paste("[", range.Rum.Purch[1], ", ", range.Rum.Purch[2], "]", sep = ""))
Tab2.Medians[8,] <- c("Distance from which ruminants were purchased", median.Dist.Purch, nrow(dist.purch), paste("[", range.Dist.Purch[1], ", ", range.Dist.Purch[2], "]", sep = ""))
Tab2.Medians[9,] <- c("Number of pans ruminants can access", median.num.pans, nrow(pans), paste("[", range.num.pans[1], ", ", range.num.pans[2], "]", sep = ""))
Tab2.Medians[10,] <- c("Number of water sources rumiantns can access", median.water.source , nrow(nat), paste("[", range.water.source[1], ", ", range.water.source[2], "]", sep = ""))
Tab2.Medians[11,] <- c("Farm size Km2", median.farm.size, nrow(mydata), paste("[", range.farm.size[1], ", ", range.farm.size[2], "]", sep = ""))
Tab2.Medians[12,] <- c("Percent of days with precipitation between January 15-March 15 2010 (%)", round(median.prct.rain*100,1), nrow(mydata), paste("[", round(range.prct.rain[1]*100,1), ", ", round(range.prct.rain[2]*100,1), "]", sep = ""))
Tab2.Medians[13,] <- c("2-Month cumulative precipitation between January 15-March 15 2010 (mm)", round(median.cum.rain.2010,1) , nrow(mydata), paste("[", round(range.cum.rain.2010 [1],1), ", ", round(range.cum.rain.2010 [2],1), "]", sep = ""))
Tab2.Medians[14,] <- c("2-Month cumulative precipitation between October 15-December 15 2009 (mm)", round(median.Ave.cum.rain,1) , nrow(mydata), paste("[", round(range.Ave.cum.rain[1],1), ", ", round(range.Ave.cum.rain[2],1), "]", sep = ""))

#Rename columns
names(Tab2.Medians) <- c("Variable", "Median", "Number of Farms (N = 120)", "Range")

#write file
write.csv(Tab2.Medians, "Publication_Figures/Table 1 Descriptive Stats - medians.csv", row.names = FALSE)

#Table 1 percents and CIs
Tab4.Percents <- data.frame(matrix(ncol = 4, nrow = 9))
Tab4.Percents[1,] <- c("Wildlife Mix with Livestock", round(Wildlife.prct.ci[1]*100,1), nrow(mydata), paste("[", round(attr(Wildlife.prct.ci, "ci")[[1]]*100,1), ", ", round(attr(Wildlife.prct.ci, "ci")[[2]]*100,1), "]", sep = ""))
Tab4.Percents[2,] <- c("Vaccinated Before RVF Outbreak", round(Vax.Before.prct.ci[1]*100,1), nrow(mydata), paste("[", round(attr(Vax.Before.prct.ci, "ci")[[1]]*100,1), ", ", round(attr(Vax.Before.prct.ci, "ci")[[2]]*100,1), "]", sep = ""))
Tab4.Percents[3,] <- c("Vaccinated During RVF Outbreak", round(Vax.During.prct.ci[1]*100,1), nrow(mydata), paste("[", round(attr(Vax.During.prct.ci, "ci")[[1]]*100,1), ", ", round(attr(Vax.During.prct.ci, "ci")[[2]]*100,1), "]", sep = ""))
Tab4.Percents[4,] <- c("Affected by RVF Outbreak", round(RVFOutbreak.ci[1]*100,1), nrow(mydata), paste("[", round(attr(RVFOutbreak.ci, "ci")[[1]]*100,1), ", ", round(attr(RVFOutbreak.ci, "ci")[[2]]*100,1), "]", sep = ""))
Tab4.Percents[5,] <- c("Any Reported Deaths from RVF", round(AnyDeath.ci[1]*100,1), nrow(mydata), paste("[", round(attr(AnyDeath.ci, "ci")[[1]]*100,1), ", ", round(attr(AnyDeath.ci, "ci")[[2]]*100,1), "]", sep = ""))
Tab4.Percents[6,] <- c("Any Reported Abortions from RVF", round(AnyAbort.ci[1]*100,1), nrow(mydata), paste("[", round(attr(AnyAbort.ci, "ci")[[1]]*100,1), ", ", round(attr(AnyAbort.ci, "ci")[[2]]*100,1), "]", sep = ""))
Tab4.Percents[7,] <- c("Small Holder Production", round(farm.type.ci[1]*100,1), nrow(mydata), paste("[", round(farm.type.ci[4]*100,1), ", ", round(farm.type.ci[7]*100,1), "]", sep = ""))
Tab4.Percents[8,] <- c("Semi-Commercial Production", round(farm.type.ci[3]*100,1), nrow(mydata), paste("[", round(farm.type.ci[6]*100,1), ", ", round(farm.type.ci[9]*100,1), "]", sep = ""))
Tab4.Percents[9,] <- c("Commercial Production", round(farm.type.ci[2]*100,1), nrow(mydata), paste("[", round(farm.type.ci[5]*100,1), ", ", round(farm.type.ci[8]*100,1), "]", sep = ""))

names(Tab4.Percents) <- c("Variable", "Percent", "Number of Farms (N = 120)", "95% Confidence Interval")

write.csv(Tab4.Percents, "Publication_Figures/Table 2 Descriptive Stats - percents.csv", row.names = FALSE)

#Percent death and abortions
#Figure

#Gather into long dataset
Liv <- gather(mydata, "Effect", "Number", X4b.All.Cattle.Die:PrctSheep.Die.2010, X4a.Cows.Abort.2010:X4b.GoatKids.Die.2010)#X4b.All.Sheep.Die, X4b.All.Goats.Die, , PrctCattle.Abort.2010:) #,

#Group by species and age
Liv <- Liv%>%
  select(Farm_ID,Effect, Number, X2a.Total.Cattle.2010, X2b.Total.Sheep.2010, X2c.Total.Goats.2010)%>%
  mutate(Species = if_else(str_detect(Effect, "Goat") | str_detect(Effect, "Doe"), "Goats", if_else(str_detect(Effect, "Cattle") | str_detect(Effect, "Cow")| str_detect(Effect, "Calv"), "Cattle", "Sheep")))%>%
  mutate(Calc = if_else(str_detect(Effect, "Prct"), "Percent", "Number" ))%>%
  mutate(Age = if_else(str_detect(Effect, "Lamb")| str_detect(Effect, "Kid")| str_detect(Effect, "Calves"), "Juvenile", "Adult" ))%>%
  mutate(Morb = if_else(str_detect(Effect, "Die"), "Death", "Abort" ))

#Include only percents on farms that had more than one animal of that species
Liv.Prct <- Liv%>%
  filter(Calc == "Percent")%>%
  filter(Number >0)

#Set factor levels for figure
Liv.Prct$Species <- factor(Liv.Prct$Species, levels = c("Sheep", "Cattle", "Goats"))

##################################
##################################
#Figure 1

#Count number of farms with deaths and abortions
Fig2dat <- filter(Liv, !Calc == "Percent")
Fig2dat <- filter(Fig2dat, Effect == "X4b.All.Cattle.Die" | Effect == "X4b.All.Sheep.Die" | Effect == "X4b.All.Goats.Die" | Effect == "X4a.Cows.Abort.2010" | Effect == "X4a.Ewes.Abort.2010" | Effect == "X4a.Does.Abort.2010" )
Fig2dat <- Fig2dat%>%
  select(Farm_ID, Species, Number, Morb, X2a.Total.Cattle.2010, X2b.Total.Sheep.2010, X2c.Total.Goats.2010)%>%
  mutate(OnFarm = if_else(Species == "Cattle" & X2a.Total.Cattle.2010 == 0, "NA", if_else(Species == "Sheep" & X2b.Total.Sheep.2010 == 0, "NA", if_else(Species == "Goats" & X2c.Total.Goats.2010 == 0, "NA", "Yes"))))%>%
  mutate(OnFarm = replace(OnFarm, OnFarm == "NA", NA))
#Make dataframe wide
Fig2dat <- pivot_wider(Fig2dat, id_cols = c(Farm_ID, Species), names_from=Morb, values_from = c(Number, OnFarm))

#Rename columns
Fig2dat <- Fig2dat%>%
  rename(Death = Number_Death)%>%
  rename(Abort = Number_Abort)%>%
  rename(OnFarm = OnFarm_Death)%>%
  select(-OnFarm_Abort)

#Calculate slope of lines for abortions and deaths
predicted.dat <- predicted.dat%>%
  mutate(prct.predicted.a = Median.Predicted.Abortions/TotalPopDomR_2010)%>%
  mutate(prct.predicted.d = Median.Predicted.Deaths/TotalPopDomR_2010)

#Join data
Fig2dat <- full_join(predicted.dat, Fig2dat)

#Columns for Fig 1a
Fig2dat <- Fig2dat%>%
  mutate(Predicted = if_else(Median.Predicted.Abortions >0 | Median.Predicted.Deaths > 0, "Yes", "No"))%>%
  mutate(Reported = if_else(Abort >0 | Death > 0, "Yes", "No"))

#Set species for facets
Fig2dat$Species <- factor(Fig2dat$Species, levels = c("Sheep", "Cattle", "Goats"))

#Figure 1b
Fig2.Num.DA <- ggplot(Fig2dat, aes(x=Abort, y=Death)) +
  geom_point() +
  facet_wrap(~Species)+
  stat_smooth(aes(x=Median.Predicted.Abortions, y = Median.Predicted.Deaths), method = "lm", formula = y~x)+
  labs(x = "Number of Abortions", y = "Number of Deaths")+
  theme(aspect.ratio = 1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=16, colour = "black"),
        axis.text = element_text(colour = "black"))

#Figure 1a
Fig2dat1 <- gather(Fig2dat, "RVFonFarm", "PredictedorReported", Predicted, Reported)

#For farms with RVF
Fig2dat1 <- Fig2dat1%>%
  filter(!is.na(OnFarm))

Fig2dat1$RVFonFarm <- factor(Fig2dat1$RVFonFarm, levels = c("Reported", "Predicted"))

#Set colors
my.cols <- c( "#66A61E", "#E7298A")
my.cols2 <- c("#66A61E", "#E6AB02")

#Figure 1a
Fig2.RVF.YN <- ggplot(Fig2dat1, aes(x=RVFonFarm, fill= PredictedorReported)) +
  geom_bar()+
  scale_fill_manual(values = my.cols) + #palette = "Dark2", values = c(.5)) + #color = c("#E7298A", "#66A61E")
  facet_wrap(~ Species)+
  labs(x="", y = "Number of Farms", fill = "RVF Present")+ #x = "Estimated By",
  theme(aspect.ratio = 1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=16, colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.text.x= element_text(angle = 45, hjust = 1))

#Combine into Figure 1
Fig5_2 <- ggarrange(Fig2.RVF.YN, Fig2.Num.DA, ncol = 1, nrow = 2, labels = c("(a)", "(b)"))

#Name file
fil.name5.2 <- "Publication_Figures/Fig 1 Number of animals died and aborted reported and predicted.png"

#export
ggexport(Fig5_2, filename = fil.name5.2, width = 613, height = 450)
