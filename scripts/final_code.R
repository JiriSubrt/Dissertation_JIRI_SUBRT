##### Dissertation QHI plants ######
# Author: Jiri Subrt          
# Date: 30/2/2022
# Some code snippets were taken from other people, I indicate that in the code
####################################

##### Libraries #####
library(rgdal)
library(raster)
library(sp)
library(rgbif)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(stringr)
library(devtools)
library(countrycode)
library(rnaturalearthdata)
library(CoordinateCleaner)
library(stargazer)
library(purrr)
library(bbplot)
library(MCMCglmm)

##### Functions #####

# Calculate the number of words in a string
nwords <- function(string, pseudo=F){
  ifelse( pseudo, 
          pattern <- "\\S+", 
          pattern <- "[[:alpha:]]+" 
  )
  str_count(string, pattern)
}

# Function for custom ggplot theme
theme_jiri <- function(){
  theme_bw()+
    theme(axis.text.x = element_text(size = 30),    
          axis.text.y = element_text(size = 28),
          plot.title=element_text(size= 30),
          axis.title = element_text(size = 30, face = "plain"),
          #text = element_text(family = "Arial"),
          panel.grid = element_blank(),                                              
          plot.margin = unit(c(1,1,1,1), units = , "cm"),                
          legend.text = element_text(size = 14),  
          legend.title = element_text(size = 14),
          legend.position = c(0.9, 0.2)) 
}

# Function for clean MCMC code outputs
# Taken from Myers-Smith, 2019
clean.MCMC <- function(x) {
  sols <- summary(x)$solutions  # pull out relevant info from model summary
  Gcovs <- summary(x)$Gcovariances
  Rcovs <- summary(x)$Rcovariances
  
  fixed <- data.frame(row.names(sols), sols, row.names = NULL)  # convert to dataframes with the row.names as the first col
  random <- data.frame(row.names(Gcovs), Gcovs, row.names = NULL)
  residual <- data.frame(row.names(Rcovs), Rcovs, row.names = NULL)
  
  names(fixed)[names(fixed) == "row.names.sols."] <- "variable"  # change the columns names to variable, so they all match
  names(random)[names(random) == "row.names.Gcovs."] <- "variable"
  names(residual)[names(residual) == "row.names.Rcovs."] <- "variable"
  
  fixed$effect <- "fixed"  # add ID column for type of effect (fixed, random, residual)
  random$effect <- "random"
  residual$effect <- "residual"
  
  modelTerms <- as.data.frame(bind_rows(fixed, random, residual))  # merge it all together
}

convert.proportional.warm.niche.cti <- function(location, climate){
  return_df <- left_join(location, climate, by = c("Species" = "species_name"))
  relcover_climate_warm <- return_df %>% 
    na.omit()
  relcover_climate_warm$TempWarmProp <- relcover_climate_warm$cover * relcover_climate_warm$temp_mean_warmest_q /100
  proportional_warm_niche <- relcover_climate_warm[, c("PLOT", "year", "Species", "TempWarmProp", "sub_name", "site_name")] %>%
    group_by(site_name, sub_name, PLOT, year) %>%
    dplyr::summarise(CTI = sum(TempWarmProp))
  return(proportional_warm_niche)
}

# Function adding info about models
# how many plots are positive and negative
# which are sig and nonsig as characters
add.additional.info <- function(subsite_model){
  w_added_data <- subsite_model %>% 
    mutate(Sign =
             case_when(estimate > 0 ~ "Positive", 
                       estimate < 0 ~ "Negative")) %>% 
    mutate(Significant =
             case_when(p.value > 0.05 ~ "Non-significant",
                       p.value < 0.05 ~ "Significant"))
  return(w_added_data)
}

# Function that returns the minimul and maximum CTI value
slice.min.max <- function(subsite_CTI){
  slice_minimum <- subsite_CTI %>% 
    group_by(Plot) %>% 
    slice_max(n = 1, year)
  
  slice_maximum <- subsite_CTI %>% 
    group_by(Plot) %>% 
    slice_min(n = 1, year)
  
  slice_min_max <- bind_rows(slice_minimum, slice_maximum)
  return(slice_min_max)
}

# Specifying prior for MCMC models
# Taken from Myers-Smith et al, 2019
prior2 <- list(R = list(V = 1, nu = 0.002), 
               G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.v = 10000), 
                        G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.v = 10000)))


##### PART 1: Vegetation change at QHI #####

# Read data
QHI_new_cover <- read.csv("data/qhi/QHI_cover_1999_2019_sas.csv")

# Initial species list
initial_species <- as.data.frame(c(unique(QHI_new_cover$Species)), c(length(unique(QHI_new_cover$Species))))
write.csv(initial_species, "outputs/species_lists/initial_species_list.csv")

# Data wrangling

# Changing species names to same format (Code from Mariana Garcia Criado)
QHI_new_cover$Species[QHI_new_cover$Species == "ERIVAG"] <- "Eriophorum vaginatum"
QHI_new_cover$Species[QHI_new_cover$Species == "SALPHL"] <- "Salix phlebophylla"
QHI_new_cover$Species[QHI_new_cover$Species == "SALRET"] <- "Salix reticulata"
QHI_new_cover$Species[QHI_new_cover$Species == "THASUB"] <- "Thamnolia subuliformis"
QHI_new_cover$Species[QHI_new_cover$Species == "SALPUL"] <- "Salix pulchra"
QHI_new_cover$Species[QHI_new_cover$Species == "CETCUC"] <- "Cetraria cucculata"
QHI_new_cover$Species[QHI_new_cover$Species == "SENATR"] <- "Senecio astropurpureus"
QHI_new_cover$Species[QHI_new_cover$Species == "Senecio atropurpureus"] <- "Senecio astropurpureus"
QHI_new_cover$Species[QHI_new_cover$Species == "CETISL"] <- "Cetraria islandica"
QHI_new_cover$Species[QHI_new_cover$Species == "SAUANG"] <- "Saussurea angustifolia"
QHI_new_cover$Species[QHI_new_cover$Species == "CLAMIT"] <- "Cladina mitis"
QHI_new_cover$Species[QHI_new_cover$Species == "POLBIS"] <- "Polygonum bistorta"
QHI_new_cover$Species[QHI_new_cover$Species == "FESBAF"] <- "Festuca baffinensis"
QHI_new_cover$Species[QHI_new_cover$Species == "PEDLAN"] <- "Pedicularis lanata"
QHI_new_cover$Species[QHI_new_cover$Species == "DRYINT"] <- "Dryas integrifolia"
QHI_new_cover$Species[QHI_new_cover$Species == "DACARC"] <- "Dactylina arctica"
QHI_new_cover$Species[QHI_new_cover$Species == "POAARC"] <- "Poa arctica"
QHI_new_cover$Species[QHI_new_cover$Species == "HIEALP"] <- "Hierochloe alpine"
QHI_new_cover$Species[QHI_new_cover$Species == "Hierochole alpine"] <- "Hierochloe alpine"
QHI_new_cover$Species[QHI_new_cover$Species == "PEDCAP"] <- "Pedicularis capitata"
QHI_new_cover$Species[QHI_new_cover$Species == "ARCLAT"] <- "Arctagrostis latifolia"
QHI_new_cover$Species[QHI_new_cover$Species == "PYRGRA"] <- "Pyrola grandiflora"
QHI_new_cover$Species[QHI_new_cover$Species == "SALARC"] <- "Salix arctica"
QHI_new_cover$Species[QHI_new_cover$Species == "LUPARC"] <- "Lupinus arcticus"
QHI_new_cover$Species[QHI_new_cover$Species == "POLVIV"] <- "Polygonum viviparum"
QHI_new_cover$Species[QHI_new_cover$Species == "PEDSUD"] <- "Pedicularis sudetica"
QHI_new_cover$Species[QHI_new_cover$Species == "PEDVER"] <- "Pedicularis vertisilata"
QHI_new_cover$Species[QHI_new_cover$Species == "CASTET"] <- "Cassiope tetragona"
QHI_new_cover$Species[QHI_new_cover$Species == "STELON"] <- "Stellaria longipes"
QHI_new_cover$Species[QHI_new_cover$Species == "ERIANG"] <- "Eriophorum angustifolium"
QHI_new_cover$Species[QHI_new_cover$Species == "OXYMAY"] <- "Oxytropis maydalliana"
QHI_new_cover$Species[QHI_new_cover$Species == "OXYCAM"] <- "Oxytropis campestris"
QHI_new_cover$Species[QHI_new_cover$Species == "ASTUMB"] <- "Astralagus umbelletus"
QHI_new_cover$Species[QHI_new_cover$Species == "PARNUD"] <- "Parrya nudicalis"
QHI_new_cover$Species[QHI_new_cover$Species == "ALEOCH"] <- "Alectoria ochroleuca"
QHI_new_cover$Species[QHI_new_cover$Species == "LAGGLA"] <- "Lagostis glauca"
QHI_new_cover$Species[QHI_new_cover$Species == "CARDIG"] <- "Cardamine digitalis"
QHI_new_cover$Species[QHI_new_cover$Species == "SAXNEL"] <- "Saxifraga nelsoniana"
QHI_new_cover$Species[QHI_new_cover$Species == "ALOALP"] <- "Alopecurus alpinus"
QHI_new_cover$Species[QHI_new_cover$Species == "LUZARC"] <- "Luzula arctica"
QHI_new_cover$Species[QHI_new_cover$Species == "Luzula arctica "] <- "Luzula arctica"
QHI_new_cover$Species[QHI_new_cover$Species == "OXYNIG"] <- "Oxytropis nigrescens"
QHI_new_cover$Species[QHI_new_cover$Species == "POAALP"] <- "Poa alpina"
QHI_new_cover$Species[QHI_new_cover$Species == "CETNIV"] <- "Cetraria nivalis"
QHI_new_cover$Species[QHI_new_cover$Species == "CARBEL"] <- "Cardamine bellidifolia"
QHI_new_cover$Species[QHI_new_cover$Species == "PAPRAD"] <- "Papaver radicatum"
QHI_new_cover$Species[QHI_new_cover$Species == "Poa ?"] <- "XXXPOA"
QHI_new_cover$Species[QHI_new_cover$Species == "Cladina (brown)"] <- "XXXBROWNCLADINA"
QHI_new_cover$Species[QHI_new_cover$Species == "Cetraria Species that is brown"] <- "XXXBROWNCETRARIA"

# Remove all species that have XXX prefix or are not at species level
# Removing: 1) XXXspecies 2) Non-complete species 3) non-identifiable species (e.g., "graminoid")
QHI_new_cover <- QHI_new_cover %>% 
  dplyr::filter(nwords(QHI_new_cover$Species) == 2 | nwords(QHI_new_cover$Species) == 3) %>% 
  filter(Species != "Cetraria spp")

# Species list clean
clean_species <- as.data.frame(c(unique(QHI_new_cover$Species)), c(length(unique(QHI_new_cover$Species))))
write.csv(clean_species, "outputs/species_lists/clean_species_list.csv")

# Filling in zeros (Code snippet adapted from Myers-Smith, 2019)
QHI_new_cover_zeros <- QHI_new_cover %>% 
  tidyr::complete(PLOT, sub_name, year, Species, fill = list(cover = 0)) %>% 
  group_by(Species, year, sub_name, PLOT) %>%
  summarise(cover = mean(cover)) %>% 
  ungroup()

## Wrangling of Herschel vegetation type ##
# Changing str() of variables
QHI_new_cover_HE <- subset(QHI_new_cover_zeros, sub_name == "QHI:HE")
QHI_new_cover_HE$PLOT <- as.factor(QHI_new_cover_HE$PLOT)
QHI_new_cover_HE$Species <- as.factor(QHI_new_cover_HE$Species)

# Removing species that are not in HE that are in KO
QHI_only_HE <- QHI_new_cover %>% 
  filter(sub_name == "QHI:HE")

# Write species only in HE
species_HE_only_list <- as.data.frame(c(unique(QHI_only_HE$Species)), c(length(unique(QHI_only_HE$Species))))
write.csv(species_HE_only_list, "outputs/species_lists/HE_species_list.csv")

## Wrangling of Komakuk vegetation type ##
# Changing str() of variables
QHI_new_cover_KO <- subset(QHI_new_cover_zeros, sub_name == "QHI:KO")
QHI_new_cover_KO$Species <- as.factor(QHI_new_cover_KO$Species)
QHI_new_cover_KO$cover <- as.numeric(QHI_new_cover_KO$cover)
QHI_new_cover_KO$PLOT <- as.factor(QHI_new_cover_KO$PLOT)
unique(QHI_new_cover_KO$year)

# Removing species that are not in KO that are in HE
QHI_only_KO <- QHI_new_cover %>% 
  filter(sub_name == "QHI:KO")

# Write species only in KO 
species_KO_only_list <- as.data.frame(c(unique(QHI_only_KO$Species)), c(length(unique(QHI_only_KO$Species))))
write.csv(species_KO_only_list, "outputs/species_lists/KO_species_list.csv")

##### Models for Vegetation change #####
# Code for MCMC models and calculating model predictions adapted from Myers-Smith et al., 2019

# Herschel
HE_plot_linear <- MCMCglmm(cover ~ year * Species - 1, 
                               random = ~ year + PLOT, data = QHI_new_cover_HE, 
                               family = "gaussian", pr = TRUE, nitt = 100000, burnin = 10000,thin = 30,
                               prior = prior2)

# Summary and checking trace plots
summary(HE_plot_linear) 
plot(HE_plot_linear$Sol)
plot(HE_plot_linear$VCV)

# Significant species
# *** Salix pulchra, Poa arctica, Eriophorum vaginatum, Arctagrostis latifolia   
# * Salix reticulata, Cetraria cucculata

# MCMC for HE with significant species

# Herschel Salix pulchra
HE_plots_sp <- subset(QHI_new_cover_HE, Species == "Salix pulchra")
HE_plots_sp <- as.data.frame(HE_plots_sp)

# Column for the trial hits and misses
HE_plots_sp$failures <- round(100 - HE_plots_sp$cover)

HE_plot_sp <- MCMCglmm(cbind(round(cover), failures) ~ I(year - 1998), 
                       random = ~ year + PLOT, data = HE_plots_sp, 
                       family = "multinomial2", pr = TRUE, nitt = 100000, 
                       burnin = 20000, prior = prior2)

# Predictions Sal pulchra
nyears <- 21
niter <- length(HE_plot_sp$Sol[,"(Intercept)"])

HE_plot_preds_salpul <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    HE_plot_preds_salpul[i,j] <- (HE_plot_sp$Sol[i,"(Intercept)"] + HE_plot_sp$Sol[i,"I(year - 1998)"]*j)
  }
}

HE_plot_preds_df_salpul <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  HE_plot_preds_df_salpul[i,] <- quantile(HE_plot_preds_salpul[,i], c(0.025, 0.5, 0.975))
}

HE_plot_preds_df_salpul <- cbind.data.frame(lower = HE_plot_preds_df_salpul[,1], 
                                            mean = HE_plot_preds_df_salpul[,2], upper = HE_plot_preds_df_salpul[,3], year = seq(1:21))

# Approximate percentage change = 17.2% 
(plogis(((-1.5867173) + 0.0457232*21))*100) - (plogis(((-1.5867173) + 0.0457232*1))*100) 

# Herschel Poa arctica
HE_plots_pa <- subset(QHI_new_cover_HE, Species == "Poa arctica")
HE_plots_pa <- as.data.frame(HE_plots_pa)

# Column for the trial hits and misses
HE_plots_pa$failures <- round(100 - HE_plots_pa$cover)

HE_plot_pa <- MCMCglmm(cbind(round(cover), failures) ~ I(year - 1998), 
                       random = ~ year + PLOT, data = HE_plots_pa, 
                       family = "multinomial2", pr = TRUE, nitt = 100000, 
                       burnin = 20000, prior = prior2)

# Predictions Poa arctica
nyears <- 21
niter <- length(HE_plot_pa$Sol[,"(Intercept)"])

HE_plot_preds_poar <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    HE_plot_preds_poar[i,j] <- 100*plogis(HE_plot_pa$Sol[i,"(Intercept)"] + HE_plot_pa$Sol[i,"I(year - 1998)"]*j)
  }
}

HE_plot_preds_df_poar <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  HE_plot_preds_df_poar[i,] <- quantile(HE_plot_preds_poar[,i], c(0.025, 0.5, 0.975))
}

HE_plot_preds_df_poar <- cbind.data.frame(lower = HE_plot_preds_df_poar[,1], 
                                          mean = HE_plot_preds_df_poar[,2], upper = HE_plot_preds_df_poar[,3], year = seq(1:21))


# Approximate percentage change = 5.3
(plogis(((-5.23434) + 0.11804*21))*100) - (plogis(((-5.23434) + 0.11804*1))*100) 

# Herschel Eriophorum vaginatum
HE_plots_ev <- subset(QHI_new_cover_HE, Species == "Eriophorum vaginatum")
HE_plots_ev <- as.data.frame(HE_plots_ev)

# Column for the trial hits and misses
HE_plots_ev$failures <- round(100 - HE_plots_ev$cover)

HE_plot_ev <- MCMCglmm(cbind(round(cover), failures) ~ I(year - 1998), 
                       random = ~ year + PLOT, data = HE_plots_ev, 
                       family = "multinomial2", pr = TRUE, nitt = 100000, 
                       burnin = 20000, prior = prior2)

summary(HE_plot_ev)

# Predictions E vaginatum
nyears <- 21
niter <- length(HE_plot_ev$Sol[,"(Intercept)"])

HE_plot_preds_ev <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    HE_plot_preds_ev[i,j] <- 100*plogis(HE_plot_ev$Sol[i,"(Intercept)"] + HE_plot_ev$Sol[i,"I(year - 1998)"]*j)
  }
}

HE_plot_preds_df_ev <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  HE_plot_preds_df_ev[i,] <- quantile(HE_plot_preds_ev[,i], c(0.025, 0.5, 0.975))
}

HE_plot_preds_df_ev <- cbind.data.frame(lower = HE_plot_preds_df_ev[,1], 
                                        mean = HE_plot_preds_df_ev[,2], upper = HE_plot_preds_df_ev[,3], year = seq(1:21))


# Approximate percentage change 
(plogis(((0.15800) + 0.07312*21))*100) - (plogis(((0.15800) + 0.07312*1))*100) 

# Herschel Arctagrostis latifolia
HE_plots_al <- subset(QHI_new_cover_HE, Species == "Arctagrostis latifolia")
HE_plots_al <- as.data.frame(HE_plots_al)

# Column for the trial hits and misses
HE_plots_al$failures <- round(100 - HE_plots_al$cover)

HE_plot_al <- MCMCglmm(cbind(round(cover), failures) ~ I(year - 1998), 
                       random = ~ year + PLOT, data = HE_plots_al, 
                       family = "multinomial2", pr = TRUE, nitt = 100000, 
                       burnin = 20000, prior = prior2)

# Predictions A latifolia
nyears <- 21
niter <- length(HE_plot_al$Sol[,"(Intercept)"])

HE_plot_preds_al <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    HE_plot_preds_al[i,j] <- 100*plogis(HE_plot_al$Sol[i,"(Intercept)"] + HE_plot_al$Sol[i,"I(year - 1998)"]*j)
  }
}

HE_plot_preds_df_al <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  HE_plot_preds_df_al[i,] <- quantile(HE_plot_preds_al[,i], c(0.025, 0.5, 0.975))
}

HE_plot_preds_df_al <- cbind.data.frame(lower = HE_plot_preds_df_al[,1], 
                                        mean = HE_plot_preds_df_al[,2], upper = HE_plot_preds_df_al[,3], year = seq(1:21))


# Approximate percentage change = 15.6%
(plogis(((-4.43682) + 0.13570*21))*100) -(plogis(((-4.43682) + 0.13570*1))*100) 


# Herschel Salix reticulata
HE_plots_sr <- subset(QHI_new_cover_HE, Species == "Salix reticulata")
HE_plots_sr <- as.data.frame(HE_plots_sr)

# Column for the trial hits and misses
HE_plots_sr$failures <- round(100 - HE_plots_sr$cover)

HE_plot_sr <- MCMCglmm(cbind(round(cover), failures) ~ I(year - 1998), 
                       random = ~ year + PLOT, data = HE_plots_sr, 
                       family = "multinomial2", pr = TRUE, nitt = 100000, 
                       burnin = 20000, prior = prior2)

# Predictions Salix reticulata
nyears <- 21
niter <- length(HE_plot_sr$Sol[,"(Intercept)"])

HE_plot_preds_sr <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    HE_plot_preds_sr[i,j] <- 100*plogis(HE_plot_sr$Sol[i,"(Intercept)"] + HE_plot_sr$Sol[i,"I(year - 1998)"]*j)
  }
}

HE_plot_preds_df_sr <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  HE_plot_preds_df_sr[i,] <- quantile(HE_plot_preds_sr[,i], c(0.025, 0.5, 0.975))
}

HE_plot_preds_df_sr <- cbind.data.frame(lower = HE_plot_preds_df_sr[,1], 
                                        mean = HE_plot_preds_df_sr[,2], upper = HE_plot_preds_df_sr[,3], year = seq(1:21))

# Approximate percentage change = 3.8%
(plogis(((-3.30660) + 0.03801*21))*100) -(plogis(((-3.30660) + 0.03801*1))*100) 

# Herschel Cetraria cucculata
HE_plots_cc <- subset(QHI_new_cover_HE, Species == "Cetraria cucculata")
HE_plots_cc <- as.data.frame(HE_plots_cc)

# Creating a column for the trial hits and misses
HE_plots_cc$failures <- round(100 - HE_plots_cc$cover)

HE_plot_cc <- MCMCglmm(cbind(round(cover), failures) ~ I(year - 1998), 
                       random = ~ year + PLOT, data = HE_plots_cc, 
                       family = "multinomial2", pr = TRUE, nitt = 100000, 
                       burnin = 20000, prior = prior2)

# Predictions Cetraria cucculata
nyears <- 21
niter <- length(HE_plot_cc$Sol[,"(Intercept)"])

HE_plot_preds_cc <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    HE_plot_preds_cc[i,j] <- 100*plogis(HE_plot_cc$Sol[i,"(Intercept)"] + HE_plot_cc$Sol[i,"I(year - 1998)"]*j)
  }
}

HE_plot_preds_df_cc <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  HE_plot_preds_df_cc[i,] <- quantile(HE_plot_preds_cc[,i], c(0.025, 0.5, 0.975))
}

HE_plot_preds_df_cc <- cbind.data.frame(lower = HE_plot_preds_df_cc[,1], 
                                        mean = HE_plot_preds_df_cc[,2], upper = HE_plot_preds_df_cc[,3], year = seq(1:21))

# Approximate percentage change = -3.5
(plogis(((-3.0684172) + (-0.1069125)*21))*100) -(plogis(((-3.0684172) + (-0.1069125)*1))*100) 


#### Herschel significant vegetation plots ####

(herschel_veg_change_plot_salpul <- ggplot() +
    geom_point(data = subset(QHI_new_cover_HE, Species == "Salix pulchra"), 
               aes(x = year, y = cover, colour = factor(Species)), alpha = 0.8, size = 4, show.legend = FALSE) +
    scale_color_manual(values = c("#ffcd44"), name = "") +
    scale_fill_manual(values = c("#ffcd44"), name = "") +
    scale_x_continuous(breaks = c(1999, 2004, 2009, 2014, 2019)) +
    geom_ribbon(data = HE_plot_preds_df_salpul, aes(x = year + 1998, ymin = lower, ymax = upper), fill = "#ffcd44", alpha = 0.2) +
    geom_line(data = HE_plot_preds_df_salpul, aes(x = year+ 1998, y = mean), colour = "#ffcd44") +
   ggtitle("(a) Salix pulchra") +
    theme_jiri()+
    theme(legend.position = c(0.275, 0.95), 
          axis.line.x = element_line(color = "black", size = 0.5),
          axis.line.y = element_line(color = "black", size = 0.5)) +
    labs(x = "", y = "Species cover (% in plot)\n"))

(herschel_veg_change_plot_poar <- ggplot() +
    geom_point(data = subset(QHI_new_cover_HE, Species == "Poa arctica"), 
               aes(x = year, y = cover, colour = factor(Species)), alpha = 0.8, size = 4, show.legend = FALSE) +
    scale_color_manual(values = c("#EE00EE"), name = "") +
    scale_fill_manual(values = c("#EE00EE"), name = "") +
    scale_x_continuous(breaks = c(1999, 2004, 2009, 2014, 2019)) +
    geom_ribbon(data = HE_plot_preds_df_poar, aes(x = year + 1998, ymin = lower, ymax = upper), fill = "#EE00EE", alpha = 0.2) +
    geom_line(data = HE_plot_preds_df_poar, aes(x = year+ 1998, y = mean), colour = "#EE00EE") +
    ggtitle("(b) Poa arctica") +
    theme_jiri()+
    theme(legend.position = c(0.275, 0.95), 
          axis.line.x = element_line(color = "black", size = 0.5),
          axis.line.y = element_line(color = "black", size = 0.5)) +
    labs(x = "", y = "Species cover (% in plot)\n"))

(herschel_veg_change_plot_ev <- ggplot() +
    geom_point(data = subset(QHI_new_cover_HE, Species == "Eriophorum vaginatum"), 
               aes(x = year, y = cover, colour = factor(Species)), alpha = 0.8, size = 4, show.legend = FALSE) +
    scale_color_manual(values = c("#FF595F"), name = "") +
    scale_fill_manual(values = c("#FF595F"), name = "") +
    scale_x_continuous(breaks = c(1999, 2004, 2009, 2014, 2019)) +
    geom_ribbon(data = HE_plot_preds_df_ev, aes(x = year + 1998, ymin = lower, ymax = upper), fill = "#FF595F", alpha = 0.2) +
    geom_line(data = HE_plot_preds_df_ev, aes(x = year+ 1998, y = mean), colour = "#FF595F") +
    ggtitle("(c) Eriophorum vaginatum") +
    theme_jiri()+
    theme(legend.position = c(0.275, 0.95), 
          axis.line.x = element_line(color = "black", size = 0.5),
          axis.line.y = element_line(color = "black", size = 0.5)) +
    labs(x = "", y = "Species cover (% in plot)\n"))

(herschel_veg_change_plot_al <- ggplot() +
    geom_point(data = subset(QHI_new_cover_HE, Species == "Arctagrostis latifolia"), 
               aes(x = year, y = cover, colour = factor(Species)), alpha = 0.8, size = 4, show.legend = FALSE) +
    scale_color_manual(values = c("#B59E1F"), name = "") +
    scale_fill_manual(values = c("#B59E1F"), name = "") +
    scale_x_continuous(breaks = c(1999, 2004, 2009, 2014, 2019)) +
    geom_ribbon(data = HE_plot_preds_df_al, aes(x = year + 1998, ymin = lower, ymax = upper), fill = "#B59E1F", alpha = 0.2) +
    geom_line(data = HE_plot_preds_df_al, aes(x = year+ 1998, y = mean), colour = "#B59E1F") +
    ggtitle("(d) Arctagrostis latifolia") +
    theme_jiri()+
    theme(legend.position = c(0.275, 0.95), 
          axis.line.x = element_line(color = "black", size = 0.5),
          axis.line.y = element_line(color = "black", size = 0.5)) +
    labs(x = "", y = "Species cover (% in plot)\n"))

(herschel_veg_change_plot_sr <- ggplot() +
    geom_point(data = subset(QHI_new_cover_HE, Species == "Salix reticulata"), 
               aes(x = year, y = cover, colour = factor(Species)), alpha = 0.8, size = 4, show.legend = FALSE) +
    scale_color_manual(values = c("#1FB592"), name = "") +
    scale_fill_manual(values = c("#1FB592"), name = "") +
    scale_x_continuous(breaks = c(1999, 2004, 2009, 2014, 2019)) +
    geom_ribbon(data = HE_plot_preds_df_sr, aes(x = year + 1998, ymin = lower, ymax = upper), fill = "#1FB592", alpha = 0.2) +
    geom_line(data = HE_plot_preds_df_sr, aes(x = year+ 1998, y = mean), colour = "#1FB592") +
    theme_jiri()+
    ggtitle("(e) Salix reticulata") +
    theme(legend.position = c(0.275, 0.95), 
          axis.line.x = element_line(color = "black", size = 0.5),
          axis.line.y = element_line(color = "black", size = 0.5)) +
    labs(x = "", y = "Species cover (% in plot)\n"))

(herschel_veg_change_plot_cc <- ggplot() +
    geom_point(data = subset(QHI_new_cover_HE, Species == "Cetraria cucculata"), 
               aes(x = year, y = cover, colour = factor(Species)), alpha = 0.8, size = 4, show.legend = FALSE) +
    scale_color_manual(values = c("#2BD6D0"), name = "") +
    scale_fill_manual(values = c("#2BD6D0"), name = "") +
    scale_x_continuous(breaks = c(1999, 2004, 2009, 2014, 2019)) +
    geom_ribbon(data = HE_plot_preds_df_cc, aes(x = year + 1998, ymin = lower, ymax = upper), fill = "#2BD6D0", alpha = 0.2) +
    geom_line(data = HE_plot_preds_df_cc, aes(x = year+ 1998, y = mean), colour = "#2BD6D0") +
    theme_jiri()+
    ggtitle("(f) Cetraria cucculata") +
    theme(legend.position = c(0.275, 0.95), 
          axis.line.x = element_line(color = "black", size = 0.5),
          axis.line.y = element_line(color = "black", size = 0.5),
          axis.text=element_text(size=30)) +
    labs(x = "", y = "Species cover (% in plot)\n"))

# Connect them together in a grid
HE_veg_change_grid <- gridExtra::grid.arrange(herschel_veg_change_plot_salpul,
                                              herschel_veg_change_plot_poar,
                                              herschel_veg_change_plot_ev,
                                              herschel_veg_change_plot_al,
                                              herschel_veg_change_plot_sr, 
                                              herschel_veg_change_plot_cc, ncol = 2)

ggsave("figures_outputs/Figure3_HE_cover_significant_species.pdf", 
       HE_veg_change_grid, width = 40, height = 55, units = "cm")

##### Komakuk bayes veg model ####

KO_plot_linear <- MCMCglmm(cover ~ year * Species - 1, 
                               random = ~ year + PLOT, data = QHI_new_cover_KO, 
                               family = "gaussian", pr = TRUE, nitt = 100000, burnin = 10000,thin = 30,
                               prior = prior2)

# Summary and checking trace plots
summary(KO_plot_linear)
plot(KO_plot_linear$Sol)
plot(KO_plot_linear$VCV)
# Significant species are Alopecurus alpinus, Salix arctica and Arctagrostis latifolia

# Komakuk A. alpinus
KO_plots_aa <- subset(QHI_new_cover_KO, Species == "Alopecurus alpinus")
KO_plots_aa <- as.data.frame(KO_plots_aa)

# Creating a column for the trial hits and misses
KO_plots_aa$failures <- round(100 - KO_plots_aa$cover)

KO_plot_aa <- MCMCglmm(cbind(round(cover), failures) ~ I(year - 1998), 
                       random = ~ year + PLOT, data = KO_plots_aa, 
                       family = "multinomial2", pr = TRUE, nitt = 100000, 
                       burnin = 20000, prior = prior2)

# Predictions A.alpinus
nyears <- 21
niter <- length(KO_plot_aa$Sol[,"(Intercept)"])

KO_plot_preds_aloalp <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    KO_plot_preds_aloalp[i,j] <- 100*plogis(KO_plot_aa$Sol[i,"(Intercept)"] + KO_plot_aa$Sol[i,"I(year - 1998)"]*j)
  }
}

KO_plot_preds_df_aloalp <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  KO_plot_preds_df_aloalp[i,] <- quantile(KO_plot_preds_aloalp[,i], c(0.025, 0.5, 0.975))
}

KO_plot_preds_df_aloalp <- cbind.data.frame(lower = KO_plot_preds_df_aloalp[,1], 
                                            mean = KO_plot_preds_df_aloalp[,2], upper = KO_plot_preds_df_aloalp[,3], year = seq(1:21))

# Approximate percentage change = 2.9%
(plogis(((-26.4236) + 1.0925*21))*100) -(plogis(((-26.4236) + 1.0925*1))*100) 

# Komakuk S.arctica
KO_plots_sa <- subset(QHI_new_cover_KO, Species == "Salix arctica")
KO_plots_sa <- as.data.frame(KO_plots_sa)

# Column for the trial hits and misses
KO_plots_sa$failures <- round(100 - KO_plots_sa$cover)

KO_plot_sa <- MCMCglmm(cbind(round(cover), failures) ~ I(year - 1998), 
                       random = ~ year + PLOT, data = KO_plots_sa, 
                       family = "multinomial2", pr = TRUE, nitt = 100000, 
                       burnin = 20000, prior = prior2)

# Predictions Salix arctica
nyears <- 21
niter <- length(KO_plot_sa$Sol[,"(Intercept)"])

KO_plot_preds_salarc <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    KO_plot_preds_salarc[i,j] <- 100*plogis(KO_plot_sa$Sol[i,"(Intercept)"] + KO_plot_sa$Sol[i,"I(year - 1998)"]*j)
  }
}

KO_plot_preds_df_salarc <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  KO_plot_preds_df_salarc[i,] <- quantile(KO_plot_preds_salarc[,i], c(0.025, 0.5, 0.975))
}

KO_plot_preds_df_salarc <- cbind.data.frame(lower = KO_plot_preds_df_salarc[,1], 
                                            mean = KO_plot_preds_df_salarc[,2], upper = KO_plot_preds_df_salarc[,3], year = seq(1:21))

# Approximate percentage change = 7.6%
(plogis(((-1.79163) + 0.02583*21))*100) -(plogis(((-1.79163) + 0.02583*1))*100) 

# Komakuk Arctagrostis latifolia
KO_plots_al <- subset(QHI_new_cover_KO, Species == "Arctagrostis latifolia")
KO_plots_al <- as.data.frame(KO_plots_al)

# Column for the trial hits and misses
KO_plots_al$failures <- round(100 - KO_plots_al$cover)

KO_plot_al <- MCMCglmm(cbind(round(cover), failures) ~ I(year - 1998), 
                       random = ~ year + PLOT, data = KO_plots_al, 
                       family = "multinomial2", pr = TRUE, nitt = 100000, 
                       burnin = 20000, prior = prior2)


# Predictions Arctagrostis latifolia
nyears <- 21
niter <- length(KO_plot_al$Sol[,"(Intercept)"])

KO_plot_preds_arclat <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    KO_plot_preds_arclat[i,j] <- 100*plogis(KO_plot_al$Sol[i,"(Intercept)"] + KO_plot_al$Sol[i,"I(year - 1998)"]*j)
  }
}

KO_plot_preds_df_arclat <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  KO_plot_preds_df_arclat[i,] <- quantile(KO_plot_preds_arclat[,i], c(0.025, 0.5, 0.975))
}

KO_plot_preds_df_arclat <- cbind.data.frame(lower = KO_plot_preds_df_arclat[,1], 
                                            mean = KO_plot_preds_df_arclat[,2], upper = KO_plot_preds_df_arclat[,3], year = seq(1:21))

# Approximate percentage change = 2%
(plogis(((-7.57814) + 0.17825*21))*100) -(plogis(((-7.57814) + 0.17825*1))*100) 


##### Komakuk significant vegtation plots #####

(komakuk_veg_change_plot_al <- ggplot() +
    geom_point(data = subset(QHI_new_cover_KO, Species == "Arctagrostis latifolia"), 
               aes(x = year, y = cover, colour = factor(Species)), alpha = 0.8, size = 4, show.legend = FALSE) +
    scale_color_manual(values = c("#B59E1F"), name = "") +
    scale_fill_manual(values = c("#B59E1F"), name = "") +
    scale_x_continuous(breaks = c(1999, 2004, 2009, 2014, 2019)) +
    geom_ribbon(data = KO_plot_preds_df_arclat, aes(x = year + 1998, ymin = lower, ymax = upper), fill = "#B59E1F", alpha = 0.2) +
    geom_line(data = KO_plot_preds_df_arclat, aes(x = year+ 1998, y = mean), colour = "#B59E1F") +
    theme_jiri()+
    theme(legend.position = c(0.275, 0.95), 
          axis.line.x = element_line(color = "black", size = 0.5),
          axis.line.y = element_line(color = "black", size = 0.5)) +
    labs(x = "", y = "Species cover (% in plot)"))

(komakuk_veg_change_plot_aloalp <- ggplot() +
    geom_point(data = subset(filter(QHI_new_cover_KO, cover != 0), Species == "Alopecurus alpinus"), 
               aes(x = year, y = cover, colour = factor(Species)), alpha = 0.8, size = 4, show.legend = FALSE) +
    scale_color_manual(values = c("#2D95E0FA"), name = "") +
    scale_fill_manual(values = c("#2D95E0FA"), name = "") +
    scale_x_continuous(breaks = c(1999, 2004, 2009, 2014, 2019)) +
    geom_ribbon(data = KO_plot_preds_df_aloalp, aes(x = year + 1998, ymin = lower, ymax = upper), fill = "#2D95E0FA", alpha = 0.2) +
    geom_line(data = KO_plot_preds_df_aloalp, aes(x = year+ 1998, y = mean), colour = "#2D95E0FA") +
    theme_jiri()+
    theme(legend.position = c(0.275, 0.95), 
          axis.line.x = element_line(color = "black", size = 0.5),
          axis.line.y = element_line(color = "black", size = 0.5)) +
    labs(x = "", y = "Species cover (% in plot)"))

(komakuk_veg_change_plot_salarc <- ggplot() +
    geom_point(data = subset(filter(QHI_new_cover_KO, cover != 0), Species == "Salix arctica"), 
               aes(x = year, y = cover, colour = factor(Species)), alpha = 0.8, size = 4, show.legend = FALSE) +
    scale_color_manual(values = c("#71B51F"), name = "") +
    scale_fill_manual(values = c("#71B51F"), name = "") +
    scale_x_continuous(breaks = c(1999, 2004, 2009, 2014, 2019)) +
    geom_ribbon(data = KO_plot_preds_df_salarc, aes(x = year + 1998, ymin = lower, ymax = upper), fill = "#71B51F", alpha = 0.2) +
    geom_line(data = KO_plot_preds_df_salarc, aes(x = year+ 1998, y = mean), colour = "#71B51F") +
    theme_jiri()+
    theme(legend.position = c(0.275, 0.95), 
          axis.line.x = element_line(color = "black", size = 0.5),
          axis.line.y = element_line(color = "black", size = 0.5)) +
    labs(x = "", y = "Species cover (% in plot)"))

# Connect plots together to a grid
KO_veg_change_grid <- gridExtra::grid.arrange(komakuk_veg_change_plot_al + ggtitle("(a) Arcogrostis latifolia"),
                                              komakuk_veg_change_plot_aloalp + ggtitle("(b) Alopecurus alpinus"),
                                              komakuk_veg_change_plot_salarc + ggtitle("(c) Salix arctica"), ncol = 1)


ggsave("figures_outputs/Figure4_KO_cover_significant_species.pdf", 
       KO_veg_change_grid, width = 20, height = 50, units = "cm")

# Model outputs
# Tables for all KO and HE species
HE_plot_all_linear_outputs <- clean.MCMC(HE_plot_linear)
KO_plot_all_linear_outputs <- clean.MCMC(KO_plot_linear)
HE_plot_sp_output <- clean.MCMC(HE_plot_sp)
HE_plot_al_output <- clean.MCMC(HE_plot_al)
HE_plot_pa_output <- clean.MCMC(HE_plot_pa)
HE_plot_ev_output <- clean.MCMC(HE_plot_ev)
HE_plot_sr_output <- clean.MCMC(HE_plot_sr)
HE_plot_cc_output <- clean.MCMC(HE_plot_cc)
KO_plot_aa_output <- clean.MCMC(KO_plot_aa)
KO_plot_sa_output <- clean.MCMC(KO_plot_sa)
KO_plot_al_output <- clean.MCMC(KO_plot_al)

# Herschel species table
stargazer(HE_plot_linear, type = "text", summary = FALSE, digits = 2)

# Komakuk species table
stargazer(KO_plot_linear, type = "html", summary = FALSE, digits = 2)

# Write csv
write.csv(HE_plot_all_linear_outputs, file = "model_outputs/HE_plot_all_linear_outputs.csv")
write.csv(KO_plot_all_linear_outputs, file = "model_outputs/KO_plot_all_linear_outputs.csv")

# Herschel species tables
stargazer(HE_plot_all_linear_outputs, type = "html", summary = FALSE, digits = 2)
stargazer(HE_plot_sp_output, type = "html", summary = FALSE, digits = 2)
stargazer(HE_plot_al_output, type = "html", summary = FALSE, digits = 2)
stargazer(HE_plot_pa_output, type = "html", summary = FALSE, digits = 2)
stargazer(HE_plot_ev_output, type = "html", summary = FALSE, digits = 2)
stargazer(HE_plot_sr_output, type = "html", summary = FALSE, digits = 2)
stargazer(HE_plot_cc_output, type = "html", summary = FALSE, digits = 2)

# Komakuk species tables
stargazer(KO_plot_all_linear_outputs, type = "html", summary = FALSE, digits = 2)
stargazer(KO_plot_aa_output, type = "html", summary = FALSE, digits = 2)
stargazer(KO_plot_sa_output, type = "html", summary = FALSE, digits = 2)
stargazer(KO_plot_al_output, type = "html", summary = FALSE, digits = 2)

##### Species Richness at QHI ####
# Code snippets adapted from Myers-Smith et al., 2019

# Add unique plots
QHI_new_cover$plot_unique <- paste(QHI_new_cover$sub_name,QHI_new_cover$PLOT,QHI_new_cover$year,sep="")

# Species richness model

#Site alpha diversity
alpha_site <- ddply(QHI_new_cover,.(sub_name, year), summarise,
                    richness = length(unique(Species)))
alpha_site$plot_unique <- paste(alpha_site$sub_name, alpha_site$PLOT, sep = "")

alpha_site_HE <- alpha_site %>% 
  filter(sub_name == "QHI:HE")

alpha_site_KO <- alpha_site %>% 
  filter(sub_name == "QHI:KO")

#### Richness HE model ####
richness.HE.model <- MCMCglmm(richness ~ I(year - 1998), data = dplyr::filter(alpha_site, Vegetation == "QHI:HE"), 
                          family = "gaussian", pr = TRUE, nitt = 100000, burnin = 20000)

# Trace plots examination
summary(richness.HE.model)
plot(richness.HE.model$Sol)
plot(richness.HE.model$VCV)

# Model predictions
nyears <- 21
niter <- length(richness.HE.model$Sol[,"(Intercept)"])

richness_HE_preds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    richness_HE_preds[i,j] <- richness.HE.model$Sol[i,"(Intercept)"] + richness.HE.model$Sol[i,"I(year - 1998)"]*j
  }
}

richness_HE_preds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  richness_HE_preds_df[i,] <- quantile(richness_HE_preds[,i], c(0.025, 0.5, 0.975))
}

richness_HE_preds_df <- cbind.data.frame(lower = richness_HE_preds_df[,1], 
                                         mean = richness_HE_preds_df[,2], upper = richness_HE_preds_df[,3], year = seq(1:21))


#### Komakuk richness model ####
richness_KO_m <- MCMCglmm(richness ~ I(year - 1998), data = dplyr::filter(alpha_site, Vegetation == "QHI:KO"), 
                          family = "gaussian", pr = TRUE, nitt = 100000, burnin = 20000)

# Trace plots
summary(richness_KO_m)
plot(richness_KO_m$Sol)
plot(richness_KO_m$VCV)

# Predictions
nyears <- 21
niter <- length(richness_KO_m$Sol[,"(Intercept)"])

richness_KO_preds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    richness_KO_preds[i,j] <- richness_KO_m$Sol[i,"(Intercept)"] + richness_KO_m$Sol[i,"I(year - 1998)"]*j
  }
}

richness_KO_preds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  richness_KO_preds_df[i,] <- quantile(richness_KO_preds[,i], c(0.025, 0.5, 0.975))
}

richness_KO_preds_df <- cbind.data.frame(lower = richness_KO_preds_df[,1], 
                                         mean = richness_KO_preds_df[,2], upper = richness_KO_preds_df[,3], year = seq(1:21))

# Approximate percentage change = 8 species
((17.7462) + 0.4192*21) - ((17.7462) + 0.4192*1)


##### Plot credible intervals richness #####

(richness.plot <- ggplot() +
    geom_point(data = alpha_site, aes(x = year, y = richness, colour = factor(sub_name)), 
               alpha = 0.8, size = 8, position = position_jitter(height = 0.3, width = 0.3)) +
    scale_color_manual(values = c("#ffa544", "#2b299b"), name = "", labels = c("Herschel (ns)", "Komakuk (*)")) +
    scale_fill_manual(values = c("#ffa544","#2b299b", labels = c("Herschel (ns)", "Komakuk (*)"))) +
    scale_x_continuous(breaks = c(1999, 2004, 2009, 2014, 2019)) +
    geom_ribbon(data = richness_HE_preds_df, aes(x = year + 1998, ymin = lower, ymax = upper), 
                fill = "#ffa544", alpha = 0.2) +
    geom_line(data = richness_HE_preds_df, aes(x = year + 1998, y = mean), colour = "#ffa544") +
    geom_ribbon(data = richness_KO_preds_df, aes(x = year + 1998, ymin = lower, ymax = upper), 
                fill = "#2b299b", alpha = 0.2) +
    geom_line(data = richness_KO_preds_df, aes(x = year + 1998, y = mean), colour = "#2b299b") +
    theme_jiri() +
    theme(legend.text=element_text(size=28),
          legend.position = c(0.8, 0.2)) +
    coord_cartesian(ylim = c(15, 30), xlim = c(1999, 2019)) +
    labs(x = "", y = "Species richness\n"))

ggsave("figures_outputs/Figure5_species_richness.pdf", 
       richness.plot, width = 35, height = 25, units = "cm")


##### PART 2: Thermophilisation analysis #####

#### Species list ####
# Species names are already clean from the first part of the code
all_sites <- QHI_new_cover

#### GBIF download ####
# checking how many species from GBIF are actually in GBIF.
# The output should give me at least one species if present in GBIF
species_list_prelim <- unique(sort(all_sites$Species)) # gives me 47 species 
species_gbif_prelim <- occ_data(scientificName = species_list_prelim, hasCoordinate = TRUE, decimalLatitude = '63,90', limit = 1)  

# These are species that GBIF did not recognised
# Alopecurus alpinus 
all_sites$Species <- ifelse(all_sites$Species == "Alopecurus alpinus", "Alopecurus alpinus Sm.", all_sites$Species)
# Polygonum viviparum 
all_sites$Species <- ifelse(all_sites$Species == "Polygonum viviparum", "Polygonum viviparum L.", all_sites$Species)
# Oxytropis maydalliana
all_sites$Species <- ifelse(all_sites$Species == "Oxytropis maydalliana", "Oxytropis maydelliana", all_sites$Species)
# Astralagus umbelletus
all_sites$Species <- ifelse(all_sites$Species == "Astralagus umbelletus", "Astragalus umbellatus", all_sites$Species)
# Pedicularis vertisilata
all_sites$Species <- ifelse(all_sites$Species == "Pedicularis vertisilata", "Pedicularis verticillata", all_sites$Species)
# Cardamine digitalis 
all_sites$Species <- ifelse(all_sites$Species == "Cardamine digitalis", "Cardamine hyperborea", all_sites$Species)
# Lagostis glauca
all_sites$Species <- ifelse(all_sites$Species == "Lagostis glauca", "Lagotis glauca", all_sites$Species)
# Kobresia myotosoides
all_sites$Species <- ifelse(all_sites$Species == "Kobresia myotosoides", "Kobresia myosuroides", all_sites$Species)
# Hierochloe alpine
all_sites$Species <- ifelse(all_sites$Species == "Hierochloe alpine", "Hierochloe alpina", all_sites$Species)
# Pedicularis longsdorfi
all_sites$Species <- ifelse(all_sites$Species == "Pedicularis longsdorfi", "Pedicularis langsdorffii", all_sites$Species)
# Senecio astropurpureus
all_sites$Species <- ifelse(all_sites$Species == "Senecio astropurpureus", "Senecio atropurpureus", all_sites$Species)
# Parrya nudicalis
all_sites$Species <- ifelse(all_sites$Species == "Parrya nudicalis", "Parrya nudicaulis.", all_sites$Species)
# Cladina mitis 
all_sites$Species <- ifelse(all_sites$Species == "Cladina mitis", "Cladina mitis (Sandst.) Hustich", all_sites$Species)
# Cetraria cucculata
all_sites$Species <- ifelse(all_sites$Species == "Cetraria cucculata", "Cetraria cucullata (Bellardi) Ach.", all_sites$Species)
# Lupinus arcticus
all_sites$Species <- ifelse(all_sites$Species == "Lupinus arcticus ", "Lupinus arcticus", all_sites$Species)
# Lupinus arcticus
all_sites$Species <- ifelse(all_sites$Species == "Luzula arctica ", "Luzula arctica", all_sites$Species)

# Write changed species to be downloaded from GBIF
species_list_gbif <- as.data.frame(c(unique(all_sites$Species)), c(length(unique(all_sites$Species))))
write.csv(species_list_gbif, "outputs/species_lists/gbif_species_list.csv")

# GBIF occurrences download official
species_list <- c(unique(all_sites$Species)) 
species_gbif <- occ_data(scientificName = species_list, hasCoordinate = TRUE, decimalLatitude = '63,90', limit = 100000)  
# All species have at least some returned occurrences!

##### Species coordinates ####
# Code taken from https://data-blog.gbif.org/post/downloading-long-species-lists-on-gbif/
coords_gbif <- vector("list", length(species_list))
names(coords_gbif) <- species_list
for (s in species_list) {
  coords <- species_gbif[[s]]$data[ , c("decimalLongitude", "decimalLatitude",
                                        "occurrenceStatus", "coordinateUncertaintyInMeters",
                                        "institutionCode", "year", "gbifID")]
  coords_gbif[[s]] <- data.frame(species = s, coords)
}
lapply(coords_gbif, head)
# collapse the list into a data frame:
coords_gbif_df <- as.data.frame(do.call(rbind, coords_gbif), row.names = 1:sum(sapply(coords_gbif, nrow)))
head(coords_gbif_df)
tail(coords_gbif_df)

# Save GBIF file
write.csv(coords_gbif_df, "data/coords_gbif_final.csv")

#### GBIF coordinates cleaning ####
#  Code from https://ropensci.github.io/CoordinateCleaner/articles/Cleaning_GBIF_data_with_CoordinateCleaner.html
# and from https://ourcodingclub.github.io/tutorials

# Flag problematic 
load("data/buffland_1deg.rda")
flags <- clean_coordinates(x = coords_gbif_df,
                           lon = "decimalLongitude",
                           lat = "decimalLatitude",
                           species = "species",
                           tests = c("capitals", "centroids", "equal", "gbif", "institutions",
                                     "zeros", "seas"),
                           outliers_method = "distance", outliers_td = 5000,
                           seas_ref = buffland) # these for now

# Exclude problematic records
coords_gbif_clean <- coords_gbif_df[flags$.summary,]

# Age of records
table(coords_gbif_clean$year)
coords_gbif_clean <- coords_gbif_clean%>%
  filter(year > 1850) # remove records from before 1850

# Save final file
write.csv(coords_gbif_clean, "data/coords_gbif_clean.csv")

#### CHELSA data ####
# Code adapted from Joseph Everest
# CHELSA climatologies (.tif) are stored on the Team Shrub Drive

# Defining the filepath to the files
folderpath.chelsa.HD <- ("D:/CHELSA") # Hard drive
filenames.chelsa <- list.files(folderpath.chelsa.HD, pattern = "*.tif")
filepath.chelsa = paste0(folderpath.chelsa.HD, "/", filenames.chelsa)

# Read in the raster file
chelsa_raster_10_10 <- raster("data/CHELSA/CHELSA_bio10_10.tif") # Mean warmest quarter

# Create raster stack
chelsa.stack_10_10 <- stack(chelsa_raster_10_10)

# Splitting my data to coords and species
species_only <- c(coords_gbif_clean$species)
coords_only <- coords_gbif_clean %>%
  dplyr::select(decimalLongitude, decimalLatitude)

# Extract variables values for each pair of coordinates
points_species <- SpatialPoints(coords_only)
chelsa.extract_10_10 <- raster::extract(chelsa.stack_10_10, points_species, df = TRUE)

# Combine dataframes
# Convert the SpatialPoints (sp) object into a dataframe 
coord.df <- as.data.frame(points_species)

# Reassign the 'ID' to the ITEX coordinates dataframe
coord.df$ID <- row.names(coord.df)
coord.df$ID <- as.numeric(coord.df$ID) # Make numeric

# Merge the two dataframes: extracted CHELSA variables and the ITEX coordinates
coord.chelsa.combo.final <- left_join(chelsa.extract_10_10, coord.df, by = c("ID" = "ID"))

coord.chelsa.combo.final <- coord.chelsa.combo.final %>% 
  mutate(CHELSA11_temp_mean_warmest_q = CHELSA_bio10_10/10) %>% # Divide by 10 to get to degC
  dplyr::select(-CHELSA_bio10_10) %>%
  mutate(species_name = species_only) %>%
  na.omit(CHELSA11_temp_mean_warmest_q)

# Climatic mean calcualtion
climatic_mean <- coord.chelsa.combo.final %>%
  dplyr::group_by(species_name) %>%
  dplyr::summarise(temp_mean_warmest_q = mean(CHELSA11_temp_mean_warmest_q))

# Save
write.csv(climatic_mean, "outputs/thermal_niches/thermal_niches_all_species.csv")

# Changing species back to their original name from ITEX for consistency
# Alopecurus alpinus 
climatic_mean$species_name <- ifelse(climatic_mean$species_name == "Alopecurus alpinus Sm.", "Alopecurus alpinus", climatic_mean$species_name)
# Polygonum viviparum 
climatic_mean$species_name <- ifelse(climatic_mean$species_name ==  "Polygonum viviparum L.", "Polygonum viviparum", climatic_mean$species_name)
# Oxytropis maydalliana
climatic_mean$species_name <- ifelse(climatic_mean$species_name == "Oxytropis maydelliana", "Oxytropis maydalliana", climatic_mean$species_name)
# Astralagus umbelletus
climatic_mean$species_name <- ifelse(climatic_mean$species_name == "Astragalus umbellatus", "Astralagus umbelletus", climatic_mean$species_name)
# Pedicularis vertisilata
climatic_mean$species_name <- ifelse(climatic_mean$species_name == "Pedicularis verticillata", "Pedicularis vertisilata", climatic_mean$species_name)
# Cardamine digitalis 
climatic_mean$species_name <- ifelse(climatic_mean$species_name == "Cardamine hyperborea", "Cardamine digitalis", climatic_mean$species_name)
# Lagostis glauca
climatic_mean$species_name <- ifelse(climatic_mean$species_name == "Lagotis glauca", "Lagostis glauca", climatic_mean$species_name)
# Kobresia myotosoides
climatic_mean$species_name <- ifelse(climatic_mean$species_name == "Kobresia myosuroides", "Kobresia myotosoides", climatic_mean$species_name)
# Hierochloe alpine
climatic_mean$species_name <- ifelse(climatic_mean$species_name == "Hierochloe alpina", "Hierochloe alpine", climatic_mean$species_name)
# Pedicularis longsdorfi
climatic_mean$species_name <- ifelse(climatic_mean$species_name == "Pedicularis langsdorffii", "Pedicularis longsdorfi", climatic_mean$species_name)
# Senecio astropurpureus
climatic_mean$species_name <- ifelse(climatic_mean$species_name == "Senecio atropurpureus", "Senecio astropurpureus", climatic_mean$species_name)
# Parrya nudicalis
climatic_mean$species_name <- ifelse(climatic_mean$species_name =="Parrya nudicaulis.", "Parrya nudicalis", climatic_mean$species_name)
# Cladina mitis 
climatic_mean$species_name <- ifelse(climatic_mean$species_name == "Cladina mitis (Sandst.) Hustich", "Cladina mitis", climatic_mean$species_name)
# Cetraria cucculata
climatic_mean$species_name <- ifelse(climatic_mean$species_name == "Cetraria cucullata (Bellardi) Ach.", "Cetraria cucculata", climatic_mean$species_name)

#### Calculating CTI for every location ####
# Herschel vegetation 
# using the "QHI_new_cover" data, which is a dataframe form the beginning,
# it does not have zeros filled in for absent species, but it is cleaned for names

QHI_HE <- QHI_new_cover %>% 
  filter(sub_name == "QHI:HE") %>%
  dplyr::select(year, PLOT, Species, cover, sub_name, site_name) 

HE_proportional_warm_niche <- convert.proportional.warm.niche.cti(QHI_HE, climatic_mean)
HE_proportional_warm_niche$PLOT <- as.factor(HE_proportional_warm_niche$PLOT)

##### HE Thermophilisation plots and model ####
# Again, MCMC code for credible intervals calculation adapted from Myers-Smith et al., 2018

# Herschel CTI
HE_cti_linear <- MCMCglmm(CTI ~ I(year-1998), 
                          random = ~ year + Plot, data = HE_proportional_warm_niche, 
                          family = "gaussian", pr = TRUE, nitt = 100000, burnin = 10000,thin = 30)

# Check trace plots
summary(HE_cti_linear)
plot(HE_cti_linear$VCV)
plot(HE_cti_linear$Sol)

# Predictions
nyears <- 21
niter <- length(HE_cti_linear$Sol[,"(Intercept)"])

HE_cti_preds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    HE_cti_preds[i,j] <- HE_cti_linear$Sol[i,"(Intercept)"] + HE_cti_linear$Sol[i,"I(year - 1998)"]*j
  }
}

HE_cti_preds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  HE_cti_preds_df[i,] <- quantile(HE_cti_preds[,i], c(0.025, 0.5, 0.975))
}

HE_cti_preds_df <- cbind.data.frame(lower = HE_cti_preds_df[,1], 
                                    mean = HE_cti_preds_df[,2], upper = HE_cti_preds_df[,3], year = seq(1:21))

# Approximate temperature change = 8.1˚C
((12.5151) + 0.4084*21) - ((12.5151) + 0.4084*1)

# Preds model CTI Herschel
(cti.he.preds.plot <- ggplot() +
    geom_point(data = HE_proportional_warm_niche, aes(x = year, y = CTI, colour = Plot), 
               alpha = 0.8, size = 7, position = position_jitter(height = 0.3, width = 0.3)) +
    scale_x_continuous(breaks = c(1999, 2004, 2009, 2014, 2019)) +
    scale_colour_viridis_d(option = "plasma") +
    geom_ribbon(data = HE_cti_preds_df, aes(x = year + 1998, ymin = lower, ymax = upper), 
                fill = "#ffa544", alpha = 0.2) +
    geom_line(data = HE_cti_preds_df, aes(x = year + 1998, y = mean), colour = "#ffa544") +
    theme_jiri() +
    theme(legend.text=element_text(size=20),
          legend.title = element_text(size = 20)) +
    xlab("\nYear\n") + 
    ylab(bquote(''*~"Community Temperature Index" ~ "°C")))

# Thermophilization rate of each plot (TRplot)
# The rate of change in CTI between the initial and final census thermal migration rate TMR 
HE_mean_thermoph_per_plot <- slice.min.max(HE_proportional_warm_niche) %>% 
  arrange(year)

HE_mean_thermophil_per_plot <- HE_mean_thermoph_per_plot %>% 
  mutate(minmaxyear = case_when(year < 2000 ~ "Minimum", 
                                year > 2010 ~ "Maximum")) %>% 
  group_by(site_name, sub_name, Plot) %>% 
  summarise(meanCTI = CTI[minmaxyear=='Maximum'] -
              CTI[minmaxyear=='Minimum'])  

HE_mean_thermophil_per_plot <- HE_mean_thermophil_per_plot %>% 
  mutate(Plot = 1:nrow(HE_mean_thermophil_per_plot)) %>% 
  mutate(sub_name = "HE") %>% 
  mutate(site_name = "QHI")

# Plot the thermophilisation migration rate for Herschel
(HE_mean_thermophilization_per_plot <- ggplot(HE_mean_thermoph_per_plot, aes(x = Plot, y = CTI)) +
    geom_point(aes(y = CTI, colour = year > 2015), size = 5, alpha = 0.5, show.legend = FALSE) +
    geom_line(aes(group = Plot), arrow = arrow(), size = 1.5, alpha = 0.9, color = "#ffa544") +
    scale_colour_viridis_d(option = "plasma") +
    xlab("\nPlot") +
    ylab(bquote(''*~"Community Temperature Index" ~ "°C (1999-2019)")) +
    theme_jiri())

# Komakuk CTI
QHI_KO <- QHI_new_cover %>% 
  filter(sub_name == "QHI:KO") %>%
  dplyr::select(year, PLOT, Species, cover, sub_name, site_name) 

KO_proportional_warm_niche <- convert.proportional.warm.niche.cti(QHI_KO, climatic_mean)
KO_proportional_warm_niche$PLOT <- as.factor(KO_proportional_warm_niche$PLOT)

# Thermophilisation plots and model
KO.cti.model.bayes <- MCMCglmm(CTI ~ I(year-1998), 
                               random = ~ year + Plot, data = KO_proportional_warm_niche, 
                               family = "gaussian", pr = TRUE, nitt = 100000, burnin = 10000,thin = 30)
# Check trace plots
summary(KO.cti.model.bayes)
plot(KO.cti.model.bayes$VCV)
plot(KO.cti.model.bayes$Sol)

# Predictions
nyears <- 21
niter <- length(KO.cti.model.bayes$Sol[,"(Intercept)"])

KO_cti_preds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    KO_cti_preds[i,j] <- KO.cti.model.bayes$Sol[i,"(Intercept)"] + KO.cti.model.bayes$Sol[i,"I(year - 1998)"]*j
  }
}

KO_cti_preds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  KO_cti_preds_df[i,] <- quantile(KO_cti_preds[,i], c(0.025, 0.5, 0.975))
}

KO_cti_preds_df <- cbind.data.frame(lower = KO_cti_preds_df[,1], 
                                    mean = KO_cti_preds_df[,2], upper = KO_cti_preds_df[,3], year = seq(1:21))

# Approximate temperature change = 5.328 ˚C
((7.8399) + 0.2664*21) - ((7.8399) + 0.2664*1)

# Preds model CTI KO plot
(cti.ko.preds.plot <- ggplot() +
    geom_point(data = KO_proportional_warm_niche, aes(x = year, y = CTI, colour = Plot), 
               alpha = 0.8, size = 7, position = position_jitter(height = 0.3, width = 0.3)) +
    scale_x_continuous(breaks = c(1999, 2004, 2009, 2014, 2019)) +
    scale_colour_viridis_d(option = "viridis") +
    geom_ribbon(data = KO_cti_preds_df, aes(x = year + 1998, ymin = lower, ymax = upper), 
                fill = "#2b299b", alpha = 0.2) +
    geom_line(data = KO_cti_preds_df, aes(x = year + 1998, y = mean), colour = "#2b299b") +
    #geom_ribbon(data = richness_KO_preds_df, aes(x = year + 1998, ymin = lower, ymax = upper), 
    #            fill = "#2b299b", alpha = 0.2) +
    #geom_line(data = richness_KO_preds_df, aes(x = year + 1998, y = mean), colour = "#2b299b") +
    theme_jiri() +
    theme(axis.title.y = element_blank(),
          legend.text=element_text(size=20),
          legend.title = element_text(size = 20)) +
    #coord_cartesian(ylim = c(15, 35), xlim = c(1999, 2019)) +
    xlab("\nYear\n") + 
    ylab(bquote(''*~"Community Temperature Index" ~ "°C\n")))

# Panel thermophilization CTI both HE and KO
panel.cti <- gridExtra::grid.arrange(cti.he.preds.plot + ggtitle("(a) Herschel Vegetation (***)"),
                                     cti.ko.preds.plot + ggtitle("(b) Komakuk Vegetation (***)"), ncol=2)


ggsave("figures_outputs/Figure6_CTI.pdf", 
       panel.cti, width = 50, height = 30, units = "cm")

# Komakuk Thermophilization rate of each plot (TRplot)
KO_mean_thermoph_per_plot <- slice.min.max(KO_proportional_warm_niche) %>% 
  arrange(year)

KO_mean_thermophil_per_plot <- KO_mean_thermoph_per_plot %>% 
  mutate(minmaxyear = case_when(year < 2000 ~ "Minimum", 
                                year > 2010 ~ "Maximum")) %>% 
  group_by(site_name,sub_name,Plot) %>% 
  summarise(meanCTI = CTI[minmaxyear=='Maximum'] -
              CTI[minmaxyear=='Minimum'])  

KO_mean_thermophil_per_plot <- KO_mean_thermophil_per_plot %>% 
  mutate(Plot = 1:nrow(KO_mean_thermophil_per_plot)) %>% 
  mutate(sub_name = "KO") %>% 
  mutate(site_name = "QHI")

# Plot TMR for Komakuk
(KO_mean_thermophilization_per_plot <- ggplot(KO_mean_thermoph_per_plot, aes(x = Plot, y = CTI)) +
    geom_point(aes(y = CTI, colour = year > 2015), size = 5, alpha = 0.5, show.legend = FALSE) +
    geom_line(aes(group = Plot), arrow = arrow(), size = 1.5, alpha = 0.9, color = "#2b299b") +
    scale_colour_viridis_d(option = "viridis") +
    xlab("\nPlot") +
    #ylab(bquote(''*~"Community Temperature Index" ~ "°C (1999-2019)")) +
    theme_jiri() +
    theme(axis.title.y = element_blank()))

(QHI_mean_thermophil_plot_HE <- ggplot(QHI_thermophilization, aes(x = Plot, y = CTI)) +
    geom_point(aes(y = CTI, colour = year > 2015), size = 10, alpha = 0.5, show.legend = FALSE) +
    geom_line(aes(group = Plot), arrow = arrow(), size = 3, alpha = 0.9, color = "#ffa544") +
    scale_colour_viridis_d(option = "viridis") +
    facet_wrap(~ sub_name) +
    geom_point(aes(y = CTI)) +
    ylab(bquote(''*~"Community Temperature Index" ~ "°C")) +
    theme_jiri() +
    theme(strip.text = element_text(size = 35)))

(QHI_mean_thermophil_plot_KO <- ggplot(QHI_thermophilization, aes(x = Plot, y = CTI)) +
    geom_point(aes(y = CTI, colour = year > 2015), size =10, alpha = 0.5, show.legend = FALSE) +
    geom_line(aes(group = Plot), arrow = arrow(), size = 3, alpha = 0.9, color = "#2b299b") +
    scale_colour_viridis_d(option = "viridis") +
    facet_wrap(~ sub_name) +
    geom_point(aes(y = CTI)) +
    ylab(bquote(''*~"Community Temperature Index" ~ "°C")) +
    theme_jiri() +
    theme(strip.text = element_text(size = 35)))

ggsave("figures_outputs/Figure6_TMR_HE.pdf", 
       QHI_mean_thermophil_plot_HE, width = 50, height = 30, units = "cm")

ggsave("figures_outputs/Figure6_TMR_KO.pdf", 
       QHI_mean_thermophil_plot_KO, width = 50, height = 30, units = "cm")


## QHI thermophilization per plot
thermophil_panel <- gridExtra::grid.arrange(HE_mean_thermophilization_per_plot + ggtitle("(a) Herschel Vegetation"),
                                            KO_mean_thermophilization_per_plot + ggtitle("(b) Komakuk Vegetation"), ncol=2)

ggsave("figures_outputs/Figure7_TMR.pdf", 
       thermophil_panel, width = 50, height = 30, units = "cm")

### Plotting for both HE and KO together
QHI_mean_thermophil <- bind_rows(HE_mean_thermophil_per_plot, KO_mean_thermophil_per_plot)
QHI_mean_thermophil$sub_name <- as.factor(QHI_mean_thermophil$sub_name)
QHI_mean_thermophil$sub_name[QHI_mean_thermophil$sub_name == "QHI:KO"] <- "Komakuk"
QHI_mean_thermophil$sub_name[QHI_mean_thermophil$sub_name == "QHI:HE"] <- "Herschel"

# Mean change in decade
(QHI_mean_thermophilization_per_plot <- ggplot(QHI_mean_thermophil, aes(x = sub_name, y = (CTI/20), fill= CTI)) +
    geom_boxplot(aes(fill = sub_name, alpha = 0.7), show.legend = FALSE) +
    #scale_color_manual(values = c("#ffa544", "#2b299b")) +
    scale_fill_manual(values = c("#ffa544", "#2b299b")) +
    xlab("\nVegetation Type") +
    ylab(bquote(''*~"Mean change of CTI" ~ "°C (1999-2019)")) +
    theme_jiri() +
    theme(axis.title.x = element_blank()))

ggsave("figures_outputs/Figure8_mean_change_CTI.pdf", 
       QHI_mean_thermophilization_per_plot, width = 50, height = 30, units = "cm")

##### RQ4 #####
# Thermal Niche analysis QHI

# I have effect sizes from vegetation change models (winner and loser info)
slopes_of_change_qhi # Excel file from "data/qhi/slopes_of_change.xlsx" - Choose third tab
slopes_qhi <- slopes_of_change_qhi 

species_strip <- str_sub(slopes_qhi$variable, start = 13)
slopes_qhi <- slopes_qhi %>% 
  mutate(species_strip) %>% 
  dplyr::select(-variable) 

slopes_HE <- slopes_qhi %>% 
  dplyr::filter(subsite == "QHI:HE") %>% 
  dplyr::select(estimate, species_strip) %>% 
  dplyr::rename(species_name = species_strip)

slopes_KO <- slopes_qhi %>% 
  dplyr::filter(subsite == "QHI:KO") %>% 
  dplyr::select(estimate, species_strip)  %>% 
  dplyr::rename(species_name = species_strip)

#Thermal niche calculation for all species
climatic_mean
# Species that are only in HE
QHI_only_HE <- QHI_new_cover %>% 
  filter(sub_name == "QHI:HE")
species_unique_HE <- c(sort(unique(QHI_only_HE$Species)))

climatic_mean_HE <- climatic_mean %>% 
  filter(species_name %in% species_unique_HE) 

# Species that are only in KO
QHI_only_KO <- QHI_new_cover %>% 
  filter(sub_name == "QHI:KO")
species_unique_KO <- c(unique(QHI_only_KO$Species))

climatic_mean_KO <- climatic_mean %>% 
  filter(species_name %in% species_unique_KO) 

# bind them
thermal_niche_KO <- climatic_mean_KO %>% 
  mutate(slope_of_change = slopes_KO$estimate)  %>% 
  dplyr::mutate(Significance = "Non-significant")

thermal_niche_KO$Significance[thermal_niche_KO$species_name == "Alopecurus alpinus"] <- "Significant"
thermal_niche_KO$Significance[thermal_niche_KO$species_name == "Arctagrostis latifolia"] <- "Significant"
thermal_niche_KO$Significance[thermal_niche_KO$species_name == "Salix arctica"] <- "Significant"
thermal_niche_KO$Significance <- as.factor(thermal_niche_KO$Significance)
thermal_niche_KO$species_name <- as.factor(thermal_niche_KO$species_name)

thermal_niche_HE <- climatic_mean_HE %>% 
  dplyr::mutate(slope_of_change = slopes_HE$estimate) %>% 
  dplyr::mutate(Significance = "Non-significant")
thermal_niche_HE$Significance[thermal_niche_HE$species_name == "Salix pulchra"] <- "Significant"
thermal_niche_HE$Significance[thermal_niche_HE$species_name == "Poa arctica"] <- "Significant"
thermal_niche_HE$Significance[thermal_niche_HE$species_name == "Eriophorum vaginatum"] <- "Significant"
thermal_niche_HE$Significance[thermal_niche_HE$species_name == "Arctagrostis latifolia"] <- "Significant"
thermal_niche_HE$Significance[thermal_niche_HE$species_name == "Salix reticulata"] <- "Significant"
thermal_niche_HE$Significance[thermal_niche_HE$species_name == "Cetraria cucculata"] <- "Significant"
thermal_niche_HE$Significance <- as.factor(thermal_niche_HE$Significance)
thermal_niche_HE$species_name <- as.factor(thermal_niche_HE$species_name)

# Save thermal niches for every veg type
thermal_niche_only_HE <- thermal_niche_HE %>% 
  select(species_name, temp_mean_warmest_q)
write.csv(thermal_niche_only_HE, "outputs/thermal_niches/thhermal_niches_HE.csv")

thermal_niche_only_KO <- thermal_niche_KO %>% 
  select(species_name, temp_mean_warmest_q)
write.csv(thermal_niche_only_KO, "outputs/thermal_niches/thermal_niches_KO.csv")


#lm(thermal niche ~ slopes of change over time) - 1 model per subsite
# incorporate whether the change in abundance over time was statistically significant or not - but I might just include that visually like with the alpha of the dots in the final plot. 

# HE
HE.thermal.niche.bayes <- MCMCglmm(temp_mean_warmest_q ~ slope_of_change, 
                                   data = thermal_niche_HE, 
                                   family = "gaussian", pr = TRUE, nitt = 100000, burnin = 10000,thin = 30)
summary(HE.thermal.niche.bayes)
plot(HE.thermal.niche.bayes$VCV)
plot(HE.thermal.niche.bayes$Sol)

# Approximate temperature change 
((9.24183) + 0.10328*48) - ((9.24183) + 0.10328*1)

# Species change significance:
# *** Salix pulchra, Poa arctica, Eriophorum vaginatum, Arctagrostis latifolia   
# * Salix reticulata, Cetratia cucculata


(HE.thermal.niche.plot <- ggplot(aes(temp_mean_warmest_q, slope_of_change),
                                 data = thermal_niche_HE) +
    geom_point(aes(color = Significance), size = 5, alpha = 0.7, show.legend = FALSE) +
    geom_smooth(method = "lm", color = "#ffa544", , fill = "#ffa544") +
    scale_x_continuous(breaks = c(6, 10, 14)) +
    scale_colour_viridis_d(option = "viridis") +
    ylab("Change in % cover over 20 years") +
    xlab(bquote(''*~"Mean warm thermal niche" ~ "(°C)")) +
    geom_text(aes(label=ifelse(Significance == "Significant",as.character(species_name),'')),size = 8.2, hjust=0.8,vjust=1.7, show.legend = FALSE) +
    theme_jiri())


# KO
# Species change significance:
# *** Alopecurus alpinus, Arctagrostis latifolia    
# * Salix arctica

# Bayesian KO
KO.thermal.niche.bayes <- MCMCglmm(temp_mean_warmest_q ~ slope_of_change, 
                                   data = thermal_niche_KO, 
                                   family = "gaussian", pr = TRUE, nitt = 100000, burnin = 10000,thin = 30)
summary(KO.thermal.niche.bayes)
plot(KO.thermal.niche.bayes$VCV)
plot(KO.thermal.niche.bayes$Sol)

plasma_pal <- c("green", viridis::viridis(n = 4))

(KO.thermal.niche.plot <- ggplot(aes(temp_mean_warmest_q, slope_of_change),
                                 data = thermal_niche_KO) +
    geom_point(aes(color = Significance), size = 5, alpha = 0.7) +
    geom_smooth(method = "lm", color = "#2b299b", fill = "#2b299b") +
    scale_x_continuous(breaks = c(6, 10, 14)) +
    #scale_color_manual(values = plasma_pal) +
    scale_colour_viridis_d(option = "viridis") +
    ylab("") +
    xlab(bquote(''*~"Mean warm thermal niche" ~ "(°C)")) +
    geom_text(aes(label=ifelse(Significance == "Significant",as.character(species_name),'')),size = 8.2, hjust=-0.1,vjust=2) +
    theme_jiri() +
    theme(legend.position = c(0.8,0.8), legend.box.background = element_rect(colour = "black"),
          legend.text=element_text(size=28),
          legend.title = element_text(size = 28)))


panel_winnerLoser <- gridExtra::grid.arrange(HE.thermal.niche.plot + ggtitle("(a) Herschel Vegetation (*)"),
                                             KO.thermal.niche.plot + ggtitle("(b) Komakuk Vegetation (ns)"), ncol=2)


ggsave("outputs/figures_outputs/Figure9_winner_loser.pdf", 
       panel_winnerLoser, width = 50, height = 30, units = "cm")


