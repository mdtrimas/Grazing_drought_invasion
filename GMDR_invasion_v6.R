########## How does multi-year, multi-intensity drought and grazing impact invasive annual bromes? ########
#Frost - dissertation chapter

#for PC
setwd("C:/Users/mdtrimas/Box Sync/R work/DxG_spcomp")
search()

#for Mac
setwd("~/Box Sync/R work/DxG_spcomp")
search()

#loading packages
install.packages("vegan")
install.packages("tidyverse")
install.packages("codyn")
install.packages("ggplot2")
install.packages("reshape2")
install.packages("plotly")
install.packages("nmle")
install.packages("lme4")
install.packages("olsrr")
install.packages("car")
install.packages("patchwork")
install.packages("lmerTest")
install.packages("piecewiseSEM")
install.packages("multcomp")
install.packages("MuMIn")

library(vegan)
library(tidyverse)
library(codyn)
library(ggplot2)
library(reshape2)
library(plotly)
library(nlme)
library(lme4)
library(olsrr)
library(car)
library(patchwork)
library(lmerTest)
library(piecewiseSEM)
library(multcomp)
library(MuMIn)


##################################################
##### Read in Data ##############

#cleaned up data is available in GitHub - will need to read in as annotated throughout

#Set ggplot2 theme to black and white
theme_set(theme_bw())
#Update ggplot2 theme - make box around the x-axis title size 30, vertically justify x-axis title to 0.35, 
#Place a margin of 15 around the x-axis title.  
#Make the x-axis title size 30. For y-axis title, make the box size 30, put the writing at a 90 degree angle, and vertically justify the title to 0.5.  
#Add a margin of 15 and make the y-axis text size 25. Make the plot title size 30 and vertically justify it to 2.  Do not add any grid lines.  
#Do not add a legend title, and make the legend size 20
theme_update(axis.title.x=element_text(size=20, vjust=-0.35, margin=margin(t=12)),
             axis.text.x=element_text(size=20), axis.title.y=element_text(size=20, angle=90, vjust=0.5,
                                                                          margin=margin(r=15)), axis.text.y=element_text(size=20), plot.title =
               element_text(size=20, vjust=2), panel.grid.major=element_blank(),
             panel.grid.minor=element_blank(),
             legend.text=element_text(size=15))

#set colorblind friendly color palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


#########################################################################
########### Creating functions #########################

#repeated measures anova with 4 independent variables
anova_t3 <- function(IndVars=IndVars, DepVar=DepVar, RndForm=RndForm, Data=Data){
  anova_out <- {}
  IndVarMatrix <- matrix(nrow=length(IndVars),ncol=length(IndVars))
  IndVars2x <- c(IndVars,IndVars)
  
  for(REORDER in 1:length(IndVars)){
    IndVarMatrix[REORDER,] <- IndVars2x[REORDER:(length(IndVars)+(REORDER-1))]
  }
  rm(IndVars2x)
  
  for(RUN in 1:length(IndVars)){
    model_formula_temp <- paste0(DepVar,"~", paste0(IndVarMatrix[RUN,], collapse="*"))
    model_temp <- lme(as.formula(model_formula_temp)
                      , data=Data
                      , random = as.formula(RndForm)
                      , correlation=corCompSymm(form = as.formula(RndForm))
                      , control=lmeControl(returnObject=TRUE)
                      , na.action = na.omit)
    anova_out_temp <- anova(model_temp)
    
    if(length(IndVars)==4){ ## Pulls model output for variable that is last
      if(RUN==1){ ### This currently works for 
        anova_partial_temp <- rbind(anova_out_temp[(length(IndVars)+1),],
                                    anova_out_temp[7:8,]
        )
      }
      if(RUN %in% 2:length(IndVars)){
        anova_partial_temp <- rbind(anova_out_temp[(length(IndVars)+1),],
                                    anova_out_temp[7,]
        )
      }
    } # End if vars==4 statement
    
    if(length(IndVars)==3){ ## Pulls model output for variable that is last
      if(RUN==1){ ### This currently works for 
        anova_partial_temp <- rbind(anova_out_temp[(length(IndVars)+1),],
                                    anova_out_temp[7:8,]
        )
      }
      if(RUN %in% 2:length(IndVars)){
        anova_partial_temp <- rbind(anova_out_temp[(length(IndVars)+1),],
                                    anova_out_temp[7,]
        )
      }
    } # End if vars==3 statement
    
    if(length(IndVars)==2){ ## Pulls model output for variable that is last
      if(RUN==1){ ### This currently works for 
        anova_partial_temp <- rbind(anova_out_temp[(length(IndVars)+1),],
                                    anova_out_temp[4,]
        )
      }
      if(RUN %in% 2:length(IndVars)){
        anova_partial_temp <- anova_out_temp[(length(IndVars)+1),]
        
      }
    } # End if vars == 2 statement
    anova_out <- rbind(anova_out, anova_partial_temp)
  } # End reorder loop
  return(anova_out)
} # End function


#repeated measures anova with 3 independent variables
anova_t3_3ind <- function(IndVars=IndVars, DepVar=DepVar, RndForm=RndForm, Data=Data){
  anova_out <- {}
  IndVarMatrix <- matrix(nrow=length(IndVars),ncol=length(IndVars))
  IndVars2x <- c(IndVars,IndVars)
  
  for(REORDER in 1:length(IndVars)){
    IndVarMatrix[REORDER,] <- IndVars2x[REORDER:(length(IndVars)+(REORDER-1))]
  }
  rm(IndVars2x)
  
  for(RUN in 1:length(IndVars)){
    model_formula_temp <- paste0(DepVar,"~", paste0(IndVarMatrix[RUN,], collapse="*"))
    model_temp <- lme(as.formula(model_formula_temp)
                      , data=Data
                      , random = as.formula(RndForm)
                      , correlation=corCompSymm(form = as.formula(RndForm))
                      , control=lmeControl(returnObject=TRUE)
                      , na.action = na.omit)
    anova_out_temp <- anova(model_temp)
    
    if(length(IndVars)==3){ ## Pulls model output for variable that is last
      if(RUN==1){ ### This currently works for 
        anova_partial_temp <- rbind(anova_out_temp[(length(IndVars)+1),],
                                    anova_out_temp[7:8,]
        )
      }
      if(RUN %in% 2:length(IndVars)){
        anova_partial_temp <- rbind(anova_out_temp[(length(IndVars)+1),],
                                    anova_out_temp[7,]
        )
      }
    } # End if vars==3 statement
    
    if(length(IndVars)==2){ ## Pulls model output for variable that is last
      if(RUN==1){ ### This currently works for 
        anova_partial_temp <- rbind(anova_out_temp[(length(IndVars)+1),],
                                    anova_out_temp[4,]
        )
      }
      if(RUN %in% 2:length(IndVars)){
        anova_partial_temp <- anova_out_temp[(length(IndVars)+1),]
        
      }
    } # End if vars == 2 statement
    anova_out <- rbind(anova_out, anova_partial_temp)
  } # End reorder loop
  return(anova_out)
} # End function


#repeated measures anova with 2 independent variables
anova_t3_2ind <- function(IndVars=IndVars, DepVar=DepVar, RndForm=RndForm, Data=Data){
  anova_out <- {}
  IndVarMatrix <- matrix(nrow=length(IndVars),ncol=length(IndVars))
  IndVars2x <- c(IndVars,IndVars)
  
  for(REORDER in 1:length(IndVars)){
    IndVarMatrix[REORDER,] <- IndVars2x[REORDER:(length(IndVars)+(REORDER-1))]
  }
  rm(IndVars2x)
  
  for(RUN in 1:length(IndVars)){
    model_formula_temp <- paste0(DepVar,"~", paste0(IndVarMatrix[RUN,], collapse="*"))
    model_temp <- lme(as.formula(model_formula_temp)
                      , data=Data
                      , random = as.formula(RndForm)
                      , correlation=corCompSymm(form = as.formula(RndForm))
                      , control=lmeControl(returnObject=TRUE)
                      , na.action = na.omit)
    anova_out_temp <- anova(model_temp)
    
    if(length(IndVars)==2){ ## Pulls model output for variable that is last
      if(RUN==1){ ### This currently works for 
        anova_partial_temp <- rbind(anova_out_temp[(length(IndVars)+1),],
                                    anova_out_temp[4,]
        )
      }
      if(RUN %in% 2:length(IndVars)){
        anova_partial_temp <- anova_out_temp[(length(IndVars)+1),]
        
      }
    } # End if vars == 2 statement
    anova_out <- rbind(anova_out, anova_partial_temp)
  } # End reorder loop
  return(anova_out)
} # End function


#########################################################
##### Absolute ANPP - Data manip ####

#Fort Keogh

#Read in this one
write.csv(FKANPP, file = "FKANPP.csv", row.names = FALSE)

#Thunder Basin
#Read in this one
write.csv(TBANPP, file = "TBANPP.csv", row.names = FALSE)

#checking on normality of residuals - FK
resid15 <- lm(data = subset(FKANPP, year == 2018), avg_annual_kg_h ~ drought*grazing_treatment)
ols_plot_resid_hist(resid15) #fairly normal
ols_test_normality(resid15) #Shapiro-wilk p 0.0149, K-S 0.621, C-M 0.00, A-D 0.069 - not normal

resid16 <- lm(data = subset(FKANPP, year == 2019), avg_annual_kg_h ~ drought*grazing_treatment)
ols_plot_resid_hist(resid16) #right skew
ols_test_normality(resid16) #failed normality tests

resid17 <- lm(data = subset(FKANPP, year == 2020), avg_annual_kg_h ~ drought*grazing_treatment)
ols_plot_resid_hist(resid17) #right skew
ols_test_normality(resid17) #failed normality tests

resid18 <- lm(data = subset(FKANPP, year == 2021), avg_annual_kg_h ~ drought*grazing_treatment)
ols_plot_resid_hist(resid18) #right skew
ols_test_normality(resid18) #failed normality tests

#checking on normality of residuals - TB
resida.1 <- lm(data = subset(TBANPP, year == 2018), avg_brome_kg_h ~ drought*grazing_treatment)
ols_plot_resid_hist(resida.1) #right skew
ols_test_normality(resida.1) #failed tests

resida.2 <- lm(data = subset(TBANPP, year == 2019), avg_brome_kg_h ~ drought*grazing_treatment)
ols_plot_resid_hist(resida.2) #right skew
ols_test_normality(resida.2) #failed normality tests

resida.3 <- lm(data = subset(TBANPP, year == 2020), avg_brome_kg_h ~ drought*grazing_treatment)
ols_plot_resid_hist(resida.3) #right skew
ols_test_normality(resida.3) #failed normality tests

resida.4 <- lm(data = subset(TBANPP, year == 2021), avg_brome_kg_h ~ drought*grazing_treatment)
ols_plot_resid_hist(resida.4) #right skew
ols_test_normality(resida.4) #failed normality tests


#########################################################
##### Absolute ANPP - Graphs ####


#Making figures of abs ANPP
#scatter plots with means and se 
FKANPP_avg_drought <- FKANPP %>%
  group_by(year, drought) %>%
  summarise(avg_ANPP = mean(avg_annual_kg_h), se_ANPP = sd(avg_annual_kg_h)/sqrt(length(avg_annual_kg_h))) %>%
  ungroup()

TBANPP_avg_drought <- TBANPP %>%
  group_by(year, drought) %>%
  summarise(avg_ANPP = mean(avg_brome_kg_h), se_ANPP = sd(avg_brome_kg_h)/sqrt(length(avg_brome_kg_h))) %>%
  ungroup()

#remake figs 2/17/22
FKANPP_avg_drought2 <- FKANPP_avg_drought %>%
  mutate(year = as.factor(year)) %>%
  filter(year != "2018")
TBANPP_avg_drought2 <- TBANPP_avg_drought %>%
  mutate(year = as.factor(year)) %>%
  filter(year != "2018")

FKdr_ANPP2 <- ggplot(data = FKANPP_avg_drought2, aes(x = drought, y = avg_ANPP, color = year, shape = year)) +
  geom_point(size = 3, position = position_dodge(0.5)) + 
  ylim(0, 1100) +
  geom_smooth(data = subset(FKANPP_avg_drought2, year == "2019"), aes(group = year), method = "lm", se = FALSE) +   
  geom_smooth(data = subset(FKANPP_avg_drought2, year == "2020"), aes(group = year), method = "lm", se = FALSE, linetype = "dashed") +   
  geom_errorbar(aes(ymin = avg_ANPP - se_ANPP, ymax = avg_ANPP + se_ANPP), width = 7, position = position_dodge(.5)) +
  scale_colour_manual(values = cbPalette) +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16)) +
  labs(x = "Precipitation Reduction (%)", y = "Annual Grass ANPP (kg ha-1)", color = "Year", shape = "Year", title = "MT") 

TBdr_ANPP2 <- ggplot(data = TBANPP_avg_drought2, aes(x = drought, y = avg_ANPP, color = year, shape = year)) +
  geom_point(size = 3, position = position_dodge(0.5)) + 
  geom_smooth(data = subset(TBANPP_avg_drought2, year != "2019"), aes(group = year), method = "lm", se = FALSE, linetype = "dashed") +   
  geom_errorbar(aes(ymin = avg_ANPP - se_ANPP, ymax = avg_ANPP + se_ANPP), width = 7, position = position_dodge(.5)) +
  scale_colour_manual(values = cbPalette) +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16)) +
  labs(x = "Precipitation Reduction (%)", y = "Brome ANPP (kg ha-1)", color = "Year", shape = "Year", title = "WY") 

FKdr_ANPP2 / TBdr_ANPP2




#graph graze vs abs ANPP with SE bars
#remake figs 2/17/22
#need to regroup and take out 2018-2019
FKANPP_avg_graze2 <- FKANPP %>%
  mutate(year = as.factor(year)) %>%
  filter(year %in% c("2020", "2021")) %>%
  mutate(realg = ifelse(year == "2020" & livestock_util_2019 == "50", "Stable", 
                ifelse(year == "2020" & livestock_util_2019 == "70", "Heavy", 
                ifelse(year == "2021" & livestock_util_2020 == "30", "Destock", 
                ifelse(year == "2021" & livestock_util_2020 == "50", "Stable",
                ifelse(year == "2021" & livestock_util_2020 == "70", "Heavy", "0")))))) %>%
  group_by(year, realg) %>%
  summarise(avg_ANPP = mean(avg_annual_kg_h), se_ANPP = sd(avg_annual_kg_h)/sqrt(length(avg_annual_kg_h))) %>%
  ungroup()

TBANPP_avg_graze2 <- TBANPP %>%
  mutate(year = as.factor(year)) %>%
  filter(year %in% c("2020", "2021")) %>%
  mutate(realg = ifelse(year == "2020" & livestock_util_2019 == "50", "Stable", 
                ifelse(year == "2020" & livestock_util_2019 == "70", "Heavy", 
                ifelse(year == "2021" & livestock_util_2020 == "30", "Destock", 
                ifelse(year == "2021" & livestock_util_2020 == "50", "Stable",
                ifelse(year == "2021" & livestock_util_2020 == "70", "Heavy", "0")))))) %>%
  group_by(year, realg) %>%
  summarise(avg_ANPP = mean(avg_brome_kg_h), se_ANPP = sd(avg_brome_kg_h)/sqrt(length(avg_brome_kg_h))) %>%
  ungroup()

#change order of x axis
FKANPP_avg_graze2$realg <- factor(FKANPP_avg_graze2$realg, levels = c("Destock", "Stable", "Heavy")) 
TBANPP_avg_graze2$realg <- factor(TBANPP_avg_graze2$realg, levels = c("Destock", "Stable", "Heavy")) 


FKgr_ANPP2 <- ggplot(data = FKANPP_avg_graze2, aes(x = realg, y = avg_ANPP, color = year, shape = year)) +
  geom_point(size = 3, position = position_dodge(0.5)) + 
  geom_errorbar(aes(ymin = avg_ANPP - se_ANPP, ymax = avg_ANPP + se_ANPP), width = 0.2, position = position_dodge(.5)) +
  scale_colour_manual(values = cbPalette) +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16)) +
  labs(x = "Grazing Treatment", y = "Annual Grass ANPP (kg ha-1)", color = "Year", shape = "Year", title = "MT") 

TBgr_ANPP2 <- ggplot(data = TBANPP_avg_graze2, aes(x = realg, y = avg_ANPP, color = year, shape = year)) +
  geom_point(size = 3, position = position_dodge(0.5)) + 
  geom_errorbar(aes(ymin = avg_ANPP - se_ANPP, ymax = avg_ANPP + se_ANPP), width = 0.2, position = position_dodge(.5)) +
  scale_colour_manual(values = cbPalette) +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16)) +
  labs(x = "Grazing Treatment", y = "Brome ANPP (kg ha-1)", color = "Year", shape = "Year", title = "WY") 

FKgr_ANPP2 / TBgr_ANPP2

#FK dxg was significant in 2021, so make that figure
FKANPP_dxg2021 <- FKANPP %>%
  mutate(year = as.factor(year)) %>%
  filter(year == "2021") %>%
  mutate(realg = ifelse(year == "2021" & livestock_util_2020 == "30", "Destock", 
                ifelse(year == "2021" & livestock_util_2020 == "50", "Stable",
                ifelse(year == "2021" & livestock_util_2020 == "70", "Heavy", "0")))) %>%
  group_by(year, realg, drought) %>%
  summarise(avg_ANPP = mean(avg_annual_kg_h), se_ANPP = sd(avg_annual_kg_h)/sqrt(length(avg_annual_kg_h))) %>%
  ungroup()

FKANPP_dxg2021$realg <- factor(FKANPP_dxg2021$realg, levels = c("Destock", "Stable", "Heavy")) 

#med g p = 0.03383, hi g p = 0.05611
FKANPP_dxg2021_2 <- ggplot(data = FKANPP_dxg2021, aes(x = factor(drought), y = avg_ANPP, color = realg, shape = realg)) +
  geom_point(size = 3, position = position_dodge(0.4)) + 
  ylim(0, 325) + 
  geom_smooth(data = subset(FKANPP_dxg2021, realg == "Stable"), aes(group = realg), method = "lm", se = FALSE) +   
  geom_smooth(data = subset(FKANPP_dxg2021, realg == "Heavy"), aes(group = realg), method = "lm", se = FALSE, linetype = "dashed") +   
  geom_errorbar(aes(ymin = avg_ANPP - se_ANPP, ymax = avg_ANPP + se_ANPP), width = .2, position = position_dodge(0.4)) +
  scale_colour_manual(values = cbPalette) +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16)) +
  labs(x = "Precipitation Reduction (%)", y = "Annual Grass ANPP (kg ha-1)", color = "Grazing Treatment", shape = "Grazing Treatment", title = "MT 2021") 



#make fig with d vs ANPP FK, interaction to the right, and TB below
FKdr_ANPP2 + FKANPP_dxg2021_2 + TBdr_ANPP2 + plot_layout(ncol = 2)

FKANPP_dxg2021_2 + plot_spacer()+ FKdr_ANPP2 + TBdr_ANPP2 + plot_layout(ncol = 2)




##for ESA 2022
FKdr_ANPP2 <- ggplot(data = FKANPP_avg_drought2, aes(x = drought, y = avg_ANPP, color = year, shape = year)) +
  geom_point(size = 3, position = position_dodge(0.5)) + 
  ylim(0, 1100) +
  geom_smooth(data = subset(FKANPP_avg_drought2, year == "2019"), aes(group = year), method = "lm", se = FALSE) +   
  geom_smooth(data = subset(FKANPP_avg_drought2, year == "2020"), aes(group = year), method = "lm", se = FALSE, linetype = "dashed") +   
  geom_errorbar(aes(ymin = avg_ANPP - se_ANPP, ymax = avg_ANPP + se_ANPP), width = 7, position = position_dodge(.5)) +
  scale_colour_manual(values = cbPalette) +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16), legend.position = "none") +
  labs(x = "Precipitation Reduction (%)", y = "ANPP (kg ha-1)", color = "Year", shape = "Year", title = "MT") 

TBdr_ANPP2 <- ggplot(data = TBANPP_avg_drought2, aes(x = drought, y = avg_ANPP, color = year, shape = year)) +
  geom_point(size = 3, position = position_dodge(0.5)) + 
  geom_smooth(data = subset(TBANPP_avg_drought2, year != "2019"), aes(group = year), method = "lm", se = FALSE, linetype = "dashed") +   
  geom_errorbar(aes(ymin = avg_ANPP - se_ANPP, ymax = avg_ANPP + se_ANPP), width = 7, position = position_dodge(.5)) +
  scale_colour_manual(values = cbPalette) +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16), axis.title.y = element_blank()) +
  labs(x = "Precipitation Reduction (%)", y = "ANPP (kg ha-1)", color = "Year", shape = "Year", title = "WY") 

FKdr_ANPP2 + TBdr_ANPP2

FKANPP_dxg2021_2 <- ggplot(data = FKANPP_dxg2021, aes(x = factor(drought), y = avg_ANPP, color = realg, shape = realg)) +
  geom_point(size = 3, position = position_dodge(0.4)) + 
  ylim(0, 325) + 
  geom_smooth(data = subset(FKANPP_dxg2021, realg == "Stable"), aes(group = realg), method = "lm", se = FALSE) +   
  geom_smooth(data = subset(FKANPP_dxg2021, realg == "Heavy"), aes(group = realg), method = "lm", se = FALSE, linetype = "dashed") +   
  geom_errorbar(aes(ymin = avg_ANPP - se_ANPP, ymax = avg_ANPP + se_ANPP), width = .2, position = position_dodge(0.4)) +
  scale_colour_manual(values = cbPalette) +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16)) +
  labs(x = "Precipitation Reduction (%)", y = "ANPP (kg ha-1)", color = "Grazing Treatment", shape = "Grazing Treatment", title = "MT 2021") 

(FKdr_ANPP2 + TBdr_ANPP2) / (FKANPP_dxg2021_2 + plot_spacer()) 

##


#########################################################
##### Absolute ANPP - Stats ####

#statistics on ANPP
FKANPP2018 <- FKANPP %>%
  filter(year == 2018)
FKANPP2019 <- FKANPP %>%
  filter(year == 2019)
FKANPP2020 <- FKANPP %>%
  filter(year == 2020)
FKANPP2021 <- FKANPP %>%
  filter(year == 2021)
TBANPP2018 <- TBANPP %>%
  filter(year == 2018)
TBANPP2019 <- TBANPP %>%
  filter(year == 2019)
TBANPP2020 <- TBANPP %>%
  filter(year == 2020)
TBANPP2021 <- TBANPP %>%
  filter(year == 2021)

#Stats splitting ANPP by year with transformations
#FK 2018 and 2021 normalish but try square root, or data ^0.75
#FK 2019-2020; TB 2018-2020 try square root, cube root, quad root, and so on. Then try ln(data+0.1) and log(data+0.1)

#Fort Keogh

#FK 2018 - sq root, data^0.75
FKANPP2018_sqrt <- FKANPP2018 %>%
  mutate(ANPP_sqrt = sqrt(avg_annual_kg_h)) %>%
  mutate(ANPP_75 = avg_annual_kg_h^0.75)

resid1.1 <- lm(data = FKANPP2018_sqrt, ANPP_sqrt ~ drought*grazing_treatment)
ols_plot_resid_hist(resid1.1) #fairly normal
ols_test_normality(resid1.1) #passed all tests except Cramer (p = 0.000)

#rerun stats with sqrt data and graphs
FKANPP2018_sqrt_stats <- lmerTest::lmer(data = FKANPP2018_sqrt, ANPP_sqrt ~ drought*grazing_treatment +
                                          (1|block) + (1|block:paddock))

anova(FKANPP2018_sqrt_stats, type = 3) #no significance


#FK 2019 - sq root, cube root, etc
FKANPP2019_sqrt <- FKANPP2019 %>%
  mutate(ANPP_sqrt = sqrt(avg_annual_kg_h)) %>%
  mutate(ANPP_cb = avg_annual_kg_h^1/3) %>%
  mutate(ANPP_qd = avg_annual_kg_h^1/4) %>%
  mutate(ANPP_ln = log(avg_annual_kg_h + 0.1)) %>%
  mutate(ANPP_log = log10(avg_annual_kg_h + 0.1)) %>%
  mutate(ANPP_5 = avg_annual_kg_h^1/5)

resid1.3 <- lm(data = FKANPP2019_sqrt, ANPP_ln ~ drought*grazing_treatment)
ols_plot_resid_hist(resid1.3) #fairly normal
ols_test_normality(resid1.3) #passed K-S, S-W, A-D but failed Cramer

#revisit stats - don't include grazing, because ANPP was collected before grazing occurred; try 2018 data as covariate
FKANPP2019_sqrt_stats2 <- lmerTest::lmer(data = FKANPP2019_sqrt, ANPP_ln ~ drought +
                                           (1|block) + (1|block:paddock))

anova(FKANPP2019_sqrt_stats2, type = 3) #drought p = 6.972e-05
AIC(FKANPP2019_sqrt_stats2) #AIC = 96.2456
AICc(FKANPP2019_sqrt_stats2) #AICc = 97.4956
rsquared(FKANPP2019_sqrt_stats2) #marg = 0.2120074


#2018 covariate
FKANPP2018_cov <- FKANPP2018 %>%
  dplyr:: select(c(plot, avg_annual_kg_h)) %>%
  mutate(absANPP2018 = avg_annual_kg_h) %>%
  dplyr:: select(-avg_annual_kg_h)

FKANPP2019_sqrt2 <- FKANPP2019_sqrt %>%
  left_join(FKANPP2018_cov)


FKANPP2019_sqrt_stats3 <- lmerTest::lmer(data = FKANPP2019_sqrt2, ANPP_ln ~ drought + absANPP2018 +
                                           (1|block) + (1|block:paddock))

anova(FKANPP2019_sqrt_stats3, type = 3) #drought p = 0.0001171, 2018ANPP not sig
AIC(FKANPP2019_sqrt_stats3) #AIC = 112.0906
AICc(FKANPP2019_sqrt_stats3) #AICc = 113.8778

#FK 2020 - sq root, cube root, etc
FKANPP2020_sqrt <- FKANPP2020 %>%
  mutate(ANPP_sqrt = sqrt(avg_annual_kg_h)) %>%
  mutate(ANPP_cb = avg_annual_kg_h^1/3) %>%
  mutate(ANPP_qd = avg_annual_kg_h^1/4) %>%
  mutate(ANPP_ln = log(avg_annual_kg_h + 0.1)) %>%
  mutate(ANPP_log = log10(avg_annual_kg_h + 0.1)) %>%
  mutate(ANPP_5 = avg_annual_kg_h^1/5)

resid1.5 <- lm(data = FKANPP2019_sqrt, ANPP_ln ~ drought*grazing_treatment)
ols_plot_resid_hist(resid1.3) #fairly normal
ols_test_normality(resid1.3) #passed K-S, S-W, A-D but failed Cramer

#rerun stats with ln data and graphs
FKANPP2020_sqrt_stats <- lmerTest::lmer(data = FKANPP2020_sqrt, ANPP_ln ~ drought*grazing_treatment +
                                          (1|block) + (1|block:paddock))

anova(FKANPP2020_sqrt_stats, type = 3) #no sig


#revisit stats - include actual grazing treatments from 2019, because ANPP was collected before grazing occurred; then try 2018 data as covariate
FKANPP2020_sqrt2 <- FKANPP2020_sqrt %>%
  left_join(FKANPP2018_cov) %>%
  mutate(livestock_util_2019 = as.factor(livestock_util_2019))

FKANPP2020_sqrt_stats2 <- lmerTest::lmer(data = FKANPP2020_sqrt2, ANPP_ln ~ drought*livestock_util_2019 +
                                           (1|block) + (1|block:paddock))

anova(FKANPP2020_sqrt_stats2, type = 3) #drought p = 0.08396, no other sig
AIC(FKANPP2020_sqrt_stats2) #AIC = 135.7612
AICc(FKANPP2020_sqrt_stats2) #AICc = 138.196
rsquared(FKANPP2020_sqrt_stats2) #marg = 0.05889012

#covariate
FKANPP2020_sqrt_stats3 <- lmerTest::lmer(data = FKANPP2020_sqrt2, ANPP_ln ~ drought*livestock_util_2019 + absANPP2018 +
                                           (1|block) + (1|block:paddock))

anova(FKANPP2020_sqrt_stats3, type = 3) #nothing significant
AIC(FKANPP2020_sqrt_stats3) #AIC = 150.9269
AICc(FKANPP2020_sqrt_stats3) #AICc = 154.1269

#covariate - trying as random effect
FKANPP2020_sqrt_stats3b <- lmerTest::lmer(data = FKANPP2020_sqrt2, ANPP_ln ~ drought*livestock_util_2019 + 
                                            (1|absANPP2018) + (1|block) + (1|block:paddock))

anova(FKANPP2020_sqrt_stats3b, type = 3) #drought p = 0.0684
AIC(FKANPP2020_sqrt_stats3b) #AIC = 137.6396
AICc(FKANPP2020_sqrt_stats3b) #AICc = 140.8396

#d + g, no interaction
FKANPP2020_sqrt_stats2b <- lmerTest::lmer(data = FKANPP2020_sqrt2, ANPP_ln ~ drought + livestock_util_2019 +
                                            (1|block) + (1|block:paddock))

anova(FKANPP2020_sqrt_stats2b, type = 3) #drought p = 0.06956, no other sig
AIC(FKANPP2020_sqrt_stats2b) #AIC = 124.8508
AICc(FKANPP2020_sqrt_stats2b) #AICc = 126.6381



#best model selection - try grazing treatment as numeric
#grazing as numeric
FKANPP2020_sqrt$livestock_util_2019 <- as.numeric(FKANPP2020_sqrt$livestock_util_2019)

FKANPP2020_sqrt_stats4 <- lmerTest::lmer(data = FKANPP2020_sqrt, ANPP_ln ~ drought*livestock_util_2019 +
                                           (1|block) + (1|block:paddock))

anova(FKANPP2020_sqrt_stats4, type = 3) #no significance
AIC(FKANPP2020_sqrt_stats4) #AIC = 147.7442
AICc(FKANPP2020_sqrt_stats4) #AICc = 150.1789

#recode graze values
FKANPP2020_sqrt <- FKANPP2020_sqrt %>%
  mutate(g2019 = ifelse(livestock_util_2019 == "50", "-1", 
                        ifelse(livestock_util_2019 == "70", "1", livestock_util_2019))) %>%
  mutate(g2019 = as.numeric(g2019))

FKANPP2020_sqrt_stats5 <- lmerTest::lmer(data = FKANPP2020_sqrt, ANPP_ln ~ drought*g2019 +
                                           (1|block) + (1|block:paddock))

anova(FKANPP2020_sqrt_stats5, type = 3) #drought p = 0.08396
AIC(FKANPP2020_sqrt_stats5) #AIC = 138.5338
AICc(FKANPP2020_sqrt_stats5) #AICc = 140.9686

#try drought + grazing, no interaction
FKANPP2020_sqrt_stats4b <- lmerTest::lmer(data = FKANPP2020_sqrt, ANPP_ln ~ drought + livestock_util_2019 +
                                            (1|block) + (1|block:paddock))

anova(FKANPP2020_sqrt_stats4b, type = 3) #drought p = 0.06956
AIC(FKANPP2020_sqrt_stats4b) #AIC = 130.8423
AICc(FKANPP2020_sqrt_stats4b) #AICc = 132.6295


FKANPP2020_sqrt_stats5b <- lmerTest::lmer(data = FKANPP2020_sqrt, ANPP_ln ~ drought + g2019 +
                                            (1|block) + (1|block:paddock))

anova(FKANPP2020_sqrt_stats5b, type = 3) #drought p = 0.06956
AIC(FKANPP2020_sqrt_stats5b) #AIC = 126.2371
AICc(FKANPP2020_sqrt_stats5b) #AICc = 128.0244


#just drought
FKANPP2020_sqrt_stats6 <- lmerTest::lmer(data = FKANPP2020_sqrt, ANPP_ln ~ drought +
                                           (1|block) + (1|block:paddock))

anova(FKANPP2020_sqrt_stats6, type = 3) #drought p = 0.06956
AIC(FKANPP2020_sqrt_stats6) #AIC = 121.9916
AICc(FKANPP2020_sqrt_stats6) #AICc = 123.2416




#FK 2021 - sq root, data^0.75
FKANPP2021_sqrt <- FKANPP2021 %>%
  mutate(ANPP_sqrt = sqrt(avg_annual_kg_h)) %>%
  mutate(ANPP_75 = avg_annual_kg_h^0.75)

resid1.6 <- lm(data = FKANPP2021_sqrt, ANPP_sqrt ~ drought*grazing_treatment)
ols_plot_resid_hist(resid1.6) #fairly normal
ols_test_normality(resid1.6) #passed all tests except Cramer (p = 0.000)

#rerun stats with sqrt data and graphs
FKANPP2021_sqrt_stats <- lmerTest::lmer(data = FKANPP2021_sqrt, ANPP_sqrt ~ drought*grazing_treatment +
                                          (1|block) + (1|block:paddock))

anova(FKANPP2021_sqrt_stats, type = 3) #drought*graze p = 0.02868

#revisit stats - include actual grazing treatments from 2019 + 2020, because ANPP was collected before grazing occurred; then try 2018 data as covariate
FKANPP2021_sqrt2 <- FKANPP2021_sqrt %>%
  left_join(FKANPP2018_cov) %>%
  mutate(graze2021 = ifelse(grazing_category == "HHMMM", "HH", 
                            ifelse(grazing_category == "MLLMM", "ML", 
                                   ifelse(grazing_category == "MMMMM", "MM", grazing_category)))) %>%
  mutate(graze2021 = as.factor(graze2021))


FKANPP2021_sqrt_stats2 <- lmerTest::lmer(data = FKANPP2021_sqrt2, ANPP_sqrt ~ drought*graze2021 +
                                           (1|block) + (1|block:paddock))

anova(FKANPP2021_sqrt_stats2, type = 3) #drought*graze p = 0.02868
AIC(FKANPP2021_sqrt_stats2) #AIC = 331.6801
AICc(FKANPP2021_sqrt_stats2) #AICc = 335.771

#trying to find which combos from interaction are significantly different
pair1 <- ls_means(FKANPP2021_sqrt_stats2)
contrast(pair1, "pairwise")

difflsmeans(FKANPP2021_sqrt_stats2)
show_tests(ls_means(FKANPP2021_sqrt_stats2))

TukeyHSD(FKANPP2021_sqrt_stats2, "graze2021")

summary(glht(FKANPP2021_sqrt_stats2, linfct = mcp(graze2021 = "Tukey")), test = adjusted(type = "BH")) #this seems like best bet; can only get factored comparisons
#no significance from this


#instead separate by grazing treatment and assess significance
FKANPP2021_sqrt_low <- FKANPP2021_sqrt2 %>%
  filter(graze2021 == "ML")
FKANPP2021_sqrt_med <- FKANPP2021_sqrt2 %>%
  filter(graze2021 == "MM")
FKANPP2021_sqrt_hi <- FKANPP2021_sqrt2 %>%
  filter(graze2021 == "HH")

FKANPP2021_sqrt_statsl <- lmerTest::lmer(data = FKANPP2021_sqrt_low, ANPP_sqrt ~ drought +
                                           (1|block) + (1|block:paddock))

anova(FKANPP2021_sqrt_statsl, type = 3) #no sig

FKANPP2021_sqrt_statsm <- lmerTest::lmer(data = FKANPP2021_sqrt_med, ANPP_sqrt ~ drought +
                                           (1|block) + (1|block:paddock))

anova(FKANPP2021_sqrt_statsm, type = 3) #p = 0.03383
rsquared(FKANPP2021_sqrt_statsm) #marg R = 0.1343245, cond R = 0.5872144

FKANPP2021_sqrt_statsh <- lmerTest::lmer(data = FKANPP2021_sqrt_hi, ANPP_sqrt ~ drought +
                                           (1|block) + (1|block:paddock))

anova(FKANPP2021_sqrt_statsh, type = 3) #p = 0.05611
rsquared(FKANPP2021_sqrt_statsh) #marg R = 0.107563, cond R = 0.5784025





#d+g, no interaction, graze as factor
FKANPP2021_sqrt_stats2b <- lmerTest::lmer(data = FKANPP2021_sqrt2, ANPP_sqrt ~ drought + graze2021 +
                                            (1|block) + (1|block:paddock))

anova(FKANPP2021_sqrt_stats2b, type = 3) #no significance
AIC(FKANPP2021_sqrt_stats2b) #AIC = 325.0392
AICc(FKANPP2021_sqrt_stats2b) #AICc = 327.4739


#covariate
FKANPP2021_sqrt_stats3 <- lmerTest::lmer(data = FKANPP2021_sqrt2, ANPP_sqrt ~ drought*graze2021 + absANPP2018 +
                                           (1|block) + (1|block:paddock))

anova(FKANPP2021_sqrt_stats3, type = 3) #drought*graze p = 0.03246
AIC(FKANPP2021_sqrt_stats3) #AIC = 343.9655
AICc(FKANPP2021_sqrt_stats3) #AICc = 349.0818

#covariate as random effect - got lots of warnings about failing to converge
FKANPP2021_sqrt_stats3b <- lmerTest::lmer(data = FKANPP2021_sqrt2, ANPP_sqrt ~ drought*graze2021 + 
                                            + (1|absANPP2018) + (1|block) + (1|block:paddock))

anova(FKANPP2021_sqrt_stats3b, type = 3) #drought*graze p = 0.02065
AIC(FKANPP2021_sqrt_stats3b) #AIC = 332.2134
AICc(FKANPP2021_sqrt_stats3b) #AICc = 337.3297


#trying additional models
#grazing as numeric, interaction of d*g
FKANPP2021_sqrt2$livestock_util_2020 <- as.numeric(FKANPP2021_sqrt2$livestock_util_2020)
FKANPP2021_sqrt_stats4 <- lmerTest::lmer(data = FKANPP2021_sqrt2, ANPP_sqrt ~ drought*livestock_util_2020 +
                                           (1|block) + (1|block:paddock))

anova(FKANPP2021_sqrt_stats4, type = 3) #drought p = 0.06251, drought*graze p = 0.07039
AIC(FKANPP2021_sqrt_stats4) #AIC = 344.9457
AICc(FKANPP2021_sqrt_stats4) #AICc = 347.3804

#grazing as numeric, d+g instead
FKANPP2021_sqrt_stats4b <- lmerTest::lmer(data = FKANPP2021_sqrt2, ANPP_sqrt ~ drought + livestock_util_2020 +
                                            (1|block) + (1|block:paddock))

anova(FKANPP2021_sqrt_stats4b, type = 3) #no significance
AIC(FKANPP2021_sqrt_stats4b) #AIC = 334.1401
AICc(FKANPP2021_sqrt_stats4b) #AICc = 335.9273

#recode grazing, interaction of d*g
FKANPP2021_sqrt3 <- FKANPP2021_sqrt2 %>%
  mutate(g2020 = ifelse(livestock_util_2020 == "30", "-1", 
                        ifelse(livestock_util_2020 == "50", "0", 
                               ifelse(livestock_util_2020 == "70", "1", livestock_util_2020))))

FKANPP2021_sqrt_stats5 <- lmerTest::lmer(data = FKANPP2021_sqrt3, ANPP_sqrt ~ drought*g2020 +
                                           (1|block) + (1|block:paddock))

anova(FKANPP2021_sqrt_stats5, type = 3) #drought*graze p = 0.02868
AIC(FKANPP2021_sqrt_stats5) #AIC = 331.6801
AICc(FKANPP2021_sqrt_stats5) #AICc = 335.771

#recode grazing, d+g no interaction
FKANPP2021_sqrt_stats5b <- lmerTest::lmer(data = FKANPP2021_sqrt3, ANPP_sqrt ~ drought + g2020 +
                                            (1|block) + (1|block:paddock))

anova(FKANPP2021_sqrt_stats5b, type = 3) #no significance
AIC(FKANPP2021_sqrt_stats5b) #AIC = 325.0392
AICc(FKANPP2021_sqrt_stats5b) #AICc = 327.4739

#just drought, no grazing
FKANPP2021_sqrt_stats6 <- lmerTest::lmer(data = FKANPP2021_sqrt3, ANPP_sqrt ~ drought +
                                           (1|block) + (1|block:paddock))

anova(FKANPP2021_sqrt_stats6, type = 3) #no significance
AIC(FKANPP2021_sqrt_stats6) #AIC = 328.3599
AICc(FKANPP2021_sqrt_stats6) #AICc = 329.6099




#Thunder Basin

#TB 2018 - sq root, cube root, etc
TBANPP2018_sqrt <- TBANPP2018 %>%
  mutate(ANPP_sqrt = sqrt(avg_brome_kg_h)) %>%
  mutate(ANPP_cb = avg_brome_kg_h^1/3) %>%
  mutate(ANPP_qd = avg_brome_kg_h^1/4) %>%
  mutate(ANPP_ln = log(avg_brome_kg_h + 0.1)) %>%
  mutate(ANPP_log = log10(avg_brome_kg_h + 0.1)) %>%
  mutate(ANPP_5 = avg_brome_kg_h^1/5)

resid1.7 <- lm(data = TBANPP2018_sqrt, ANPP_sqrt ~ drought*grazing_treatment)
ols_plot_resid_hist(resid1.7) #right skew
ols_test_normality(resid1.7) #passed K-S but failed other tests

#rerun stats with sqrt data
TBANPP2018_sqrt_stats <- lmerTest::lmer(data = TBANPP2018_sqrt, ANPP_sqrt ~ drought*grazing_treatment +
                                          (1|block) + (1|block:paddock))

anova(TBANPP2018_sqrt_stats, type = 3) #no sig


#TB 2019 - sq root, cube root, etc
TBANPP2019_sqrt <- TBANPP2019 %>%
  mutate(ANPP_sqrt = sqrt(avg_brome_kg_h)) %>%
  mutate(ANPP_cb = avg_brome_kg_h^1/3) %>%
  mutate(ANPP_qd = avg_brome_kg_h^1/4) %>%
  mutate(ANPP_ln = log(avg_brome_kg_h + 0.1)) %>%
  mutate(ANPP_log = log10(avg_brome_kg_h + 0.1)) %>%
  mutate(ANPP_5 = avg_brome_kg_h^1/5)

resid1.14 <- lm(data = TBANPP2019_sqrt, ANPP_ln ~ drought*grazing_treatment)
ols_plot_resid_hist(resid1.14) #left skew
ols_test_normality(resid1.14) #passed K-S, but failed other tests

#rerun stats with ln data
TBANPP2019_sqrt_stats <- lmerTest::lmer(data = TBANPP2019_sqrt, ANPP_ln ~ drought*grazing_treatment +
                                          (1|block) + (1|block:paddock))

anova(TBANPP2019_sqrt_stats, type = 3) #no sig

#revisit stats with just drought for 2019
TBANPP2019_sqrt_stats2 <- lmerTest::lmer(data = TBANPP2019_sqrt, ANPP_ln ~ drought +
                                           (1|block) + (1|block:paddock))

anova(TBANPP2019_sqrt_stats2, type = 3) #no sig


#TB 2020 - sq root, cube root, etc
TBANPP2020_sqrt <- TBANPP2020 %>%
  mutate(ANPP_sqrt = sqrt(avg_brome_kg_h)) %>%
  mutate(ANPP_cb = avg_brome_kg_h^1/3) %>%
  mutate(ANPP_qd = avg_brome_kg_h^1/4) %>%
  mutate(ANPP_ln = log(avg_brome_kg_h + 0.1)) %>%
  mutate(ANPP_log = log10(avg_brome_kg_h + 0.1)) %>%
  mutate(ANPP_5 = avg_brome_kg_h^1/5)

resid1.20 <- lm(data = TBANPP2020_sqrt, ANPP_ln ~ drought*grazing_treatment)
ols_plot_resid_hist(resid1.20) #normal ish looking
ols_test_normality(resid1.20) #passed K-S, but failed other tests (p value 0.052 S-W and 0.0415 A-D)


#rerun stats with ln data
TBANPP2020_sqrt_stats <- lmerTest::lmer(data = TBANPP2020_sqrt, ANPP_ln ~ drought*grazing_treatment +
                                          (1|block) + (1|block:paddock))

anova(TBANPP2020_sqrt_stats, type = 3) #no sig

#2020 stats with drought and 2019 grazing treatment
TBANPP2020_sqrt$livestock_util_2019 <- as.factor(TBANPP2020_sqrt$livestock_util_2019)
TBANPP2020_sqrt_stats2 <- lmerTest::lmer(data = TBANPP2020_sqrt, ANPP_ln ~ drought*livestock_util_2019 +
                                           (1|block) + (1|block:paddock))

anova(TBANPP2020_sqrt_stats2, type = 3) #drought p = 0.05396
rsquared(TBANPP2020_sqrt_stats2) #marg = 0.02888395


#TB 2021 - sq root, cube root, etc - none worked, so stick with untransformed
TBANPP2021_sqrt <- TBANPP2021 %>%
  mutate(ANPP_sqrt = sqrt(avg_brome_kg_h)) %>%
  mutate(ANPP_cb = avg_brome_kg_h^1/3) %>%
  mutate(ANPP_qd = avg_brome_kg_h^1/4) %>%
  mutate(ANPP_ln = log(avg_brome_kg_h + 0.1)) %>%
  mutate(ANPP_log = log10(avg_brome_kg_h + 0.1)) %>%
  mutate(ANPP_5 = avg_brome_kg_h^1/5) %>%
  mutate(ANPP_sq = avg_brome_kg_h^2) %>%
  mutate(ANPP_cr = avg_brome_kg_h^0.75)

resid1.25 <- lm(data = TBANPP2021_sqrt, ANPP_sqrt ~ drought*grazing_treatment)
ols_plot_resid_hist(resid1.25) #right skew
ols_test_normality(resid1.25) #failed all tests 

resid1.26 <- lm(data = TBANPP2021_sqrt, ANPP_ln ~ drought*grazing_treatment)
ols_plot_resid_hist(resid1.26) #left skew
ols_test_normality(resid1.26) #failed all tests 

resid1.27 <- lm(data = TBANPP2021_sqrt, ANPP_sq ~ drought*grazing_treatment)
ols_plot_resid_hist(resid1.27) #right skew
ols_test_normality(resid1.27) #failed all tests 

resid1.28 <- lm(data = TBANPP2021_sqrt, ANPP_cb ~ drought*grazing_treatment)
ols_plot_resid_hist(resid1.28) #right skew
ols_test_normality(resid1.28) #failed all tests 

resid1.29 <- lm(data = TBANPP2021_sqrt, ANPP_qd ~ drought*grazing_treatment)
ols_plot_resid_hist(resid1.29) #right skew
ols_test_normality(resid1.29) #failed all tests 

resid1.30 <- lm(data = TBANPP2021_sqrt, ANPP_log ~ drought*grazing_treatment)
ols_plot_resid_hist(resid1.30) #right skew
ols_test_normality(resid1.30) #failed all tests 

resid1.31 <- lm(data = TBANPP2021_sqrt, ANPP_5 ~ drought*grazing_treatment)
ols_plot_resid_hist(resid1.31) #right skew
ols_test_normality(resid1.31) #failed all tests 

resid1.32 <- lm(data = TBANPP2021_sqrt, ANPP_cr ~ drought*grazing_treatment)
ols_plot_resid_hist(resid1.32) #right skew
ols_test_normality(resid1.32) #failed all tests 

#no transformation worked so using untransformed data
TBbromeANPP2021stats <- lmerTest::lmer(data = TBANPP2021, avg_brome_kg_h ~ drought*grazing_treatment +
                                         (1|block) + (1|block:paddock))

anova(TBbromeANPP2021stats, type = 3) #no sig

#revisist stats with drought and 2020 grazing
TBANPP2021_2 <- TBANPP2021 %>%
  mutate(graze2021 = ifelse(grazing_category == "HHMMM", "HH", 
                            ifelse(grazing_category == "MLLMM", "ML", 
                                   ifelse(grazing_category == "MMMMM", "MM", grazing_category)))) %>%
  mutate(graze2021 = as.factor(graze2021))


TBbromeANPP2021stats2 <- lmerTest::lmer(data = TBANPP2021_2, avg_brome_kg_h ~ drought*graze2021 +
                                          (1|block) + (1|block:paddock))

anova(TBbromeANPP2021stats2, type = 3) #drought p = 0.07127
rsquared(TBbromeANPP2021stats2) #marg = 0.08106277


####################################################################3







####################################################################################
###### Phenology Data - through time - Data manip #############

#Read in this one
write.csv(TBphen_clean, file = "TBphen_clean.csv", row.names = FALSE)

#Read in this one
write.csv(FKphen_clean, file = "FKphen_clean.csv", row.names = FALSE)


#FK drought
FKphen_2 <- FKphen_clean %>%
  filter(month %in% c(5, 6, 7, 8)) %>%
  drop_na(pct_green) %>%
  group_by(year, month, date, drought) %>%
  summarise(avg_pct_green = mean(pct_green)) %>%
  ungroup() %>%
  mutate(month_day = ifelse(date == "2019-05-15", "May 15", 
                    ifelse(date == "2019-05-28", "May 28", 
                    ifelse(date == "2019-06-10", "June 10", 
                    ifelse(date == "2019-06-24", "June 24", 
                    ifelse(date == "2019-07-08", "July 8",
                    ifelse(date == "2019-07-19", "July 19", 
                    ifelse(date == "2019-08-05", "August 5", 
                    ifelse(date == "2019-08-27", "August 27", 
                    ifelse(date == "2020-05-06", "May 6", 
                    ifelse(date == "2020-05-18", "May 18", 
                    ifelse(date == "2020-06-03", "June 3", 
                    ifelse(date == "2020-06-16", "June 16", 
                    ifelse(date == "2020-07-01", "July 1", 
                    ifelse(date == "2020-07-13", "July 13", 
                    ifelse(date == "2020-07-26", "July 26", 
                    ifelse(date == "2020-08-12", "August 12", 
                    ifelse(date == "2020-08-24", "August 24", 
                    ifelse(date == "2021-05-07", "May 7", 
                    ifelse(date == "2021-05-17", "May 17", 
                    ifelse(date == "2021-06-03", "June 3", 
                    ifelse(date == "2021-06-16", "June 16", 
                    ifelse(date == "2021-07-01", "July 1", 
                    ifelse(date == "2021-07-14", "July 14", 
                    ifelse(date == "2021-08-05", "August 5", 
                    ifelse(date == "2021-08-18", "August 18", month))))))))))))))))))))))))))

#FK graze
FKphen_3 <- FKphen_clean %>%
  filter(month %in% c(5, 6, 7, 8)) %>%
  drop_na(pct_green) %>%
  group_by(year, month, date, grazing_treatment) %>%
  summarise(avg_pct_green = mean(pct_green)) %>%
  ungroup() %>%
  mutate(month_day = ifelse(date == "2019-05-15", "May 15", 
                    ifelse(date == "2019-05-28", "May 28", 
                    ifelse(date == "2019-06-10", "June 10", 
                    ifelse(date == "2019-06-24", "June 24", 
                  ifelse(date == "2019-07-08", "July 8",
                  ifelse(date == "2019-07-19", "July 19", 
                  ifelse(date == "2019-08-05", "August 5", 
                  ifelse(date == "2019-08-27", "August 27", 
                  ifelse(date == "2020-05-06", "May 6", 
                  ifelse(date == "2020-05-18", "May 18", 
                  ifelse(date == "2020-06-03", "June 3", 
                  ifelse(date == "2020-06-16", "June 16", 
                  ifelse(date == "2020-07-01", "July 1", 
                  ifelse(date == "2020-07-13", "July 13", 
                  ifelse(date == "2020-07-26", "July 26", 
                  ifelse(date == "2020-08-12", "August 12", 
                  ifelse(date == "2020-08-24", "August 24", 
                  ifelse(date == "2021-05-07", "May 7", 
                  ifelse(date == "2021-05-17", "May 17", 
                  ifelse(date == "2021-06-03", "June 3", 
                  ifelse(date == "2021-06-16", "June 16", 
                  ifelse(date == "2021-07-01", "July 1", 
                  ifelse(date == "2021-07-14", "July 14", 
                  ifelse(date == "2021-08-05", "August 5", 
                  ifelse(date == "2021-08-18", "August 18", month))))))))))))))))))))))))))

FKphen_2$month_day <- factor(FKphen_2$month_day, levels = c("May 6", "May 7", "May 15", "May 17", "May 18", "May 28", 
                                                            "June 3", "June 10", "June 16", "June 24",
                                                            "July 1", "July 8", "July 13", "July 14", "July 19", "July 26", 
                                                            "August 5", "August 12", "August 18", "August 24", "August 27"))

FKphen_3$month_day <- factor(FKphen_3$month_day, levels = c("May 6", "May 7", "May 15", "May 17", "May 18", "May 28", 
                                                            "June 3", "June 10", "June 16", "June 24",
                                                            "July 1", "July 8", "July 13", "July 14", "July 19", "July 26", 
                                                            "August 5", "August 12", "August 18", "August 24", "August 27"))


#TB drought
TBphen_2 <- TBphen_clean %>%
  filter(month %in% c(5, 6, 7, 8)) %>%
  #filter(species == "BRAR_pct_green") %>%
  drop_na(pct_green) %>%
  group_by(year, month, date, drought) %>%
  summarise(avg_pct_green = mean(pct_green)) %>%
  ungroup() %>%
  mutate(month_day = ifelse(date == "2019-05-02", "May 2", 
                    ifelse(date == "2019-05-19", "May 19", 
                    ifelse(date == "2019-06-03", "June 3", 
                    ifelse(date == "2019-06-14", "June 14", 
                    ifelse(date == "2019-06-24", "June 24",
                    ifelse(date == "2019-07-08", "July 8", 
                    ifelse(date == "2019-07-23", "July 23", 
                    ifelse(date == "2019-08-09", "August 9", 
                    ifelse(date == "2019-08-18", "August 18", 
                    ifelse(date == "2019-08-27", "August 27", 
                    ifelse(date == "2020-05-01", "May 1", 
                    ifelse(date == "2020-05-06", "May 6", 
                    ifelse(date == "2020-05-19", "May 19", 
                    ifelse(date == "2020-06-01", "June 1", 
                    ifelse(date == "2020-06-16", "June 16", 
                    ifelse(date == "2020-06-18", "June 18", 
                    ifelse(date == "2020-06-29", "June 29", 
                    ifelse(date == "2020-06-30", "June 30", 
                    ifelse(date == "2020-07-16", "July 16", 
                    ifelse(date == "2020-07-30", "July 30", 
                    ifelse(date == "2020-08-11", "August 11", 
                    ifelse(date == "2020-08-27", "August 27", 
                    ifelse(date == "2021-05-06", "May 6", 
                    ifelse(date == "2021-05-19", "May 19", 
                    ifelse(date == "2021-06-02", "June 2", 
                    ifelse(date == "2021-06-16", "June 16", 
                    ifelse(date == "2021-06-29", "June 29", 
                    ifelse(date == "2021-07-14", "July 14", 
                    ifelse(date == "2021-07-27", "July 27", 
                    ifelse(date == "2021-08-11", "August 11", 
                    ifelse(date == "2021-08-25", "August 25", month))))))))))))))))))))))))))))))))

#TB graze
TBphen_3 <- TBphen_clean %>%
  filter(month %in% c(5, 6, 7, 8)) %>%
  #filter(species == "BRAR_pct_green") %>%
  drop_na(pct_green) %>%
  group_by(year, month, date, grazing_treatment) %>%
  summarise(avg_pct_green = mean(pct_green)) %>%
  ungroup() %>%
  mutate(month_day = ifelse(date == "2019-05-02", "May 2", 
                            ifelse(date == "2019-05-19", "May 19", 
                            ifelse(date == "2019-06-03", "June 3", 
                            ifelse(date == "2019-06-14", "June 14", 
                            ifelse(date == "2019-06-24", "June 24",
                            ifelse(date == "2019-07-08", "July 8", 
                            ifelse(date == "2019-07-23", "July 23", 
                            ifelse(date == "2019-08-09", "August 9", 
                            ifelse(date == "2019-08-18", "August 18", 
                            ifelse(date == "2019-08-27", "August 27", 
                            ifelse(date == "2020-05-01", "May 1", 
                            ifelse(date == "2020-05-06", "May 6", 
                            ifelse(date == "2020-05-19", "May 19", 
                            ifelse(date == "2020-06-01", "June 1", 
                            ifelse(date == "2020-06-16", "June 16", 
                            ifelse(date == "2020-06-18", "June 18", 
                            ifelse(date == "2020-06-29", "June 29", 
                            ifelse(date == "2020-06-30", "June 30", 
                            ifelse(date == "2020-07-16", "July 16", 
                            ifelse(date == "2020-07-30", "July 30", 
                            ifelse(date == "2020-08-11", "August 11", 
                            ifelse(date == "2020-08-27", "August 27", 
                            ifelse(date == "2021-05-06", "May 6", 
                            ifelse(date == "2021-05-19", "May 19", 
                            ifelse(date == "2021-06-02", "June 2", 
                            ifelse(date == "2021-06-16", "June 16", 
                            ifelse(date == "2021-06-29", "June 29", 
                            ifelse(date == "2021-07-14", "July 14", 
                            ifelse(date == "2021-07-27", "July 27", 
                            ifelse(date == "2021-08-11", "August 11", 
                            ifelse(date == "2021-08-25", "August 25", month))))))))))))))))))))))))))))))))

TBphen_2$month_day <- factor(TBphen_2$month_day, levels = c("May 1", "May 2", "May 6", "May 19", 
                                                            "June 1", "June 2", "June 3", "June 14", "June 16", "June 18", "June 24", "June 29", "June 30", 
                                                            "July 8", "July 16", "July 14", "July 23", "July 27", "July 30",
                                                            "August 9", "August 11", "August 18", "August 25", "August 27"))

TBphen_3$month_day <- factor(TBphen_3$month_day, levels = c("May 1", "May 2", "May 6", "May 19", 
                                                            "June 1", "June 2", "June 3", "June 14", "June 16", "June 18", "June 24", "June 29", "June 30", 
                                                            "July 8", "July 16", "July 14", "July 23", "July 27", "July 30",
                                                            "August 9", "August 11", "August 18", "August 25", "August 27"))

#redo again but take out august
FKphen_4 <- FKphen_2 %>%
  filter(month != 8)                                                                                                                                                                        

FKphen_5 <- FKphen_3 %>%
  filter(month != 8)

TBphen_4 <- TBphen_2 %>%
  filter(month != 8)

TBphen_5 <- TBphen_3 %>%
  filter(month != 8)

FKphen_4$month_day <- factor(FKphen_4$month_day, levels = c("May 6", "May 7", "May 15", "May 17", "May 18", "May 28", 
                                                            "June 3", "June 10", "June 16", "June 24",
                                                            "July 1", "July 8", "July 13", "July 14", "July 19", "July 26"))

FKphen_5$month_day <- factor(FKphen_5$month_day, levels = c("May 6", "May 7", "May 15", "May 17", "May 18", "May 28", 
                                                            "June 3", "June 10", "June 16", "June 24",
                                                            "July 1", "July 8", "July 13", "July 14", "July 19", "July 26"))

TBphen_4$month_day <- factor(TBphen_4$month_day, levels = c("May 1", "May 2", "May 6", "May 19", 
                                                            "June 1", "June 2", "June 3", "June 14", "June 16", "June 18", "June 24", "June 29", "June 30", 
                                                            "July 8", "July 16", "July 14", "July 23", "July 27", "July 30"))

TBphen_5$month_day <- factor(TBphen_5$month_day, levels = c("May 1", "May 2", "May 6", "May 19", 
                                                            "June 1", "June 2", "June 3", "June 14", "June 16", "June 18", "June 24", "June 29", "June 30", 
                                                            "July 8", "July 16", "July 14", "July 23", "July 27", "July 30"))



#try making each graph separately and replace each date with a numbered day of that year
#ex jan 1 = 1

FKphen_4_2019 <- FKphen_4 %>%
  filter(year == 2019) %>%
  mutate(num_day = ifelse(month_day == "May 15", "134", 
                  ifelse(month_day == "May 28", "147", 
                  ifelse(month_day == "June 10", "160", 
                  ifelse(month_day == "June 24", "174", 
                  ifelse(month_day == "July 8", "188", 
                  ifelse(month_day == "July 19", "198", month_day)))))))

FKphen_4_2020 <- FKphen_4 %>%
  filter(year == 2020) %>%
  mutate(num_day = ifelse(month_day == "May 6", "126", 
                  ifelse(month_day == "May 18", "138", 
                  ifelse(month_day == "June 3", "154", 
                  ifelse(month_day == "June 16", "167", 
                  ifelse(month_day == "July 1", "182", 
                  ifelse(month_day == "July 13", "194", 
                  ifelse(month_day == "July 26", "207", month_day))))))))

FKphen_4_2021 <- FKphen_4 %>%
  filter(year == 2021) %>%
  mutate(num_day = ifelse(month_day == "May 7", "126", 
                  ifelse(month_day == "May 17", "136", 
                  ifelse(month_day == "June 3", "153", 
                  ifelse(month_day == "June 16", "166", 
                  ifelse(month_day == "July 1", "181", 
                  ifelse(month_day == "July 14", "194", month_day)))))))

TBphen_4_2019 <- TBphen_4 %>%
  filter(year == 2019) %>%
  mutate(num_day = ifelse(month_day == "May 2", "121", 
                  ifelse(month_day == "May 19", "138", 
                  ifelse(month_day == "June 3", "153", 
                  ifelse(month_day == "June 14", "164",
                  ifelse(month_day == "June 24", "174",
                  ifelse(month_day == "July 8", "188", 
                  ifelse(month_day == "July 23", "203", month_day))))))))

TBphen_4_2020 <- TBphen_4 %>%
  filter(year == 2020) %>%
  mutate(num_day = ifelse(month_day == "May 1", "121", 
                  ifelse(month_day == "May 6", "126", 
                  ifelse(month_day == "May 19", "139", 
                  ifelse(month_day == "June 1", "152", 
                  ifelse(month_day == "June 16", "167", 
                  ifelse(month_day == "June 18", "169",
                  ifelse(month_day == "June 29", "180",
                  ifelse(month_day == "June 30", "181",
                  ifelse(month_day == "July 16", "197", 
                  ifelse(month_day == "July 30", "211", month_day)))))))))))

TBphen_4_2021 <- TBphen_4 %>%
  filter(year == 2021) %>%
  mutate(num_day = ifelse(month_day == "May 6", "125", 
                  ifelse(month_day == "May 19", "138", 
                  ifelse(month_day == "June 2", "152", 
                  ifelse(month_day == "June 16", "166", 
                  ifelse(month_day == "June 29", "179", 
                  ifelse(month_day == "July 14", "194", 
                  ifelse(month_day == "July 27", "207", month_day))))))))



####################################################################################
###### Phenology Data - through time - Graphs #############
#rearrange color order
cbPalette2 <- c("#0072B2", "#56B4E9", "#009E73", "#F0E442",  "#E69F00")

### through time standardizing by control
#FK
#need to make new df with control data then merge with other data, subtract (treatment - control) and use that diff in analyses and figs
FKphen_con <- FKphen_clean %>%
  filter(drought == 0) %>%
  group_by(date, year, block, paddock) %>%
  summarise(avg_0 = mean(pct_green)) %>%   #average 2 control plots together
  ungroup() %>%
  drop_na(avg_0)

FKphen_clean2 <- FKphen_clean %>%
  full_join(FKphen_con) %>%
  drop_na(pct_green) %>%
  filter(drought != 0) %>%  #take out control plots because we dont want to overinflate 0's
  mutate(pct_diff = pct_green - avg_0)


FKphen_std <- FKphen_clean2 %>%
  filter(month %in% c(5, 6, 7)) %>%
  group_by(year, month, date, drought) %>%
  summarise(avg_pct_green = mean(pct_diff)) %>%
  ungroup() %>%
  mutate(month_day = ifelse(date == "2019-05-15", "May 15", 
        ifelse(date == "2019-05-28", "May 28", 
        ifelse(date == "2019-06-10", "June 10", 
        ifelse(date == "2019-06-24", "June 24", 
        ifelse(date == "2019-07-08", "July 8",
        ifelse(date == "2019-07-19", "July 19", 
        ifelse(date == "2020-05-06", "May 6", 
        ifelse(date == "2020-05-18", "May 18", 
        ifelse(date == "2020-06-03", "June 3", 
        ifelse(date == "2020-06-16", "June 16", 
        ifelse(date == "2020-07-01", "July 1", 
        ifelse(date == "2020-07-13", "July 13", 
        ifelse(date == "2020-07-26", "July 26", 
        ifelse(date == "2021-05-07", "May 7", 
        ifelse(date == "2021-05-17", "May 17", 
        ifelse(date == "2021-06-03", "June 3", 
        ifelse(date == "2021-06-16", "June 16", 
        ifelse(date == "2021-07-01", "July 1", 
        ifelse(date == "2021-07-14", "July 14", month))))))))))))))))))))

#add in julian day for dates
FKphen_con_2019 <- FKphen_std %>%
  filter(year == 2019) %>%
  mutate(num_day = ifelse(month_day == "May 15", "134", 
                  ifelse(month_day == "May 28", "147", 
                  ifelse(month_day == "June 10", "160", 
                  ifelse(month_day == "June 24", "174", 
                  ifelse(month_day == "July 8", "188", 
                  ifelse(month_day == "July 19", "198", month_day)))))))

FKphen_con_2020 <- FKphen_std %>%
  filter(year == 2020) %>%
  mutate(num_day = ifelse(month_day == "May 6", "126", 
                  ifelse(month_day == "May 18", "138", 
                  ifelse(month_day == "June 3", "154", 
                  ifelse(month_day == "June 16", "167", 
                  ifelse(month_day == "July 1", "182", 
                  ifelse(month_day == "July 13", "194", 
                  ifelse(month_day == "July 26", "207", month_day))))))))

FKphen_con_2021 <- FKphen_std %>%
  filter(year == 2021) %>%
  mutate(num_day = ifelse(month_day == "May 7", "126", 
                  ifelse(month_day == "May 17", "136", 
                  ifelse(month_day == "June 3", "153", 
                  ifelse(month_day == "June 16", "166", 
                  ifelse(month_day == "July 1", "181", 
                  ifelse(month_day == "July 14", "194", month_day)))))))


#TB
TBphen_edit2 <- TBphen_edit %>%
  filter(month %in% c(5, 6, 7)) %>%
  full_join(metadata) %>%
  drop_na(pct_green) %>%
  filter(plot != 24 | date != "2019-06-24") %>%  #filter out dates that switched BRAR/BRTE
  filter(plot != 24 | date != "2019-07-08") %>%
  filter(plot != 24 | date != "2019-07-23") %>%
  filter(plot != 35 | date != "2019-07-08") %>%
  filter(plot != 35 | date != "2019-07-23")

TBphen_edit_con <- TBphen_edit2 %>%
  filter(drought %in% c(1, 2)) %>% #filter just control plots
  mutate(pct_con = pct_green) %>%
  dplyr::select(-c(pct_green, species)) %>% 
  spread(key = species2, value = pct_con) %>% #spread data to make 2 control columns
  rename(BRAR_con = BRAR, BRTE_con = BRTE, control = drought, plot_con = plot)

TBphen_edit_trt <- TBphen_edit2 %>%
  filter(drought != 1) %>%   #filter out control plots
  filter(drought != 2) %>%
  mutate(pct_trt = pct_green) %>%
  dplyr::select(-c(pct_green, species)) %>% 
  spread(key = species2, value = pct_trt) %>% #spread data to make 2 treatment columns
  rename(BRAR_trt = BRAR, BRTE_trt = BRTE)

TBphen_edit_combo <- TBphen_edit_con %>%
  left_join(TBphen_edit_trt)

TBphen_edit_combo2 <- TBphen_edit_combo %>%
  mutate(BRAR = BRAR_trt - BRAR_con, BRTE = BRTE_trt - BRTE_con) %>%
  group_by(date, month, year, block, paddock, grazing_category, grazing_treatment, livestock_util_2019, livestock_util_2020, livestock_util_2021, plot, drought) %>%
  gather(species, green_diff, 20:21) %>%
  ungroup() %>%
  drop_na(green_diff)

#need to average across date and plot so there arent duplicates of % diff when both controls had BRAR or BRTE
TBphen_edit_combo2.1 <- TBphen_edit_combo2 %>%
  group_by(date, month, year, block, paddock, grazing_category, grazing_treatment, livestock_util_2019, livestock_util_2020, livestock_util_2021, plot, drought) %>%
  summarise(green_diff = mean(green_diff)) %>%
  ungroup()

TBphen_edit_combo3 <- TBphen_edit_combo2.1 %>%
  group_by(year, month, date, drought) %>%
  summarise(avg_pct_diff = mean(green_diff)) %>%
  ungroup() %>%
  mutate(month_day = ifelse(date == "2019-05-02", "May 2", 
                    ifelse(date == "2019-05-19", "May 19", 
                    ifelse(date == "2019-06-03", "June 3", 
                    ifelse(date == "2019-06-14", "June 14", 
                    ifelse(date == "2019-06-24", "June 24",
                    ifelse(date == "2019-07-08", "July 8", 
                    ifelse(date == "2019-07-23", "July 23", 
                    ifelse(date == "2020-05-01", "May 1", 
                    ifelse(date == "2020-05-06", "May 6", 
                    ifelse(date == "2020-05-19", "May 19", 
                    ifelse(date == "2020-06-01", "June 1", 
                    ifelse(date == "2020-06-16", "June 16", 
                    ifelse(date == "2020-06-18", "June 18", 
                    ifelse(date == "2020-06-29", "June 29", 
                    ifelse(date == "2020-06-30", "June 30", 
                    ifelse(date == "2020-07-16", "July 16", 
                    ifelse(date == "2020-07-30", "July 30",
                    ifelse(date == "2021-05-06", "May 6", 
                    ifelse(date == "2021-05-19", "May 19", 
                    ifelse(date == "2021-06-02", "June 2", 
                    ifelse(date == "2021-06-16", "June 16", 
                    ifelse(date == "2021-06-29", "June 29", 
                    ifelse(date == "2021-07-14", "July 14", 
                    ifelse(date == "2021-07-27", "July 27",  month)))))))))))))))))))))))))


TBphen_edit_combo_2019 <- TBphen_edit_combo3 %>%
  filter(year == 2019) %>%
  mutate(num_day = ifelse(month_day == "May 2", "121", 
                ifelse(month_day == "May 19", "138", 
                ifelse(month_day == "June 3", "153", 
                ifelse(month_day == "June 14", "164",
                ifelse(month_day == "June 24", "174",
                ifelse(month_day == "July 8", "188", 
                ifelse(month_day == "July 23", "203", month_day))))))))

TBphen_edit_combo_2020 <- TBphen_edit_combo3 %>%
  filter(year == 2020) %>%
  mutate(num_day = ifelse(month_day == "May 1", "121", 
                ifelse(month_day == "May 6", "126", 
                ifelse(month_day == "May 19", "139", 
                ifelse(month_day == "June 1", "152", 
                ifelse(month_day == "June 16", "167", 
                ifelse(month_day == "June 18", "169",
                ifelse(month_day == "June 29", "180",
                ifelse(month_day == "June 30", "181",
                ifelse(month_day == "July 16", "197", 
                ifelse(month_day == "July 30", "211", month_day)))))))))))

TBphen_edit_combo_2021 <- TBphen_edit_combo3 %>%
  filter(year == 2021) %>%
  mutate(num_day = ifelse(month_day == "May 6", "125", 
                ifelse(month_day == "May 19", "138", 
                ifelse(month_day == "June 2", "152", 
                ifelse(month_day == "June 16", "166", 
                ifelse(month_day == "June 29", "179", 
                ifelse(month_day == "July 14", "194", 
                ifelse(month_day == "July 27", "207", month_day))))))))





#figs
#rearrange color order
cbPalette2 <- c("#0072B2", "#56B4E9", "#009E73", "#F0E442",  "#E69F00")

FKphen_std_2019 <- ggplot(FKphen_con_2019, aes(x = as.numeric(num_day), y = avg_pct_green, color = as.factor(drought), shape = as.factor(drought))) + 
  geom_point(size = 3) + 
  geom_line(aes(group = as.factor(drought))) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10)) +
  scale_colour_manual(values = cbPalette2) +
  ylim(-35, 40) +
  xlim(120, 200) +
  labs(x = "Day of Year", y = "Greenness (% Difference)", color = "Precipitation Reduction (%)", shape = "Precipitation Reduction (%)", title = "MT 2019") + 
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16), legend.position = "none") 

FKphen_std_2020 <- ggplot(FKphen_con_2020, aes(x = as.numeric(num_day), y = avg_pct_green, color = as.factor(drought), shape = as.factor(drought))) + 
  geom_point(size = 3) + 
  geom_line(aes(group = as.factor(drought))) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10)) +
  scale_colour_manual(values = cbPalette2) +
  ylim(-35, 40) +
  xlim(120, 210) +
  labs(x = "Day of Year", y = "Greenness (% Difference)", color = "Precipitation Reduction (%)", shape = "Precipitation Reduction (%)", title = "MT 2020") + 
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16), legend.position = "none")

FKphen_std_2021 <- ggplot(FKphen_con_2021, aes(x = as.numeric(num_day), y = avg_pct_green, color = as.factor(drought), shape = as.factor(drought))) + 
  geom_point(size = 3) + 
  geom_line(aes(group = as.factor(drought))) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10)) +
  scale_colour_manual(values = cbPalette2) +
  ylim(-35, 40) +
  xlim(120, 200) +
  labs(x = "Day of Year", y = "Greenness (% Difference)", color = "Precipitation Reduction (%)", shape = "Precipitation Reduction (%)", title = "MT 2021") + 
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16))

TBphen_std_2019 <- ggplot(TBphen_edit_combo_2019, aes(x = as.numeric(num_day), y = avg_pct_diff, color = as.factor(drought), shape = as.factor(drought))) + 
  geom_point(size = 3) + 
  geom_line(aes(group = as.factor(drought))) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10)) +
  scale_colour_manual(values = cbPalette2) +
  ylim(-35, 40) +
  xlim(120, 205) +
  labs(x = "Day of Year", y = "Greenness (% Difference)", color = "Precipitation Reduction (%)", shape = "Precipitation Reduction (%)", title = "WY 2019") + 
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16), legend.position = "none") 

TBphen_std_2020 <- ggplot(TBphen_edit_combo_2020, aes(x = as.numeric(num_day), y = avg_pct_diff, color = as.factor(drought), shape = as.factor(drought))) + 
  geom_point(size = 3) + 
  geom_line(aes(group = as.factor(drought))) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10)) +
  scale_colour_manual(values = cbPalette2) +
  ylim(-35, 40) +
  xlim(120, 215) +
  labs(x = "Day of Year", y = "Greenness (% Difference)", color = "Precipitation Reduction (%)", shape = "Precipitation Reduction (%)", title = "WY 2020") + 
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16), legend.position = "none") 

TBphen_std_2021 <- ggplot(TBphen_edit_combo_2021, aes(x = as.numeric(num_day), y = avg_pct_diff, color = as.factor(drought), shape = as.factor(drought))) + 
  geom_point(size = 3) + 
  geom_line(aes(group = as.factor(drought))) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10)) +
  scale_colour_manual(values = cbPalette2) +
  ylim(-35, 40) +
  xlim(120, 210) +
  labs(x = "Day of Year", y = "Greenness (% Difference)", color = "Precipitation Reduction (%)", shape = "Precipitation Reduction (%)", title = "WY 2021") + 
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16)) 

(FKphen_std_2019 + FKphen_std_2020 + FKphen_std_2021) / (TBphen_std_2019 + TBphen_std_2020 + TBphen_std_2021)
#FKphenstd 1700x500


####################################################################################
###### Phenology Data - through time - Statistics #############

#through time statistics - repeated measures Anovas
#stats on may-july phenology 
FKphen_stats <- FKphen_clean %>%
  filter(month %in% c(5, 6, 7))




#redo stats with standardized pct green
#FK
#stats on may-july phenology 
FKphen_stats2 <- FKphen_clean2 %>%
  filter(month %in% c(5, 6, 7))

#stats for each year with real drought and grazing treatments
#2019 just do drought
anova_t3_2ind(IndVars = c("date", "drought"), 
              DepVar = "pct_diff", 
              RndForm = "~1 | block/paddock/plot", 
              Data = subset(FKphen_stats2, year == 2019)) #drought p < 0.0001; date p = 0.5187; date*drought p = 0.6271

#2020 drought + 2019 grazing
FKphen_stats3 <- FKphen_stats2 %>%
  mutate(livestock_util_2019 = as.factor(livestock_util_2019))

anova_t3_3ind(IndVars = c("date", "drought", "livestock_util_2019"), 
              DepVar = "pct_diff", 
              RndForm = "~1 | block/paddock/plot", 
              Data = subset(FKphen_stats3, year == 2020)) #drought p = 0.0182

#2021 drought + 2019/2020 grazing
FKphen_stats4 <- FKphen_stats2 %>%
  mutate(graze2021 = ifelse(grazing_category == "HHMMM", "HH", 
                    ifelse(grazing_category == "MLLMM", "ML", 
                    ifelse(grazing_category == "MMMMM", "MM", grazing_category)))) %>%
  mutate(graze2021 = as.factor(graze2021))

anova_t3_3ind(IndVars = c("date", "drought", "graze2021"), 
              DepVar = "pct_diff", 
              RndForm = "~1 | block/paddock/plot", 
              Data = subset(FKphen_stats4, year == 2021)) #drought p = 0.0038


#TB
#2019 just do drought
anova_t3_2ind(IndVars = c("date", "drought"), 
              DepVar = "green_diff", 
              RndForm = "~1 | block/paddock/plot", 
              Data = subset(TBphen_edit_combo2, year == 2019)) #no sig

#2020 drought + 2019 grazing
TBphen_stats <- TBphen_edit_combo2 %>%
  mutate(livestock_util_2019 = as.factor(livestock_util_2019))

anova_t3_3ind(IndVars = c("date", "drought", "livestock_util_2019"), 
              DepVar = "green_diff", 
              RndForm = "~1 | block/paddock/plot", 
              Data = subset(TBphen_stats, year == 2020)) #no sig

#2021 drought + 2019/2020 grazing
TBphen_stats2 <- TBphen_edit_combo2 %>%
  mutate(graze2021 = ifelse(grazing_category == "HHMMM", "HH", 
                    ifelse(grazing_category == "MLLMM", "ML", 
                    ifelse(grazing_category == "MMMMM", "MM", grazing_category)))) %>%
  mutate(graze2021 = as.factor(graze2021))

anova_t3_3ind(IndVars = c("date", "drought", "graze2021"), 
              DepVar = "green_diff", 
              RndForm = "~1 | block/paddock/plot", 
              Data = subset(TBphen_stats2, year == 2021)) #drought p = 0.0855



###############################################################################

            
                
                
                
##############################################################################
######### Phenology single time point - Data manip ############
#single time point graphs
#group and average for just 1 date

FKphen2019_june24 <- FKphen_clean %>%
  filter(date == "2019-06-24") %>%
  group_by(site, year, drought) %>%
  summarise(avg = mean(pct_green), se = sd(pct_green)/sqrt(length(pct_green))) %>%
  ungroup()

FKphen2020_june16 <- FKphen_clean %>%
  filter(date == "2020-06-16") %>%
  group_by(site, year, drought) %>%
  summarise(avg = mean(pct_green), se = sd(pct_green)/sqrt(length(pct_green))) %>%
  ungroup()

FKphen2021_june16 <- FKphen_clean %>%
  filter(date == "2021-06-16") %>%
  group_by(site, year, drought) %>%
  summarise(avg = mean(pct_green), se = sd(pct_green)/sqrt(length(pct_green))) %>%
  ungroup()


TBphen2019_july8 <- TBphen_clean %>%
  filter(date == "2019-07-08") %>%
  drop_na(pct_green) %>%
  group_by(site, year, drought) %>%
  summarise(avg = mean(pct_green), se = sd(pct_green)/sqrt(length(pct_green))) %>%
  ungroup()

TBphen2020_june18 <- TBphen_clean %>%
  filter(date == "2020-06-18") %>%
  #filter(species == "BRAR_pct_green") %>%
  drop_na(pct_green) %>%
  group_by(site, year, drought) %>%
  summarise(avg = mean(pct_green), se = sd(pct_green)/sqrt(length(pct_green))) %>%
  ungroup()

TBphen2021_june16 <- TBphen_clean %>%
  filter(date == "2021-06-16") %>%
  #filter(species == "BRAR_pct_green") %>%
  drop_na(pct_green) %>%
  group_by(site, year, drought) %>%
  summarise(avg = mean(pct_green), se = sd(pct_green)/sqrt(length(pct_green))) %>%
  ungroup()






##use standardized diff
#FK
FKphen2019_june24_std <- FKphen_clean2 %>%
  filter(date == "2019-06-24") %>%
  group_by(site, year, drought) %>%
  summarise(avg = mean(pct_diff), se = sd(pct_diff)/sqrt(length(pct_diff))) %>%
  ungroup()

FKphen2020_june16_std <- FKphen_clean2 %>%
  filter(date == "2020-06-16") %>%
  group_by(site, year, drought) %>%
  summarise(avg = mean(pct_diff), se = sd(pct_diff)/sqrt(length(pct_diff))) %>%
  ungroup()

FKphen2021_june16_std <- FKphen_clean2 %>%
  filter(date == "2021-06-16") %>%
  group_by(site, year, drought) %>%
  summarise(avg = mean(pct_diff), se = sd(pct_diff)/sqrt(length(pct_diff))) %>%
  ungroup()


  #join datasets
FK1pt_std <- FKphen2019_june24_std %>%
  full_join(FKphen2020_june16_std) %>%
  full_join(FKphen2021_june16_std) %>%
  mutate(year = as.factor(year))

#TB
TBphen2019_july8_std <- TBphen_edit_combo2.1 %>%
  filter(date == "2019-07-08") %>%
  group_by(year, drought) %>%
  summarise(avg = mean(green_diff), se = sd(green_diff)/sqrt(length(green_diff))) %>%
  ungroup()

TBphen2020_june18_std <- TBphen_edit_combo2.1 %>%
  filter(date == "2020-06-18") %>%
  group_by(year, drought) %>%
  summarise(avg = mean(green_diff), se = sd(green_diff)/sqrt(length(green_diff))) %>%
  ungroup()

TBphen2021_june16_std <- TBphen_edit_combo2.1 %>%
  filter(date == "2021-06-16") %>%
  group_by(year, drought) %>%
  summarise(avg = mean(green_diff), se = sd(green_diff)/sqrt(length(green_diff))) %>%
  ungroup()

TB1pt_std <- TBphen2019_july8_std %>%
  full_join(TBphen2020_june18_std) %>%
  full_join(TBphen2021_june16_std) %>%
  mutate(year = as.factor(year))


##############################################################################
######### Phenology single time point - Graphs ############

#standardized fig
FKdr_1pt_std <- ggplot(data = FK1pt_std, aes(x = drought, y = avg, color = year, shape = year)) +
  geom_point(size = 3, position = position_dodge(0.5)) + 
  ylim(-40, 40) +
  geom_smooth(data = subset(FK1pt_std, year != "2020"), aes(group = year), method = "lm", se = FALSE) +
  geom_smooth(data = subset(FK1pt_std, year == "2020"), aes(group = year), method = "lm", se = FALSE, linetype = "dashed") +
  geom_errorbar(aes(ymin = avg - se, ymax = avg + se), width = 7, position = position_dodge(.5)) +
  scale_colour_manual(values = cbPalette) +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16)) +
  labs(x = "Precipitation Reduction (%)", y = "Greenness (% Difference)", color = "Year", shape = "Year", title = "MT") 

TBdr_1pt_std <- ggplot(data = TB1pt_std, aes(x = drought, y = avg, color = year, shape = year)) +
  geom_point(size = 3, position = position_dodge(0.5)) + 
  ylim(-50, 50) +
  geom_smooth(data = subset(TB1pt_std, year == "2019"), aes(group = year), method = "lm", se = FALSE) +
  geom_errorbar(aes(ymin = avg - se, ymax = avg + se), width = 7, position = position_dodge(.5)) +
  scale_colour_manual(values = cbPalette) +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16)) +
  labs(x = "Precipitation Reduction (%)", y = "Greenness (% Difference)", color = "Year", shape = "Year", title = "WY") 

#FKstd1 900x700
FKdr_1pt_std / TBdr_1pt_std



##############################################################################
######### Phenology single time point - Statistics ############


##standardized stats
#FK
#FK 2019
FKphen2019_june24_std2 <- FKphen_clean2 %>%
  filter(date == "2019-06-24")

FKphen2019_june24_std2stats <- lmerTest::lmer(data = FKphen2019_june24_std2, pct_diff ~ drought + 
                                              (1|block) + (1|block:paddock))

anova(FKphen2019_june24_std2stats, type = 3) #drought p = 3.773e-5
rsquared(FKphen2019_june24_std2stats) #marg = 0.2895493

resid7.5 <- lm(data = FKphen2019_june24_std2, pct_diff ~ drought)
ols_plot_resid_hist(resid7.5) #fairly normal
ols_test_normality(resid7.5) #passed normality tests 


#FK 2020
FKphen2020_june16_std2 <- FKphen_clean2 %>%
  filter(date == "2020-06-16")

FKphen2020_june16_std3 <- FKphen2020_june16_std2 %>%
  mutate(livestock_util_2019 = as.factor(livestock_util_2019))

FKphen2020_june16_std2stats <- lmerTest::lmer(data = FKphen2020_june16_std3, pct_diff ~ drought*livestock_util_2019 + 
                                              (1|block) + (1|block:paddock))

anova(FKphen2020_june16_std2stats, type = 3) #drought p = 0.06321
rsquared(FKphen2020_june16_std2stats) #marg = 0.05214125

resid7.5 <- lm(data = FKphen2020_june16_std3, pct_diff ~ drought*livestock_util_2019)
ols_plot_resid_hist(resid7.5) #left skew
ols_test_normality(resid7.5) #failed normality tests 

#try trans
FKphen2020_june16_std3_trans <- FKphen2020_june16_std3 %>%
  #mutate(pct_sqrt = sqrt(pct_diff)) %>%
  mutate(pct_sq = pct_diff^2) %>%
  #mutate(pct_ln = log(pct_diff + 10)) %>%
  #mutate(pct_log = log10(pct_diff + 10)) %>%
  mutate(pct_cb = pct_diff^1/3) %>%
  mutate(pct_qd = pct_diff^1/4) %>%
  mutate(pct_cr = pct_diff^0.75)

resid7.5 <- lm(data = FKphen2020_june16_std3_trans, pct_sq ~ drought*livestock_util_2019)  #cb and sq look equally bad ; stick with untransformed
ols_plot_resid_hist(resid7.5) #left skew
ols_test_normality(resid7.5) #failed normality tests 




#FK 2021
FKphen2021_june16_std2 <- FKphen_clean2 %>%
  filter(date == "2021-06-16")

FKphen2021_june16_std3 <- FKphen2021_june16_std2 %>%
  mutate(graze2021 = ifelse(grazing_category == "HHMMM", "HH", 
                    ifelse(grazing_category == "MLLMM", "ML", 
                    ifelse(grazing_category == "MMMMM", "MM", grazing_category)))) %>%
  mutate(graze2021 = as.factor(graze2021))

FKphen2021_june16_std2stats <- lmerTest::lmer(data = FKphen2021_june16_std3, pct_diff ~ drought*graze2021 + 
                                              (1|block) + (1|block:paddock))

anova(FKphen2021_june16_std2stats, type = 3) #drought p = 0.0002044
rsquared(FKphen2021_june16_std2stats) #marg = 0.3110462

resid7.5 <- lm(data = FKphen2021_june16_std3, pct_diff ~ drought*graze2021)
ols_plot_resid_hist(resid7.5) #left skew
ols_test_normality(resid7.5) #failed normality tests 


#TB
#TB 2019
TBphen2019_july8_std2 <- TBphen_edit_combo2 %>%
  filter(date == "2019-07-08")

TBphen2019_july8_std2stats <- lmerTest::lmer(data = TBphen2019_july8_std2, green_diff ~ drought + 
                                                (1|block) + (1|block:paddock))

anova(TBphen2019_july8_std2stats, type = 3) #drought p = 0.03581
rsquared(TBphen2019_july8_std2stats) #marg = 0.1119971

resid7.5 <- lm(data = TBphen2019_july8_std2, green_diff ~ drought)
ols_plot_resid_hist(resid7.5) #fairly normal
ols_test_normality(resid7.5) #passed normality tests 


#TB 2020
TBphen2020_june18_std2 <- TBphen_edit_combo2 %>%
  filter(date == "2020-06-18")

TBphen2020_june18_std3 <- TBphen2020_june18_std2 %>%
  mutate(livestock_util_2019 = as.factor(livestock_util_2019))

TBphen2020_june18_std2stats <- lmerTest::lmer(data = TBphen2020_june18_std3, green_diff ~ drought*livestock_util_2019 + 
                                                (1|block) + (1|block:paddock))

anova(TBphen2020_june18_std2stats, type = 3) #graze p = 0.02831, drought*graze p = 0.04709

resid7.5 <- lm(data = TBphen2020_june18_std3, green_diff ~ drought*livestock_util_2019)
ols_plot_resid_hist(resid7.5) #normal
ols_test_normality(resid7.5) #pass normality tests 


#TB 2021
TBphen2021_june16_std2 <- TBphen_edit_combo2 %>%
  filter(date == "2021-06-16")

TBphen2021_june16_std3 <- TBphen2021_june16_std2 %>%
  mutate(graze2021 = ifelse(grazing_category == "HHMMM", "HH", 
                    ifelse(grazing_category == "MLLMM", "ML", 
                    ifelse(grazing_category == "MMMMM", "MM", grazing_category)))) %>%
  mutate(graze2021 = as.factor(graze2021))

TBphen2021_june16_std2stats <- lmerTest::lmer(data = TBphen2021_june16_std3, green_diff ~ drought*graze2021 + 
                                                (1|block) + (1|block:paddock))

anova(TBphen2021_june16_std2stats, type = 3) #no sig

resid7.5 <- lm(data = TBphen2021_june16_std3, green_diff ~ drought*graze2021)
ols_plot_resid_hist(resid7.5) #normal
ols_test_normality(resid7.5) #pass normality tests 

########################################################################





##########################################################################
###### Cleaning Sp Comp Data ##### 

#Read in this one

write.csv(bromerelcov, file = "bromerelcov.csv", row.names = FALSE)

#dropping plots
bromerelcov_dropped <- bromerelcov %>%   #dropping plots that have 0's for BRAR and BRTE for all 4 years
  filter(!(site == "TB" & plot %in% c("26", "28"))) %>% #here, both 26 and 28 had neither BRAR or BRTE, so dropping those plots totally
  filter(!(site == "TB" & symbol == "BRTE" & plot %in% c("20", "21", "23", "27", "30", "36", "45", "48")))



#######################################################################
########## Absolute sp comp - data manip and graphs ###########

#calculate with means and se 
cov_avg_drought <- bromerelcov_dropped %>%
  group_by(site, year, symbol, drought) %>%
  summarise(avg_abs = mean(aerial_cover), se_abs = sd(aerial_cover)/sqrt(length(aerial_cover)), avg_rel = mean(rel_brome_cov), se_rel = sd(rel_brome_cov)/sqrt(length(rel_brome_cov))) %>%
  ungroup()

cov_avg_graze <- bromerelcov_dropped %>%
  group_by(site, year, symbol, grazing_treatment) %>%
  summarise(avg_abs = mean(aerial_cover), se_abs = sd(aerial_cover)/sqrt(length(aerial_cover)), avg_rel = mean(rel_brome_cov), se_rel = sd(rel_brome_cov)/sqrt(length(rel_brome_cov))) %>%
  ungroup()


#separate graphs by site and species - drought
FKBRAR_avgdr <- cov_avg_drought %>%
  filter(site == "FK") %>%
  filter(symbol == "BRAR")

FKBRTE_avgdr <- cov_avg_drought %>%
  filter(site == "FK") %>%
  filter(symbol == "BRTE")

TBBRAR_avgdr <- cov_avg_drought %>%
  filter(site == "TB") %>%
  filter(symbol == "BRAR")

TBBRTE_avgdr <- cov_avg_drought %>%
  filter(site == "TB") %>%
  filter(symbol == "BRTE")

#Figures
#remake figs 2/18/22
FKBRAR_avgdr2 <- FKBRAR_avgdr %>%
  mutate(year = as.factor(year))
FKBRTE_avgdr2 <- FKBRTE_avgdr %>%
  mutate(year = as.factor(year))
TBBRAR_avgdr2 <- TBBRAR_avgdr %>%
  mutate(year = as.factor(year))
TBBRTE_avgdr2 <- TBBRTE_avgdr %>%
  mutate(year = as.factor(year))

FKBRAR_scatter_dr_abs2 <- ggplot(data = subset(FKBRAR_avgdr2, year != "2018"), aes(x = drought, y = avg_abs, color = year, shape = year)) +
  geom_point(size = 3, position = position_dodge(1)) + 
  scale_colour_manual(values = cbPalette) +
  geom_errorbar(aes(ymin = avg_abs - se_abs, ymax = avg_abs + se_abs), width = 7, position = position_dodge(1)) +
  labs(x = "Precipitation Reduction (%)", y = "Percent Cover (%)", title = "MT BRAR", color = "Year", shape = "Year") +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16))

FKBRTE_scatter_dr_abs2 <- ggplot(data = subset(FKBRTE_avgdr2, year != "2018"), aes(x = drought, y = avg_abs, color = year, shape = year)) +
  geom_point(size = 3, position = position_dodge(1)) + 
  ylim(-0.1, 20) +
  scale_colour_manual(values = cbPalette) +
  geom_errorbar(aes(ymin = avg_abs - se_abs, ymax = avg_abs + se_abs), width = 7, position = position_dodge(1)) +
  labs(x = "Precipitation Reduction (%)", y = "Percent Cover (%)", title = "MT BRTE", color = "Year", shape = "Year") +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16))

TBBRAR_scatter_dr_abs2 <- ggplot(data = subset(TBBRAR_avgdr2, year != "2018"), aes(x = drought, y = avg_abs, color = year, shape = year)) +
  geom_point(size = 3, position = position_dodge(1)) + 
  scale_colour_manual(values = cbPalette) +
  geom_errorbar(aes(ymin = avg_abs - se_abs, ymax = avg_abs + se_abs), width = 7, position = position_dodge(1)) +
  labs(x = "Precipitation Reduction (%)", y = "Percent Cover (%)", title = "WY BRAR", color = "Year", shape = "Year") +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16))

TBBRTE_scatter_dr_abs2 <- ggplot(data = subset(TBBRTE_avgdr2, year != "2018"), aes(x = drought, y = avg_abs, color = year, shape = year)) +
  geom_point(size = 3, position = position_dodge(1)) + 
  ylim(-0.1, 15) +
  scale_colour_manual(values = cbPalette) +
  geom_errorbar(aes(ymin = avg_abs - se_abs, ymax = avg_abs + se_abs), width = 7, position = position_dodge(1)) +
  labs(x = "Precipitation Reduction (%)", y = "Percent Cover (%)", title = "WY BRTE", color = "Year", shape = "Year") +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16))

(FKBRAR_scatter_dr_abs2 + FKBRTE_scatter_dr_abs2) / (TBBRAR_scatter_dr_abs2 + TBBRTE_scatter_dr_abs2)


#separate graphs by site and species - grazing

#remake figs 2/17/22
#need to regroup and take out 2018-2019
FKBRAR_covgr <- bromerelcov_dropped %>%
  mutate(year = as.factor(year)) %>%
  filter(site == "FK" & symbol == "BRAR") %>%
  filter(year %in% c("2020", "2021")) %>%
  mutate(realg = ifelse(year == "2020" & livestock_util_2019 == "50", "Stable", 
                ifelse(year == "2020" & livestock_util_2019 == "70", "Heavy", 
                ifelse(year == "2021" & livestock_util_2020 == "30", "Destock", 
                ifelse(year == "2021" & livestock_util_2020 == "50", "Stable",
                ifelse(year == "2021" & livestock_util_2020 == "70", "Heavy", "0")))))) %>%
  group_by(year, realg) %>% #calculating rel and abs at same time
  summarise(avg_rel = mean(rel_brome_cov), se_rel = sd(rel_brome_cov)/sqrt(length(rel_brome_cov)), avg_abs = mean(aerial_cover), se_abs = sd(aerial_cover)/sqrt(length(aerial_cover))) %>%
  ungroup()

FKBRTE_covgr <- bromerelcov_dropped %>%
  mutate(year = as.factor(year)) %>%
  filter(site == "FK" & symbol == "BRTE") %>%
  filter(year %in% c("2020", "2021")) %>%
  mutate(realg = ifelse(year == "2020" & livestock_util_2019 == "50", "Stable", 
                ifelse(year == "2020" & livestock_util_2019 == "70", "Heavy", 
                ifelse(year == "2021" & livestock_util_2020 == "30", "Destock", 
                ifelse(year == "2021" & livestock_util_2020 == "50", "Stable",
                ifelse(year == "2021" & livestock_util_2020 == "70", "Heavy", "0")))))) %>%
  group_by(year, realg) %>% #calculating rel and abs at same time
  summarise(avg_rel = mean(rel_brome_cov), se_rel = sd(rel_brome_cov)/sqrt(length(rel_brome_cov)), avg_abs = mean(aerial_cover), se_abs = sd(aerial_cover)/sqrt(length(aerial_cover))) %>%
  ungroup()

TBBRAR_covgr <- bromerelcov_dropped %>%
  mutate(year = as.factor(year)) %>%
  filter(site == "TB" & symbol == "BRAR") %>%
  filter(year %in% c("2020", "2021")) %>%
  mutate(realg = ifelse(year == "2020" & livestock_util_2019 == "50", "Stable", 
                 ifelse(year == "2020" & livestock_util_2019 == "70", "Heavy", 
                 ifelse(year == "2021" & livestock_util_2020 == "30", "Destock", 
                 ifelse(year == "2021" & livestock_util_2020 == "50", "Stable",
                 ifelse(year == "2021" & livestock_util_2020 == "70", "Heavy", "0")))))) %>%
  group_by(year, realg) %>% #calculating rel and abs at same time
  summarise(avg_rel = mean(rel_brome_cov), se_rel = sd(rel_brome_cov)/sqrt(length(rel_brome_cov)), avg_abs = mean(aerial_cover), se_abs = sd(aerial_cover)/sqrt(length(aerial_cover))) %>%
  ungroup()

TBBRTE_covgr <- bromerelcov_dropped %>%
  mutate(year = as.factor(year)) %>%
  filter(site == "TB" & symbol == "BRTE") %>%
  filter(year %in% c("2020", "2021")) %>%
  mutate(realg = ifelse(year == "2020" & livestock_util_2019 == "50", "Stable", 
                ifelse(year == "2020" & livestock_util_2019 == "70", "Heavy", 
                ifelse(year == "2021" & livestock_util_2020 == "30", "Destock", 
                ifelse(year == "2021" & livestock_util_2020 == "50", "Stable",
                ifelse(year == "2021" & livestock_util_2020 == "70", "Heavy", "0")))))) %>%
  group_by(year, realg) %>% #calculating rel and abs at same time
  summarise(avg_rel = mean(rel_brome_cov), se_rel = sd(rel_brome_cov)/sqrt(length(rel_brome_cov)), avg_abs = mean(aerial_cover), se_abs = sd(aerial_cover)/sqrt(length(aerial_cover))) %>%
  ungroup()


#change order of x axis
FKBRAR_covgr$realg <- factor(FKBRAR_covgr$realg, levels = c("Destock", "Stable", "Heavy")) 
FKBRTE_covgr$realg <- factor(FKBRTE_covgr$realg, levels = c("Destock", "Stable", "Heavy")) 
TBBRAR_covgr$realg <- factor(TBBRAR_covgr$realg, levels = c("Destock", "Stable", "Heavy")) 
TBBRTE_covgr$realg <- factor(TBBRTE_covgr$realg, levels = c("Destock", "Stable", "Heavy")) 

FKBRAR_covgr3 <- ggplot(data = FKBRAR_covgr, aes(x = realg, y = avg_abs, color = year, shape = year)) +
  geom_point(size = 3, position = position_dodge(.5)) +
  geom_errorbar(aes(ymin = avg_abs - se_abs, ymax = avg_abs + se_abs), width = .2, position = position_dodge(.5)) +
  scale_colour_manual(values = cbPalette) +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16)) +
  labs(x = "Grazing Treatment", y = "Percent Cover (%)", color = "Year", shape = "Year", title = "MT BRAR") 
FKBRTE_covgr3 <- ggplot(data = FKBRTE_covgr, aes(x = realg, y = avg_abs, color = year, shape = year)) +
  geom_point(size = 3, position = position_dodge(.5)) +
  geom_errorbar(aes(ymin = avg_abs - se_abs, ymax = avg_abs + se_abs), width = .2, position = position_dodge(.5)) +
  scale_colour_manual(values = cbPalette) +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16)) +
  labs(x = "Grazing Treatment", y = "Percent Cover (%)", color = "Year", shape = "Year", title = "MT BRTE") 
TBBRAR_covgr3 <- ggplot(data = TBBRAR_covgr, aes(x = realg, y = avg_abs, color = year, shape = year)) +
  geom_point(size = 3, position = position_dodge(.5)) +
  geom_errorbar(aes(ymin = avg_abs - se_abs, ymax = avg_abs + se_abs), width = .2, position = position_dodge(.5)) +
  scale_colour_manual(values = cbPalette) +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16)) +
  labs(x = "Grazing Treatment", y = "Percent Cover (%)", color = "Year", shape = "Year", title = "WY BRAR") 
TBBRTE_covgr3 <- ggplot(data = TBBRTE_covgr, aes(x = realg, y = avg_abs, color = year, shape = year)) +
  geom_point(size = 3, position = position_dodge(.5)) +
  geom_errorbar(aes(ymin = avg_abs - se_abs, ymax = avg_abs + se_abs), width = .2, position = position_dodge(.5)) +
  scale_colour_manual(values = cbPalette) +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16)) +
  labs(x = "Grazing Treatment", y = "Percent Cover (%)", color = "Year", shape = "Year", title = "WY BRTE") 

(FKBRAR_covgr3 + FKBRTE_covgr3) / (TBBRAR_covgr3 + TBBRTE_covgr3)






#######################################################################
########## Absolute sp comp - Stats ###########

#statistics
#Test for normality of residuals
FKBRAR_covdrop <- bromerelcov_dropped %>%
  filter(site == "FK" & symbol == "BRAR")

FKBRTE_covdrop <- bromerelcov_dropped %>%
  filter(site == "FK" & symbol == "BRTE")

TBBRAR_covdrop <- bromerelcov_dropped %>%
  filter(site == "TB" & symbol == "BRAR")

TBBRTE_covdrop <- bromerelcov_dropped %>%
  filter(site == "TB" & symbol == "BRTE")

#FK BRAR 2018
resid4.0 <- lm(data = subset(FKBRAR_covdrop, year == 2018), aerial_cover ~ drought*grazing_treatment)
ols_plot_resid_hist(resid4.0) #normalish
ols_test_normality(resid4.0) #passed normality tests

#FK BRAR 2019
resid4.1 <- lm(data = subset(FKBRAR_covdrop, year == 2019), aerial_cover ~ drought*grazing_treatment)
ols_plot_resid_hist(resid4.1) #normalish
ols_test_normality(resid4.1) #passed normality tests

#FK BRAR 2020
resid4.2 <- lm(data = subset(FKBRAR_covdrop, year == 2020), aerial_cover ~ drought*grazing_treatment)
ols_plot_resid_hist(resid4.2) #normalish
ols_test_normality(resid4.2) #passed normality tests

#FK BRAR 2021
resid4.3 <- lm(data = subset(FKBRAR_covdrop, year == 2021), aerial_cover ~ drought*grazing_treatment)
ols_plot_resid_hist(resid4.3) #right skew
ols_test_normality(resid4.3) #failed normality tests


#FK BRTE 2018
resid4.4 <- lm(data = subset(FKBRTE_covdrop, year == 2018), aerial_cover ~ drought*grazing_treatment)
ols_plot_resid_hist(resid4.4) #right skew
ols_test_normality(resid4.4) #failed normality tests

#FK BRTE 2019
resid4.5 <- lm(data = subset(FKBRTE_covdrop, year == 2019), aerial_cover ~ drought*grazing_treatment)
ols_plot_resid_hist(resid4.5) #right skew
ols_test_normality(resid4.5) #failed normality tests

#FK BRTE 2020
resid4.6 <- lm(data = subset(FKBRTE_covdrop, year == 2020), aerial_cover ~ drought*grazing_treatment)
ols_plot_resid_hist(resid4.6) #right skew
ols_test_normality(resid4.6) #failed normality tests

#FK BRTE 2021
resid4.7 <- lm(data = subset(FKBRTE_covdrop, year == 2021), aerial_cover ~ drought*grazing_treatment)
ols_plot_resid_hist(resid4.7) #right skew
ols_test_normality(resid4.7) #failed normality tests


#TB BRAR 2018
resid4.8 <- lm(data = subset(TBBRAR_covdrop, year == 2018), aerial_cover ~ drought*grazing_treatment)
ols_plot_resid_hist(resid4.8) #right skew
ols_test_normality(resid4.8) #failed normality tests

#TB BRAR 2019
resid4.9 <- lm(data = subset(TBBRAR_covdrop, year == 2019), aerial_cover ~ drought*grazing_treatment)
ols_plot_resid_hist(resid4.9) #right skew
ols_test_normality(resid4.9) #failed normality tests

#TB BRAR 2020
resid4.10 <- lm(data = subset(TBBRAR_covdrop, year == 2020), aerial_cover ~ drought*grazing_treatment)
ols_plot_resid_hist(resid4.10) #right skew
ols_test_normality(resid4.10) #failed normality tests

#TB BRAR 2021
resid4.11 <- lm(data = subset(TBBRAR_covdrop, year == 2021), aerial_cover ~ drought*grazing_treatment)
ols_plot_resid_hist(resid4.11) #right skew
ols_test_normality(resid4.11) #failed normality tests


#TB BRTE 2018
resid4.12 <- lm(data = subset(TBBRTE_covdrop, year == 2018), aerial_cover ~ drought*grazing_treatment)
ols_plot_resid_hist(resid4.12) #peak in center
ols_test_normality(resid4.12) #failed normality tests

#TB BRTE 2019
resid4.13 <- lm(data = subset(TBBRTE_covdrop, year == 2019), aerial_cover ~ drought*grazing_treatment)
ols_plot_resid_hist(resid4.13) #right skew
ols_test_normality(resid4.13) #failed normality tests

#TB BRTE 2020
resid4.14 <- lm(data = subset(TBBRTE_covdrop, year == 2020), aerial_cover ~ drought*grazing_treatment)
ols_plot_resid_hist(resid4.14) #right skew
ols_test_normality(resid4.14) #failed normality tests

#TB BRTE 2021
resid4.15 <- lm(data = subset(TBBRTE_covdrop, year == 2021), aerial_cover ~ drought*grazing_treatment)
ols_plot_resid_hist(resid4.15) #right skew
ols_test_normality(resid4.15) #failed normality tests


#FK BRAR 2018 - no trans needed
FKBRAR_covdrop2018 <- FKBRAR_covdrop %>%
  filter(year == 2018)

FKBRAR_abs_2018_stats <- lmerTest::lmer(data = FKBRAR_covdrop2018, aerial_cover ~ drought*grazing_treatment + 
                                          (1|block) + (1|block:paddock)) 
anova(FKBRAR_abs_2018_stats, type = 3) #nothing significant


#FK BRAR 2019 - no trans needed
FKBRAR_covdrop2019 <- FKBRAR_covdrop %>%
  filter(year == 2019)

FKBRAR_abs_2019_stats <- lmerTest::lmer(data = FKBRAR_covdrop2019, aerial_cover ~ drought*grazing_treatment + 
                                          (1|block) + (1|block:paddock))
anova(FKBRAR_abs_2019_stats, type = 3) #nothing significant


#revisit stats - don't include grazing, because ANPP was collected before grazing occurred; try 2018 data as covariate
FKBRAR_abs_2019_stats2 <- lmerTest::lmer(data = FKBRAR_covdrop2019, aerial_cover ~ drought +
                                           (1|block) + (1|block:paddock))

anova(FKBRAR_abs_2019_stats2, type = 3) #no sig
AIC(FKBRAR_abs_2019_stats2) #AIC = 449.4569
AICc(FKBRAR_abs_2019_stats2) #AICc = 450.7069

#2018 covariate
FKBRAR_covdrop2018_2 <- FKBRAR_covdrop2018 %>%
  dplyr:: select(c(plot, aerial_cover)) %>%
  mutate(abscov2018 = aerial_cover) %>%
  dplyr:: select(-aerial_cover)

FKBRAR_covdrop20192 <- FKBRAR_covdrop2019 %>%
  left_join(FKBRAR_covdrop2018_2)


FKBRAR_abs_2019_stats3 <- lmerTest::lmer(data = FKBRAR_covdrop20192, aerial_cover ~ drought + abscov2018 +
                                           (1|block) + (1|block:paddock))

anova(FKBRAR_abs_2019_stats3, type = 3) #abscov2018 p = 3.092e-6, no other sig
AIC(FKBRAR_abs_2019_stats3) #AIC = 433.7783
AICc(FKBRAR_abs_2019_stats3) #AICc = 435.5655




#FK BRAR 2020 - no trans needed
FKBRAR_covdrop2020 <- FKBRAR_covdrop %>%
  filter(year == 2020)

FKBRAR_abs_2020_stats <- lmerTest::lmer(data = FKBRAR_covdrop2020, aerial_cover ~ drought*grazing_treatment + 
                                          (1|block) + (1|block:paddock))
anova(FKBRAR_abs_2020_stats, type = 3) #nothing significant


#revisit stats - include actual grazing treatments from 2019, because ANPP was collected before grazing occurred; then try 2018 data as covariate
FKBRAR_covdrop2020_2 <- FKBRAR_covdrop2020 %>%
  left_join(FKBRAR_covdrop2018_2) %>%
  mutate(livestock_util_2019 = as.factor(livestock_util_2019))

FKBRAR_abs_2020_stats2 <- lmerTest::lmer(data = FKBRAR_covdrop2020_2, aerial_cover ~ drought*livestock_util_2019 +
                                           (1|block) + (1|block:paddock))

anova(FKBRAR_abs_2020_stats2, type = 3) #no sig
AIC(FKBRAR_abs_2020_stats2) #AIC = 438.8341
AICc(FKBRAR_abs_2020_stats2) #AICc = 441.2688

#d+g w graze as factor
FKBRAR_abs_2020_stats2b <- lmerTest::lmer(data = FKBRAR_covdrop2020_2, aerial_cover ~ drought + livestock_util_2019 +
                                            (1|block) + (1|block:paddock))

anova(FKBRAR_abs_2020_stats2b, type = 3) #no sig
AIC(FKBRAR_abs_2020_stats2b) #AIC = 434.6809
AICc(FKBRAR_abs_2020_stats2b) #AICc = 436.4681


#covariate
FKBRAR_abs_2020_stats3 <- lmerTest::lmer(data = FKBRAR_covdrop2020_2, aerial_cover ~ drought*livestock_util_2019 + abscov2018 +
                                           (1|block) + (1|block:paddock))

anova(FKBRAR_abs_2020_stats3, type = 3) #relANPP2018 p = 2.211e-5, no other sig
AIC(FKBRAR_abs_2020_stats3) #AIC = 429.3751
AICc(FKBRAR_abs_2020_stats3) #AICc = 432.5751

#random effect of 2018
FKBRAR_abs_2020_stats3b <- lmerTest::lmer(data = FKBRAR_covdrop2020_2, aerial_cover ~ drought*livestock_util_2019 + 
                                            (1|abscov2018) + (1|block) + (1|block:paddock))

anova(FKBRAR_abs_2020_stats3b, type = 3) #no sig
AIC(FKBRAR_abs_2020_stats3b) #AIC = 438.9308
AICc(FKBRAR_abs_2020_stats3b) #AICc = 442.1308

#numeric graze, d*g
FKBRAR_covdrop2020$livestock_util_2019 <- as.numeric(FKBRAR_covdrop2020$livestock_util_2019)
FKBRAR_abs_2020_stats4 <- lmerTest::lmer(data = FKBRAR_covdrop2020, aerial_cover ~ drought*livestock_util_2019 +
                                           (1|block) + (1|block:paddock))

anova(FKBRAR_abs_2020_stats4, type = 3) #no sig
AIC(FKBRAR_abs_2020_stats4) #AIC = 450.817
AICc(FKBRAR_abs_2020_stats4) #AICc = 453.2518

#numeric graze, d + g
FKBRAR_abs_2020_stats5 <- lmerTest::lmer(data = FKBRAR_covdrop2020, aerial_cover ~ drought + livestock_util_2019 +
                                           (1|block) + (1|block:paddock))

anova(FKBRAR_abs_2020_stats5, type = 3) #no sig
AIC(FKBRAR_abs_2020_stats5) #AIC = 440.6724
AICc(FKBRAR_abs_2020_stats5) #AICc = 442.4596

#recode grazing, d*g
FKBRAR_covdrop2020b <- FKBRAR_covdrop2020 %>%
  mutate(g2019 = ifelse(livestock_util_2019 == "50", "-1", 
                        ifelse(livestock_util_2019 == "70", "1", livestock_util_2019))) %>%
  mutate(g2019 = as.numeric(g2019))
FKBRAR_abs_2020_stats6 <- lmerTest::lmer(data = FKBRAR_covdrop2020b, aerial_cover ~ drought*g2019 +
                                           (1|block) + (1|block:paddock))

anova(FKBRAR_abs_2020_stats6, type = 3) #no sig
AIC(FKBRAR_abs_2020_stats6) #AIC = 441.6066
AICc(FKBRAR_abs_2020_stats6) #AICc = 444.0414

#recode graze, d+g
FKBRAR_abs_2020_stats7 <- lmerTest::lmer(data = FKBRAR_covdrop2020b, aerial_cover ~ drought + g2019 +
                                           (1|block) + (1|block:paddock))

anova(FKBRAR_abs_2020_stats7, type = 3) #no sig
AIC(FKBRAR_abs_2020_stats7) #AIC = 436.0672
AICc(FKBRAR_abs_2020_stats7) #AICc = 437.8544

#just drought
FKBRAR_abs_2020_stats8 <- lmerTest::lmer(data = FKBRAR_covdrop2020b, aerial_cover ~ drought +
                                           (1|block) + (1|block:paddock))

anova(FKBRAR_abs_2020_stats8, type = 3) #no sig
AIC(FKBRAR_abs_2020_stats8) #AIC = 439.6363
AICc(FKBRAR_abs_2020_stats8) #AICc = 440.8863






#FK BRAR 2021 - try transformations
FKBRAR_covdrop2021 <- FKBRAR_covdrop %>%
  filter(year == 2021)

FKBRAR_covdrop2021_trans <- FKBRAR_covdrop2021 %>%
  mutate(cov_sqrt = sqrt(aerial_cover)) %>%
  mutate(cov_sq = aerial_cover^2) %>%
  mutate(cov_cb = aerial_cover^1/3) %>%
  mutate(cov_qd = aerial_cover^1/4) %>%
  mutate(cov_ln = log(aerial_cover + 0.1)) %>%
  mutate(cov_log = log10(aerial_cover + 0.1)) %>%
  mutate(cov_cr = aerial_cover^0.75)

resid4.18 <- lm(data = FKBRAR_covdrop2021_trans, cov_ln ~ drought*grazing_treatment)
ols_plot_resid_hist(resid4.18) #normal
ols_test_normality(resid4.18) #pass normality tests except cramer

#stats with ln
FKBRAR_abs_2021_stats <- lmerTest::lmer(data = FKBRAR_covdrop2021_trans, cov_ln ~ drought*grazing_treatment +
                                          (1|block) + (1|block:paddock))
anova(FKBRAR_abs_2021_stats, type = 3) #nothing significant


#revisit stats - include actual grazing treatments from 2019 + 2020, because ANPP was collected before grazing occurred; then try 2018 data as covariate
FKBRAR_covdrop2021_trans2 <- FKBRAR_covdrop2021_trans %>%
  left_join(FKBRAR_covdrop2018_2) %>%
  mutate(graze2021 = ifelse(grazing_category == "HHMMM", "HH", 
                            ifelse(grazing_category == "MLLMM", "ML", 
                                   ifelse(grazing_category == "MMMMM", "MM", grazing_category)))) %>%
  mutate(graze2021 = as.factor(graze2021))


FKBRAR_abs_2021_stats2 <- lmerTest::lmer(data = FKBRAR_covdrop2021_trans2, cov_ln ~ drought*graze2021 +
                                           (1|block) + (1|block:paddock))

anova(FKBRAR_abs_2021_stats2, type = 3) #no sig
AIC(FKBRAR_abs_2021_stats2) #AIC = 210.5236
AICc(FKBRAR_abs_2021_stats2) #AICc = 214.6146

#graze factor, d+g
FKBRAR_abs_2021_stats2b <- lmerTest::lmer(data = FKBRAR_covdrop2021_trans2, cov_ln ~ drought + graze2021 +
                                            (1|block) + (1|block:paddock))

anova(FKBRAR_abs_2021_stats2b, type = 3) #no sig
AIC(FKBRAR_abs_2021_stats2b) #AIC = 193.0515
AICc(FKBRAR_abs_2021_stats2b) #AICc = 195.4862


#covariate
FKBRAR_abs_2021_stats3 <- lmerTest::lmer(data = FKBRAR_covdrop2021_trans2, cov_ln ~ drought*graze2021 + abscov2018 +
                                           (1|block) + (1|block:paddock))

anova(FKBRAR_abs_2021_stats3, type = 3) #abscov2018 p = 0.000411, no other sig
AIC(FKBRAR_abs_2021_stats3) #AIC = 207.1463
AICc(FKBRAR_abs_2021_stats3) #AICc = 212.2625

#2018 as random
FKBRAR_abs_2021_stats3b <- lmerTest::lmer(data = FKBRAR_covdrop2021_trans2, cov_ln ~ drought*graze2021 + 
                                            (1|abscov2018) + (1|block) + (1|block:paddock))

anova(FKBRAR_abs_2021_stats3b, type = 3) #graze p = 0.07725
AIC(FKBRAR_abs_2021_stats3b) #AIC = 210.6212
AICc(FKBRAR_abs_2021_stats3b) #AICc = 215.7375

#numeric graze, d*g
FKBRAR_covdrop2021_trans2$livestock_util_2020 <- as.numeric(FKBRAR_covdrop2021_trans2$livestock_util_2020)
FKBRAR_abs_2021_stats4 <- lmerTest::lmer(data = FKBRAR_covdrop2021_trans2, cov_ln ~ drought*livestock_util_2020 +
                                           (1|block) + (1|block:paddock))

anova(FKBRAR_abs_2021_stats4, type = 3) #no sig
AIC(FKBRAR_abs_2021_stats4) #AIC = 218.575
AICc(FKBRAR_abs_2021_stats4) #AICc = 221.0098

#numeric graze, d + g
FKBRAR_abs_2021_stats5 <- lmerTest::lmer(data = FKBRAR_covdrop2021_trans2, cov_ln ~ drought + livestock_util_2020 +
                                           (1|block) + (1|block:paddock))

anova(FKBRAR_abs_2021_stats5, type = 3) #no sig
AIC(FKBRAR_abs_2021_stats5) #AIC = 202.1662
AICc(FKBRAR_abs_2021_stats5) #AICc = 203.9534

#recode graze, d*g
FKBRAR_covdrop2021_trans3 <- FKBRAR_covdrop2021_trans2 %>%
  mutate(g2020 = ifelse(livestock_util_2020 == "30", "-1", 
                        ifelse(livestock_util_2020 == "50", "0", 
                               ifelse(livestock_util_2020 == "70", "1", livestock_util_2020)))) %>%
  mutate(g2020 = as.numeric(g2020))
FKBRAR_abs_2021_stats6 <- lmerTest::lmer(data = FKBRAR_covdrop2021_trans3, cov_ln ~ drought*g2020 +
                                           (1|block) + (1|block:paddock))

anova(FKBRAR_abs_2021_stats6, type = 3) #no sig
AIC(FKBRAR_abs_2021_stats6) #AIC = 206.5921
AICc(FKBRAR_abs_2021_stats6) #AICc = 209.0269

#recode graze, d+g
FKBRAR_abs_2021_stats7 <- lmerTest::lmer(data = FKBRAR_covdrop2021_trans3, cov_ln ~ drought + g2020 +
                                           (1|block) + (1|block:paddock))

anova(FKBRAR_abs_2021_stats7, type = 3) #no sig
AIC(FKBRAR_abs_2021_stats7) #AIC = 196.1747
AICc(FKBRAR_abs_2021_stats7) #AICc = 197.9619

#just drought
FKBRAR_abs_2021_stats8 <- lmerTest::lmer(data = FKBRAR_covdrop2021_trans3, cov_ln ~ drought +
                                           (1|block) + (1|block:paddock))

anova(FKBRAR_abs_2021_stats8, type = 3) #no sig
AIC(FKBRAR_abs_2021_stats8) #AIC = 193.2239
AICc(FKBRAR_abs_2021_stats8) #AICc = 194.4739







#FK BRTE 2018 - try transformations
FKBRTE_covdrop2018 <- FKBRTE_covdrop %>%
  filter(year == 2018)

#stats with raw data bc no transformations worked
FKBRTE_abs_2018_stats <- lmerTest::lmer(data = FKBRTE_covdrop2018, aerial_cover ~ drought*grazing_treatment + 
                                          (1|block) + (1|block:paddock))
anova(FKBRTE_abs_2018_stats, type = 3) #nothing significant


#FK BRTE 2019 - try transformations
FKBRTE_covdrop2019 <- FKBRTE_covdrop %>%
  filter(year == 2019)

#stats with raw data bc no transformations worked
FKBRTE_abs_2019_stats <- lmerTest::lmer(data = FKBRTE_covdrop2019, aerial_cover ~ drought*grazing_treatment + 
                                          (1|block) + (1|block:paddock))
anova(FKBRTE_abs_2019_stats, type = 3) #nothing significant


#revisit stats - don't include grazing, because ANPP was collected before grazing occurred; try 2018 data as covariate
FKBRTE_abs_2019_stats2 <- lmerTest::lmer(data = FKBRTE_covdrop2019, aerial_cover ~ drought +
                                           (1|block) + (1|block:paddock))

anova(FKBRTE_abs_2019_stats2, type = 3) #no sig
AIC(FKBRTE_abs_2019_stats2) #AIC = 408.0123
AICc(FKBRTE_abs_2019_stats2) #AICc = 409.2623

#2018 covariate
FKBRTE_covdrop2018_2 <- FKBRTE_covdrop2018 %>%
  dplyr:: select(c(plot, aerial_cover)) %>%
  mutate(abscov2018 = aerial_cover) %>%
  dplyr:: select(-aerial_cover)

FKBRTE_covdrop20192 <- FKBRTE_covdrop2019 %>%
  left_join(FKBRTE_covdrop2018_2)


FKBRTE_abs_2019_stats3 <- lmerTest::lmer(data = FKBRTE_covdrop20192, aerial_cover ~ drought + abscov2018 +
                                           (1|block) + (1|block:paddock))

anova(FKBRTE_abs_2019_stats3, type = 3) #abscov2018 p = 2.276e-8, no other sig
AIC(FKBRTE_abs_2019_stats3) #AIC = 382.9139
AICc(FKBRTE_abs_2019_stats3) #AICc = 384.7012




#FK BRTE 2020 - try transformations
FKBRTE_covdrop2020 <- FKBRTE_covdrop %>%
  filter(year == 2020)

#stats with raw data bc no transformations worked
FKBRTE_abs_2020_stats <- lmerTest::lmer(data = FKBRTE_covdrop2020, aerial_cover ~ drought*grazing_treatment + 
                                          (1|block) + (1|block:paddock))
anova(FKBRTE_abs_2020_stats, type = 3) #nothing significant


#revisit stats - include actual grazing treatments from 2019, because ANPP was collected before grazing occurred; then try 2018 data as covariate
FKBRTE_covdrop2020_2 <- FKBRTE_covdrop2020 %>%
  left_join(FKBRTE_covdrop2018_2) %>%
  mutate(livestock_util_2019 = as.factor(livestock_util_2019))

FKBRTE_abs_2020_stats2 <- lmerTest::lmer(data = FKBRTE_covdrop2020_2, aerial_cover ~ drought*livestock_util_2019 +
                                           (1|block) + (1|block:paddock))

anova(FKBRTE_abs_2020_stats2, type = 3) #no sig
AIC(FKBRTE_abs_2020_stats2) #AIC = 445.4724
AICc(FKBRTE_abs_2020_stats2) #AICc = 447.9072

#d+g, factor graze
FKBRTE_abs_2020_stats2b <- lmerTest::lmer(data = FKBRTE_covdrop2020_2, aerial_cover ~ drought + livestock_util_2019 +
                                            (1|block) + (1|block:paddock))

anova(FKBRTE_abs_2020_stats2b, type = 3) #no sig
AIC(FKBRTE_abs_2020_stats2b) #AIC = 440.8222
AICc(FKBRTE_abs_2020_stats2b) #AICc = 442.6095


#covariate
FKBRTE_abs_2020_stats3 <- lmerTest::lmer(data = FKBRTE_covdrop2020_2, aerial_cover ~ drought*livestock_util_2019 + abscov2018 +
                                           (1|block) + (1|block:paddock))

anova(FKBRTE_abs_2020_stats3, type = 3) #relANPP2018 p = 3.113e-12, no other sig
AIC(FKBRTE_abs_2020_stats3) #AIC = 404.2688
AICc(FKBRTE_abs_2020_stats3) #AICc = 407.4688

#2018 as random
FKBRTE_abs_2020_stats3b <- lmerTest::lmer(data = FKBRTE_covdrop2020_2, aerial_cover ~ drought*livestock_util_2019 + 
                                            (1|abscov2018) + (1|block) + (1|block:paddock))

anova(FKBRTE_abs_2020_stats3b, type = 3) #relANPP2018 p = 3.113e-12, no other sig
AIC(FKBRTE_abs_2020_stats3b) #AIC = 385.2542
AICc(FKBRTE_abs_2020_stats3b) #AICc = 388.4542

#numeric g, d*g
FKBRTE_covdrop2020$livestock_util_2019 <- as.numeric(FKBRTE_covdrop2020$livestock_util_2019)
FKBRTE_abs_2020_stats4 <- lmerTest::lmer(data = FKBRTE_covdrop2020, aerial_cover ~ drought*livestock_util_2019 +
                                           (1|block) + (1|block:paddock))

anova(FKBRTE_abs_2020_stats4, type = 3) #no sig
AIC(FKBRTE_abs_2020_stats4) #AIC = 457.4554
AICc(FKBRTE_abs_2020_stats4) #AICc = 459.8901

#numeric g, d + g
FKBRTE_abs_2020_stats5 <- lmerTest::lmer(data = FKBRTE_covdrop2020, aerial_cover ~ drought + livestock_util_2019 +
                                           (1|block) + (1|block:paddock))

anova(FKBRTE_abs_2020_stats5, type = 3) #no sig
AIC(FKBRTE_abs_2020_stats5) #AIC = 446.8137
AICc(FKBRTE_abs_2020_stats5) #AICc = 448.6009

#recode g, d*g
FKBRTE_covdrop2020b <- FKBRTE_covdrop2020 %>%
  mutate(g2019 = ifelse(livestock_util_2019 == "50", "-1", 
                        ifelse(livestock_util_2019 == "70", "1", livestock_util_2019))) %>%
  mutate(g2019 = as.numeric(g2019))
FKBRTE_abs_2020_stats6 <- lmerTest::lmer(data = FKBRTE_covdrop2020b, aerial_cover ~ drought*g2019 +
                                           (1|block) + (1|block:paddock))

anova(FKBRTE_abs_2020_stats6, type = 3) #no sig
AIC(FKBRTE_abs_2020_stats6) #AIC = 448.245
AICc(FKBRTE_abs_2020_stats6) #AICc = 450.6798

#recode g, d+g
FKBRTE_abs_2020_stats7 <- lmerTest::lmer(data = FKBRTE_covdrop2020b, aerial_cover ~ drought + g2019 +
                                           (1|block) + (1|block:paddock))

anova(FKBRTE_abs_2020_stats7, type = 3) #no sig
AIC(FKBRTE_abs_2020_stats7) #AIC = 442.2085
AICc(FKBRTE_abs_2020_stats7) #AICc = 443.9957

#just d
FKBRTE_abs_2020_stats8 <- lmerTest::lmer(data = FKBRTE_covdrop2020b, aerial_cover ~ drought +
                                           (1|block) + (1|block:paddock))

anova(FKBRTE_abs_2020_stats8, type = 3) #no sig
AIC(FKBRTE_abs_2020_stats8) #AIC = 443.3925
AICc(FKBRTE_abs_2020_stats8) #AICc = 444.6425




#FK BRTE 2021 - try transformations
FKBRTE_covdrop2021 <- FKBRTE_covdrop %>%
  filter(year == 2021)

#stats with raw data bc no transformations worked
FKBRTE_abs_2021_stats <- lmerTest::lmer(data = FKBRTE_covdrop2021, aerial_cover ~ drought*grazing_treatment + 
                                          (1|block) + (1|block:paddock))
anova(FKBRTE_abs_2021_stats, type = 3) #nothing significant


#revisit stats - include actual grazing treatments from 2019 + 2020, because ANPP was collected before grazing occurred; then try 2018 data as covariate
FKBRTE_covdrop2021_2 <- FKBRTE_covdrop2021 %>%
  left_join(FKBRTE_covdrop2018_2) %>%
  mutate(graze2021 = ifelse(grazing_category == "HHMMM", "HH", 
                            ifelse(grazing_category == "MLLMM", "ML", 
                                   ifelse(grazing_category == "MMMMM", "MM", grazing_category)))) %>%
  mutate(graze2021 = as.factor(graze2021))


FKBRTE_abs_2021_stats2 <- lmerTest::lmer(data = FKBRTE_covdrop2021_2, aerial_cover ~ drought*graze2021 +
                                           (1|block) + (1|block:paddock))

anova(FKBRTE_abs_2021_stats2, type = 3) #no sig
AIC(FKBRTE_abs_2021_stats2) #AIC = 452.83
AICc(FKBRTE_abs_2021_stats2) #AICc = 456.9209

#g as factor, d+g
FKBRTE_abs_2021_stats2b <- lmerTest::lmer(data = FKBRTE_covdrop2021_2, aerial_cover ~ drought + graze2021 +
                                            (1|block) + (1|block:paddock))

anova(FKBRTE_abs_2021_stats2b, type = 3) #no sig
AIC(FKBRTE_abs_2021_stats2b) #AIC = 444.9165
AICc(FKBRTE_abs_2021_stats2b) #AICc = 447.3513


#covariate
FKBRTE_abs_2021_stats3 <- lmerTest::lmer(data = FKBRTE_covdrop2021_2, aerial_cover ~ drought*graze2021 + abscov2018 +
                                           (1|block) + (1|block:paddock))

anova(FKBRTE_abs_2021_stats3, type = 3) #abscov2018 p = 1.375e-10; graze p = 0.01373
AIC(FKBRTE_abs_2021_stats3) #AIC = 416.7079
AICc(FKBRTE_abs_2021_stats3) #AICc = 421.8242

#2018 as random
FKBRTE_abs_2021_stats3b <- lmerTest::lmer(data = FKBRTE_covdrop2021_2, aerial_cover ~ drought*graze2021 + 
                                            (1|abscov2018) + (1|block) + (1|block:paddock))

anova(FKBRTE_abs_2021_stats3b, type = 3) #graze p = 0.06338
AIC(FKBRTE_abs_2021_stats3b) #AIC = 419.2207
AICc(FKBRTE_abs_2021_stats3b) #AICc = 424.337

#numeric g, d*g
FKBRTE_covdrop2021_2$livestock_util_2020 <- as.numeric(FKBRTE_covdrop2021_2$livestock_util_2020)
FKBRTE_abs_2021_stats4 <- lmerTest::lmer(data = FKBRTE_covdrop2021_2, aerial_cover ~ drought*livestock_util_2020 +
                                           (1|block) + (1|block:paddock))

anova(FKBRTE_abs_2021_stats4, type = 3) #no sig
AIC(FKBRTE_abs_2021_stats4) #AIC = 467.1596
AICc(FKBRTE_abs_2021_stats4) #AICc = 469.5944

#numeric g, d+g
FKBRTE_abs_2021_stats5 <- lmerTest::lmer(data = FKBRTE_covdrop2021_2, aerial_cover ~ drought + livestock_util_2020 +
                                           (1|block) + (1|block:paddock))

anova(FKBRTE_abs_2021_stats5, type = 3) #no sig
AIC(FKBRTE_abs_2021_stats5) #AIC = 455.6338
AICc(FKBRTE_abs_2021_stats5) #AICc = 457.421

#rescale g, d*g
FKBRTE_covdrop2021_3 <- FKBRTE_covdrop2021_2 %>%
  mutate(g2020 = ifelse(livestock_util_2020 == "30", "-1", 
                        ifelse(livestock_util_2020 == "50", "0", 
                               ifelse(livestock_util_2020 == "70", "1", livestock_util_2020)))) %>%
  mutate(g2020 = as.numeric(g2020))

FKBRTE_abs_2021_stats6 <- lmerTest::lmer(data = FKBRTE_covdrop2021_3, aerial_cover ~ drought*g2020 +
                                           (1|block) + (1|block:paddock))

anova(FKBRTE_abs_2021_stats6, type = 3) #no sig
AIC(FKBRTE_abs_2021_stats6) #AIC = 455.1767
AICc(FKBRTE_abs_2021_stats6) #AICc = 457.6115

#rescale g, d+g
FKBRTE_abs_2021_stats7 <- lmerTest::lmer(data = FKBRTE_covdrop2021_3, aerial_cover ~ drought + g2020 +
                                           (1|block) + (1|block:paddock))

anova(FKBRTE_abs_2021_stats7, type = 3) #no sig
AIC(FKBRTE_abs_2021_stats7) #AIC = 449.6423
AICc(FKBRTE_abs_2021_stats7) #AICc = 451.4295

#d only
FKBRTE_abs_2021_stats8 <- lmerTest::lmer(data = FKBRTE_covdrop2021_3, aerial_cover ~ drought +
                                           (1|block) + (1|block:paddock))

anova(FKBRTE_abs_2021_stats8, type = 3) #no sig
AIC(FKBRTE_abs_2021_stats8) #AIC = 452.0575
AICc(FKBRTE_abs_2021_stats8) #AICc = 453.3075






#TB BRAR 2018 - try transformations
TBBRAR_covdrop2018 <- TBBRAR_covdrop %>%
  filter(year == 2018)

TBBRAR_covdrop2018_trans <- TBBRAR_covdrop2018 %>%
  mutate(cov_sqrt = sqrt(aerial_cover)) %>%
  mutate(cov_sq = aerial_cover^2) %>%
  mutate(cov_cb = aerial_cover^1/3) %>%
  mutate(cov_qd = aerial_cover^1/4) %>%
  mutate(cov_ln = log(aerial_cover + 0.1)) %>%
  mutate(cov_log = log10(aerial_cover + 0.1)) %>%
  mutate(cov_cr = aerial_cover^0.75)

resid4.48 <- lm(data = TBBRAR_covdrop2018_trans, cov_ln ~ drought*grazing_treatment)
ols_plot_resid_hist(resid4.48) #normalish
ols_test_normality(resid4.48) #passed normality tests

#stats with ln
TBBRAR_abs_2018_stats <- lmerTest::lmer(data = TBBRAR_covdrop2018_trans, cov_ln ~ drought*grazing_treatment + 
                                          (1|block) + (1|block:paddock))
anova(TBBRAR_abs_2018_stats, type = 3) #nothing significant


#TB BRAR 2019 - try transformations
TBBRAR_covdrop2019 <- TBBRAR_covdrop %>%
  filter(year == 2019)

#stats with raw data bc no transformations worked
TBBRAR_abs_2019_stats <- lmerTest::lmer(data = TBBRAR_covdrop2019, aerial_cover ~ drought*grazing_treatment + 
                                          (1|block) + (1|block:paddock))
anova(TBBRAR_abs_2019_stats, type = 3) #nothing significant

#revisit with just drought
TBBRAR_abs_2019_stats2 <- lmerTest::lmer(data = TBBRAR_covdrop2019, aerial_cover ~ drought + 
                                           (1|block) + (1|block:paddock))
anova(TBBRAR_abs_2019_stats2, type = 3) #nothing significant


#TB BRAR 2020 - try transformations
TBBRAR_covdrop2020 <- TBBRAR_covdrop %>%
  filter(year == 2020)

TBBRAR_covdrop2020_trans <- TBBRAR_covdrop2020 %>%
  mutate(cov_sqrt = sqrt(aerial_cover)) %>%
  mutate(cov_sq = aerial_cover^2) %>%
  mutate(cov_cb = aerial_cover^1/3) %>%
  mutate(cov_qd = aerial_cover^1/4) %>%
  mutate(cov_ln = log(aerial_cover + 0.1)) %>%
  mutate(cov_log = log10(aerial_cover + 0.1)) %>%
  mutate(cov_cr = aerial_cover^0.75)

resid4.57 <- lm(data = TBBRAR_covdrop2020_trans, cov_ln ~ drought*grazing_treatment)
ols_plot_resid_hist(resid4.57) #normalish
ols_test_normality(resid4.57) #passed normality tests except cramer

#stats with ln
TBBRAR_abs_2020_stats <- lmerTest::lmer(data = TBBRAR_covdrop2020_trans, cov_ln ~ drought*grazing_treatment + 
                                          (1|block) + (1|block:paddock))
anova(TBBRAR_abs_2020_stats, type = 3) #nothing significant

#revisit with 2019 grazing
TBBRAR_covdrop2020_trans$livestock_util_2019 <- as.factor(TBBRAR_covdrop2020_trans$livestock_util_2019)
TBBRAR_abs_2020_stats2 <- lmerTest::lmer(data = TBBRAR_covdrop2020_trans, cov_ln ~ drought*livestock_util_2019 + 
                                           (1|block) + (1|block:paddock))
anova(TBBRAR_abs_2020_stats2, type = 3) #nothing significant



#TB BRAR 2021 - try transformations
TBBRAR_covdrop2021 <- TBBRAR_covdrop %>%
  filter(year == 2021)

TBBRAR_covdrop2021_trans <- TBBRAR_covdrop2021 %>%
  mutate(cov_sqrt = sqrt(aerial_cover)) %>%
  mutate(cov_sq = aerial_cover^2) %>%
  mutate(cov_cb = aerial_cover^1/3) %>%
  mutate(cov_qd = aerial_cover^1/4) %>%
  mutate(cov_ln = log(aerial_cover + 0.1)) %>%
  mutate(cov_log = log10(aerial_cover + 0.1)) %>%
  mutate(cov_cr = aerial_cover^0.75)

resid4.58 <- lm(data = TBBRAR_covdrop2021_trans, cov_ln ~ drought*grazing_treatment)
ols_plot_resid_hist(resid4.58) #normalish
ols_test_normality(resid4.58) #pass normality tests

#stats with ln
TBBRAR_abs_2021_stats <- lmerTest::lmer(data = TBBRAR_covdrop2021_trans, cov_ln ~ drought*grazing_treatment + 
                                          (1|block) + (1|block:paddock))
anova(TBBRAR_abs_2021_stats, type = 3) #nothing significant

#revisit with 2020 grazing
TBBRAR_covdrop2021_trans2 <- TBBRAR_covdrop2021_trans %>%
  mutate(graze2021 = ifelse(grazing_category == "HHMMM", "HH", 
                            ifelse(grazing_category == "MLLMM", "ML", 
                                   ifelse(grazing_category == "MMMMM", "MM", grazing_category)))) %>%
  mutate(graze2021 = as.factor(graze2021))


TBBRAR_abs_2021_stats2 <- lmerTest::lmer(data = TBBRAR_covdrop2021_trans2, cov_ln ~ drought*graze2021 + 
                                           (1|block) + (1|block:paddock))
anova(TBBRAR_abs_2021_stats2, type = 3) #nothing significant



#TB BRTE 2018 - try transformations
TBBRTE_covdrop2018 <- TBBRTE_covdrop %>%
  filter(year == 2018)

TBBRTE_covdrop2018_trans <- TBBRTE_covdrop2018 %>%
  mutate(cov_sqrt = sqrt(aerial_cover)) %>%
  mutate(cov_sq = aerial_cover^2) %>%
  mutate(cov_cb = aerial_cover^1/3) %>%
  mutate(cov_qd = aerial_cover^1/4) %>%
  mutate(cov_ln = log(aerial_cover + 0.1)) %>%
  mutate(cov_log = log10(aerial_cover + 0.1)) %>%
  mutate(cov_cr = aerial_cover^0.75)

resid4.59 <- lm(data = TBBRTE_covdrop2018_trans, cov_ln ~ drought*grazing_treatment)
ols_plot_resid_hist(resid4.59) #normalish
ols_test_normality(resid4.59) #passed normality tests except A-D and cramer

#stats with ln
TBBRTE_abs_2018_stats <- lmerTest::lmer(data = TBBRTE_covdrop2018_trans, cov_ln ~ drought*grazing_treatment + 
                                          (1|block) + (1|block:paddock))
anova(TBBRTE_abs_2018_stats, type = 3) #grazing p = 0.01686


#TB BRTE 2019 - try transformations
TBBRTE_covdrop2019 <- TBBRTE_covdrop %>%
  filter(year == 2019)

#stats with raw data bc no transformations worked
TBBRTE_abs_2019_stats <- lmerTest::lmer(data = TBBRTE_covdrop2019, aerial_cover ~ drought*grazing_treatment + 
                                          (1|block) + (1|block:paddock))
anova(TBBRTE_abs_2019_stats, type = 3) #nothing significant

#revisit with just drought
TBBRTE_abs_2019_stats2 <- lmerTest::lmer(data = TBBRTE_covdrop2019, aerial_cover ~ drought + 
                                           (1|block) + (1|block:paddock))
anova(TBBRTE_abs_2019_stats2, type = 3) #nothing significant


#TB BRTE 2020 - try transformations
TBBRTE_covdrop2020 <- TBBRTE_covdrop %>%
  filter(year == 2020)

#stats with raw data bc no transformations worked
TBBRTE_abs_2020_stats <- lmerTest::lmer(data = TBBRTE_covdrop2020, aerial_cover ~ drought*grazing_treatment + 
                                          (1|block) + (1|block:paddock))
anova(TBBRTE_abs_2020_stats, type = 3) #nothing significant

#revisit with 2019 graze
TBBRTE_covdrop2020$livestock_util_2019 <- as.factor(TBBRTE_covdrop2020$livestock_util_2019)
TBBRTE_abs_2020_stats2 <- lmerTest::lmer(data = TBBRTE_covdrop2020, aerial_cover ~ drought*livestock_util_2019 + 
                                           (1|block) + (1|block:paddock))
anova(TBBRTE_abs_2020_stats2, type = 3) #nothing significant


#TB BRTE 2021 - try transformations
TBBRTE_covdrop2021 <- TBBRTE_covdrop %>%
  filter(year == 2021)

#stats with raw data bc no transformations worked
TBBRTE_abs_2021_stats <- lmerTest::lmer(data = TBBRTE_covdrop2021, aerial_cover ~ drought*grazing_treatment + 
                                          (1|block) + (1|block:paddock))
anova(TBBRTE_abs_2021_stats, type = 3) #nothing significant

#revisit with grazing from 2020
TBBRTE_covdrop2021_2 <- TBBRTE_covdrop2021 %>%
  mutate(graze2021 = ifelse(grazing_category == "HHMMM", "HH", 
                            ifelse(grazing_category == "MLLMM", "ML", 
                                   ifelse(grazing_category == "MMMMM", "MM", grazing_category)))) %>%
  mutate(graze2021 = as.factor(graze2021))


TBBRTE_abs_2021_stats2 <- lmerTest::lmer(data = TBBRTE_covdrop2021_2, aerial_cover ~ drought*graze2021 + 
                                           (1|block) + (1|block:paddock))
anova(TBBRTE_abs_2021_stats2, type = 3) #nothing significant





#add to supplemental stats with not dropped - only need to do TB bc FK didnt have any plots without bromes ever
TBBRAR_cov2019 <- bromerelcov %>%
  filter(site == "TB" & symbol == "BRAR") %>%
  filter(year == "2019")

TBBRAR_cov2020 <- bromerelcov %>%
  filter(site == "TB" & symbol == "BRAR") %>%
  filter(year == "2020")

TBBRAR_cov2021 <- bromerelcov %>%
  filter(site == "TB" & symbol == "BRAR") %>%
  filter(year == "2021")

TBBRTE_cov2019 <- bromerelcov %>%
  filter(site == "TB" & symbol == "BRTE") %>%
  filter(year == "2019")

TBBRTE_cov2020 <- bromerelcov %>%
  filter(site == "TB" & symbol == "BRTE") %>%
  filter(year == "2020")

TBBRTE_cov2021 <- bromerelcov %>%
  filter(site == "TB" & symbol == "BRTE") %>%
  filter(year == "2021")

#TB 2019 BRAR - no transformation
TBBRAR_abs_2019_stats3 <- lmerTest::lmer(data = TBBRAR_cov2019, aerial_cover ~ drought + 
                                           (1|block) + (1|block:paddock))
anova(TBBRAR_abs_2019_stats3, type = 3) #nothing significant

#TB 2020 BRAR - ln data+0.1
TBBRAR_cov2020_trans <- TBBRAR_cov2020 %>%
  mutate(cov_ln = log(aerial_cover + 0.1)) 

TBBRAR_cov2020_trans$livestock_util_2019 <- as.factor(TBBRAR_cov2020_trans$livestock_util_2019)
TBBRAR_abs_2020_stats3 <- lmerTest::lmer(data = TBBRAR_cov2020_trans, cov_ln ~ drought*livestock_util_2019 + 
                                           (1|block) + (1|block:paddock))
anova(TBBRAR_abs_2020_stats3, type = 3) #nothing significant

#TB 2021 BRAR - ln data+0.1
TBBRAR_cov2021_trans <- TBBRAR_cov2021 %>%
  mutate(cov_ln = log(aerial_cover + 0.1)) 

TBBRAR_cov2021_trans2 <- TBBRAR_cov2021_trans %>%
  mutate(graze2021 = ifelse(grazing_category == "HHMMM", "HH", 
                            ifelse(grazing_category == "MLLMM", "ML", 
                                   ifelse(grazing_category == "MMMMM", "MM", grazing_category)))) %>%
  mutate(graze2021 = as.factor(graze2021))


TBBRAR_abs_2021_stats3 <- lmerTest::lmer(data = TBBRAR_cov2021_trans2, cov_ln ~ drought*graze2021 + 
                                           (1|block) + (1|block:paddock))
anova(TBBRAR_abs_2021_stats3, type = 3) #nothing significant


#TB 2019 BRTE - no transformation
TBBRTE_abs_2019_stats3 <- lmerTest::lmer(data = TBBRTE_cov2019, aerial_cover ~ drought + 
                                           (1|block) + (1|block:paddock))
anova(TBBRTE_abs_2019_stats3, type = 3) #nothing significant

#TB 2020 BRTE - no trans
TBBRTE_cov2020$livestock_util_2019 <- as.factor(TBBRTE_cov2020$livestock_util_2019)
TBBRTE_abs_2020_stats3 <- lmerTest::lmer(data = TBBRTE_cov2020, aerial_cover ~ drought*livestock_util_2019 + 
                                           (1|block) + (1|block:paddock))
anova(TBBRTE_abs_2020_stats3, type = 3) #nothing significant

#TB 2021 BRTE - no trans

#revisit with 2020 grazing
TBBRTE_cov2021_2 <- TBBRTE_cov2021 %>%
  mutate(graze2021 = ifelse(grazing_category == "HHMMM", "HH", 
                            ifelse(grazing_category == "MLLMM", "ML", 
                                   ifelse(grazing_category == "MMMMM", "MM", grazing_category)))) %>%
  mutate(graze2021 = as.factor(graze2021))


TBBRTE_abs_2021_stats3 <- lmerTest::lmer(data = TBBRTE_cov2021_2, aerial_cover ~ drought*graze2021 + 
                                           (1|block) + (1|block:paddock))
anova(TBBRTE_abs_2021_stats3, type = 3) #nothing significant

##########################################################################




######################################################
######## Soil moisture figures ########

#Read in this one
write.csv(FKsoilmoist2, file = "FKsoilmoist2.csv", row.names = FALSE)

#Read in this one
write.csv(TBsoilmoist2, file = "TBsoilmoist2.csv", row.names = FALSE)


#FK 2019
FK2019SM <- FKsoilmoist %>%
  filter(year == "2019") %>%
  drop_na(soil_moist) %>%
  dplyr:: select(c(date, plot, soil_moist)) 

FK2019SM2 <- FK2019SM %>%
  complete(date = seq.Date(min(date), max(date), by = "day"), plot) %>%
  group_by(plot) %>%
  fill(soil_moist) %>%
  ungroup() %>%
  full_join(metadata) %>%
  drop_na(soil_moist) %>%
  filter(site == "FK") %>%
  mutate(drought = ifelse(drought == 1, 0, drought)) %>%
  mutate(drought = ifelse(drought == 2, 0, drought)) %>%
  mutate(drought = as.factor(drought))

FK2019SM3 <- FK2019SM2 %>%
  group_by(drought) %>%
  summarise(avg = mean(soil_moist), se = sd(soil_moist)/sqrt(length(soil_moist))) %>%
  ungroup()

#FK 2020
FK2020SM <- FKsoilmoist %>%
  filter(year == "2020") %>%
  drop_na(soil_moist) %>%
  dplyr:: select(c(date, plot, soil_moist)) 

FK2020SM2 <- FK2020SM %>%
  complete(date = seq.Date(min(date), max(date), by = "day"), plot) %>%
  group_by(plot) %>%
  fill(soil_moist) %>%
  ungroup() %>%
  full_join(metadata) %>%
  drop_na(soil_moist) %>%
  filter(site == "FK") %>%
  mutate(drought = ifelse(drought == 1, 0, drought)) %>%
  mutate(drought = ifelse(drought == 2, 0, drought)) %>%
  mutate(drought = as.factor(drought))

FK2020SM3 <- FK2020SM2 %>%
  group_by(drought) %>%
  summarise(avg = mean(soil_moist), se = sd(soil_moist)/sqrt(length(soil_moist))) %>%
  ungroup()

#FK 2021
FK2021SM <- FKsoilmoist %>%
  filter(year == "2021") %>%
  drop_na(soil_moist) %>%
  dplyr:: select(c(date, plot, soil_moist)) 

FK2021SM2 <- FK2021SM %>%
  complete(date = seq.Date(min(date), max(date), by = "day"), plot) %>%
  group_by(plot) %>%
  fill(soil_moist) %>%
  ungroup() %>%
  full_join(metadata) %>%
  drop_na(soil_moist) %>%
  filter(site == "FK") %>%
  mutate(drought = ifelse(drought == 1, 0, drought)) %>%
  mutate(drought = ifelse(drought == 2, 0, drought)) %>%
  mutate(drought = as.factor(drought))

FK2021SM3 <- FK2021SM2 %>%
  group_by(drought) %>%
  summarise(avg = mean(soil_moist), se = sd(soil_moist)/sqrt(length(soil_moist))) %>%
  ungroup()


#TB 2019
TB2019SM <- TBsoilmoist %>%
  filter(year == "2019") %>%
  drop_na(soil_moist) %>%
  dplyr:: select(c(date, plot, soil_moist)) 

TB2019SM2 <- TB2019SM %>%
  complete(date = seq.Date(min(date), max(date), by = "day"), plot) %>%
  group_by(plot) %>%
  fill(soil_moist) %>%
  ungroup() %>%
  full_join(metadata) %>%
  drop_na(soil_moist) %>%
  filter(site == "TB") %>%
  mutate(drought = ifelse(drought == 1, 0, drought)) %>%
  mutate(drought = ifelse(drought == 2, 0, drought)) %>%
  mutate(drought = as.factor(drought))

TB2019SM3 <- TB2019SM2 %>%
  group_by(drought) %>%
  summarise(avg = mean(soil_moist), se = sd(soil_moist)/sqrt(length(soil_moist))) %>%
  ungroup()

#TB 2020
TB2020SM <- TBsoilmoist %>%
  filter(year == "2020") %>%
  drop_na(soil_moist) %>%
  dplyr:: select(c(date, plot, soil_moist)) 

TB2020SM2 <- TB2020SM %>%
  complete(date = seq.Date(min(date), max(date), by = "day"), plot) %>%
  group_by(plot) %>%
  fill(soil_moist) %>%
  ungroup() %>%
  full_join(metadata) %>%
  drop_na(soil_moist) %>%
  filter(site == "TB") %>%
  mutate(drought = ifelse(drought == 1, 0, drought)) %>%
  mutate(drought = ifelse(drought == 2, 0, drought)) %>%
  mutate(drought = as.factor(drought))

TB2020SM3 <- TB2020SM2 %>%
  group_by(drought) %>%
  summarise(avg = mean(soil_moist), se = sd(soil_moist)/sqrt(length(soil_moist))) %>%
  ungroup()

#TB 2021
TB2021SM <- TBsoilmoist %>%
  filter(year == "2021") %>%
  drop_na(soil_moist) %>%
  dplyr:: select(c(date, plot, soil_moist)) 

TB2021SM2 <- TB2021SM %>%
  complete(date = seq.Date(min(date), max(date), by = "day"), plot) %>%
  group_by(plot) %>%
  fill(soil_moist) %>%
  ungroup() %>%
  full_join(metadata) %>%
  drop_na(soil_moist) %>%
  filter(site == "TB") %>%
  mutate(drought = ifelse(drought == 1, 0, drought)) %>%
  mutate(drought = ifelse(drought == 2, 0, drought)) %>%
  mutate(drought = as.factor(drought))

TB2021SM3 <- TB2021SM2 %>%
  group_by(drought) %>%
  summarise(avg = mean(soil_moist), se = sd(soil_moist)/sqrt(length(soil_moist))) %>%
  ungroup()

################################################
###### Soil moisture - statistics ########


##### Soil moisture analyses #####
#this is avergaing across all time points - get most robust SM picture
#MT 2019
FK2019SM2_2 <- FK2019SM2 %>%
  group_by(plot, site, block, paddock, drought, grazing_treatment, grazing_category, livestock_util_2019, livestock_util_2020, livestock_util_2021) %>%
  summarise(avgSM = mean(soil_moist)) %>%
  ungroup()
FK2019SM2_3 <- FK2019SM2_2 %>%
  group_by(drought, site) %>%
  summarise(avg_SM = mean(avgSM), se_SM = sd(avgSM)/sqrt(length(avgSM))) %>%
  ungroup()

FKsm19_stats <- lmerTest::lmer(data = FK2019SM2_2, avgSM ~ drought + (1|block) + (1|block:paddock))
anova(FKsm19_stats) #p = 7.307e-13
summary(glht(FKsm19_stats, linfct = mcp(drought = "Tukey"), test = adjusted(type = "BH"))) #most levels are very significantly different from one another

FK2019SM4_2 <- FK2019SM2_3 %>%
  mutate(sig = ifelse(drought == "0", "a", 
              ifelse(drought == "25", "b*", 
              ifelse(drought == "50", "b", 
              ifelse(drought == "75", "c", 
              ifelse(drought == "99", "c", "0"))))))

#MT 2020
FK2020SM2_2 <- FK2020SM2 %>%
  group_by(plot, site, block, paddock, drought, grazing_treatment, grazing_category, livestock_util_2019, livestock_util_2020, livestock_util_2021) %>%
  summarise(avgSM = mean(soil_moist)) %>%
  ungroup()
FK2020SM2_3 <- FK2020SM2_2 %>%
  group_by(drought, site) %>%
  summarise(avg_SM = mean(avgSM), se_SM = sd(avgSM)/sqrt(length(avgSM))) %>%
  ungroup()

FKsm20_stats <- lmerTest::lmer(data = FK2020SM2_2, avgSM ~ drought + (1|block) + (1|block:paddock))
anova(FKsm20_stats) #p = 1.009e-11
summary(glht(FKsm20_stats, linfct = mcp(drought = "Tukey"), test = adjusted(type = "BH"))) #most levels are very significantly different from one another

FK2020SM4_2 <- FK2020SM2_3 %>%
  mutate(sig = ifelse(drought == "0", "a", 
              ifelse(drought == "25", "b", 
              ifelse(drought == "50", "b", 
              ifelse(drought == "75", "c*", 
              ifelse(drought == "99", "c", "0"))))))

#MT 2021
FK2021SM2_2 <- FK2021SM2 %>%
  group_by(plot, site, block, paddock, drought, grazing_treatment, grazing_category, livestock_util_2019, livestock_util_2020, livestock_util_2021) %>%
  summarise(avgSM = mean(soil_moist)) %>%
  ungroup()
FK2021SM2_3 <- FK2021SM2_2 %>%
  group_by(drought, site) %>%
  summarise(avg_SM = mean(avgSM), se_SM = sd(avgSM)/sqrt(length(avgSM))) %>%
  ungroup()

FKsm21_stats <- lmerTest::lmer(data = FK2021SM2_2, avgSM ~ drought + (1|block) + (1|block:paddock))
anova(FKsm21_stats) #p = 0.2442


#WY 2019
TB2019SM2_2 <- TB2019SM2 %>%
  group_by(plot, site, block, paddock, drought, grazing_treatment, grazing_category, livestock_util_2019, livestock_util_2020, livestock_util_2021) %>%
  summarise(avgSM = mean(soil_moist)) %>%
  ungroup()
TB2019SM2_3 <- TB2019SM2_2 %>%
  group_by(drought, site) %>%
  summarise(avg_SM = mean(avgSM), se_SM = sd(avgSM)/sqrt(length(avgSM))) %>%
  ungroup()

TBsm19_stats <- lmerTest::lmer(data = TB2019SM2_2, avgSM ~ drought + (1|block) + (1|block:paddock))
anova(TBsm19_stats) #p = 0.0002875
summary(glht(TBsm19_stats, linfct = mcp(drought = "Tukey"), test = adjusted(type = "BH"))) #most levels are very significantly different from one another

TB2019SM4_2 <- TB2019SM2_3 %>%
  mutate(sig = ifelse(drought == "0", "a", 
              ifelse(drought == "25", "a, b", 
              ifelse(drought == "50", "a, b", 
              ifelse(drought == "75", "b", 
              ifelse(drought == "99", "b, c", "0"))))))

#WY 2020
TB2020SM2_2 <- TB2020SM2 %>%
  group_by(plot, site, block, paddock, drought, grazing_treatment, grazing_category, livestock_util_2019, livestock_util_2020, livestock_util_2021) %>%
  summarise(avgSM = mean(soil_moist)) %>%
  ungroup()
TB2020SM2_3 <- TB2020SM2_2 %>%
  group_by(drought, site) %>%
  summarise(avg_SM = mean(avgSM), se_SM = sd(avgSM)/sqrt(length(avgSM))) %>%
  ungroup()

TBsm20_stats <- lmerTest::lmer(data = TB2020SM2_2, avgSM ~ drought + (1|block) + (1|block:paddock))
anova(TBsm20_stats) #p = 0.003725
summary(glht(TBsm20_stats, linfct = mcp(drought = "Tukey"), test = adjusted(type = "BH"))) #most levels are very significantly different from one another

TB2020SM4_2 <- TB2020SM2_3 %>%
  mutate(sig = ifelse(drought == "0", "a", 
              ifelse(drought == "25", "a", 
              ifelse(drought == "50", "a", 
              ifelse(drought == "75", "b", 
              ifelse(drought == "99", "b", "0"))))))

#WY 2021
TB2021SM2_2 <- TB2021SM2 %>%
  group_by(plot, site, block, paddock, drought, grazing_treatment, grazing_category, livestock_util_2019, livestock_util_2020, livestock_util_2021) %>%
  summarise(avgSM = mean(soil_moist)) %>%
  ungroup()
TB2021SM2_3 <- TB2021SM2_2 %>%
  group_by(drought, site) %>%
  summarise(avg_SM = mean(avgSM), se_SM = sd(avgSM)/sqrt(length(avgSM))) %>%
  ungroup()

TBsm21_stats <- lmerTest::lmer(data = TB2021SM2_2, avgSM ~ drought + (1|block) + (1|block:paddock))
anova(TBsm21_stats) #p = 0.1285



#graphs
FKsm19fill2 <- ggplot(data = FK2019SM4_2, aes(x = drought, y = avg_SM)) +
  ylim(0, 30) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = avg_SM - se_SM, ymax = avg_SM + se_SM), width = .2, position = position_dodge(.05)) +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16)) +
  labs(x = "Precipitation Reduction (%)", y = "Soil Moisture (VWC%)", title = "MT 2019") +
  geom_text(aes(label = sig, y = avg_SM + se_SM), vjust = -0.5, size = 6)

FKsm20fill2 <- ggplot(data = FK2020SM4_2, aes(x = drought, y = avg_SM)) +
  ylim(0, 30) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = avg_SM - se_SM, ymax = avg_SM + se_SM), width = .2, position = position_dodge(.05)) +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16)) +
  labs(x = "Precipitation Reduction (%)", y = "Soil Moisture (VWC%)", title = "MT 2020") +
  geom_text(aes(label = sig, y = avg_SM + se_SM), vjust = -0.5, size = 6)

FKsm21fill2 <- ggplot(data = FK2021SM2_3, aes(x = drought, y = avg_SM)) +
  ylim(0, 30) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = avg_SM - se_SM, ymax = avg_SM + se_SM), width = .2, position = position_dodge(.05)) +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16)) +
  labs(x = "Precipitation Reduction (%)", y = "Soil Moisture (VWC%)", title = "MT 2021") 

TBsm19fill2 <- ggplot(data = TB2019SM4_2, aes(x = drought, y = avg_SM)) +
  ylim(0, 30) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = avg_SM - se_SM, ymax = avg_SM + se_SM), width = .2, position = position_dodge(.05)) +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16)) +
  labs(x = "Precipitation Reduction (%)", y = "Soil Moisture (VWC%)", title = "WY 2019") +
  geom_text(aes(label = sig, y = avg_SM + se_SM), vjust = -0.5, size = 6)

TBsm20fill2 <- ggplot(data = TB2020SM4_2, aes(x = drought, y = avg_SM)) +
  ylim(0, 30) + 
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = avg_SM - se_SM, ymax = avg_SM + se_SM), width = .2, position = position_dodge(.05)) +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16)) +
  labs(x = "Precipitation Reduction (%)", y = "Soil Moisture (VWC%)", title = "WY 2020") +
  geom_text(aes(label = sig, y = avg_SM + se_SM), vjust = -0.5, size = 6)

TBsm21fill2 <- ggplot(data = TB2021SM2_3, aes(x = drought, y = avg_SM)) +
  ylim(0, 30) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = avg_SM - se_SM, ymax = avg_SM + se_SM), width = .2, position = position_dodge(.05)) +
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16)) +
  labs(x = "Precipitation Reduction (%)", y = "Soil Moisture (VWC%)", title = "WY 2021") 

FKsm19fill2 + TBsm19fill2 + FKsm20fill2 + TBsm20fill2 + FKsm21fill2 + TBsm21fill2 + plot_layout(ncol = 2)




##################################################################



##########################################################
######### Precip figure ##########
#comes from NOAA data
FKprecip$month_year <- factor(FKprecip$month_year, levels = c("Jan_2019", "Feb_2019", "Mar_2019", "Apr_2019", "May_2019", 
                                                          "Jun_2019", "Jul_2019", "Aug_2019", "Sep_2019", "Oct_2019", 
                                                          "Nov_2019", "Dec_2019",
                                                          "Jan_2020", "Feb_2020", "Mar_2020", "Apr_2020", "May_2020", 
                                                          "Jun_2020", "Jul_2020", "Aug_2020", "Sep_2020", "Oct_2020", 
                                                          "Nov_2020", "Dec_2020", 
                                                          "Jan_2021", "Feb_2021", "Mar_2021", "Apr_2021", "May_2021", 
                                                          "Jun_2021", "Jul_2021", "Aug_2021", "Sep_2021", "Oct_2021", 
                                                          "Nov_2021", "Dec_2021"))


TBprecip$month_year <- factor(TBprecip$month_year, levels = c("Jan_2019", "Feb_2019", "Mar_2019", "Apr_2019", "May_2019", 
                                                              "Jun_2019", "Jul_2019", "Aug_2019", "Sep_2019", "Oct_2019", 
                                                              "Nov_2019", "Dec_2019",
                                                              "Jan_2020", "Feb_2020", "Mar_2020", "Apr_2020", "May_2020", 
                                                              "Jun_2020", "Jul_2020", "Aug_2020", "Sep_2020", "Oct_2020", 
                                                              "Nov_2020", "Dec_2020", 
                                                              "Jan_2021", "Feb_2021", "Mar_2021", "Apr_2021", "May_2021", 
                                                              "Jun_2021", "Jul_2021", "Aug_2021", "Sep_2021", "Oct_2021", 
                                                              "Nov_2021", "Dec_2021"))


FKprecip_fig <- ggplot(data = FKprecip) +
  geom_bar(aes(x = month_year, y = precip_mm), stat = "identity") +
  geom_point(aes(x = month_year, y = cumm)) + 
  geom_line(aes(x = month_year, y = cumm), group = 1) + 
  labs(x = "", y = "Precipitation (mm)", title = "MT") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16)) +
  geom_hline(yintercept = 345.44, linetype = "dashed", color = "red") + 
  scale_x_discrete(labels=c("Jan_2019" = "Jan", "Feb_2019" = "Feb", "Mar_2019" = "Mar", "Apr_2019" = "Apr",
                            "May_2019" = "May", "Jun_2019" = "Jun", "Jul_2019" = "Jul",
                            "Aug_2019" = "Aug", "Sep_2019" = "Sep", "Oct_2019" = "Oct", "Nov_2019" = "Nov", "Dec_2019" = "Dec", 
                            "Jan_2020" = "Jan", "Feb_2020" = "Feb", "Mar_2020" = "Mar", "Apr_2020" = "Apr", "May_2020" = "May", 
                            "Jun_2020" = "Jun", "Jul_2020" = "Jul", "Aug_2020" = "Aug", "Sep_2020" = "Sep", "Oct_2020" = "Oct", 
                            "Nov_2020" = "Nov", "Dec_2020" = "Dec", "Jan_2021" = "Jan", "Feb_2021" = "Feb", "Mar_2021" = "Mar", 
                            "Apr_2021" = "Apr", "May_2021" = "May", "Jun_2021" = "Jun", "Jul_2021" = "Jul", "Aug_2021" = "Aug", "Sep_2021" = "Sep", 
                            "Oct_2021" = "Oct", "Nov_2021" = "Nov", "Dec_2021" = "Dec"))  


TBprecip_fig <- ggplot(data = TBprecip) +
  geom_bar(aes(x = month_year, y = precip_mm), stat = "identity") +
  geom_point(aes(x = month_year, y = cumm)) + 
  geom_line(aes(x = month_year, y = cumm), group = 1) + 
  labs(x = "Month", y = "Precipitation (mm)", title = "WY") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  theme(plot.title = element_text(vjust = - 7, hjust = 0.07, size = 16)) +
  geom_hline(yintercept = 366.268, linetype = "dashed", color = "red") + 
  scale_x_discrete(labels=c("Jan_2019" = "Jan", "Feb_2019" = "Feb", "Mar_2019" = "Mar", "Apr_2019" = "Apr",
                            "May_2019" = "May", "Jun_2019" = "Jun", "Jul_2019" = "Jul",
                            "Aug_2019" = "Aug", "Sep_2019" = "Sep", "Oct_2019" = "Oct", "Nov_2019" = "Nov", "Dec_2019" = "Dec", 
                            "Jan_2020" = "Jan", "Feb_2020" = "Feb", "Mar_2020" = "Mar", "Apr_2020" = "Apr", "May_2020" = "May", 
                            "Jun_2020" = "Jun", "Jul_2020" = "Jul", "Aug_2020" = "Aug", "Sep_2020" = "Sep", "Oct_2020" = "Oct", 
                            "Nov_2020" = "Nov", "Dec_2020" = "Dec", "Jan_2021" = "Jan", "Feb_2021" = "Feb", "Mar_2021" = "Mar", 
                            "Apr_2021" = "Apr", "May_2021" = "May", "Jun_2021" = "Jun", "Jul_2021" = "Jul", "Aug_2021" = "Aug", "Sep_2021" = "Sep", 
                            "Oct_2021" = "Oct", "Nov_2021" = "Nov", "Dec_2021" = "Dec"))  


FKprecip_fig / TBprecip_fig

#rain 1800 x 1400

