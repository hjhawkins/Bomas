# Plots to examine effect size of kraaling on bare ground and basal cover
# (Other effects include LSUdays, wattle, slope, fire and rain).

#Heidi Hawkins, 2021
#***********************************************************************************************
if(!require(pacman)){install.packages("pacman", dependencies=TRUE); library(pacman)}
p_load(Rmisc, gfcanalysis, dplyr, car, stlplus, tidyverse, tidymodels, Metrics, lubridate, multcomp, multcompView, lsmeans, ggplot2, gtable, ggthemes, data.table, ggpubr)

setwd("D:/KRAALS")
#setwd("C:/Users/01423355/Downloads/OneDrive")

#Plot theme panel.background = element_rect(fill = "transparent", colour = NA),
theme_plotCI_big <- theme(axis.text = element_text(size = 16, color="#5C5C61")) +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "#5C5C61", size=1))+ theme(panel.background = element_rect(fill = "white"))+ theme(axis.title = element_text(size = 16, colour = "#5C5C61"))
theme_plotCI_med <- theme(axis.text = element_text(size = 14, color="#5C5C61")) +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "#5C5C61", size=1))+ theme(panel.background = element_rect(fill = "white"))+ theme(axis.title = element_text(size = 14, colour = "#5C5C61"))
theme_plotCI_small <- theme(axis.text = element_text(size = 12, color="#5C5C61")) +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "#5C5C61", size=1))+ theme(panel.background = element_rect(fill = "white"))+ theme(axis.title = element_text(size = 12, colour = "#5C5C61"))
theme_plotCI_smaller <- theme(axis.text = element_text(size = 11, color="#5C5C61")) +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "#5C5C61", size=1))+ theme(panel.background = element_rect(fill = "white"))+ theme(axis.title = element_text(size = 12, colour = "#5C5C61"))

#load data and check formats and values
data <- read.csv("DataPaired.csv")
names(data)
str(data)
as_tibble(data)
summary(data)
view(data)


#CORRELATION MATRIX: using ggcorrplot for MIXED factor/numerical df 
library(correlation)
library(ggcorrplot)
names(data)
data_select <- data[c(1,3,4:6,14,15,18,21,24,27,30,33,36,39)]
names(data_select)
corr <- correlation(data_select, method = "auto") %>%
  ggcorrplot(hc.order = TRUE, #hierarchical clustering
             type = "lower",
             #method = "circle",
             insig = "blank",
             outline.color = "white",
             lab = TRUE)

corr
cor_test(data_select, "Rain7End","INFIL.KC")


#Vizualize the effect size of kraaling on bare ground and grass basal cover
bare <- data %>%
  ggplot(aes(Bare.C, Bare.KC, color = Wattle)) + 
  geom_point(size=2) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), size = 1, color = "#0193D7", fill = "#7ECBEF") +
  geom_hline(yintercept=0, linetype = "dashed", color = "red", size = 1) +
  ylab("Effect size of kraaling (Kraal-Control)") +
  xlab("Initial bare ground (%)") +
  theme_plotCI_med
bare <- bare +
  theme(legend.title=element_blank(),legend.text=element_text(size=12, color="#5C5C61"))
bare

grass <- data %>%
  ggplot(aes(Grass.C, Grass.KC, color = Wattle)) + 
  geom_point(size=2) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), size = 1, color = "#0193D7", fill = "#7ECBEF") +
  geom_hline(yintercept=0, linetype = "dashed", color = "red", size = 1) +
  ylab("") +
  xlab("Initial grass basal cover (%)") +
  theme_plotCI_med
grass <- grass +
  theme(legend.title=element_blank(),legend.text=element_text(size=12, color="#5C5C61"))
grass

#What about forbs, biomass (dpm) and soil [elements]/SOM?
forb <- data %>%
  ggplot(aes(Forbs.C, Forbs.KC, color = Wattle)) + 
  geom_point(size=2) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), size = 1, color = "#0193D7", fill = "#7ECBEF") +
  geom_hline(yintercept=0, linetype = "dashed", color = "red", size = 1) +
  ylab("Effect size of kraaling (Kraal-Control)") +
  xlab("Initial basal cover of forbs (%)") +
  theme_plotCI_med
forb <- forb +
  theme(legend.title=element_blank(),legend.text=element_text(size=12, color="#5C5C61"))
forb

bio <- data %>%
  ggplot(aes(BIO.C, BIO.KC, color = Wattle)) + 
  geom_point(size=2) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), size = 1, color = "#0193D7", fill = "#7ECBEF") +
  geom_hline(yintercept=0, linetype = "dashed", color = "red", size = 1) +
  ylab("") +
  xlab(expression(paste("Initial herbaceous biomass (kg dry matter", ~ha^-1, ")"))) +
  theme_plotCI_med
bio <- bio +
  theme(legend.title=element_blank(),legend.text=element_text(size=12, color="#5C5C61"))
bio

soilP <- data %>%
  ggplot(aes(soilP.C, soilP.KC, color = Wattle)) + 
  geom_point(size=2) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), size = 1, color = "#0193D7", fill = "#7ECBEF") +
  geom_hline(yintercept=0, linetype = "dashed", color = "red", size = 1) +
  ylab("Effect size of kraaling (Kraal-Control)") +
  xlab("Initial soil P (%)") +
  theme_plotCI_med
soilP <- soilP +
  theme(legend.title=element_blank(),legend.text=element_text(size=12, color="#5C5C61"))
soilP

soilS <- data %>%
  ggplot(aes(soilS.C, soilS.KC, color = Wattle)) + 
  geom_point(size=2) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), size = 1, color = "#0193D7", fill = "#7ECBEF") +
  geom_hline(yintercept=0, linetype = "dashed", color = "red", size = 1) +
  ylab("Effect size of kraaling (Kraal-Control)") +
  xlab("Initial soil S (%)") +
  theme_plotCI_med
soilS <- soilS +
  theme(legend.title=element_blank(),legend.text=element_text(size=12, color="#5C5C61"))
soilS

SOM <- data %>%
  ggplot(aes(SOM.C, SOM.KC, color = Wattle)) + 
  geom_point(size=2) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), size = 1, color = "#0193D7", fill = "#7ECBEF") +
  geom_hline(yintercept=0, linetype = "dashed", color = "red", size = 1) +
  ylab("Effect size of kraaling (Kraal-Control)") +
  xlab("Initial soil organic matter (%)") +
  theme_plotCI_med
SOM <- SOM +
  theme(legend.title=element_blank(),legend.text=element_text(size=12, color="#5C5C61"))
SOM


#Arrange plots together
effect_veg <- ggarrange(bare, grass, forb, bio, labels = c("A", "B", "C", "D"),
                      common.legend = TRUE, legend = "bottom")
effect_veg

effect_soil <- ggarrange(soilP, soilS, SOM, ncol = 1, labels = c("A", "B", "C"),
                        common.legend = TRUE, legend = "bottom")
effect_soil


#***************************************************************************************
#PRINT graphs
ggsave("effect_veg.pdf", plot = effect_veg, width = 25, height = 22, units = "cm",
       dpi = 1000, limitsize = TRUE)
ggsave("effect_veg.png", plot = effect_veg, width = 25, height = 22, units = "cm",
       dpi = 1000, limitsize = TRUE)
ggsave("effect_soil.pdf", plot = effect_soil, width = 14, height = 34, units = "cm",
       dpi = 1000, limitsize = TRUE)
ggsave("effect_soil.png", plot = effect_soil, width = 14, height = 34, units = "cm",
       dpi = 1000, limitsize = TRUE)

dev.off()
