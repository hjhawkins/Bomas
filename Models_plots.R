#Plots and Linear mixed effects models for soil P and S with kraaling

#Heidi Hawkins, 2021
#***********************************************************************************************
if(!require(pacman)){install.packages("pacman", dependencies=TRUE); library(pacman)}
p_load(Hmisc, PerformanceAnalytics, corrplot, ggcorrplot, psych, Rmisc, gfcanalysis, dplyr, ggplot2, gtable, ggthemes, data.table, multcomp, multcompView, ggpubr, car, tidyverse, lme4, lmerTest)

#Hmisc for straight Pearson or Spearman correlation matric
#PerformanceAnalytics for correlation charts
#corrplot for correlation plots
#ggcorrplot for mixed factor/numerical data

setwd("D:/KRAALS")
#setwd("C:/Users/01423355/Downloads/OneDrive")

#Plot theme panel.background = element_rect(fill = "transparent", colour = NA),
theme_plotCI_big <- theme(axis.text = element_text(size = 16, color="#5C5C61")) +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "#5C5C61", size=1))+ theme(panel.background = element_rect(fill = "white"))+ theme(axis.title = element_text(size = 16, colour = "#5C5C61"))
theme_plotCI_med <- theme(axis.text = element_text(size = 14, color="#5C5C61")) +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "#5C5C61", size=1))+ theme(panel.background = element_rect(fill = "white"))+ theme(axis.title = element_text(size = 14, colour = "#5C5C61"))
theme_plotCI_small <- theme(axis.text = element_text(size = 12, color="#5C5C61")) +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "#5C5C61", size=1))+ theme(panel.background = element_rect(fill = "white"))+ theme(axis.title = element_text(size = 12, colour = "#5C5C61"))

#********************
#DATA EXPLORATION
#********************
#load data and check formats and values
data <- read.csv("DataLong.csv")
names(data)
str(data)
head(data)

#DATA DISTRIBUTION for zero inflated (LSU days) and bounded data (soil element%s): 
hist(data$LSU.days) #OR
hist(data$Calcium) #etc use logit transformed data for soil elements

#Logit transform bounded data and fix some column formats 
data$LSU.days <- as.numeric(data$LSU.days)
data$Treatment <- as.factor(data$Treatment)
data$Site.history <- as.factor(data$Site.history)
data$Phosphorus <- log(data$SoilP) #does not improve
data$Sulphur <- logit(data$SoilS) #is already normal
data$Potassium <- logit(data$SoilK) #improves
data$Nitrogen <- logit(data$SoilN) #similar
data$Calcium <- logit(data$SoilCa) #improves
#LSUdays is zero inflated (Poisson distribution?). Log or other transformations to satisfy model requirements perform poorly (see https://besjournals.onlinelibrary.wiley.com/doi/10.1111/j.2041-210X.2010.00021.x)
# while quasi-Poisson and negative binomial models perform well. For R see https://cran.r-project.org/web/packages/pscl/vignettes/countreg.pdf

#BUT checkING various model (Poisson, quasi-Poisson, Negative binomial regression..)
fm_pois <- glm(LSU.days ~ ., data = data, family = poisson)
fm_qpois <- glm(LSU.days ~ ., data = data, family = quasipoisson)
fm_nbin <- MASS::glm.nb(LSU.days ~ ., data = data)
fm_zinb0 <- zeroinfl(LSU.days ~ ., data = data, dist = "negbin")

summary(fm_pois)
summary(fm_qpois)
summary(fm_nbin)
summary(fm_zinb0)

#FINDS None of these models performs well, ie data is not truely Poisson/binomial etc distribution, 
#cube root can result in a more normal distribution and we can try both in models
data$Kraaling <- data$LSU.days^(1/3) #cube root reduces zero inflation most
data$LSU.days_sqrt <- sqrt(data$LSU.days) #square root
hist(data$Kraaling)
hist(data$LSU.days_sqrt)
view(data)
str(data)

#CORRELATION MATRIX: using ggcorrplot for MIXED factor/numerical df 
library(correlation)
library(ggcorrplot)
#reduce autocorrelating factors and format matrix (after initial full correlation)

#Prep dataframe
names(data)
data_select <- data[c(1:4,13,14,15,16,18,19,21,22,23,25,31:36)] #using transformed data
names(data_select)

corr <- correlation(data_select, method = "auto") %>%
  ggcorrplot(hc.order = TRUE, #hierarchical clustering
             type = "lower",
             #method = "circle",
             insig = "blank",
             outline.color = "white",
             lab = TRUE)

corr
#here we see largest effect of kraaling (factor kraal/LSUdays/cube root LSUdays) is on Soil P and S

#(for numerical comparisons using correlation plot or heat map see: http://www.sthda.com/english/wiki/correlation-matrix-a-quick-start-guide-to-analyze-format-and-visualize-a-correlation-matrix-using-r-software )
#Using psych package
library(psych)
pairs.panels(data_select, 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
)

#a few odd anovas
names(data_select)
aov1 <- aov(Bare ~ Site.history, data = data_select)
summary(aov1)
res <- aov1$residuals
hist(res, xlab="Residuals") #ok
#posthoc TukeyHSD
tukey<-TukeyHSD(aov1)
tukey


#********************
#PLOTS
#********************
#PLOTS of soil P and S per treatment and site history
#First assign letters for significant comparisons (Mike's prac code) from data means and tukey test, for plot
#SOIL P
names(data)
(data_mean <- data.frame(data[,c(12,14,27)] %>% group_by(Treatment, Site.history) %>% summarise_each(list(mean))))
view(data_mean)

#2-WAY ANOVA and Tukey posthoc tests for plots, Treatment signif: P = 0.0426 *, Site history not: P= 0.2223
#P: Transformations for soil P did not improve normality, #P not signif with K Wallis test but is with anova
#S: no transformation needed so used ANOVA (also N)
#K, Ca improved by logit transform so used ANOVA

#Soil P
k <- kruskal.test(SoilP ~ Treatment, data = data) 
k
aov <- aov(SoilP ~ Treatment + Site.history, data = data)
summary(aov)
tukey <- glht(aov, linfct = mcp(Treatment = "Tukey")) #make sure var is a factor (not char or other)
summary(tukey)
tukey1 <- glht(aov, linfct = mcp(Site.history = "Tukey")) #make sure var is a factor (not char or other)
summary(tukey1)

#Assign letters for significant comparisons from data means and tukey test, for plot
names(data_mean)
names(tukey)
names(tukey1)

#CLD: compact letter display
tukey.cls <- cld(summary(tukey), level = 0.05, decreasing = FALSE)
#Copy the letters into a data.frame
Kraal_tukey1 <- data.frame(rownames(as.data.frame(tukey.cls$mcletters$Letters)), as.data.frame(tukey.cls$mcletters$Letters))
view(Kraal_tukey1)
#Give names to data frame
names(Kraal_tukey1) <- c("Treatment", "Letter")
view(Kraal_tukey1)
#Combine with the data.frame that has the letters in it
Kraal_tukeyP <- merge(Kraal_tukey1, data_mean)
view(Kraal_tukeyP)
#Make sure that the order of the treatment is preserved for plotting
Kraal_tukeyP$Treatment <- factor(Kraal_tukeyP$Treatment, levels = c("control", "kraal", ordered = TRUE))
Kraal_tukeyP$Site.history <- factor(Kraal_tukeyP$Site.history, levels = c("No wattle or fire", "Wattle no fire", "Wattle plus fire", ordered = TRUE))

#SOIL Ca
names(data)
(data_meanCa <- data.frame(data[,c(12,14,29)] %>% group_by(Treatment, Site.history) %>% summarise_each(list(mean))))
view(data_meanCa)
data_SEsummaryCa1 <- summarySE(data, measurevar=c("SoilCa"), groupvars=c("Treatment"))
view(data_SEsummaryCa1)
data_SEsummaryCa2 <- summarySE(data, measurevar=c("SoilCa"), groupvars=c("Site.history"))
view(data_SEsummaryCa2)


#2-WAY ANOVA and Tukey posthoc tests for plots, P = 0.0507 for Treatment, P = 0.0433 for Site history
aovCa <- aov(SoilCa ~ Treatment + Site.history, data = data)
summary(aovCa)
tukeyCa <- glht(aovCa, linfct = mcp(Site.history = "Tukey")) #make sure var is a factor (not char or other)
summary(tukeyCa)

#SOIL S
names(data)
(data_meanS <- data.frame(data[,c(12,14,30)] %>% group_by(Treatment, Site.history) %>% summarise_each(list(mean))))
view(data_meanS)

#2-WAY ANOVA and Tukey posthoc tests for plots, P = 0.0304 for Treatment, P = 1 for Site history
aovS <- aov(SoilS ~ Treatment + Site.history, data = data)
summary(aovS)
tukeyS <- glht(aovS, linfct = mcp(Treatment = "Tukey")) #make sure var is a factor (not char or other)
summary(tukeyS)
tukeyS1 <- glht(aovS, linfct = mcp(Site.history = "Tukey")) #make sure var is a factor (not char or other)
summary(tukeyS1)

#Assign letters for significant comparisons from data means and tukey test, for plot
names(data_meanS)
names(tukeyS)
names(tukeyS1)

#CLD: compact letter display
tukey.cls <- cld(summary(tukeyS), level = 0.05, decreasing = FALSE)
#Copy the letters into a data.frame
Kraal_tukeyS <- data.frame(rownames(as.data.frame(tukey.cls$mcletters$Letters)), as.data.frame(tukey.cls$mcletters$Letters))
view(Kraal_tukeyS)
#Give names to data frame
names(Kraal_tukeyS) <- c("Treatment", "Letter")
view(Kraal_tukeyS)
#Combine with the data.frame that has the letters in it
Kraal_tukeyS <- merge(Kraal_tukeyS, data_meanS)
view(Kraal_tukeyS)
#Make sure that the order of the Land tenure types is preserved for plotting
Kraal_tukeyS$Treatment <- factor(Kraal_tukeyP$Treatment, levels = c("control", "kraal", ordered = TRUE))
Kraal_tukeyS$Site.history <- factor(Kraal_tukeyP$Site.history, levels = c("No wattle or fire", "Wattle no fire", "Wattle plus fire", ordered = TRUE))

#SOM (marginally signif for site history, P = 0.0548)
names(data)
str(data)
(data_meanSOM <- data.frame(data[,c(12,14,25)] %>% group_by(Treatment, Site.history) %>% summarise_each(list(mean))))
view(data_meanSOM)
#means and se for site history
(data_meanSOM1 <- data.frame(data[,c(12,14,25)] %>% group_by(Site.history) %>% summarise_each(list(mean))))
view(data_meanSOM1)
data_SEsummarySOM1 <- summarySE(data, measurevar=c("SOM"), groupvars=c("Site.history"))
view(data_SEsummarySOM1)

#2-WAY ANOVA and Tukey posthoc tests for plots, P = 0.851 for Treatment P = 0.0548 for Site history
aovSOM <- aov(SOM ~ Treatment + Site.history, data = data)
summary(aovSOM)

tukeySOM <- glht(aovSOM, linfct = mcp(Site.history = "Tukey")) #make sure var is a factor (not char or other)
summary(tukeySOM)
tukeySOM1 <- glht(aovSOM, linfct = mcp(Treatment = "Tukey")) #make sure var is a factor (not char or other)
summary(tukeySOM1)

#Assign letters for significant comparisons from data means and tukey test, for plot
names(data_meanSOM)
names(tukeySOM)
names(tukeySOM1)

#CLD: compact letter display
tukey.cls <- cld(summary(tukeySOM), level = 0.05, decreasing = FALSE)
#Copy the letters into a data.frame
Kraal_tukeySOM <- data.frame(rownames(as.data.frame(tukey.cls$mcletters$Letters)), as.data.frame(tukey.cls$mcletters$Letters))
view(Kraal_tukeySOM)
#Give names to data frame
names(Kraal_tukeySOM) <- c("Site.history", "Letter")
view(Kraal_tukeySOM)
#Combine with the data.frame that has the letters in it
Kraal_tukeySOM <- merge(Kraal_tukeySOM, data_meanSOM)
view(Kraal_tukeySOM)
#Make sure that the order of the Land tenure types is preserved for plotting
Kraal_tukeySOM$Treatment <- factor(Kraal_tukeyP$Treatment, levels = c("control", "kraal", ordered = TRUE))
Kraal_tukeySOM$Site.history <- factor(Kraal_tukeyP$Site.history, levels = c("No wattle or fire", "Wattle no fire", "Wattle plus fire", ordered = TRUE))

#********************
#PLOT DATA with stats
#********************
names(data)

#Boxplot with letters
P <- data %>%
  ggplot(aes(x=interaction(Treatment, Site.history), SoilP)) +
  geom_boxplot(aes(fill = Treatment), size = .3, outlier.size = .6) +
  stat_summary(data = data, aes(x=interaction(Treatment, Site.history), SoilP), fun = mean, geom = "point", size=8, shape=21, color="black", fill="white") +
  xlab("Site history") +
  ylab("Total soil P (%)") +
  #Now put the letters on the graph from the post-hoc Tukey tests
  geom_text(data=Kraal_tukeyP, hjust=0.5, vjust=0.25, aes(x=interaction(Treatment,Site.history), y=SoilP, label = Letter), size=5) +
  scale_x_discrete(labels = c("No wattle or fire", "", "Wattle no fire", "", "Wattle plus fire", "")) +
  theme(axis.text.x = element_text(hjust=0.1)) +
  theme_plotCI_big +
  theme(legend.text=element_text(size=12, color="#5C5C61"), legend.title = element_text(size=12, color="#5C5C61")) +
  labs(x="") +
  coord_cartesian(ylim = c(0,0.08))
P

S <- data %>%
  ggplot(aes(x=interaction(Treatment, Site.history), SoilS)) +
  geom_boxplot(aes(fill = Treatment), size = .3, outlier.size = .6) +
  stat_summary(data = data, aes(x=interaction(Treatment, Site.history), SoilS), fun = mean, geom = "point", size=8, shape=21, color="black", fill="white") +
  xlab("Site history") +
  ylab("Total soil S (%)") +
  #Now put the letters on the graph from the post-hoc Tukey tests
  geom_text(data=Kraal_tukeyS, hjust=0.5, vjust=0.25, aes(x=interaction(Treatment,Site.history), y=SoilS, label = Letter), size=5) +
  scale_x_discrete(labels = c("No wattle or fire", "", "Wattle no fire", "", "Wattle plus fire", "")) +
  theme(axis.text.x = element_text(hjust=0.1)) +
  theme_plotCI_big +
  theme(legend.text=element_text(size=12, color="#5C5C61"), legend.title = element_text(size=12, color="#5C5C61")) +
  coord_cartesian(ylim = c(0,0.045))
S

#Arrange plots together
Soilbox <- ggarrange(P, S, ncol= 1, labels = c("A", "B"),
                        common.legend = TRUE, legend = "bottom")
Soilbox


#********************
#LMER MODELS for soil
#********************
#predictors based on correlation matrix above

#PHOSPHORUS
#testing significance, null model without fixed effect of interest, note REML must be false
kraalsmodel.null1 = lmer(SoilP ~ Seasons + Site.history + (1|Site), data = data, REML=FALSE)
kraalsmodel.null2 = lmer(SoilP ~ Site.history + (1|Site), data = data, REML=FALSE)
kraalsmodel.null3 = lmer(SoilP ~ Seasons + (1|Site), data = data, REML=FALSE)

#testing significance, full model with fixed effect of interest, note REML must be false
kraalsmodel.model1 = lmer(SoilP ~ Treatment + LSU.days_cr + Seasons + Site.history + (1|Site), data = data, REML=FALSE)
kraalsmodel.model2 = lmer(SoilP ~ Treatment + LSU.days_cr + Site.history + (1|Site), data = data, REML=FALSE)
kraalsmodel.model3 = lmer(SoilP ~ Treatment + LSU.days_cr + Seasons + (1|Site), data = data, REML=FALSE)

#Test signif. using Likelihood ratio
anova(kraalsmodel.null1,kraalsmodel.model1)
anova(kraalsmodel.null2,kraalsmodel.model2)
anova(kraalsmodel.null3,kraalsmodel.model3)

#LMER with lowest AIC value (complexity vs predictive power) is Phosphorus ~ Treatment + LSUdays_cr + Site.history + (1|Site)
#with AIC of -155 and p = 0.04905 *

#Lets check assumptions of linearity, homoskadacity, normality
plot(fitted(kraalsmodel.model2), residuals(kraalsmodel.model2))
hist(residuals(kraalsmodel.model2))
qqnorm(residuals(kraalsmodel.model2))
#looks good 

#SULPHUR
#testing significance, null model without fixed effect of interest, note REML must be false
kraalsmodel.null1 = lmer(SoilS ~ Slope + Site.history + SoilP + (1|Site), data = data, REML=FALSE)
kraalsmodel.null2 = lmer(SoilS ~ Site.history + SoilP + (1|Site), data = data, REML=FALSE)
kraalsmodel.null3 = lmer(SoilS ~ SoilP + (1|Site), data = data, REML=FALSE)
kraalsmodel.null4 = lmer(SoilS ~ (1|Site), data = data, REML=FALSE)

#testing significance, full model with fixed effect of interest, note REML must be false
kraalsmodel.model1 = lmer(SoilS ~ Treatment + LSU.days_cr + Slope + Site.history + SoilP + (1|Site), data = data, REML=FALSE)
kraalsmodel.model2 = lmer(SoilS ~ Treatment + LSU.days_cr + Site.history + SoilP + (1|Site), data = data, REML=FALSE)
kraalsmodel.model3 = lmer(SoilS ~ Treatment + LSU.days_cr + SoilP + (1|Site), data = data, REML=FALSE)
kraalsmodel.model4 = lmer(SoilS ~ Treatment + LSU.days_cr + (1|Site), data = data, REML=FALSE)

#Test signif. using Likelihood ratio
anova(kraalsmodel.null1,kraalsmodel.model1)
anova(kraalsmodel.null2,kraalsmodel.model2)
anova(kraalsmodel.null3,kraalsmodel.model3)
anova(kraalsmodel.null4,kraalsmodel.model4)

#LMER with lowest AIC value (complexity vs predictive power) is model 4 SoilS ~ Treatment + LSUdays_cr + (1 | Site)
#with AIC of -177.19 and p =  0.009458 **

#Lets check assumptions of linearity, homoskadacity, normality
plot(fitted(kraalsmodel.model4), residuals(kraalsmodel.model4))
hist(residuals(kraalsmodel.model4))
qqnorm(residuals(kraalsmodel.model4))
#looks good 

#********************
#PRINT PLOTS (1000 dpi for production, otherwise 300dpi ok)
#********************
ggsave("correlation.tiff", plot = corr, width = 22, height = 22, units = "cm",
       dpi = 400, limitsize = TRUE)
ggsave("correlation.png", plot = corr, width = 22, height = 22, units = "cm",
       dpi = 400, limitsize = TRUE)
ggsave("Soilbox.tiff", plot = Soilbox, width = 16, height = 22, units = "cm",
       dpi = 400, limitsize = TRUE)
ggsave("Soilbox.png", plot = Soilbox, width = 16, height = 22, units = "cm",
       dpi = 400, limitsize = TRUE)

dev.off()
