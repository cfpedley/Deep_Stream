# STEP 3

# FILES REQUIRED: Rarefied fungal richness tables (.csv files) produced via R script '2_Filtering_and_rarefying', files with soil data (these were obtained by from my supervisor, with some edits)

# This code aims to explore soil parameters to make sure correlated parameters aren't used as predictors in the same models. We'll also check the distribution of each parameter, make plots of soil parameters vs fire history, and perform statistical analyses to see if there is a relationship between soil parameters and fire history.


library(tidyverse)
library(ggplot2)
library(GGally)
library(plotly)
library(viridisLite)
library(lme4)
library(lmerTest)
library(pbkrtest)
library(emmeans)


# Files containing soil data
# plot level
soil_data <- read.csv("C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Soil/Soil_plot_means.csv")
# sample level
sample_summaries <- read.csv("C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Sample_summaries.csv")

# Files containing richness data
richness <- read.csv("C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Rarefied/Richness_rarefy_15000.csv")
glomero_richness <- read.csv("C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Rarefied/Glomero_richness_rarefy_15000.csv")

richness_soil <- inner_join(sample_summaries, richness, by = "Sample")
grichness_soil <- inner_join(sample_summaries, glomero_richness, by = "Sample")


## HISTOGRAMS OF EACH SOIL PARAMETER

#pH
jpeg('C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Soil/pH.jpg', width=1000, height=800, quality = 100, res=300, pointsize = 7)
hist(soil_data$pH, main="Histogram of soil pH values", xlab = "pH", ylab = "Frequency")
dev.off()
#pH logged
jpeg('C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Soil/pH_logged.jpg', width=1000, height=800, quality = 100, res=300, pointsize = 7)
hist(soil_data$pH_logged, main="Histogram of logged soil pH values", xlab = "pH logged", ylab = "Frequency")
dev.off()

#P
jpeg('C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Soil/P.jpg', width=1000, height=800, quality = 100, res=300, pointsize = 7)
hist(soil_data$P, main="Histogram of soil P values", xlab = "P", ylab = "Frequency")
dev.off()
#P logged
jpeg('C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Soil/P_logged.jpg', width=1000, height=800, quality = 100, res=300, pointsize = 7)
hist(soil_data$P_logged, main="Histogram of logged soil P values", xlab = "P logged", ylab = "Frequency")
dev.off()

#Total C
jpeg('C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Soil/total_C.jpg', width=1000, height=800, quality = 100, res=300, pointsize = 7)
hist(soil_data$total_C, main="Histogram of soil total C values", xlab = "total C (%)", ylab = "Frequency")
dev.off()
#Total C logged
jpeg('C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Soil/total_C_logged.jpg', width=1000, height=800, quality = 100, res=300, pointsize = 7)
hist(soil_data$C_logged, main="Histogram of logged soil total C values", xlab = "logged total C (%)", ylab = "Frequency")
dev.off()



## MATRIX OF RELATIONSHIPS BETWEEN PARAMETERS

#Matrix showing every abiotic soil parameter
jpeg('C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Soil/Soil_matrix.jpg', width=5000, height=4500, quality=100, res=300)
ggpairs(soil_data, columns = c(3,5,7:9,11:18),
      ggplot2::aes(color = Fire_history)) +
  scale_fill_viridis_d(option="viridis") +
  scale_colour_viridis_d(option="viridis") +
  theme_bw()
dev.off()

#Matrix showing only predictors to be used in stats models (soil carbon, pH, and phosphorus)
jpeg('C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Soil/Soil_matrix_predictors.jpg', width=900, height=700, quality=100, res=150)
ggpairs(soil_data, columns = c(3,5,9),
        ggplot2::aes(color = Fire_history)) +
  scale_fill_viridis_d(option="viridis") +
  scale_colour_viridis_d(option="viridis") +
  theme_bw() +
  theme(text = element_text(size = 20))
dev.off()



## PLOTS OF SOIL PREDICTORS BY FIRE HISTORY

# pH
jpeg('C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Soil/pH_fire_history.jpg', width=1900, height=1300, quality=100, res=500, pointsize = 10)
ggplot(sample_summaries, aes(x=Fire_history, y=pH, fill=Fire_history)) + 
  geom_point(shape = 21, size = 3) +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0.10, 0.25))) +
  scale_fill_viridis_d() +
  xlab("Fire history") + ylab("pH") +
  ggtitle("Plot of soil pH vs fire history") +
  guides(fill = "none")
dev.off()

# Total C
jpeg('C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Soil/Total_C_fire_history.jpg', width=1900, height=1300, quality=100, res=500, pointsize = 10)
ggplot(sample_summaries, aes(x=Fire_history, y=total_C, fill=Fire_history)) + 
  geom_point(shape = 21, size = 3) +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0.10, 0.25))) +
  scale_fill_viridis_d() +
  xlab("Fire history") + ylab("Total carbon (%)") +
  ggtitle("Plot of soil carbon vs fire history") +
  guides(fill = "none")
dev.off()

# P (this one is slightly different from the two figures above because two plots had the same phosphorus value, so I needed to plot the values differently)
P_summary <- sample_summaries %>%
  group_by(Plot) %>%
  summarise(P = mean(P)) %>%
  mutate(Fire_history = c("Control", "Control", "Control", "Spring", "Spring", "Spring", "Summer", "Summer", "Summer", "Unburnt", "Unburnt", "Unburnt"))
jpeg('C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Soil/P_fire_history.jpg', width=1900, height=1300, quality=100, res=500, pointsize = 10)
ggplot(P_summary, aes(x=Fire_history, y=P, fill=Fire_history)) + 
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 2) +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0.10, 0.25))) +
  scale_fill_viridis_d() +
  xlab("Fire history") + ylab("Phosphorus (mg/L)") +
  ggtitle("Plot of soil P vs fire history") +
  guides(fill = "none")
dev.off()


# LINEAR MODELS OF SOIL PARAMETERS AND FIRE HISTORY

# RESPONSE: pH logged
# PREDICTOR: fire_history
pH_regression = lm(pH_logged ~ Fire_history, data = sample_summaries)
summary(pH_regression)
anova(pH_regression)
pH_means <- emmeans(pH_regression, spec="Fire_history")
pairs(pH_means)

# RESPONSE: total C logged
# PREDICTOR: fire_history
C_regression = lm(C_logged ~ Fire_history, data = sample_summaries)
summary(C_regression)
anova(C_regression)
C_means <- emmeans(C_regression, spec="Fire_history")
pairs(C_means)

# RESPONSE: P logged
# PREDICTOR: Fire_history
P_regression = lm(P_logged ~ Fire_history, data = sample_summaries)
summary(P_regression)
anova(P_regression)
P_means <- emmeans(P_regression, spec="Fire_history")
pairs(P_means)






## ARCHIVE
  

## BAR GRAPHS OF SOIL PREDICTORS AND FIRE HISTORY
# Note that these figures were not used in my thesis because there were two few data points to bother making bar graphs, but I'm leaving this code here in case anyone else wants to make bar graphs of soil parameters in the future

## pH
# Calculating means and 95% confidence intervals
z <- 1.96
pH_means <- sample_summaries %>% 
  group_by(Fire_history) %>%
  summarise(pH = pH, mean = mean(pH), sd = sd(pH), n = n(), ci = z*sd/sqrt(n), max = mean+ci, min = mean-ci)
head(pH_means)
write.table(pH_means, file = "C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Soil/Means_pH.csv", sep=",", row.names = F)
jpeg('C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Soil/pH_fire_history.jpg', width=1900, height=1300, quality=100, res=300, pointsize = 10)
ggplot(data = pH_means, aes(x = Fire_history, y = mean, 
                             ymin = min, ymax = max, fill=Fire_history)) + 
  geom_bar(stat = "identity", width = 0.5,  position = position_dodge(0.6)) + theme_bw() + 
  geom_errorbar(width = 0.2, position=position_dodge(0.6)) +  
  geom_point(aes(y= pH), size=1, shape=21, position=position_dodge(0.6)) +
  scale_fill_viridis_d(option = "viridis", begin = 0, end = 1) +
  xlab("Fire history") + ylab("pH") +
  guides(fill = "none") +
  ggtitle("Soil pH by fire history")
dev.off()

## Total C
# Calculating means and 95% confidence intervals
z <- 1.96
C_means <- sample_summaries %>% 
  group_by(Fire_history) %>%
  summarise(total_C = total_C, mean = mean(total_C), sd = sd(total_C), n = n(), ci = z*sd/sqrt(n), max = mean+ci, min = mean-ci)
head(C_means)
write.table(C_means, file = "C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Soil/Means_total_C.csv", sep=",", row.names = F)
jpeg('C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Soil/Total_C_fire_history.jpg', width=1900, height=1300, quality=100, res=300, pointsize = 10)
ggplot(data = C_means, aes(x = Fire_history, y = mean, 
                            ymin = min, ymax = max, fill=Fire_history)) + 
  geom_bar(stat = "identity", width = 0.5,  position = position_dodge(0.6)) + theme_bw() + 
  geom_errorbar(width = 0.2, position=position_dodge(0.6)) +  
  scale_fill_viridis_d(option = "viridis", begin = 0, end = 1) +
  xlab("Fire history") + ylab("Total C") +
  geom_point(aes(y= total_C), size=1, shape=21, position=position_dodge(0.6)) +
  guides(fill = "none") +
  ggtitle("Soil total C by fire history")
dev.off()

## P
# Calculating means and 95% confidence intervals
z <- 1.96
P_means <- sample_summaries %>% 
  group_by(Fire_history) %>%
  summarise(P = P, mean = mean(P), sd = sd(P), n = n(), ci = z*sd/sqrt(n), max = mean+ci, min = mean-ci)
head(P_means)
write.table(C_means, file = "C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Soil/Means_P.csv", sep=",", row.names = F)
jpeg('C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Soil/P_fire_history.jpg', width=1900, height=1300, quality=100, res=300, pointsize = 10)
ggplot(data = P_means, aes(x = Fire_history, y = mean, 
                           ymin = min, ymax = max, fill=Fire_history)) + 
  geom_bar(stat = "identity", width = 0.5,  position = position_dodge(0.6)) + theme_bw() + 
  geom_errorbar(width = 0.2, position=position_dodge(0.6)) +  
  scale_fill_viridis_d(option = "viridis", begin = 0, end = 1) +
  xlab("Fire history") + ylab("Phosphorus") +
  geom_point(aes(y= P), size=1, shape=21, position=position_dodge(0.6)) +
  guides(fill = "none") +
  ggtitle("Soil P by fire history")
dev.off()



## PLOTS OF SOIL PREDICTORS AND TOTAL SPECIES RICHNESS
# Note: again, I didn't use these in my thesis, just leaving here for future reference

# pH
jpeg('C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Soil/pH_richness.jpg', width=1900, height=1300, quality=100, res=300, pointsize = 10)
ggplot(data = richness_soil, aes(x = pH, y = Richness)) +
  geom_point() +
  geom_smooth(method=lm, colour="#22A884FF", fill="#22A884FF") +
  xlab("Soil pH") + ylab("Fungal species richness") +
  ggtitle("Richness by soil pH") +
  theme_bw()
dev.off()

# total C
jpeg('C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Soil/Total_C_richness.jpg', width=1900, height=1300, quality=100, res=300, pointsize = 10)
ggplot(data = richness_soil, aes(x = total_C, y = Richness)) +
  geom_point() +
  geom_smooth(method=lm, colour="#22A884FF", fill="#22A884FF") +
  xlab("Soil total C") + ylab("Fungal species richness") +
  ggtitle("Richness by soil carbon") +
  theme_bw()
dev.off()

# P
jpeg('C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Soil/P_richness.jpg', width=1900, height=1300, quality=100, res=300, pointsize = 10)
ggplot(data = richness_soil, aes(x = P, y = Richness)) +
  geom_point() +
  geom_smooth(method=lm, colour="#22A884FF", fill="#22A884FF") +
  xlab("Soil P") + ylab("Fungal species richness") +
  ggtitle("Richness by soil phosphorus") +
  theme_bw()
dev.off()

## PLOTS OF SOIL PREDICTORS AND GLOMEROMYCOTA SPECIES RICHNESS

# pH
jpeg('C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Soil/pH_glomero_richness.jpg', width=1900, height=1300, quality=100, res=300, pointsize = 10)
ggplot(data = grichness_soil, aes(x = pH, y = Richness)) +
  geom_point() +
  geom_smooth(method=lm, colour="#F68F46FF", fill="#F68F46FF") +
  xlab("Soil pH") + ylab("Glomeromycota species richness") +
  ggtitle("Glomero richness by soil pH") +
  theme_bw()
dev.off()

# total C
jpeg('C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Soil/Total_C_glomero_richness.jpg', width=1900, height=1300, quality=100, res=300, pointsize = 10)
ggplot(data = grichness_soil, aes(x = total_C, y = Richness)) +
  geom_point() +
  geom_smooth(method=lm, colour="#F68F46FF", fill="#F68F46FF") +
  xlab("Soil total C") + ylab("Glomeromycota species richness") +
  ggtitle("Glomero richness by soil carbon") +
  theme_bw()
dev.off()

# P
jpeg('C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Soil/P_glomero_richness.jpg', width=1900, height=1300, quality=100, res=300, pointsize = 10)
ggplot(data = grichness_soil, aes(x = P, y = Richness)) +
  geom_point() +
  geom_smooth(method=lm, colour="#F68F46FF", fill="#F68F46FF") +
  xlab("Soil P") + ylab("Glomeromycota species richness") +
  ggtitle("Glomero richness by soil phosphorus") +
  theme_bw()
dev.off()
