# STEP 4.1

# FILES REQUIRED: Rarefied fungal richness tables (from '2_Filering_and_rarefying.R') and file with descriptions of abiotic soil parameters, fire histories, and soil origins.

# The aim of this code is to make two bar plots showing either total species richness or Glomeromycota richness grouped by fire history and soil origin. We'll also make files containing means of richness within different groups.


## RICHNESS BAR PLOTS

library(tidyverse)
library(ggplot2)
library(ggtext)
library(gridtext)

rich_data <- read.csv("C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Rarefied/Richness_rarefy_15000.csv")

glomero_data <- read.csv("C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Rarefied/Glomero_richness_rarefy_15000.csv")

sample_summaries <- read.csv("C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Sample_summaries.csv")


all_data <- inner_join(rich_data, sample_summaries, by = "Sample")
all_data$Soil_Fire <- paste(all_data$Soil_origin, all_data$Fire_history)

all_glomero <- inner_join(glomero_data, sample_summaries, by = "Sample")
all_glomero$Soil_Fire <- paste(all_glomero$Soil_origin, all_glomero$Fire_history)


## MEANS OF RICHNESS VALUES BY FIRE HISTORY AND/OR SOIL ORIGIN

#Calculating means and 95% confidence intervals
z <- 1.96

means_all <- all_data %>% 
  group_by(Soil_Fire) %>%
  summarise(sample = Sample, richness = Richness, fire_history = Fire_history, soil_origin = Soil_origin, mean = mean(Richness), sd = sd(Richness), n = n(), ci = z*sd/sqrt(n), max = mean+ci, min = mean-ci)
str(means_all)
write.table(means_all, file = "C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Richness/Means_soil_fire.csv", sep=",", row.names = F)

means_all_fire <- all_data %>% 
  group_by(Fire_history) %>%
  summarise(mean = mean(Richness), sd = sd(Richness), n = n(), ci = z*sd/sqrt(n), max = mean+ci, min = mean-ci)
head(means_all_fire)
write.table(means_all_fire, file = "C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Richness/Means_fire.csv", sep=",", row.names = F)

means_all_soil <- all_data %>% 
  group_by(Soil_origin) %>%
  summarise(mean = mean(Richness), sd = sd(Richness), n = n(), ci = z*sd/sqrt(n), max = mean+ci, min = mean-ci)
head(means_all_soil)
write.table(means_all_soil, file = "C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Richness/Means_soil.csv", sep=",", row.names = F)


# Making bar plot grouped by fire and soil

jpeg('C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Richness/Richness_soil_fire.jpg', width=1900, height=1300, quality=100, res=300, pointsize = 10)
ggplot(data = means_all, aes(x = fire_history, y = mean, 
                                      ymin = min, ymax = max, fill=soil_origin)) + 
  geom_bar(stat = "identity", width = 0.5,  position = position_dodge(0.6)) + theme_bw() + 
  geom_errorbar(width = 0.2, position=position_dodge(0.6)) +  
  scale_fill_viridis_d(option = "viridis", begin = 0.4, end = 0.7) +
  geom_point(aes(y= richness), size=1, shape=21, position=position_dodge(0.6)) +
  xlab("Fire history") + ylab("Species richness") +
  labs(fill = "Soil origin") +
  ggtitle("Soil fungal species richness by fire history and soil origin")
dev.off()



### REPEATING SAME STEPS USING GLOMEROMYCOTA RICHNESS VALUES

## BY FIRE HISTORY

#Calculating means and 95% confidence intervals
z <- 1.96

means_glomero <- all_glomero %>% 
  group_by(Soil_Fire) %>%
  summarise(sample = Sample, richness = Richness, fire_history = Fire_history, soil_origin = Soil_origin, mean = mean(Richness), sd = sd(Richness), n = n(), ci = z*sd/sqrt(n), max = mean+ci, min = mean-ci)
str(means_glomero)
write.table(means_glomero, file = "C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Richness/Glomero_means_soil_fire.csv", sep=",", row.names = F)

means_glomero_fire <- all_glomero %>% 
  group_by(Fire_history) %>%
  summarise(mean = mean(Richness), sd = sd(Richness), n = n(), ci = z*sd/sqrt(n), max = mean+ci, min = mean-ci)
head(means_glomero_fire)
write.table(means_all_fire, file = "C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Richness/Glomero_means_fire.csv", sep=",", row.names = F)

means_glomero_soil <- all_glomero %>% 
  group_by(Soil_origin) %>%
  summarise(mean = mean(Richness), sd = sd(Richness), n = n(), ci = z*sd/sqrt(n), max = mean+ci, min = mean-ci)
head(means_glomero_soil)
write.table(means_all_soil, file = "C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Richness/Glomero_means_soil.csv", sep=",", row.names = F)


#Making bar plot grouped by fire and soil

jpeg('C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Richness/Glomero_richness_soil_fire.jpg', width=1900, height=1300, quality=100, res=300, pointsize = 10)
ggplot(data = means_glomero, aes(x = fire_history, y = mean, 
                                      ymin = min, ymax = max, fill=soil_origin)) + 
  geom_bar(stat = "identity", width = 0.5,  position = position_dodge(0.6)) + theme_bw() + 
  geom_errorbar(width = 0.2, position=position_dodge(0.6)) +  
  scale_fill_viridis_d(option = "plasma", begin = 0.7, end = 0.9) +
  geom_point(aes(y= richness), size=1, shape=21, position=position_dodge(0.6)) +
  xlab("Fire history") + ylab("Glomeromycota richness") +
  labs(fill = "Soil origin") +
  ggtitle("Glomeromycota species richness by fire history and soil origin")
dev.off()

