# STEP 6

# FILES REQUIRED: Rarefied fungal richness tables (from '2_Filering_and_rarefying.R'), file with descriptions of abiotic soil parameters, fire histories, and soil origins, and file with vegetation richness data (obtained from supervisor).

# The aim of this code is to prepare a richness .csv file based on the vegetation survey document Nicola sent me. I'll then use these richness calculations to analyse the relationship between vegetation richness and soil fungal richness.

library(tidyverse)
library(ggplot2)
library(ggpubr)


veg_original <- read.csv("C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Deep Stream quadrat presence 3 times_2022-08-27.csv")

sample_summaries <- read.csv("C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Sample_summaries.csv")

fungi_richness <- read.csv("C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Rarefied/Richness_rarefy_15000.csv") 
fungi_richness <- inner_join(fungi_richness, sample_summaries, by = "Sample")

glomero_richness <- read.csv("C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Rarefied/Glomero_richness_rarefy_15000.csv") 
glomero_richness <- inner_join(glomero_richness, sample_summaries, by = "Sample")


# Calculating vegetation richness values and joining with fungi richness data frames

veg_richness <- veg_original %>%
  group_by(subplot_id) %>%
  summarise(Richness = sum(month_2))
colnames(veg_richness) <- c("Subplot", "Veg_richness")
veg_richness$Subplot <- gsub("-O-","-",as.character(veg_richness$Subplot))
write.table(veg_richness, file = "C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Vegetation/Vegetation_richness.csv", sep=",", row.names = FALSE)

veg_fun_rich <- inner_join(fungi_richness, veg_richness, by = "Subplot")
vf_int_rich <- veg_fun_rich %>%
  filter(Soil_origin == "Inter-tussock")

veg_glomero_rich <- inner_join(glomero_richness, veg_richness, by = "Subplot")
vg_int_rich <- veg_glomero_rich %>%
  filter(Soil_origin == "Inter-tussock")


# Making plot of total fungi richness and plant richness

jpeg('C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Vegetation/Richness_vf_int.jpg', width=1900, height=1300, quality=100, res=300, pointsize = 10)                 
ggplot(vf_int_rich, aes(y=Veg_richness, x=Richness, colour=Fire_history)) +
  geom_point(size=3) +
  scale_color_viridis_d(begin=0, end=1) +
  geom_smooth(aes(y=Veg_richness, x=Richness),method=lm, se=F, inherit.aes=F, colour="Black") +
  theme_bw() +
  ylab("Plant species richness") + xlab("Inter-tussock fungal species richness") +
  labs(colour = "Fire history") +
  ggtitle("Correlation between plant and fungal species richness") +
  facet_wrap(~Fire_history) +
  theme(legend.position="none")
dev.off()


# Checking distributions and correlation test

jpeg('C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Vegetation/Plant_histogram.jpg', width=1500, height=1000, quality=100, res=300, pointsize = 10)
hist(vf_int_rich$Veg_richness, xlab="Vascular plant species richness", main = "Histogram of plant richness")
dev.off()
shapiro.test(vf_int_rich$Veg_richness)

jpeg('C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Vegetation/Fungi_histogram.jpg', width=1500, height=1000, quality=100, res=300, pointsize = 10)
hist(vf_int_rich$Richness, xlab="Fungal species richness", main = "Histogram of fungal richness")
dev.off()
shapiro.test(vf_int_rich$Richness)

# CONTROL
fungi.cor.co <- vf_int_rich %>%
  filter(Fire_history == "Control")
cor.test(fungi.cor.co$Veg_richness, fungi.cor.co$Richness, method = "pearson")
# SPRING
fungi.cor.sp <- vf_int_rich %>%
  filter(Fire_history == "Spring")
cor.test(fungi.cor.sp$Veg_richness, fungi.cor.sp$Richness, method = "pearson")
# SUMMER
fungi.cor.su <- vf_int_rich %>%
  filter(Fire_history == "Summer")
cor.test(fungi.cor.su$Veg_richness, fungi.cor.su$Richness, method = "pearson")
# UNBURNT
fungi.cor.ub <- vf_int_rich %>%
  filter(Fire_history == "Unburnt")
cor.test(fungi.cor.ub$Veg_richness, fungi.cor.ub$Richness, method = "pearson")
# OVERALL
cor.test(veg_fun_rich$Veg_richness, veg_fun_rich$Richness, method = "pearson") # Pearson's test is used when both of the variables are normally distributed
(cor(veg_fun_rich$Veg_richness, veg_fun_rich$Richness, method = "pearson"))^2 #R^2


# Making plot of Glomero richness and plant richness

jpeg('C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Vegetation/Richness_vg_int.jpg', width=1900, height=1300, quality=100, res=300, pointsize = 10)                 
ggplot(vg_int_rich, aes(y=Veg_richness, x=Richness, colour=Fire_history)) +
  geom_point(size=3) +
  scale_color_viridis_d(option = "plasma", begin=0, end=1) +
  geom_smooth(aes(y=Veg_richness, x=Richness),method=lm, se=F, inherit.aes=F, colour="Black") +
  theme_bw() +
  ylab("Plant species richness") + xlab("Inter-tussock Glomeromycota species richness") +
  labs(colour = "Fire history") +
  ggtitle("Correlation between plant and Glomeromycota species richness") +
  facet_wrap(~Fire_history) +
  theme(legend.position="none")
dev.off()


# Checking distributions and correlation test

jpeg('C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Vegetation/Plant_histogram.jpg', width=1500, height=1000, quality=100, res=300, pointsize = 10)
hist(vg_int_rich$Veg_richness, xlab="Vascular plant species richness", main = "Histogram of plant richness")
dev.off()
shapiro.test(vg_int_rich$Veg_richness)

jpeg('C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Vegetation/Fungi_histogram.jpg', width=1500, height=1000, quality=100, res=300, pointsize = 10)
hist(vg_int_rich$Richness, xlab="Fungal species richness", main = "Histogram of fungal richness")
dev.off()
shapiro.test(vg_int_rich$Richness)

# CONTROL
glomero.cor.co <- vg_int_rich %>%
  filter(Fire_history == "Control")
cor.test(glomero.cor.co$Veg_richness, glomero.cor.co$Richness, method = "pearson")
# SPRING
glomero.cor.sp <- vg_int_rich %>%
  filter(Fire_history == "Spring")
cor.test(glomero.cor.sp$Veg_richness, glomero.cor.sp$Richness, method = "pearson")
# SUMMER
glomero.cor.su <- vg_int_rich %>%
  filter(Fire_history == "Summer")
cor.test(glomero.cor.su$Veg_richness, glomero.cor.su$Richness, method = "pearson")
# UNBURNT
glomero.cor.ub <- vg_int_rich %>%
  filter(Fire_history == "Unburnt")
cor.test(glomero.cor.ub$Veg_richness, glomero.cor.ub$Richness, method = "pearson")
# OVERALL
glomero.cor.test(glomero_fun_rich$Veg_richness, glomero_fun_rich$Richness, method = "pearson") # Pearson's test is used when both of the variables are normally distributed
(cor(veg_fun_rich$Veg_richness, veg_fun_rich$Richness, method = "pearson"))^2 #R^2
