# STEP 4.2

# FILES REQUIRED: Rarefied fungal richness tables (from '2_Filering_and_rarefying.R') and file with descriptions of abiotic soil parameters, fire histories, and soil origins.

# lmer makes a linear mixed effects model, which finds the relationships between your predictors and your response, while taking into account random effects (in this case the plot - we expect similar richness values for samples within the same plot). Then you can use an ANOVA to see if the model finds significant relationships between predictors and the response.
# Use emmeans package for post hoc tests. The pairs function finds differences between groups - e.g., is there a difference between control and summer burn plots? emmeans can also give you 95% CIs of these pairwise comparisons - if the range contains 0, you cannot reject the null hypothesis.
# lm can be used for linear models, where there are no random effects (e.g., fire effect vs pH, where both variables are already at a plot level)


library(lme4)
library(lmerTest)
library(pbkrtest)
library(emmeans)


richness <- read.csv("C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Rarefied/Richness_rarefy_15000.csv")
sample_summaries <- read.csv("C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Sample_summaries.csv")
data <- inner_join(sample_summaries, richness, by = "Sample")

glomero_richness <- read.csv("C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Rarefied/Glomero_richness_rarefy_15000.csv")
glomero_data <- inner_join(sample_summaries, glomero_richness, by = "Sample")


## 1. Total Richness Mixed Models

# RESPONSE: richness
# PREDICTORS: fire_history + soil_origin
# RANDOM EFFECT: plot

total_rich_model = lmer(Richness ~ Fire_history + Soil_origin + (1|Plot), data = data)
summary(total_rich_model)
anova(total_rich_model) # ANOVA to see if model finds significant effects of predictors

# Fire_history pairs
fire_emm = emmeans(total_rich_model, spec='Fire_history')  # calculate estimated population means of fire histories
fire_emm
pairs(fire_emm) # test whether statistically significant difference between fire histories
confint(pairs(fire_emm)) # shows 95% confidence intervals - if this range doesn't contain 0, we can reject the null hypothesis
jpeg('C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Stats - richness/Richness_fire_CI.jpg', width=1400, height=1100, quality=100, res=300, pointsize = 10)
plot(pairs(fire_emm), xlab = 'Difference in fungal species richness', ylab = 'Fire history pair') # plot confidence intervals
dev.off()

# Soil_origin pairs
soil_emm = emmeans(total_rich_model, spec='Soil_origin')
soil_emm
pairs(soil_emm)
confint(pairs(soil_emm))



## 2. Glomeromycota Richness Mixed Models

# RESPONSE: glomeromycota_richness
# PREDICTORS: fire_history + soil_origin
# RANDOM EFFECT: plot

glomero_rich_model = lmer(Richness ~ Fire_history + Soil_origin + (1|Plot), data = glomero_data)
summary(glomero_rich_model)
anova(glomero_rich_model)

# Fire_history pairs
glomero_fire_emm = emmeans(glomero_rich_model, spec='Fire_history')
glomero_fire_emm
pairs(glomero_fire_emm)
confint(pairs(glomero_fire_emm))
jpeg('C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Stats - richness/Glomero_fire_CI.jpg', width=1400, height=1100, quality=100, res=300, pointsize = 10)
plot(pairs(glomero_fire_emm), xlab = 'Difference in Glomeromycota richness', ylab = 'Fire history pair') # plot confidence intervals
dev.off()

# Soil_origin pairs
glomero_soil_emm = emmeans(glomero_rich_model, spec='Soil_origin')
glomero_soil_emm
pairs(glomero_soil_emm)
confint(pairs(glomero_soil_emm))


# RESPONSE: Glomeromycota richness
# PREDICTORS: P logged
# RANDOM EFFECT: plot

glomero_P_model = lmer(Richness ~ P_logged + (1|Plot), data = glomero_data)
summary(glomero_P_model)
anova(glomero_P_model)

