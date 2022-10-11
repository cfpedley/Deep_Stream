# STEP 5

# FILES REQUIRED: Rarefied fungal ASV tables (from '2_Filering_and_rarefying.R') and file with descriptions of abiotic soil parameters, fire histories, and soil origins.

# This code aims to analyse soil fungal community composition vs fire history and soil origin. We'll do a PCoA analysis to visualise relationships and then a PERMANOVA to back up any patterns with actual stats.


library(tidyverse)
library(vegan)
library(glue)


sample_summaries <- read.csv("C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Sample_summaries.csv")

# We won't just read in the ASV table, we also need to remove any ASV rows with no reads left after rarefying earlier
ASVs <- read.csv("C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Rarefied/ASV_table_rarefied.csv") %>%
  mutate(Sample = sample_summaries$Sample) %>%
  pivot_longer(-Sample) %>%
  group_by(name) %>%
  mutate(total = sum(value)) %>% 
  filter(total != 0) %>%
  ungroup() %>%
  select(-total) %>%
  pivot_wider(Sample) %>%
  as.data.frame()
# This leaves 7167/7606 ASVs

# Now we need to make the ASV table back into its nice matrix format
rownames(ASVs) <- ASVs$Sample
ASVs <- ASVs[,-1]
ASVs <- as.matrix(ASVs)



## 1. PCoA ordination, bray curtis dissimilarity

dist <- vegdist(ASVs, method="bray") #make distance matrix

eig <- cmdscale(dist, eig=T, add=T)
pcoa <- eig$points
colnames(pcoa) <- c("PC1", "PC2")

# finding % of variation explained by each axis and pulling out x and y axes
percent_exp <- round((100* eig$eig / sum(eig$eig)), digits = 1)
percent_exp[1:2]

# viewing how much adding axes increasing the amount of variation we're explaining
tibble(percent_explained = cumsum(percent_exp),
       axis = 1:length(percent_exp)) %>%
  ggplot(aes(x=axis,y=percent_explained)) +
  geom_line()

labs <- c(glue("PCO1 ({percent_exp[1]}%)"),
          glue("PCO2 ({percent_exp[2]}%)"))
labs
# Values are 13.4% and 6.7%

jpeg('C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Stats - composition/PCoA.jpg', width=2000, height=1500, quality=100, res=300, pointsize = 10)
pcoa %>%
  as_tibble(rownames="Sample") %>%
  inner_join(.,sample_summaries, by="Sample") %>%
  ggplot(aes(x=PC1, y=PC2, shape=Soil_origin, fill=Fire_history)) +
  geom_point(size=3.5) +
  scale_fill_viridis_d(option = "viridis", begin = 0, end = 1) +
  scale_shape_manual(values = c(21,24)) +
  theme_bw() +
  labs(x=labs[1], y=labs[2], fill = "Fire history", shape = "Soil origin") +
  guides(fill=guide_legend(override.aes = list(shape=21))) 
dev.off()


## 2. PERMANOVA

perm <- adonis2(dist ~ Fire_history * Soil_origin, strata = sample_summaries$Plot, data = sample_summaries)
perm


### REPEATING WITH GLOMEROMYCOTA ASVs


glomero_ASVs <- read.csv("C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Rarefied/Glomero_table_rarefied.csv") %>%
  mutate(Sample = sample_summaries$Sample) %>%
  pivot_longer(-Sample) %>%
  group_by(name) %>%
  mutate(total = sum(value)) %>%
  filter(total != 0) %>%
  ungroup() %>%
  select(-total) %>%
  pivot_wider(Sample) %>%
  as.data.frame()

glomero_ASVs <- read.csv("C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Rarefied/Glomero_table_rarefied.csv") %>%
  mutate(Sample = sample_summaries$Sample) %>%
  pivot_longer(-Sample) %>%
  group_by(name) %>%
  mutate(total = sum(value)) %>%
  filter(total != 0) %>%
  ungroup() %>%
  select(-total) %>%
  pivot_wider(Sample) %>%
  as.data.frame()

# Now we need to make the ASV table back into its nice matrix format
rownames(glomero_ASVs) <- glomero_ASVs$Sample
glomero_ASVs <- glomero_ASVs[,-1]
glomero_ASVs <- as.matrix(glomero_ASVs)



## 1. PCoA ordination, bray curtis dissimilarity

glomero_dist <- vegdist(glomero_ASVs, method="bray") #make distance matrix

glomero_eig <- cmdscale(glomero_dist, eig=T, add=T)
glomero_pcoa <- glomero_eig$points
colnames(glomero_pcoa) <- c("PC1", "PC2")

# finding % of variation explained by each axis and pulling out x and y axes
g_percent_exp <- round((100* glomero_eig$eig / sum(glomero_eig$eig)), digits = 1)
g_percent_exp[1:2]

# viewing how much adding axes increasing the amount of variation we're explaining
tibble(percent_explained = cumsum(g_percent_exp),
       axis = 1:length(g_percent_exp)) %>%
  ggplot(aes(x=axis,y=percent_explained)) +
  geom_line()

glomero_labs <- c(glue("PCO1 ({g_percent_exp[1]}%)"),
          glue("PCO2 ({g_percent_exp[2]}%)"))
glomero_labs
# Values are 15.3% and 8.4%

jpeg('C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Stats - composition/Glomero_PCoA.jpg', width=2000, height=1500, quality=100, res=300, pointsize = 10)
glomero_pcoa %>%
  as_tibble(rownames="Sample") %>%
  inner_join(.,sample_summaries, by="Sample") %>%
  ggplot(aes(x=PC1, y=PC2, shape=Soil_origin, fill=Fire_history)) +
  geom_point(size=3.5) +
  scale_fill_viridis_d(option = "plasma", begin = 0, end = 1) +
  scale_shape_manual(values = c(21,24)) +
  theme_bw() +
  labs(x=glomero_labs[1], y=glomero_labs[2], fill = "Fire history", shape = "Soil origin") +
  guides(fill=guide_legend(override.aes = list(shape=21))) +
  ggtitle("Glomeromycota PCoA")
dev.off()


## 2. PERMANOVA

glomero_perm <- adonis2(glomero_dist ~ Fire_history * Soil_origin, strata = sample_summaries$Plot, data = sample_summaries)
glomero_perm

