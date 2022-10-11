# STEP 2

# FILES REQUIRED: ASV table and taxa table (.csv files) produced via R script '1_DADA2_final_protocol'

## The aim of this code is to create a rarefaction curve to decide on a value for rarefying our data. Then to rarefy the richness and create a normalised ASV table. Then to re-do these steps for only the Glomeromycota species.

library(tidyverse)
library(vegan)
library(phyloseq)
library(ggplot2)
library(GUniFrac)

## Reading in raw file
ASV_original <- read.csv("C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/DADA2/ASV_table.csv")


## Making quick rarefaction curve to check number of reads to rarefy to based on where the curve flattens (in this case, 15000 reads)

rarecurve(ASV_original, step=20, 
          xlab="Number of reads", ylab="Number of ASVs", label=F) # Quick check



## RAREFYING THE DATA TO 15000 READS

# Inspecting number of reads per sample
#Read data from csv file and group by Sample
read_count <- ASV_original %>%
  mutate(Sample=row.names(ASV_original))
read_count <- pivot_longer(read_count, !Sample, names_to = "ASV", values_to = "Count")
read_count <- group_by(read_count,Sample)

num_reads <- read_count %>% 
  summarise(total_reads = sum(Count))
reads_vector <- num_reads$total_reads


# Save histogram of reads per sample
jpeg('C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Rarefied/Histogram_reads_per_sample.jpg', width=1800, height=1200, quality = 100, res=300, pointsize = 10)
ggplot() + 
  theme(legend.position="NULL") +
  geom_histogram(aes(x=reads_vector, fill="Sample", alpha=0.2), binwidth = 4000) +
  scale_fill_viridis_d(option="viridis", begin=0.7, end=0.4) +
  ylab("Number of samples") + xlab("Total reads") + 
  guides(alpha = "none") +
  guides(fill = "none") +
  ggtitle("Histogram of reads per sample") +
  theme_bw()
dev.off()


# Removing samples with <15000 reads
ASV_filt <- ASV_original %>%
  mutate(Sample=rownames(ASV_original)) %>%
  pivot_longer(., !Sample, names_to = "ASV", values_to = "Count") %>%
  group_by(Sample) %>%
  mutate(Reads = sum(Count)) %>%
  pivot_wider(.,names_from = ASV, values_from = Count) %>%
  filter(Reads > 15000)
write.table(ASV_filt, file = "C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Rarefied/ASV_table_15000.csv", sep=",", row.names = F)

#Rarefy richness per sample to 15000
ASV_to_rarefy <- read.csv("C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Rarefied/ASV_table_15000.csv", row.names = 1)
ASV_to_rarefy <- subset(ASV_to_rarefy, select=c(-Reads))
rarefy_15000 <- rarefy(ASV_to_rarefy, 15000) %>%
  as_tibble(rownames = "Sample")
rarefy_15000 <- rename(rarefy_15000, Richness = value)
write.table(rarefy_15000, file = "C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Rarefied/Richness_rarefy_15000.csv", sep=",", row.names = F)

#Rarefy full ASV table
is.data.frame(ASV_to_rarefy)
ASV_rarefied <- GUniFrac::Rarefy(ASV_to_rarefy, 15000)
ASV_rarefied_df <- do.call(rbind.data.frame, ASV_rarefied)
rownames(ASV_rarefied_df) <- gsub("otu.tab.rff.","",
                                  as.character(rownames(ASV_rarefied_df)))
write.table(ASV_rarefied_df, file = "C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Rarefied/ASV_table_rarefied.csv", sep=",")



## FINAL RAREFACTION CURVE WITH LINE AT 15000

rarecurve.df <- rarecurve(ASV_original, tidy = T) #put rarefaction values into data frame

jpeg("C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Rarefied/Rarefaction_curve_fwd_15000.jpg", width=1500, height=1000, quality = 100, res=300)
ggplot(rarecurve.df, aes(x=Sample, y=Species, group=Site)) + 
  xlim(0,60000) +
  geom_line(size=0.3) +  theme_classic() + 
  geom_vline(xintercept = 15000, color = "Red") + 
  ylab("Number of species (ASVs)") + xlab("Number of reads") + 
  ggtitle("Rarefaction curve of Deep Stream 2020 forward reads")
dev.off()



## MAKE HISTOGRAM OF ASVs PER SAMPLE

# Get richness values as vector from csv
richness_hist <- read.csv("C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Rarefied/richness_rarefy_15000.csv")
richness_vector <- richness_hist$Richness

# Save histogram of ASVs per sample
jpeg('C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Rarefied/Histogram_rarefied_ASVs.jpg', width=1800, height=1200, quality = 100, res=300, pointsize = 10)
ggplot() + 
  theme(legend.position="NULL") +
  geom_histogram(aes(x=richness_vector, fill="2020", alpha=0.2), binwidth = 50) +
  scale_fill_viridis_d(option="viridis", begin=0.7, end=0.4) +
  ylab("Number of samples") + xlab("Total ASVs") + 
  guides(alpha = "none") +
  guides(fill = "none") +
  ggtitle("Histogram of ASVs per sample") +
  theme_bw()
dev.off()


## CALCULATE GLOMERO RICHNESS

taxa <- read.csv("C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/DADA2/taxa_table.csv")
Glomero_ASV_table <- ASV_rarefied_df %>%
  mutate(Sample=rownames(ASV_rarefied_df)) %>%
  pivot_longer(., !Sample, names_to = "ASV", values_to = "Count")
Glomero_ASV_table <- inner_join(Glomero_ASV_table, taxa,  by = "ASV")
Glomero_ASV_table <- filter(Glomero_ASV_table, Phylum == "Glomeromycota")
Glomero_ASV_table <- subset(Glomero_ASV_table, select=c(Sample, ASV, Count))
Glomero_ASV_table <- pivot_wider(Glomero_ASV_table, names_from = ASV, values_from = Count)
write.table(Glomero_ASV_table, file = "C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Rarefied/Glomero_table_rarefied.csv", sep=",", row.names = F)

Glomero_rarefy_15000 <- read.csv("C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Rarefied/Glomero_table_rarefied.csv", row.names = 1)
Glomero_rarefy_15000 <- rarefy(Glomero_rarefy_15000, 15000) %>%
  as_tibble(rownames = "Sample")
Glomero_rarefy_15000 <- rename(Glomero_rarefy_15000, Richness = value)
write.table(Glomero_rarefy_15000, file = "C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/Rarefied/Glomero_richness_rarefy_15000.csv", sep=",", row.names = F)
