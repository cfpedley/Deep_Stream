# STEP 1

#STARTING POINT: 
#Illumina-sequenced paired-end fastq files that have been demultiplexed by sample. Barcodes already removed. We used forward and reverse reads of the full fungal ITS region from 72 soil DNA samples. Soil samples collected at Deep Stream, Otago (144 fastq.gz files in total). DNA had been extracted using Qiagen PowerSoil Pro and DNA extracts were sequenced using Illumina MySeq.

#END PRODUCT: 
#ASV table with number of times each ASV was observed in each sample. Taxonomy assigned to output ITS sequence variants using UNITE database. Both of these tables saved as .csv files.

#NOTE: all paths need to be updated to match your device.


#Loading packages

library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")
library(tidyverse)


#Accessing files on machine

path <- "C:/Deep_Stream"
list.files(path)


#Assigning forward and reverse reads

fnFs <- sort(list.files(path, pattern = "R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "R2_001.fastq.gz", full.names = TRUE))
fnFs
fnRs


#Identifying primers

FWD <- "CTTGGTCATTTAGAGGAAGTAA"  ## ITS1F
REV <- "TCCTCCGCTTATTGATATGC"  ## ITS4


# Making set of each primer in all possible orientations

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients
REV.orients


# Make new directory with ambiguous bases filtered (reads with Ns removed)

fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = FALSE)
fnFs.filtN
fnRs.filtN


# PrimerHits table

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

#OUTPUT
#                 Forward Complement Reverse RevComp
#FWD.ForwardReads  207733          0       0       0
#FWD.ReverseReads       0          0       0      21
#REV.ForwardReads       0          0       0      63
#REV.ReverseReads  221659          0       0       0


# Cutadapt

# If version doesn't show up, move to internal C drive and make sure path is set to that location. Make sure to use two backslashes.

cutadapt <- "C:/Deep_Stream/cutadapt.exe"
system2(cutadapt, args = "--version") # Run shell commands from R

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-m", 1, "-n", 2, # -n 2 required to remove FWD and REV from reads, -m means minimum length (a minimum length of 1 gets rid of files that were just primers and would be a length of 0 after trimming)
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}


# Sanity check: making sure there are no primers left in our reads

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))


# Forward and reverse fastq filenames have the format:

cutFs <- sort(list.files(path.cut, pattern = "R1_001.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "R2_001.fastq.gz", full.names = TRUE))


# Extract sample names, assuming filenames have format:

get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
sample.names


# Quality profile plots (numbers in brackets dictate which/how many samples shown)

jpeg("C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/DADA2/Quality profiles FWD 1.jpg", width=2500, height=1800, quality = 100, res=300)
plotQualityProfile(cutFs[1:6]) #Forward reads
dev.off()

jpeg("C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/DADA2/Quality profiles FWD 2.jpg", width=2500, height=1800, quality = 100, res=300)
plotQualityProfile(cutFs[51:56])
dev.off()

plotQualityProfile(cutRs[1:6]) #Reverse reads - I don't really need these ones because I'm just using the forward reads
plotQualityProfile(cutRs[51:56])


# Filter and trim

filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, 
                    maxN = 0, maxEE = c(2, 2), truncQ = 2, minLen = 50, 
                    rm.phix = TRUE, compress = TRUE, multithread = FALSE)
head(out)


# Learn the Error Rates

errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

# Visualise estimated error rates (sanity check)

jpeg("C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/DADA2/Error rates FWD.jpg", width=2500, height=1800, quality = 100, res=300)
plotErrors(errF, nominalQ = TRUE)
dev.off()

plotErrors(errR, nominalQ = TRUE)


# Dereplicate identical reads

derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)


# Name the derep-class objects by the sample names

names(derepFs) <- sample.names
names(derepRs) <- sample.names


# Sample Inference

dadaFs <- dada(derepFs, err = errF, multithread = FALSE)
dadaRs <- dada(derepRs, err = errR, multithread = FALSE)


# Construct Sequence Table

seqtab <- makeSequenceTable(dadaFs)
dim(seqtab)


# Remove chimeras

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE, verbose=TRUE)


# Inspect distribution of sequence lengths

sequence_lengths <- nchar(getSequences(seqtab.nochim))
jpeg("C:/Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/DADA2/Histogram of sequence lengths.jpg", width=1200, height=1000, quality = 100, res=300)
hist(nchar(getSequences(seqtab.nochim)), main = paste("Histogram of sequence lengths"), xlab = "Number of bases", ylab = "Number of ASVs")
dev.off()


# Track how many reads made it through each stage of the pipeline - there should be no step where a majority of reads were lost

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "nonchim")
rownames(track) <- sample.names
track
write.table(track, "C://Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/DADA2/track.csv", sep=",")



#Assigning ASVs to seqIDs fasta file - named "ASV_[number]"

seqIDs <- paste("ASV",(1:length(colnames(seqtab.nochim))), sep = "_")
uniquesToFasta(seqtab.nochim, "C://Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/DADA2/ASV_sequences.fasta", ids = seqIDs)


#ASV table - ASVs named as ASV_number, rows are samples

otutab <- as.data.frame(seqtab.nochim)
colnames(otutab) <- seqIDs
write.table(otutab,file = "C://Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/DADA2/ASV_table.csv", sep=",")


# Assign taxonomy

taxa <- assignTaxonomy(seqtab.nochim, "C:\\Deep_Stream\\sh_general_release_dynamic_10.05.2021.fasta", minBoot= 80, outputBootstraps = TRUE, multithread = FALSE, tryRC = TRUE)

# Tidy up taxonomy assignments table

taxa.print <- as.data.frame(taxa)
taxa.print$tax.Kingdom <- gsub("k__","",as.character(taxa.print$tax.Kingdom))
taxa.print$tax.Phylum <- gsub("p__","",as.character(taxa.print$tax.Phylum))
taxa.print$tax.Class <- gsub("c__","",as.character(taxa.print$tax.Class))
taxa.print$tax.Order <- gsub("o__","",as.character(taxa.print$tax.Order))
taxa.print$tax.Family <- gsub("f__","",as.character(taxa.print$tax.Family))
taxa.print$tax.Genus <- gsub("g__","",as.character(taxa.print$tax.Genus))
taxa.print$tax.Species <- gsub("s__","",as.character(taxa.print$tax.Species))
colnames(taxa.print) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "boot.Kingdom", "boot.Phylum", "boot.Class", "boot.Order", "boot.Family", "boot.Genus", "boot.Species")
head(taxa.print)

# Create taxonomic assignments csv file

table(taxa.print$Kingdom) # Check how many ASVs were produced
taxa.print <- taxa.print %>% 
  mutate(., ASV = colnames(otutab)) # Make a new column with names of ASV
write.table(taxa.print, file = "C://Users/Charlotte/Documents/School/Vic/2022/MBIO489 - research project/Bioinformatics and analysis/DADA2/taxa_table.csv", sep=",")



## FINISHED!!

## Our final outputs should be: 
# a histogram of sequence lengths after trimming and merging
# a couple of graphs of sequence quality scores for 6 samples
# a graph of error rates by quality scores
# a .fasta file with the ASV names and their sequences
# a 'track.csv' file with details about how many reads made it through each stage of processing
# an 'ASV_table.csv' file showing the number of reads of each ASV in each sample
# and a 'taxa_table.csv' with the names of the taxa assigned to each ASV 
