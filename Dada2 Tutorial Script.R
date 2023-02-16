#Dada2 Tutorial
#Good luck me

#load the dada2 package
library(dada2); packageVersion("dada2")

#path to directory containing unzipped fastq files
path <- "C:\\Users/saman/OneDrive - The Pennsylvania State University/Ganda Lab/Spring 2023/Dada2_Tutorial/MiSeq_SOP/"
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#visualize quality profiles of forward reads
plotQualityProfile(fnFs[1:2])

#visualize quality profiles of reverse reads
plotQualityProfile(fnRs[1:2])

#save data
save.image(file = "Dada2_tutorial.RData")
