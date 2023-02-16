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
#split at first underscore and only use whats before that for the name
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#visualize quality profiles of forward reads
plotQualityProfile(fnFs[1:2])
#gray-scale is heat map of frequency of quality score at each base position
#green line shows mean quality score
#orange line shows quartiles of quality score distribution
#red line shows scaled proportion of reads that extend to that position

#visualize quality profiles of reverse reads
plotQualityProfile(fnRs[1:2])
#reverse reads typically look worse than forward; possibly due to being the second run of lasers
#not an issue for Dada2 as it factors quality info into its error model

#If you are using a less-overlapping primer set, like V1-V2 or V3-V4, your truncLen() must be large enough to maintain 20 + biological.length.variation nucleotides of overlap between them.

#assign filenames for filtered.fastq.gz files
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#use standard filtering parameters
#maxN=0 as Dada2 requires no N's
#maxEE sets max number of expected errors allowed in a read
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(out)
#if not enough reads are passing this filter relax the maxEE (esp reverse reads, maxEE=c(2,5))
#for Internal Transcribed Spacers seq, don't use trunc length


#save data
save.image(file = "Dada2_tutorial.RData")

