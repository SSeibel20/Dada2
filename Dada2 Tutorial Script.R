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

#learnErrors alternates estimation of the error rates and inference sample composition until they converge
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

#check estimated error rates by visualizing them
plotErrors(errF, nominalQ=TRUE)

#apply core sample inference algorithm
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

#inspect dada-class
#shows true sequence variants from unique sequences in the first sample
dadaFs[[1]]

#merge forward and reverse reads to make full denoised sequences
#align forward with reverse complement of reverse reads and merge into contiq sequences
#by default, merges forward and reverse sequences with 12+ base overlap
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# Inspect the merger data.frame from the first sample
head(mergers[[1]])

#construct amplicon sequence variant table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

#inspect distribution of sequence lengths
#rows correspond to samples
#columns correspond to sequence variants
table(nchar(getSequences(seqtab)))

#to remove nonspecific priming issues to get amplicons of target length
#seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 250:256] 

#remove chimeras
#dada cannot remove chimeras on its own
#can be identified by exactly reconstructing left and right segments from two or more abundant parent sequences

seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

#track reads through the pipeline 
getN <- function(x) sum(getUniques(x))
track <-cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
#if processing a single sample, remove the sapply calls: ex replace sapply(dadaFS, getN) with getN(dadaFs)

#assign taxonomy
#uses naive Bayesian classifier method to assign
#assignTaxonomy function takes a set of sequences and a training set of reference sequences as input
#assignTaxonomy function provides taxonomic assignments with at least minBoot bootstrap confidence

taxa <- assignTaxonomy(seqtab.nochim, "C:\\Users/saman/OneDrive - The Pennsylvania State University/Ganda Lab/Spring 2023/Dada2_Tutorial/MiSeq_SOP/silva_nr_v132_train_set.fa.gz", multithread=TRUE)

#to add species
#taxa <- addSpecies(taxa, "~/tax/silva_species_assignment_v132.fa.gz")

#inspect taxonomic assignments
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

#save data
save.image(file = "Dada2_tutorial.RData")

