#----- Processing 16S data with DADA2 and Phyloseq  ----#
# See https://benjjneb.github.io/dada2/tutorial_1_8.html

####Libraries####
version()
library(dada2)
library(phyloseq)
library(vegan)
library(tidyverse)
library(reshape2)
library(stringr)
library(data.table)
library(broom)
library(ape)
library(qualpalr)
library(viridis)
library(Biostrings)
library(ShortRead)
library(seqinr)

####Environment Setup####
theme_set(theme_bw())
setwd("~/project/")

setwd("~/Desktop/CG-FT_sequences_ONLY/")

####File Path Setup####
#this is so dada2 can quickly iterate through all the R1 and R2 files in your read set
path <- "~/Desktop/CG-FT_sequences_ONLY/" # CHANGE ME to the directory containing the fastq files
list.files(path)
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE)) #change the pattern to match all your R1 files
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1) #change the delimiter in quotes and the number at the end of this command to decide how to split up the file name, and which element to extract for a unique sample name

#### fastq quality plots ####
plotQualityProfile(fnFs[100:110]) #this plots the quality profiles for each sample, if you have a lot of samples, it's best to look at just a few of them, the plots take a minute or two to generate even only showing 10-20 samples.
plotQualityProfile(fnRs[1:20])

#### Identify primers ####
FWD <- "GTGYCAGCMGCCGCGGTAA"  ## CHANGE to forward primer sequence
REV <- "CCGYCAATTYMTTTRAGTTT"  ## CHANGE to reverse primer sequence
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

fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))

# The presence of ambiguous bases (Ns) in the sequencing reads makes accurate mapping of short primer sequences difficult. 
# Next we are going to “pre-filter” the sequences just to remove those with Ns, but perform no other filtering.
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE, compress = TRUE)

# Count the number of times the primers appear in the forward and reverse read, 
# while considering all possible primer orientations. Here we do it for one example sample (the number in [[]] to get an idea.
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found 
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[61]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[61]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[61]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[61]]))

#### Example output ####
#                   Forward Complement Reverse RevComp
# FWD.ForwardReads   24145          0       0       0
# FWD.ReverseReads       0          0       0       0
# REV.ForwardReads       0          0       0       0
# REV.ReverseReads   24631          0       0       0
# FWD primer is found in the forward reads in its forward orientation, and in some of the reverse reads in its reverse-complement orientation 
# (due to read-through when the ITS region is short). 
# Similarly the REV primer is found with its expected orientations

#REV <- REV.orients[["RevComp"]] #IMPORTANT!!! change orientation ONLY IF NECESSARY. see the dada2 ITS_workflow guide section "Identify Primers" for details. it is linked at the top of this guide.



#### Use cutadapt to remove primer sequences ####
# Here we use cutadapt in qiime2
# The following is done in the terminal:
cutadapt <- "/usr/local/bin/cutadapt" # CHANGE ME to the cutadapt path on your machine
system2(cutadapt, args = "--version")

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

#Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

#sanity check, should report zero for all orientations and read sets
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[121]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[121]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[121]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[121]]))

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "R1", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "R2", full.names = TRUE))

####filter and trim reads####
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

####trim & filter####
#filter and trim command. dada2 can canonically handle lots of errors, I am typically permissive in the maxEE parameter set here, in order to retain the maximum number of reads possible. error correction steps built into the dada2 pipeline have no trouble handling data with this many expected errors.
#it is best, after primer removal, to not truncate with 18s data, or with data from any region in which the length is broadly variable. you may exclude organisms that have a shorter insert than the truncation length (definitely possible, good example is giardia). defining a minimum sequence length is best.
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, truncLen=c(0,0), minLen = c(150,120),
                     maxN=c(0,0), maxEE=c(8,10), truncQ=c(2,2), rm.phix=TRUE, matchIDs=TRUE,
                     compress=TRUE, multithread=TRUE)
retained <- as.data.frame(out)
retained$percentage_retained <- retained$reads.out/retained$reads.in*100
View(retained)

####learn error rates####
#the next three sections (learn error rates, dereplication, sample inference) are the core of dada2's sequence processing pipeline. read the dada2 paper and their online documentation (linked at top of this guide) for more information on how these steps work
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE) #assess this graph. it shows the error rates observed in your dataset. strange or unexpected shapes in the plot should be considered before moving on.

####dereplication####
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names #this is just to ensure that all your R objects have the same sample names in them
names(derepFs) <- sample.names
names(derepRs) <- sample.names

####sample inference####
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

dadaFs[[1]]
dadaRs[[1]]

####OPTIONAL: remove low-sequence samples before merging####
#a "subscript out of bounds" error at the next step (merging) may indicate that you aren't merging any reads in one or more samples.
#NB, NOT getting this error doesn't necessarily mean that all of your samples end up with more than 0 merged reads, as i found out while processing a large 18s dataset. your guess is as good as mine as to why this error does or does not appear, but filtering out the samples that cause it is necessary for completion of the pipeline.
#samples_to_keep <- as.numeric(out[,"reads.out"]) > 100 #example of simple method used above after the filter and trim step. if you already did this but still got an error when merging, try the steps below
getN <- function(x) sum(getUniques(x)) #keeping track of read retention, number of unique sequences after ASV inference
track <- cbind(sapply(derepFs, getN), sapply(derepRs, getN), sapply(dadaFs, getN), sapply(dadaRs, getN))
samples_to_keep <- track[,4] > 50 #your threshold. try different ones to get the lowest one that will work. #this method accounts for dereplication/ASVs left after inference
samples_to_remove <- names(samples_to_keep)[which(samples_to_keep == FALSE)] #record names of samples you have the option of removing

####merge paired reads####
mergers <- mergePairs(dadaFs[samples_to_keep], derepFs[samples_to_keep], dadaRs[samples_to_keep], derepRs[samples_to_keep], verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

####construct sequence table####
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

####View Sequence Length Distribution Post-Merging####
#most useful with merged data. this plot will not show you much for forward reads only, which should have a uniform length distribution.
length.histogram <- as.data.frame(table(nchar(getSequences(seqtab)))) #tabulate sequence length distribution
plot(x=length.histogram[,1], y=length.histogram[,2]) #view length distribution plot

####remove low-count singleton ASVs####
#create phyloseq otu_table
otus <- otu_table(t(seqtab), taxa_are_rows = TRUE)

#some metrics from the sequence table
otu_pres_abs <- otus
otu_pres_abs[otu_pres_abs >= 1] <- 1 #creating a presence/absence table
otu_pres_abs_rowsums <- rowSums(otu_pres_abs) #counts of sample per ASV
length(otu_pres_abs_rowsums) #how many ASVs
length(which(otu_pres_abs_rowsums == 1)) #how many ASVs only present in one sample

#what are the counts of each ASV
otu_rowsums <- rowSums(otus) #raw counts per ASV
otu_singleton_rowsums <- as.data.frame(otu_rowsums[which(otu_pres_abs_rowsums == 1)]) #raw read counts in ASVs only presesnt in one sample
hist(otu_singleton_rowsums[,1], breaks=500, xlim = c(0,200), xlab="# Reads in ASV") #histogram plot of above
length(which(otu_singleton_rowsums <= 50)) #how many ASVs are there with N reads or fewer? (N=50 in example)

#IF you want to filter out rare variants (low-read-count singleton ASVs) you can use phyloseq's "transform_sample_counts" to create a relative abundance table, and then filter your ASVs by choosing a threshold of relative abundance: otus_rel_ab = transform_sample_counts(otus, function(x) x/sum(x))
dim(seqtab) # sanity check
dim(otus) # (this should be the same as last command, but the dimensions reversed)
otus_rel_ab <- transform_sample_counts(otus, function(x) x/sum(x)) #create relative abundance table
df <- as.data.frame(unclass(otus_rel_ab)) #convert to plain data frame
df[is.na(df)] <- 0 #if there are samples with no merged reads in them, and they passed the merge step (a possiblity, converting to a relative abundance table produes all NaNs for that sample. these need to be set to zero so we can do the calculations in the next steps.)
otus_rel_ab.rowsums <- rowSums(df) #compute row sums (sum of relative abundances per ASV. for those only present in one sample, this is a value we can use to filter them for relative abundance on a per-sample basis)
a <- which(as.data.frame(otu_pres_abs_rowsums) == 1) #which ASVs are only present in one sample
b <- which(otus_rel_ab.rowsums <= 0.001) #here is where you set your relative abundance threshold #which ASVs pass our filter for relative abundance
length(intersect(a,b)) #how many of our singleton ASVs fail on this filter
rows_to_remove <- intersect(a,b) #A also in B (we remove singleton ASVs that have a lower relative abundance value than our threshold)
otus_filt <- otus[-rows_to_remove,] #filter OTU table we created earlier
dim(otus_filt) #how many ASVs did you retain?
seqtab.nosingletons <- t(as.matrix(unclass(otus_filt))) #convert filtered OTU table back to a sequence table matrix to continue with dada2 pipeline

####remove chimeras####
#here we remove "bimeras" or chimeras with two sources. look at "method" to decide which type of pooling you'd like to use when judging each sequence as chimeric or non-chimeric
seqtab.nosingletons.nochim <- removeBimeraDenovo(seqtab.nosingletons, method="pooled", multithread=TRUE, verbose=TRUE) #this step can take a few minutes to a few hours, depending on the size of your dataset
#Identified 2894 bimeras out of 12903 input sequences.
dim(seqtab.nosingletons.nochim)
#[1]   478 10009
sum(seqtab.nosingletons.nochim)/sum(seqtab.nosingletons) #proportion of nonchimeras #it should be relatively high after filtering out your singletons/low-count ASVs, even if you lose a lot of ASVs, the number of reads lost should be quite low

####track read retention through steps####
getN <- function(x) sum(getUniques(x))
track <- cbind(out[samples_to_keep,], sapply(dadaFs[samples_to_keep], getN), sapply(dadaRs[samples_to_keep], getN), sapply(mergers, getN), rowSums(seqtab.nosingletons), rowSums(seqtab.nosingletons.nochim))
# If processing only a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
track <- cbind(track, 100-track[,6]/track[,5]*100, 100-track[,7]/track[,6]*100)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nosingletons", "nochimeras", "percent_singletons", "percent_chimeras")
rownames(track) <- sample.names[samples_to_keep]

####save output from sequnce table construction steps####
write.table(data.frame("row_names"=rownames(track),track),"read_retention.16s.txt", row.names=FALSE, quote=F, sep="\t")
write.table(data.frame("row_names"=rownames(seqtab.nosingletons.nochim),seqtab.nosingletons.nochim),"sequence_table.16s.merged.txt", row.names=FALSE, quote=F, sep="\t")

#if you must save your sequence table and load it back in before doing taxonomy assignments, here is how to reformat the object so that dada2 will accept it again
seqtab.nosingletons.nochim <- fread("sequence_table.16s.merged.txt", sep="\t", header=T, colClasses = c("row_names"="character"), data.table=FALSE)
row.names(seqtab.nosingletons.nochim) <- seqtab.nosingletons.nochim[,1] #set row names
seqtab.nosingletons.nochim <- seqtab.nosingletons.nochim[,-1] #remove column with the row names in it
seqtab.nosingletons.nochim <- as.matrix(seqtab.nosingletons.nochim) #cast the object as a matrix

####assign taxonomy####
#note, this takes ages if you have a large dataset. saving the sequences as a fasta file (with writeFasta) and using QIIME's taxonomy assignment command will save you time, and is only slightly less accurate than the dada2 package's taxonomy assignment function.
taxa <- assignTaxonomy(seqtab.nosingletons.nochim, "~/Desktop/lab_member_files/taxonomy_databases/silva_for_dada2/v132_for_parfreylab/16s/silva_132.16s.99_rep_set.dada2.fa.gz", multithread=TRUE)
taxa <- readRDS("taxonomy_table_fucus.RDS")

#NA taxa are hard to separate later if they have no label. apply "Unassigned" label here now.
unique(taxa[,1]) #possible labels here: eukaryotic, archaeal, bacterial, and "NA" taxa. 
NAs <- is.na(taxa[,1]) #test for NA
NAs <- which(NAs == TRUE) #get indices of NA values
taxa[NAs,1] <- "Unassigned" #apply new label to identified indices
colnames(taxa) <- c("Rank1", "Rank2", "Rank3", "Rank4", "Rank5", "Rank6", "Rank7", "Accession")

#create phyloseq object from "seqtab.nosingletons.nochim", "metadata", and "taxa"
ps.dada2 <- phyloseq(otu_table(seqtab.nosingletons.nochim, taxa_are_rows=FALSE), 
                          sample_data(metadata), 
                          tax_table(taxa))


# Update uncultured/unknown/NA taxonomy assignments
cg.ft.tax <- as.data.frame(unclass(tax_table(ps.dada2)))
cg.ft.tax$ASV <- rownames(cg.ft.tax)
cg.ft.tax <- cg.ft.tax %>% mutate_if(is.factor,as.character)
colnames(cg.ft.tax)
# "Domain"    "Phylum"    "Class"     "Order"     "Family"    "Genus"     "Species"   "Accession"

cg.ft.tax$Class <- ifelse(is.na(cg.ft.tax$Class), ifelse(str_detect(cg.ft.tax$Phylum, "unknown"), cg.ft.tax$Phylum, paste("unknown", cg.ft.tax$Phylum, sep="_")), cg.ft.tax$Class)
cg.ft.tax$Order <- ifelse(is.na(cg.ft.tax$Order), ifelse(str_detect(cg.ft.tax$Class, "unknown"), cg.ft.tax$Class, paste("unknown", cg.ft.tax$Class, sep="_")), cg.ft.tax$Order)
cg.ft.tax$Family <- ifelse(is.na(cg.ft.tax$Family), ifelse(str_detect(cg.ft.tax$Order, "unknown"), cg.ft.tax$Order, paste("unknown", cg.ft.tax$Order, sep="_")), cg.ft.tax$Family)
cg.ft.tax$Genus <- ifelse(is.na(cg.ft.tax$Genus), ifelse(str_detect(cg.ft.tax$Family, "unknown"), cg.ft.tax$Family, paste("unknown", cg.ft.tax$Family, sep="_")), cg.ft.tax$Genus)
cg.ft.tax$Species <- ifelse(is.na(cg.ft.tax$Species), ifelse(str_detect(cg.ft.tax$Genus, "unknown"), cg.ft.tax$Genus, paste("unknown", cg.ft.tax$Genus, sep="_")), cg.ft.tax$Species)

cg.ft.tax$Class <- ifelse(str_detect(cg.ft.tax$Class, "uncultured"), ifelse(str_detect(cg.ft.tax$Phylum, "unknown"), cg.ft.tax$Phylum, paste("unknown", cg.ft.tax$Phylum, sep="_")), cg.ft.tax$Class)
cg.ft.tax$Order <- ifelse(str_detect(cg.ft.tax$Order, "uncultured"), ifelse(str_detect(cg.ft.tax$Class, "unknown"), cg.ft.tax$Class, paste("unknown", cg.ft.tax$Class, sep="_")), cg.ft.tax$Order)
cg.ft.tax$Family <- ifelse(str_detect(cg.ft.tax$Family, "uncultured"), ifelse(str_detect(cg.ft.tax$Order, "unknown"), cg.ft.tax$Order, paste("unknown", cg.ft.tax$Order, sep="_")), cg.ft.tax$Family)
cg.ft.tax$Genus <- ifelse(str_detect(cg.ft.tax$Genus, "uncultured"), ifelse(str_detect(cg.ft.tax$Family, "unknown"), cg.ft.tax$Family, paste("unknown", cg.ft.tax$Family, sep="_")), cg.ft.tax$Genus)
cg.ft.tax$Species <- ifelse(str_detect(cg.ft.tax$Species, "uncultured"), ifelse(str_detect(cg.ft.tax$Genus, "unknown"), cg.ft.tax$Genus, paste("unknown", cg.ft.tax$Genus, sep="_")), cg.ft.tax$Species)

# Set rownames
rownames(cg.ft.tax) <- cg.ft.tax$ASV
cg.ft.tax.mat <- as.matrix.data.frame(cg.ft.tax, rownames.force = T)

# Re-create taxa table 
cg.ft.TAX <- tax_table(cg.ft.tax.mat)
cg.ft.OTU <- otu_table(ps.dada2)

# Read in metadata
cg.ft.meta <- read.csv(file = "~/Desktop/Desktop2020/CG_FT/Data/CG_FT_rock_fucus_combined_phyloseq_r1500_metadata_Jan2021.csv")
cg.ft.meta <- cg.ft.meta %>% select(- 'X') #remove unneeded column
cg.ft.meta <- cg.ft.meta %>% mutate_if(is.factor, as.character) # convert from factors
rownames(cg.ft.meta) <- cg.ft.meta$sampleid
cg.ft.SAM <- sample_data(cg.ft.meta)

# Make updated phyloseq object
cg.ft.fr.new <- phyloseq(cg.ft.TAX, cg.ft.OTU, cg.ft.SAM)
saveRDS(cg.ft.fr.new, file = "~/Desktop/Desktop2020/CG_FT/Data/CG_FT_combined_phyloseq_r1500_Jan2021.RDS")

