###### DADA2 tutorial ("https://benjjneb.github.io/dada2/tutorial.html")
###### This pipeline works for dataset stored in different path than the script

##################### Parameters for trim, trunc, and fastq format ##########################################
###### Sibling matched study: trimLeft = 10, fastq format: R1(forward) and R2(reverse)
###### PRJNA578223: truncLen = c(290,220), fastq format: 1(forward) and 2(reverse)
###### PRJNA589343: truncLen=c(250), fastq format: single fastq file
#############################################################################################################

################################## Install_Libraries Section ################################################
## If you don't have installed the necessary libraries (BiocManager, dada2, DECIPHER), 
## Uncomment the following three lines:
#chooseCRANmirror()
#install.packages("BiocManager")
#BiocManager::install(c("dada2", "phyloseq","DECIPHER","this.path"))
############################################################################################################

## Loadin libraries
library(DECIPHER); packageVersion("DECIPHER") #IDTaxa algorithm
library(dada2); packageVersion("dada2")
library(ggplot2); packageVersion("ggplot2")
library(stringr); packageVersion("stringr")
library(Biostrings); packageVersion("Biostrings")

## Choosing paths to get and set data
Dataset_dir <- choose.dir(default = "", caption = "Select dataset location")
IDTaxa_Directory <- choose.files(default = "", caption = "Select IDTaxa")
GetData_Path <- choose.dir(default = "", caption = "Select fastq files location")
SetData_Path <- choose.dir(default = "", caption = "Set filtered files (gzip) storage")
## Result storage and create folders 
Graphic_dir <- file.path(paste(Dataset_dir,"/results/",str_replace_all(Sys.time(),':','-'), sep = ""))
Data_dir <- file.path(paste(Graphic_dir,"/data", sep =""))
No_Chim_dir <- file.path(paste(Data_dir,"/No_chim", sep =""))
Chim_dir <- file.path(paste(Data_dir,"/Chim", sep =""))
dir.create(paste(Dataset_dir,"/results",sep=""))
dir.create(Graphic_dir)
dir.create(Data_dir)
dir.create(Chim_dir)
dir.create(No_Chim_dir)
list.files(Dataset_dir)

#Load Files, check fastq filenames format
fnFs <- sort(list.files(GetData_Path, pattern=".fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
#Save the reading order of the fastq in a csv file to compare it with the metadata and to generate labels file
write.table(fnFs, paste(Graphic_dir,"/fastqfiles.csv", sep=""))

#Quality Profile 
plotQualityProfile(fnFs[1:2])
ggsave(filename = "Quality.png", path = Graphic_dir)

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(SetData_Path, paste0(sample.names, "_filt.fastq.gz"))
names(filtFs) <- sample.names

#Actual filtering (truncLen based on quality profiles, truncate to where median quality drops below 30
out <- filterAndTrim(fnFs, filtFs, truncLen=c(250),
                     maxN=0, maxEE=c(2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(out)
write.table(out, paste(Graphic_dir,"/Filtering.csv", sep=""))

## You can plot and/or save the filtered data, if you consider it necessary,  just uncomment the following four lines
plotQualityProfile(filtFs[1:2])
ggsave(filename = "Quality_Filtered.png", path = Graphic_dir)

#Error rate learning
errF <- learnErrors(filtFs, multithread=TRUE)
write.table(errF, paste(Graphic_dir,"/Reads_Errors.csv", sep=""))

## Ploting errors
plotErrors(errF, nominalQ=TRUE)
ggsave(filename = "Reads_Errors.png", path = Graphic_dir)

#dada inference time
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)

#Inspect object
dadaFs[[1]]

## Abundance table
seqtab <- makeSequenceTable(dadaFs)
dim(seqtab)
## Save the file with all information
write.table(seqtab, paste(Graphic_dir,"/seqtab_pre_removal.csv", sep=""))
## Separate only the data
write.table(seqtab, paste(Chim_dir,"/data_0_chim.csv", sep=""), row.names=FALSE, col.names = FALSE, sep= ",")
#Feature list
Features.chim <- data.frame(Sequence = getSequences(seqtab))
#Write to file
write.table(Features.chim,paste(Chim_dir,"/features_0.csv", sep=""), row.names=FALSE, col.names = FALSE, sep= ",")

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

#chimera filter
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
#non-chimera proportion
sum(seqtab.nochim)/sum(seqtab)
## Save abundance table with chim-removal
## Save the file with all information
write.table(seqtab.nochim, paste(Graphic_dir,"/seqtab_post_removal.csv", sep=""))
## Separate only the data
write.table(seqtab.nochim, paste(No_Chim_dir,"/data_0_Nochim.csv", sep=""), row.names=FALSE, col.names = FALSE, sep= ",")
#Feature list
Features.nochim <- data.frame(Sequence = getSequences(seqtab.nochim))
#Write to file
write.table(Features.nochim,paste(No_Chim_dir,"/features_0.csv", sep=""), row.names=FALSE, col.names = FALSE, sep= ",")

#Sanity check
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(track) <- sample.names
head(track)
write.table(track, paste(Graphic_dir,"/Tracking.csv", sep=""))

## Taxonomy withchim
dna_chim <- DNAStringSet(getSequences(seqtab)) # Create a DNAStringSet from the ASVs
load(IDTaxa_Directory) # This file corresponds to "SILVA_SSU_r138_2019.RData" selected at the beginning
ids_chim <- IdTaxa(dna_chim, trainingSet, strand="top", processors=NULL, verbose=FALSE) # Use all processors
ranks_chim <- c("domain", "phylum", "class", "order", "family", "genus", "species") # Ranks of interest
## Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid_chim <- t(sapply(ids_chim, function(x) {
  m <- match(ranks_chim, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(taxid_chim) <- ranks_chim; rownames(taxid_chim) <- getSequences(seqtab)
write.table(taxid_chim, paste(Chim_dir,"/TaxID_chimera.csv", sep=""))

#taxonomy nochim
dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs
load(IDTaxa_Directory) # This file corresponds to "SILVA_SSU_r138_2019.RData" selected at the beginning
ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE) # Use all processors
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # Ranks of interest
## Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim)
write.table(taxid, paste(No_Chim_dir,"/TaxID.csv", sep=""))

## Print assignments
taxid.print <- taxid
taxid_chim.print <- taxid_chim
rownames(taxid.print) <- NULL
rownames(taxid_chim.print) <- NULL
head(taxid.print)
head(taxid_chim.print)
write.table(taxid.print, paste(Graphic_dir,"/Taxonomy_ID.csv", sep=""))
write.table(taxid_chim.print, paste(Graphic_dir,"/Taxonomy_ID_Chim.csv", sep=""))

## This is the end of the Script
######################################################################################################################