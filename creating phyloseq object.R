# what I learned. Allison's code consistantly uses the tables
  # that have acutall sequences as the headers. (seqtab.nochim, taxa)
  # siedu renames his sequences as ASV 1-x first and then
  # uses all tables named in this way. You have to pick one
  # of these methods first and then stick to it!


library("import")
library("knitr")
library("BiocStyle")
library("ggplot2")
library("gridExtra")
library("Rcpp")
library("dada2")
library("phyloseq")
library("DECIPHER")
library("ape")
library("phangorn")
library("ShortRead")
# sample names
sample.names = 193:302
head(sample.names)
sample.names = c(sample.names, 384)
sample.names
length(sample.names)

# ATTEMPT 2: loading in these files directly from their text files
asv_tab = read.table("Analysis/ASVs_counts.txt")
load("Bison_Study_Samples_seqtab.nonchim_taxa.RData")
# Load in phylogenetic tree made in Mothur
tree = read_tree("mothur\\ASVs.phylip.tre")

# Read in meta data file using 16s mockpipeline
Bison_meta = read.csv("bison_metadata_analysis.csv", header=TRUE)
row.names(Bison_meta)= sample.names
# I found my problem! Bison meta contains info on all samples in
  # the sequencing run (#1-384) but the sequence tables and other
  # tables created in the analysis part have only the samples
  # that are included in this project. So the row names are not
  # matching up between the different tables.
# NVM fixed this and still get the same error

# Check object classes
class(asv_tab) # "matrix" "array"
class(asv_tax) # "matrix" "array"
class(Bison_meta) #"data.frame"


# Create Phyloseq object using 16s mockpipeline code
ps1 = phyloseq(otu_table(seqtab.nochim,
                          taxa_are_rows = FALSE))#nochimeras

ps2 <- merge_phyloseq(ps1, sample_data(Bison_meta))

ps3 <- merge_phyloseq(ps2, tax_table(taxa))

BiocManager::install("Biostrings")

ps3_dna <- Biostrings::DNAStringSet(taxa_names(ps3))

names(ps3_dna) <- taxa_names(ps3)

ps3 <- merge_phyloseq(ps3, ps3_dna)

taxa_names(ps3) <- paste0("ASV", seq(ntaxa(ps3)))

ps3

head(refseq(ps3))

head(otu_table(ps3))

taxa_names(ps3) <- paste0("ASV_", seq(ntaxa(ps3)))
taxa_names(ps3)

row.names(Bison_meta) <- as.character(Bison_meta[, 1])                       

ps4 = merge_phyloseq(ps3, phy_tree(tree))
View(ps4)
Bison_ps = ps4
saveRDS(Bison_ps, file = "Bison_ps_attempt2.rds")
