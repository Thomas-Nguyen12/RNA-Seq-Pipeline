
args <- commandArgs()

sample_sheet <- args[1]
output_directory <- args[2]



library(DESeq2) 
library(tximportData)
library(tximport)
library(GenomicFeatures) 
library(pacman) 
library(readr) 
library(rhdf5) 



library(rtracklayer)

# setting the working directory. THis should be the output directory
setwd(output_directory)



# ----- PART 1. CREATING SAMPLE SHEET AND 
print ("Loading the Sample names...") 
sample_names <- sample_sheet$sample
print (sample_names) 
print ("Loading the conditions...") 
condition <- sample_sheet$condition  
print (condition) 
print ("Creating the sample sheet...") 


# making a dataframe with the file locations 
print ("Adding the file paths for abundance.h5") 
files <- file.path("./", paste0("results_", sample_sheet$sample_names), "abundance.h5") 
print (files)

# adding headers to the dataframe
print ("Adding headers to the dataframe...")
names(files) <- sample_names
files

# ------ PART 2. CREATING THE TxDB database
print ("Creating the txdb database") 
# A txdb, in the context of bioinformatics, is a database that stores 
# transcript annotations for a given genome assembly.


txdb <- makeTxDbFromGFF("gencode.v25.annotation.gtf")
# checking the columns 
columns(txdb) 

keytypes(txdb)
k <- keys(txdb, keytype = "TXNAME")
tx_map <- select(txdb, keys = k, columns="GENEID", keytype="TXNAME")
head(tx_map)













# ------- running tximport 

tx2gene <- tx_map

tx2gene$TXNAME <- sub("\\.\\d+$", "", tx2gene$TXNAME)

txi.kallisto <- tximport(files, type='kallisto', tx2gene=tx2gene, ignoreTxVersion=TRUE)

# ------ DESeq2 will be comparing CT27 and Stromal samples (I will then isolate the top 200 differentially expressed genes) 


# Creating a DESeq2 dataset from a txi object
ddsTxi <- DESeqDataSetFromTximport(txi.kallisto,
                                   colData = sample_sheet,
                                   design = ~ condition)
ddsTxi
# Pre-filtering the lowest expressed genes (less than 10 counts)
keep <- rowSums(counts(ddsTxi)) >= 10
dds <- ddsTxi[keep,]
dds



dds$condition <- factor(dds$condition, levels = unique(condition)) )




# differential expression analysis 
dds <- DESeq(dds)
res <- results(dds)

res 
resultsNames(dds)




