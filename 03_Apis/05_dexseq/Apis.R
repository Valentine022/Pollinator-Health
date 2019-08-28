rm( list = ls () )

basic_libraries <- c("ggplot2", "ggpubr", "gplots", "RColorBrewer", "genefilter", "lattice")

for (lib in basic_libraries) {
  if (require(package = lib, character.only = TRUE)) {
    print("Successful")
  } else {
    print("Installing")
    install.packages(lib)
    library(lib, character.only = TRUE )
  }
}

suppressPackageStartupMessages(library("DEXSeq"))
# Read data into R 
sampleInfo = read.csv("sample_information.csv", header=TRUE)
colnames(sampleInfo)[1] <- "SampleID"
#rownames(sampleInfo) <- sampleInfo$SampleID
#sampleInfo$SampleID <- NULL
sampleInfo$Caste <- as.factor(sampleInfo$Caste)
sampleInfo$Treatment <- as.factor(sampleInfo$Treatment)


#load the counts files
countfiles <- list.files("./",
                         pattern = c(".counts$"),
                         full.names = TRUE)

#load the reference genome in GFF format
flattenedfile <- list.files("./", pattern = "gtf.gff$", full.names = TRUE)

#create DEXSeq object
dxd <- DEXSeqDataSetFromHTSeq(countfiles, sampleData = sampleInfo,
                              design = ~ SampleID + exon + Treatment: exon, 
                              flattenedfile = flattenedfile)

#normalisation of exon counts by library sizes:
dxd <- estimateSizeFactors(dxd)

#estimate dispersion
dxd <- estimateDispersions(dxd)

#test for differential exon usage between treatment groups
dxd <- testForDEU(dxd)

#estimate fold changes in exon usage between castes
dxd <- estimateExonFoldChanges(dxd, fitExpToVar = "Treatment")

#list of genes we are interested in :
genes <- c("GB42644", "GB53055",  "GB43275", "GB47845", "GB53427", "GB53428")

#summarize the results in DEXSeqResults object
dexseq_object_results <- DEXSeqResults(dxd)
dexseq_object_results_df <- as.data.frame(dexseq_object_results)

## Create a directory for output:
output_directory <-  paste("results/differential_exon_raw_all/", sep = "")
dir.create(path = output_directory,
           recursive = TRUE)
## Write results to file:
write.table(dexseq_object_results_df,
            file = paste(output_directory, "dexseq_object_results_all_1.txt",
                         sep = ""),
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t",
            quote = FALSE)

## Subset columns of interest (groupID, pvalue and padj)
dexseq_object_results_pvalues <- dexseq_object_results_df[, c(1, 6, 7)]
write.table(dexseq_object_results_pvalues,
            file = paste(output_directory, "dexseq_object_results_all_Apis.txt",
                         sep = ""),
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t",
            quote = FALSE)

#save MA plot in a pdf file
pdf(paste("plot_MA_Apis.pdf", sep =""))
plotMA(dexseq_object_results,
       cex = 0.5,
       alpha = 0.05,
       ylim = c(-10, 10))
dev.off()

#examine DEU for all for genes and save plots to a pdf file
pdf(paste("plot_Apis.pdf", sep = ""))


for (gene in genes) {
  plotDEXSeq(dexseq_object_results, geneID = gene, fitExpToVar = "Treatment",
             norCounts = TRUE, displayTranscripts = TRUE, legend = TRUE, 
             cex.axis = 1.2, cex = 1.3, lwd = 2)
}
