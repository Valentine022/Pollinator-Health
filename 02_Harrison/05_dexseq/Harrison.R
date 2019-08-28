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
                              design = ~ SampleID + exon + Caste: exon, 
                              flattenedfile = flattenedfile)

#normalisation of exon counts by library sizes:
dxd <- estimateSizeFactors(dxd)

#estimate dispersion
dxd <- estimateDispersions(dxd)

#test for differential exon usage between treatment groups
dxd <- testForDEU(dxd)

#estimate fold changes in exon usage between castes
dxd <- estimateExonFoldChanges(dxd, fitExpToVar = "Development")

#list of genes we are interested in :
genes <- c("LOC100647624", "LOC100647301",  "LOC100647350", "LOC100645032", "LOC100648987", "LOC100649515", "LOC100643274", "LOC100649796", "LOC100643282", "LOC100649612")

#summarize the results in DEXSeqResults object
dexseq_object_results <- DEXSeqResults(dxd)
dexseq_object_results_df <- as.data.frame(dexseq_object_results)

## Create a directory for output:
output_directory <-  paste("results/", "life_stages_differential_exon_raw_all/", sep = "")
dir.create(path = output_directory,
           recursive = TRUE)
## Write results to file:
write.table(dexseq_object_results_df,
            file = paste(output_directory, "dexseq_object_results_all_1_.txt",
                         sep = ""),
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t",
            quote = FALSE)

## Subset columns of interest (groupID, pvalue and padj)
dexseq_object_results_pvalues <- dexseq_object_results_df[, c(1, 6, 7)]
write.table(dexseq_object_results_pvalues,
            file = paste(output_directory, "dexseq_object_results_all_1_.txt",
                         sep = ""),
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t",
            quote = FALSE)

#save MA plot in a pdf file
pdf(paste("plot_MA_1_.pdf", sep = ""))
plotMA(dexseq_object_results,
       cex = 0.5,
       alpha = 0.05,
       ylim = c(-10, 10))
dev.off()

#examine DEU for all for genes and save plots to a pdf file
pdf(paste("plot_1_.pdf", sep = ""))


for (gene in genes) {
  plotDEXSeq(dexseq_object_results, geneID = gene, fitExpToVar = "Caste",
             norCounts = TRUE, displayTranscripts = TRUE, legend = TRUE, 
             cex.axis = 1.2, cex = 1.3, lwd = 2)
}


#r-log transformation of DEXSeq object
rld <- rlog(dxd)
#rld2 <- vst(dxd)
assay_rld <- assay(rld)

#find the indexes of the four genes we are interested in in assay_rld, so we can use them to plot the heatmap
patterns <-  grep("LOC100647624+|LOC100647301+|LOC100647350+|LOC100645032+|LOC100648987+|LOC100649515+|LOC100643274+|LOC100649796+|LOC100643282+|LOC100649612+", rownames(assay_rld), perl=TRUE, value=TRUE)
genes_indexes <- c()
for(pattern in patterns){
  Indx <- which(rownames(assay_rld) == pattern)
  genes_indexes <- append(genes_indexes, Indx)
}
names_for_rld <- colData(rld)$SampleID
colnames(assay_rld) <- names_for_rld

#save PCA plots in a pdf file
pdf(paste("plot_PCA_2208192.pdf", sep = ""), paper='a4r', width = 10, height = 10)
plt1 <- DESeq2::plotPCA(rld, intgroup = "Development", ntop = genes_indexes, returnData = TRUE)
print(plt1)
dev.off()

pdf(paste("plot_PCA3.pdf", sep = ""), paper='a4r', width = 10, height = 10)
DESeq2::plotPCA(rld, intgroup = "Development", ntop = genes_indexes, returnData = TRUE)
print(plt1)
str(plt1)

