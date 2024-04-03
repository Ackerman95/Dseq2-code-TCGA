# Dseq2-code-KAS
 Differential expression analysis: DESeq2
Differential expression analysis: DESeq2

Realized in R studio V 4.1.3 

Packages needed ------
# BiocManager::install('EnhancedVolcano')
# devtools::install_github('kevinblighe/EnhancedVolcano')
# BiocManager::install("airway")
# BiocManager::install("DESeq2")
# BiocManager::install("org.Hs.eg.db")

## Load Libraries ------
library(tidyverse)
library(EnhancedVolcano)
library(airway)
library(DESeq2)
library(org.Hs.eg.db)
library(magrittr)
library(readxl)


####################### DATA GENES PRACTICAL TRUE DATA ##########################

## Define The working directory
setwd("")

## load the data matrix (countData)
mcountData <- as.matrix(read_excel("Data_genes.xlsx", sheet = 1, col_types = c("text", rep("numeric", 1082)), col_names = TRUE))

## prepare the data for the differential equation
rownames(mcountData) <- mcountData[, 1]   ## Set index rowname
mcountData <- mcountData[-1, -1]   ## remove the 1st column and line
storage.mode(mcountData) <- "integer"    ## Change from characters to numbers

## Prepare the data coldata (colData)
colData <- data.frame(condition=ifelse(grepl("Col0C", colnames(mcountData)), "control", "triggered"))   ## define analysis condition
rownames(colData) <- colnames(mcountData)  ## Attribute categories
colData$condition <- factor(colData$condition)  ## Change to factor

## Run the function (part 1) ----
dds <- DESeqDataSetFromMatrix(mcountData, colData, formula(~ condition))

## Conduct differential expression using DESeq2 in order to create 2 sets of results: -----
dds2 <- DESeq(dds, betaPrior = FALSE)
res <- results(dds2, contrast = c('condition','control','triggered'))
res2 <- lfcShrink(dds2, contrast = c('condition','control','triggered'), res = res, type = 'normal')
normCounts <- counts(dds2, normalized = T)
## Plot  basic volcano plot ------
EnhancedVolcano(res2,
                lab = rownames(res2),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Bleue versus Rouge',
                pCutoff = 10e-8,     ## can be  set to : 10e-30, 10e-10
                FCcutoff = 0.4,       ## can be set to : 0.5, 0.3
                pointSize = 3.0,
                labSize = 6.0)

