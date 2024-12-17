# DEGS

if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install('EnhancedVolcano')
BiocManager::install("devtools")
BiocManager::install("pachterlab/sleuth")


library(sleuth)
library(tidyverse)  
library(EnhancedVolcano)
library(pheatmap)

#read in sample tables - be sure to set correct path 

metadata <- read.table(file = "ExpTable_GOULD.txt", sep='\t', header=TRUE, stringsAsFactors = FALSE)

#this command sets up paths to the kallisto output that we will process in the following steps
metadata <- dplyr::mutate(metadata,
                          path = file.path('output', Run_s, 'abundance.h5'))
metadata <- dplyr::rename(metadata, sample = Run_s)

#let's check the metadata
metadata


#Read in headers for the transcripts that we aligned to with kallisto
#These will be mapped in the sleuth_prep command below

ttn<-read_delim("TTN_xl.csv", col_names = FALSE)

colnames(ttn)<-c("target_id","gene")
head(ttn)
so <- sleuth_prep(metadata, full_model = ~treat, target_mapping = ttn, 
                  extra_bootstrap_summary = TRUE, read_bootstrap_tpm = TRUE, 
                  aggregation_column = "gene")


#fit model specified above
so <- sleuth_fit(so)

#print the model
models(so)

# [  full  ]
# formula:  ~treat 
# data modeled:  obs_counts 
# transform sync'ed:  TRUE 
# coefficients:
# 	(Intercept)
#  	treatTTC

#calculate the Wald test statistic for 'beta' coefficient on every transcript 
so <- sleuth_wt(so, 'treathead')

#extract the wald test results for each transcript 
transcripts_all <- sleuth_results(so, 'treathead', show_all = FALSE, pval_aggregate = FALSE)

transcripts_allAGG <- sleuth_results(so, 'treathead', show_all = FALSE, pval_aggregate = TRUE)

#filtered by significance 
transcripts_sig <- dplyr::filter(transcripts_all, qval <= 0.05)

#filtered by significance 
transcripts_sigAGG <- dplyr::filter(transcripts_allAGG, qval <= 0.05)

transcripts_50 <- dplyr::filter(transcripts_all, qval <= 0.05) %>%
  head(50)

transcripts_50AGG <- dplyr::filter(transcripts_allAGG, qval <= 0.05) %>%
  head(50)

head(transcripts_50AGG)
head(transcripts_50)

#extract the gene symbols, qval, and b values from the Wlad test results
forVolcano<-data.frame(transcripts_all$gene, transcripts_all$qval, transcripts_all$b)

#rename the columns of the dataframe
colnames(forVolcano)<-c("gene","qval","b")

#plot
EnhancedVolcano(forVolcano,
                lab = forVolcano$gene,
                x = 'b',
                y = 'qval',
                xlab = "\u03B2",
                labSize = 3,
                legendPosition = "none")

k_table <- kallisto_table(so, normalized = TRUE)

k_DEG <- k_table %>%
  right_join(transcripts_50, "target_id")

k_DEG_select<-k_DEG %>%
  #apply log10 transformation to the tpm data
  mutate(log_tpm = log10(tpm+1)) %>%
  #select the specifc columns to plot
  dplyr::select(target_id, sample, log_tpm, gene) %>%
  #create "label" from the transcript id and gene symbol
  mutate(label = paste(target_id, gene))%>%
  #pivot data frame to a wide format
  pivot_wider(names_from = sample, values_from = log_tpm) %>%
  #drop the target_id and gene variables
  dplyr::select(!target_id & !gene) %>%
  #convert label to row name
  column_to_rownames("label") %>%
  #convert to matrix
  as.matrix(rownames.force = TRUE) 

#plot with pheatmap!
pheatmap(k_DEG_select, cexRow = 0.4, cexCol = 0.4, scale = "none")

#filter for transcripts enriched in the TTC treatment
transcripts_up <- dplyr::filter(transcripts_all, qval <= 0.05, b > 0)

up<-transcripts_up %>%
  dplyr::select(gene)

#filter for transcripts depleted in the TTC treatment
transcripts_down <- dplyr::filter(transcripts_all, qval <= 0.05, b < 0)

down<-transcripts_down %>%
  dplyr::select(gene)

#output the full transcript list
all<-transcripts_all %>%
  dplyr::select(gene)

library(clipr)
#copy to clipboard and paste into ShinyGo website
write_clip(as.character(up))

#copy to clipboard and paste into ShinyGo "background"
write_clip(as.character(all))



##### DO WITH ANNOTATED

metadata2 <- read.table(file = "ExpTable_GOULD.txt", sep='\t', header=TRUE, stringsAsFactors = FALSE)
#this command sets up paths to the kallisto output that we will process in the following steps
metadata2 <- dplyr::mutate(metadata2,
                          path = file.path('output', Run_s, 'abundance.h5'))
metadata2 <- dplyr::rename(metadata2, sample = Run_s)

metadata2
#Read in headers for the transcripts that we aligned to with kallisto
#These will be mapped in the sleuth_prep command below

ttn2<-read_delim("ENSEGOT_GN.tsv", col_names = FALSE)

colnames(ttn2)<-c("target_id","gene")

so2 <- sleuth_prep(metadata2, full_model = ~treat, target_mapping = ttn2, 
                   extra_bootstrap_summary = TRUE, read_bootstrap_tpm = TRUE, 
                   aggregation_column = "gene")

#fit model specified above
so2 <- sleuth_fit(so2)

#print the model
models(so2)

#[  full  ]
#formula:  ~treat 
#data modeled:  obs_counts 
#transform sync'ed:  TRUE 
#coefficients:
#	(Intercept)
# 	treathead

#calculate the Wald test statistic for 'beta' coefficient on every transcript 
so2 <- sleuth_wt(so2, 'treathead')

#extract the wald test results for each transcript 
transcripts_all2 <- sleuth_results(so2, 'treathead', show_all = FALSE, pval_aggregate = FALSE)
head(transcripts_all2, 10)

#filtered by significance 
transcripts_sig2 <- dplyr::filter(transcripts_all2, qval <= 0.05)
head(transcripts_sig2)

#top 50
transcripts_50_2 <- dplyr::filter(transcripts_all2, qval <= 0.05) %>%
  head(50)
head(transcripts_50_2)

genes_all <- sleuth_results(so2, 'treathead', show_all = FALSE, pval_aggregate = TRUE)
head(genes_all)
#filtered by significance 
genes_all_sig <- dplyr::filter(genes_all, qval <= 0.05)
head(genes_all_sig)

#top 50
genes_50 <- dplyr::filter(genes_all, qval <= 0.05) %>%
  head(50)
head(genes_50)
write.csv(genes_50, "genes_50.csv")

#extract the gene symbols, qval, and b values from the Wlad test results
forVolacano<-data.frame(transcripts_all2$gene, transcripts_all2$qval, transcripts_all2$b)

#rename the columns of the dataframe
colnames(forVolacano)<-c("gene","qval","b")

#plot
EnhancedVolcano(forVolacano,
                lab = forVolacano$gene,
                x = 'b',
                y = 'qval',
                xlab = "\u03B2",
                labSize = 3,
                legendPosition = "none")

k_table <- kallisto_table(so2, normalized = TRUE)

k_DEG <- k_table %>%
  right_join(transcripts_50_2, "target_id")

k_DEG_select<-k_DEG %>%
  #apply log10 transformation to the tpm data
  mutate(log_tpm = log10(tpm+1)) %>%
  #select the specifc columns to plot
  dplyr::select(target_id, sample, log_tpm, gene) %>%
  #create "label" from the transcript id and gene symbol
  mutate(label = paste(target_id, gene))%>%
  #pivot data frame to a wide format
  pivot_wider(names_from = sample, values_from = log_tpm) %>%
  #drop the target_id and gene variables
  dplyr::select(!target_id & !gene) %>%
  #convert label to row name
  column_to_rownames("label") %>%
  #convert to matrix
  as.matrix(rownames.force = TRUE) 

#plot with pheatmap!
pheatmap(k_DEG_select, cexRow = 0.4, cexCol = 0.4, scale = "none")

#filter for transcripts enriched in the TTC treatment
transcripts_up <- dplyr::filter(transcripts_all2, qval <= 0.05, b > 0)

up<-transcripts_up %>%
  dplyr::select(gene)

#filter for transcripts depleted in the TTC treatment
transcripts_down <- dplyr::filter(transcripts_all2, qval <= 0.05, b < 0)

down<-transcripts_down %>%
  dplyr::select(gene)

#output the full transcript list
all<-transcripts_all2 %>%
  dplyr::select(gene)

library(clipr)

#copy to clipboard and paste into ShinyGo website
write_clip(as.character(up))
#copy to clipboard and paste into ShinyGo "background"
write_clip(as.character(all))

write_clip(as.character(down))

library(clipr)
#copy to clipboard and paste into ShinyGo website

##### try placing genes of interest
#EnhancedVolcano(res2,
               # lab = rownames(res2),
                #x = "log2FoldChange",
               # y = "padj",
              #  selectLab = c("ENSG00000106565","ENSG00000187758"), )


#plot
EnhancedVolcano(forVolacano,
                lab = forVolacano$gene,
                x = 'b',
                y = 'qval',
                xlab = "\u03B2",
                selectLab = c("DHRS3", "SCARB2", "CYP36B1", "CYP27A1", "CYP27C1", "CYP24C1"),
                labSize = 3,
                legendPosition = "none")


