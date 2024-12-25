# Angiogenesis-paper
# load libraries

library(DESeq2)
library(tidyverse)

# Retrieving Data

untar("HIF-Project-count-matrix.tar.gz", exdir = "NB_HIF_countData")

count_data <- read.delim('NB_HIF_countData/datasets/Column_join_on_data_139,_data_137,_and_others_bdb11cef953663b0a79032f02293af7b.tabular')

count_data <- as.data.frame(count_data)

#naming the columns
new_col_names <- c("Geneid", "Scr16.Doxo.1", "Scr16.Doxo.2", "Scr16.untreated.1", "Scr16.untreated.2", 
                   "Scr13.Doxo.1", "Scr13Doxo.2", "Scr13.untreated.1", "Scr13.untreated.2", "ENSEMBL", 
                   "Genes", "GENENAME"  )

colnames(count_data) = new_col_names

#arranging Columns and deleting unneeded ones, Filtering the NA (non-annotated genes) and duplicated genes 
count_data <- count_data[,c(11,1,2,3,4,5,6,7,8,9,10,12)] %>% 
  select(-c(Geneid, ENSEMBL, GENENAME)) %>% 
  na.omit()  %>% 
  unique() 

rownames(count_data) <- count_data$Genes 

count_data <- count_data[, -1]

write.csv(count_data, file = "NB_DOX_Vs_CTRL_row_count_matrix.csv") #saving the count matrix

# load libraries
library(DESeq2)
library(tidyverse)
library(ggrepel)
library(RColorBrewer)
library(EnhancedVolcano)
library(pheatmap) 


#to read the raw count matrix 
count_data <- read.csv("NB_DOX_Vs_CTRL_row_count_matrix.csv", header = T , row.names = 1)

# create in sample info data frame
sample_name <- c("Scr16.Doxo.1", "Scr16.Doxo.2", "Scr16.untreated.1",
                 "Scr16.untreated.2", "Scr13.Doxo.1", "Scr13Doxo.2", 
                 "Scr13.untreated.1", "Scr13.untreated.2")

condition <- c("Doxo", "Doxo", "CTRL", "CTRL", 
               "Doxo", "Doxo", "CTRL", "CTRL")

sample_info <- data.frame(sample_name, condition)

rownames(sample_info) = sample_info$sample_name

write.csv(sample_info,'NB_sample_info.csv', row.names = FALSE)

# making sure the row names in sample_info matches to column names in count_data
all(colnames(count_data) %in% rownames(sample_info))

# are they in the same order?
all(colnames(count_data) == rownames(sample_info))


# construct a DESeqDataSet object

dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = sample_info,
                              design = ~ condition)

dds

rld <- rlog(dds)

#Plot PCA 

plotPCA (rld)

# set the factor level
dds$condition <- relevel(dds$condition, ref = "CTRL")


# Run DESeq ----------------------
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "Doxo", "CTRL"))

res


# MA plot
plotMA(res)

resNorm <- lfcShrink(dds = dds, res = res, type = "normal", coef = 2)

plotMA(resNorm)

# select genes that are significant

idx = which( resNorm$padj <= 0.05 & 
              abs(resNorm$log2FoldChange)  >= 0.3 & 
               resNorm$baseMean >= 20 )
sigRes = resNorm[idx, ]

# Save our significant genes and the whole results
sig_Res = cbind( genes = rownames(sigRes), sigRes )

res_Norm = cbind( genes = rownames(resNorm), resNorm )

write.csv(sig_Res, file = "NB_DOXO_HIF_sigGenes_adjp0.05_lfc0.3_bm20.csv", row.names = FALSE)

write.csv(res_Norm, file = "NB_DOXO_HIF_whole.csv", row.names = FALSE)

#ploting our DEGs

DEGs_Volcano <- EnhancedVolcano(resNorm,
                                lab = rownames(resNorm),
                                x = 'log2FoldChange',
                                y = 'pvalue',
                                xlab = bquote(~Log[2]~ 'fold change'),
                                subtitle = '',
                                caption = '',
                                pCutoff = 0.05,
                                pCutoffCol = "padj",
                                FCcutoff = 0.3,
                                pointSize = 1.5,
                                labSize = 3,
                                legendPosition = "right"
)

ggsave(DEGs_Volcano, filename = "DEGs_Volcano.png", width = 10, height = 6, dpi = 300, type = "cairo-png")

#Heatmap

mat <- assay(rld)
orderedSig <- sigRes[order(sigRes$padj), ]
top50DEG <- head(orderedSig, n= 50)
DEgenes <- mat[rownames(top50DEG),]

annotation <- as.data.frame(colData(rld)[, c("condition", "sample_name" )])
colnames(annotation) <-  c("Condition", "Sample_name" )
coul <- colorRampPalette(brewer.pal(5, "PiYG"))(25)

DEGs_heatmap <- pheatmap(DEgenes, color = coul, scale = "row", 
                         clustering_distance_rows = "correlation", 
                         annotation_col = annotation, 
                         fontsize = 15,
                         show_colnames = F, main="Top 50 Differentially Expressed Genes")

ggsave(DEGs_heatmap, filename = "DEGs_heatmap.png", width = 10, height = 15, dpi = 300, type = "cairo-png")


#loading the needed packages
library(tidyverse)
library(RColorBrewer)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(ggupset)

#load data
df <- as.data.frame(res)

# Add column in the data frame to annotate the up-regulated, down-regulated and non-significant gene expression
df <- df %>% mutate(diffexpressed = case_when(
  log2FoldChange > 0.3 & padj < 0.05 ~ 'UP',
  log2FoldChange < 0.3 & padj < 0.05 ~ 'DOWN',
  padj > 0.05 ~ 'NO'
))

# Remove NAs
df <- df[!is.na(df$diffexpressed) , ]


# Get the genes that are present in your dataframe
df$gene_symbol <- rownames(df)
genes_in_data <- df$gene_symbol

# Read in the .gmt file
Bg_Genes <- read.gmt("c2.cp.kegg_medicus.v2023.2.Hs.symbols.gmt") 
# Subset to the genes that are present in our dataset
Bg_Genes <- Bg_Genes[Bg_Genes$gene %in% genes_in_data,] 
# Save the filtered background gene set
saveRDS(Bg_Genes, 'kegg_Bg_Genes.RDS')

# Remove non-significant genes
df <- df[df$diffexpressed != 'NO', ]

# Split the dataframe into a list of sub-dataframes: upregulated, downregulated genes
DEGs_results_list <- split(df, df$diffexpressed)

# Run clusterProfiler on each sub-dataframe
results_PEA <- lapply(names(DEGs_results_list),
                      function(x) enricher(gene = DEGs_results_list[[x]]$gene_symbol,
                                           TERM2GENE = Bg_Genes))
names(results_PEA) <- names(DEGs_results_list)

#Convert the list of enrichResults for each sample_pattern to a dataframe with the pathways
Res_PEA_df <- lapply(names(results_PEA), function(x) rbind(results_PEA[[x]]@result))
names(Res_PEA_df) <- names(results_PEA)
Res_PEA_df <- do.call(rbind, Res_PEA_df)
Res_PEA_Up <- results_PEA[['UP']]@result
Res_PEA_Down <- results_PEA[['DOWN']]@result

#Save the enriched pathways with padj < 0.5 and counts > 5
target_PEs <- unique(row.names(Res_PEA_df)[Res_PEA_df$p.adjust < 0.05 & Res_PEA_df$Count > 5]) # select only target pathways have p adjusted < 0.05 and at least 6 genes
PEA_output <- Res_PEA_df[row.names(Res_PEA_df) %in% target_PEs, ]
write.csv(PEA_output, file = "PEA_output.csv")

# To shorten the pathway names
Res_PEA_Up$Description <- gsub('(H|h)iv', 'HIV', 
                              gsub('pd 1', 'PD-1',
                                   gsub('ecm', 'ECM', 
                                        gsub('(I|i)nterleukin', 'IL', 
                                             gsub('(R|r)na', 'RNA', 
                                                  gsub('(D|d)na', 'DNA',
                                                       gsub(' i ', ' I ', 
                                                            gsub('(A|a)tp ', 'ATP ', 
                                                                 gsub('(N|n)adh ', 'NADH ', 
                                                                      gsub('(N|n)ad ', 'NAD ',
                                                                           gsub('t cell', 'T cell',
                                                                                gsub('b cell', 'B cell',
                                                                                     gsub('built from .*', ' (...)',
                                                                                          gsub('mhc', 'MHC',
                                                                                               gsub('mhc class i', 'MHC I', 
                                                                                                    gsub('mhc class ii', 'MHC II', 
                                                                                                         stringr::str_to_sentence(
                                                                                                           gsub('_', ' ',  
                                                                                                                gsub('GOBP_|KEGG_MEDICUS_REFERENCE_|REACTOME_', '', Res_PEA_Up$Description)))))))))))))))))))


Res_PEA_Down$Description <- gsub('(H|h)iv', 'HIV', 
                               gsub('pd 1', 'PD-1',
                                    gsub('ecm', 'ECM', 
                                         gsub('(I|i)nterleukin', 'IL', 
                                              gsub('(R|r)na', 'RNA', 
                                                   gsub('(D|d)na', 'DNA',
                                                        gsub(' i ', ' I ', 
                                                             gsub('(A|a)tp ', 'ATP ', 
                                                                  gsub('(N|n)adh ', 'NADH ', 
                                                                       gsub('(N|n)ad ', 'NAD ',
                                                                            gsub('t cell', 'T cell',
                                                                                 gsub('b cell', 'B cell',
                                                                                      gsub('built from .*', ' (...)',
                                                                                           gsub('mhc', 'MHC',
                                                                                                gsub('mhc class i', 'MHC I', 
                                                                                                     gsub('mhc class ii', 'MHC II', 
                                                                                                          str_to_sentence(
                                                                                                            gsub('_', ' ',  
                                                                                                                 gsub('GOBP_|KEGG_MEDICUS_REFERENCE_|REACTOME_', '', Res_PEA_Down$Description)))))))))))))))))))


#Plotting
enrichres_up <- new("enrichResult",
                 result = Res_PEA_Up,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 qvalueCutoff = 0.05,
                 organism = "human",
                 ontology = "UNKNOWN",
                 gene = df$gene_symbol,
                 keytype = "UNKNOWN",
                 universe = unique(Bg_Genes$gene),
                 geneSets = Bg_Genes)

enrichres_down <- new("enrichResult",
                    result = Res_PEA_Down,
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "BH",
                    qvalueCutoff = 0.05,
                    organism = "human",
                    ontology = "UNKNOWN",
                    gene = df$gene_symbol,
                    keytype = "UNKNOWN",
                    universe = unique(Bg_Genes$gene),
                    geneSets = Bg_Genes)

bar_up <- barplot(enrichres_up, showCategory = 20, font.size= 20)  + 
  theme(legend.text=element_text(size= 20),  
        legend.title=element_text(size= 20))
cnet_up <- cnetplot(enrichres_up,  color.params = list(foldChange = TRUE, edge = FALSE, category = "#E06663", gene ="#717BA8")) 

ggsave(bar_up, filename = "bar_up.png", width = 15, height = 10, dpi = 300, type = "cairo-png")
ggsave(cnet_up, filename = "cnet_up.png", width = 10, height = 6, dpi = 300, type = "cairo-png")

bar_down <- barplot(enrichres_down, showCategory = 20, font.size= 20) + 
  theme(legend.text=element_text(size= 20),  
        legend.title=element_text(size= 20))
cnet_down <- cnetplot(enrichres_down,  color.params = list(foldChange = TRUE, edge = FALSE, category = "#E06663", gene ="#717BA8")) 

ggsave(bar_down, filename = "bar_down.png", width = 15, height = 25, dpi = 300, type = "cairo-png")
ggsave(cnet_down, filename = "cnet_down.png", width = 10, height = 8, dpi = 300, type = "cairo-png")
