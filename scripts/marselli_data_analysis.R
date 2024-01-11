# 0. Libraries and function definitions  -----------------------------
#'[' ==================================================================================================================
#'['            Load RNAseq data from Marselli (27 T2D and 58 ND) for edgeR analysis             
#'[' ==================================================================================================================
#'['
#'[' 

# Libraries
library(tidyverse) # general purpose data wrangling
library(GEOquery)  # for extracting metadata from GEO identifiers
library(edgeR)     # for DEG analysis - use the LRT not QLF
library(ggthemes)  # for ggplot mods
library(ggrepel)   # for labeling ggplots
library(RColorBrewer) # for a colorful plot
library(readxl)
library(ggtext)
library(VennDiagram)


# 1. Loading count tables ----------------------------------------------------
#'[' ==================================================================================================================
#'['      Loading bulk islet RNAseq count tables, metadata, and gene identifier cross-reference table (MANE)
#'[' ==================================================================================================================

# * 1.1 Loading MANE gene table -------------------------------------------------
#'* Load the MANE gene symbol lookup table - acquired from the ensembl MANE project https://tark.ensembl.org/web/access_mane_data/ *
#'* Files are here: https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.1/  *

# read in MANE lookup table, create cross reference tables, confirm no duplicate rows.
MANE_lookup <- read.delim(file = "data/MANE/MANE.GRCh38.v1.1.summary_clean.txt")
MANE_lookup <- MANE_lookup[!duplicated(MANE_lookup$symbol), ]
# sum(duplicated(MANE_lookup$hgnc_id)) # 43 duplications in the hgnc_id in MANE
# sum(duplicated(MANE_lookup$symbol)) # 0 duplicates in symbol
# dups_in_MANE <- MANE_lookup[duplicated(MANE_lookup$hgnc_id), ]

MANE_entrez_ensembl <- MANE_lookup %>% dplyr::select(entrezgene_id, ensembl_gene_id)
MANE_entrez_ensembl$entrezgene_id <- as.character(MANE_entrez_ensembl$entrezgene_id)

 # sum(duplicated(MANE_entrez_ensembl$entrezgene_id)) # 0 duplicates in entrez id
 # sum(duplicated(MANE_entrez_ensembl$ensembl_gene_id)) # 0 duplicates in ensembl id

MANE_ensembl_symbol <- MANE_lookup %>% dplyr::select(ensembl_gene_id, symbol)
# sum(duplicated(MANE_ensembl_symbol$ensembl_gene_id)) # 0 dups
# sum(duplicated(MANE_ensembl_symbol$symbol)) # 0 dups

MANE_ensembl_entrez_symbol <- MANE_lookup %>% dplyr::select(entrezgene_id, ensembl_gene_id, symbol)
MANE_ensembl_entrez_symbol$entrezgene_id <- as.character(MANE_ensembl_entrez_symbol$entrezgene_id)
# sum(duplicated(MANE_ensembl_entrez_symbol$symbol)) # 0 dups

# * 1.2 Load Marselli -----------------------------------------------------
#'* Load Marselli reads and metadata - GSE159984

# load counts table from GEO for GSE159984 - the script came from the GEO analyzer
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(urld, "acc=GSE159984", "file=GSE159984_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames="GeneID")

# load gene annotations - from GEO analyzer
apath <- paste(urld, "type=rnaseq_counts", "file=Human.GRCh38.p13.annot.tsv.gz", sep="&")
annot <- data.table::fread(apath, header=T, quote="", stringsAsFactors=F, data.table=F)
rownames(annot) <- annot$GeneID

# sample selection - GEO analyzer created this part to help select only the ND and T2D samples, excluding all the samples treated with glucolipotoxic stress.
# One GSM is missing its read count table, GSM4852056 - "human islets 24" from a T2D donor. 
gsms <- paste0("00000000000000000000000000011111111111111111111111",
               "11111111111111111111111111111111111XXXXXXXXXXXXXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
sml <- strsplit(gsms, split="")[[1]]

# filter out excluded samples (marked as "X") 
sel <- which(sml != "X")
sml <- sml[sel]
tbl <- tbl[ ,sel]

# group membership for samples
gs <- factor(sml)
groups <- make.names(c("T2D","ND"))
levels(gs) <- groups
sample_info <- data.frame(Group = gs, row.names = colnames(tbl))

###### The code below is to extract the read count table and genes for other purposes - call it tbl2
###### label the read count table 'tbl' with the Gene symbol from 'annot'

marselli <- as.data.frame(tbl)
marselli$GeneID <- rownames(marselli)
rownames(marselli) <- NULL
marselli <- marselli[ ,c(86,1:85)]

marselli_metadata <- read_excel("data/marselli/marselli_total_clean_metadata.xlsx") #metadata from multiple sources was manually curated from supplemental data, GEO, and files provided by Dr. Piron.

column_mapping <- setNames(marselli_metadata$SAMN_and_disease, marselli_metadata$GSM_id) # make a corresponding index to map the SAMN_disease col to the tbl

# use function below to rename the column headers from GSM id to the SAMN_disease
colnames(marselli) <- sapply(colnames(marselli), function(x){
  if (x %in% names(column_mapping)) {
    return(column_mapping[x])
  } else {
    return(x)
  }
})

# write tables as tab delimited and csv
write.table(marselli, file = "data/marselli/marselli_GSE159984_raw_read_counts.txt", sep = '\t', quote = F)
write.csv2(marselli, file = "data/marselli/marselli_GSE159984_raw_read_counts.csv", quote = F, row.names = F)

# sanity check
sum(duplicated(marselli$GeneID)) # 0 duplicated GeneIDs

# * 1.3 Mapping Marselli --------------------------------------------------
##'* Mapping gene ids for marselli *

# create first dataframe for downstream merging of all bulk data by ensembl_id
marselli_mapids <- inner_join(marselli, 
                              MANE_entrez_ensembl, 
                              by = c("GeneID" = "entrezgene_id")) # generates table with col1 = ensembl_id and col2:86 are the 85 samples

marselli_mapids <- marselli_mapids[ , c(87,2:86)]
# sum(duplicated(marselli_mapids$ensembl_gene_id)) # 0 dups

# create second dataframe with symbol in col1 for performing edgeR on this independent dataset
marselli_mapids2 <- inner_join(marselli, 
                               MANE_ensembl_entrez_symbol, 
                               by = c("GeneID" = "entrezgene_id")) # creates col1 = entrez, 2:86 samples, 87 ensembl_id, 88 symbol
marselli_mapids2 <- marselli_mapids2 %>%
  dplyr::select(-c(ensembl_gene_id, GeneID))                              # drop the excess columns

marselli_mapids2 <- marselli_mapids2[ , c(86,1:85)]                # move gene symbol to first column

# arrange columns by ND and T2D
marselli_mapids <- marselli_mapids %>%
  dplyr::select(ensembl_gene_id, contains("ND"), contains("T2D"))

marselli_mapids2 <- marselli_mapids2 %>%
  dplyr::select(symbol, contains("ND"), contains("T2D"))

# * 1.4 Reformat Marselli metadata ----------------------------------------
##'* reformat marselli metadata *

marselli_meta2 <- dplyr::select(marselli_metadata,
                                       "Human_Islet_ID" = "SAMN_and_disease", 
                                       "Human_sample_ID" = "SAMN_ID",
                                       "Site_ID" = "alt_id_1",
                                       "GEO_sample_label" = "alt_id_2",
                                       "Marselli_et_al_ID" = "alt_id_3",
                                       "Disease_state" = "donor_diabetic_type",
                                       "Sex" = "sex",
                                       "Age" = "age",
                                       "BMI" = "bmi",
                                       "Diabetes_medication" = "Anti-diabetic therapy",
                                       "COD" = "cause_death_code",
                                       "ICU_glycemia" = "ICU glycemia (mg/dl)",
                                       "Diabetes_duration_yrs" = "Diabetes duration (years)",
                                       "Cold_ischemia_time_h" = "CIT (h)")

marselli_meta2 <- mutate(marselli_meta2,
                         Ethnicity = "NA")

# arrange rows by ND and T2D, and then by Human_Islet_ID
custom_order <- c("ND", "T2D")
marselli_meta2 <- marselli_meta2 %>%
  arrange(factor(Disease_state, levels = custom_order))

# 2. edgeR on Marselli data -----------------------------------------
#'[' ==================================================================================================================
#'[' Run edgeR on Marselli data            
#'[' ==================================================================================================================

#convert marselli_mapids2 to have the genesymbol column be the rownames for edgeR to work
marselli_symbolrows <- marselli_mapids2
rownames(marselli_symbolrows) <- marselli_symbolrows$symbol
marselli_symbolrows[ , 1] <- NULL

# create data frames that describe how the samples are grouped.
marselli_sample_names <- marselli_meta2$Human_Islet_ID
marselli_conditions <- marselli_meta2$Disease_state

marselli_metadata_edger <- data.frame(
  sample_names = marselli_sample_names,
  conditions = marselli_conditions,
  stringsAsFactors = F
)

marselli_metadata_edger$conditions <- factor(marselli_metadata_edger$conditions, levels = unique(marselli_metadata_edger$conditions))
print(marselli_metadata_edger)

# create DGEList object from count matrix
marselli_dge <- DGEList(counts = marselli_symbolrows, genes = row.names(marselli_symbolrows))
marselli_dge$samples$conditions <- marselli_metadata_edger$conditions

#filtering out low expresion genes
marselli_keep <- filterByExpr(marselli_dge, group=marselli_conditions)
table(marselli_keep)
marselli_dge <- marselli_dge[marselli_keep, keep.lib.sizes=FALSE]

## TMM Normalization
marselli_dge <- calcNormFactors(marselli_dge, method = "TMM")
marselli_dge$samples

## store the normalized counts
marselli_cpm <- cpm(marselli_dge, log=F)
marselli_logcpm <- cpm(marselli_dge, log=T)

write.table(marselli_logcpm, "output/marselli_logcpm.txt", row.names = TRUE, sep = '\t', quote = F)

# spot check a few genes
# barplot(cpm["INS", ], names.arg = colnames(cpm), col = "blue", main = "INS Expression")

# create a design matrix
marselli_design <- model.matrix(~0+marselli_conditions, data = marselli_dge$samples)

marselli_design

marselli_dge <- estimateDisp(marselli_dge, marselli_design)
plotBCV(marselli_dge)

# contrasts
marselli_contrasts <- makeContrasts(T2D_v_ND = marselli_conditionsT2D - marselli_conditionsND, levels = marselli_design)

## perform LRT
marselli_fit <- glmFit(marselli_dge, marselli_design)
marselli_lrt <- glmLRT(marselli_fit, contrast = marselli_contrasts)
de_genes_lrt_marselli <- topTags(marselli_lrt, n = Inf)$table
plotMD(marselli_lrt)
summary(decideTests(marselli_lrt)) # 73 genes are significant at 0.05 FDR. 

## write tables
write.table(de_genes_lrt_marselli, file = "output/edgeR_lrt_marselli.txt", sep = '\t', quote = F, row.names = F)

# merge the de_genes result with the merged read count table - marselli_mapids2 already has symbol in col1
de_genes_lrt_marselli_counts <- full_join(de_genes_lrt_marselli, marselli_mapids2, by = c("genes" = "symbol"))
write.table(de_genes_lrt_marselli_counts, file = "output/edgeR_lrt_marselli_counts.txt", sep = '\t', quote = F, row.names = F)

# make a merge that also has the normalized counts from cpm(dge)
marselli_cpm_v2 <- as.data.frame(marselli_cpm)
marselli_cpm_v2$symbol <- rownames(marselli_cpm_v2)
rownames(marselli_cpm_v2) <- NULL
marselli_cpm_v2 <- marselli_cpm_v2[ , c(86,1:85)]

de_genes_lrt_marselli_cpm <- full_join(de_genes_lrt_marselli, marselli_cpm_v2, by = c("genes" = "symbol"))
write.table(de_genes_lrt_marselli_cpm, file = "output/edgeR_lrt_marselli_cpm.txt", sep = '\t', quote = F, row.names = F)

##'* end of Marselli edgeR analysis - further chunks are for plot creation for marselli *

# * 2.1 Marselli volcano and md plots -------------------------------------
###### volcano plot

marselli_volcano_plot <- ggplot(de_genes_lrt_marselli, aes(x = logFC, y = -log10(FDR))) +
  geom_point(color = ifelse(de_genes_lrt_marselli$FDR <= 0.05 & abs(de_genes_lrt_marselli$logFC) >= 0.58, "red", "gray"),
             alpha = 0.7,
             size = 3) +
  geom_text_repel(data = subset(de_genes_lrt_marselli, FDR <= 0.05 & abs(logFC) >= 0.58),
                  aes(label = genes), color = "black",
                  size = 4,
                  max.overlaps = 20,  # Increase if needed
                  box.padding = 0.5,  # Increase padding around each label
                  point.padding = 0.5,  # Increase padding around each point
                  nudge_x = 0.2,  # Nudge labels slightly in the x direction
                  nudge_y = 0.2) +  # Nudge labels slightly in the y direction+
  theme_minimal() +
  xlab("log<sub>2</sub> (T2D/ND)") + 
  ylab("-log<sub>10</sub> FDR") + 
  
  scale_x_continuous(breaks = seq(-5, 5, by = 1)) +  # Adjust the breaks as needed
  theme(axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.title.x = element_markdown(),
        axis.title.y = element_markdown(),
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 14),
        axis.line = element_line(color = "black", size = 1),
        panel.grid = element_blank(),
  )

# Display the plot
print(marselli_volcano_plot)

# Adjust the width of the plot
ggsave("output/marselli_edger_lrt_volcano_plot.svg", marselli_volcano_plot, width = 8.5, height = 6, device = "svg")
ggsave("output/marselli_edger_lrt_volcano_plot.pdf", marselli_volcano_plot, width = 8.5, height = 6, device = "pdf")
ggsave("output/marselli_edger_lrt_volcano_plot.tiff", marselli_volcano_plot, width = 8.5, height = 6, device = "tiff", dpi = "retina")

# Create an MD plot with gray and red
marselli_md_plot <- ggplot(de_genes_lrt_marselli %>% dplyr::arrange(FDR), aes(x = logCPM, y = logFC)) +
  geom_point(color = ifelse(de_genes_lrt_marselli$FDR <= 0.05 & abs(de_genes_lrt_marselli$logFC) >= 0.58, "red", "gray"),
             alpha = 0.7,
             size = 3) +
  #scale_color_gradient(low = "gray", high = "red") +  # Define the color gradient
  geom_text_repel(data = subset(de_genes_lrt_marselli, FDR <= 0.05 & abs(logFC) >= 0.58),
                  aes(label = genes), color = "black",
                  size = 5,
                  max.overlaps = 20,  # Increase if needed
                  box.padding = 0.5,  # Increase padding around each label
                  point.padding = 0.5,  # Increase padding around each point
                  nudge_x = 0.2,  # Nudge labels slightly in the x direction
                  nudge_y = 0.2) +
  geom_segment(data = subset(de_genes_lrt_marselli, FDR <= 0.05 & abs(logFC) >= 1),
               aes(xend = logCPM, yend = logFC),
               x = NA, y = NA, color = "black", linewidth = 0.5) +  # Add lines
  theme_minimal() +
  xlab("logCPM") + 
  ylab("log<sub>2</sub> (T2D/ND)") + 
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  theme(axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.title.x = element_markdown(),
        axis.title.y = element_markdown(),
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 14),
        legend.title = element_blank(),
        legend.position = "none",
        axis.line = element_line(color = "black", size = 1),
        panel.grid = element_blank()
  )  # Remove gridlines)  # Add black solid axis lines)

# Display the MD plot
print(marselli_md_plot)

ggsave("output/marselli_edger_lrt_MD_plot.svg", marselli_md_plot, width = 8.5, height = 6, device = "svg")
ggsave("output/marselli_edger_lrt_MD_plot.pdf", marselli_md_plot, width = 8.5, height = 6, device = "pdf")
ggsave("output/marselli_edger_lrt_MD_plot.tiff", marselli_md_plot, width = 8.5, height = 6, device = "tiff", dpi = "retina")

# Create an MD plot with FDR shown as gradient
colors <- colorRampPalette(c("gray", "red"))(200)

# Calculate the maximum -log10(FDR) value for legend scaling
max_logFDR <- max(-log10(de_genes_lrt_marselli$FDR))

marselli_md_gradient_plot <- ggplot(de_genes_lrt_marselli, aes(x = logCPM, y = logFC)) +
  geom_point(aes(color = -log10(FDR)),
             alpha = 0.7,
             size = 3) +
  geom_text_repel(data = subset(de_genes_lrt_marselli, FDR <= 0.05 & abs(logFC) >= 1),
                  aes(label = genes), color = "black",
                  size = 5,
                  max.overlaps = 20,  # Increase if needed
                  box.padding = 0.5,  # Increase padding around each label
                  point.padding = 0.5,  # Increase padding around each point
                  nudge_x = 0.2,  # Nudge labels slightly in the x direction
                  nudge_y = 0.2) +
  geom_segment(data = subset(de_genes_lrt_marselli, FDR <= 0.05 & abs(logFC) >= 1),
               aes(xend = logCPM, yend = logFC),
               x = NA, y = NA, color = "black", linewidth = 0.5) +  # Add lines
  scale_color_gradientn(colors = colors, limits = c(0, max_logFDR)) +  # Define the color gradient
  theme_minimal() +
  labs(x = "logCPM", y = "log2FC") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  theme(axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 14),
        axis.line = element_line(color = "black", size = 1),
        panel.grid = element_blank(),
        legend.position = c(1.0, 1.0),  # Adjust the legend position within the plot area
        legend.justification = c("right", "top"))  # Position the legend in the upper-right corner


# Add the colorbar legend
marselli_md_gradient_plot + guides(color = guide_colorbar())

# Display the MD plot
print(marselli_md_gradient_plot)

ggsave("output/marselli_edger_lrt_MD_gradient_plot.svg", marselli_md_gradient_plot, width = 6, height = 6, device = "svg")
ggsave("output/marselli_edger_lrt_MD_gradient_plot.pdf", marselli_md_gradient_plot, width = 6, height = 6, device = "pdf")
ggsave("output/marselli_edger_lrt_MD_gradient_plot.tiff", marselli_md_gradient_plot, width = 6, height = 6, device = "tiff", dpi = "retina")


# * 2.2 Marselli individual gene plots ------------------------------------
#'[' ==================================================================================================================
#'[' Visualizing independent gene expression from Marselli edgeR analyzed data as violin or boxplot             
#'[' ==================================================================================================================

gene_name <- "APOE" # replace this gene name with gene of interest, run script to export plots
gene_data <- marselli_cpm[gene_name, ]

# Create a data frame for the box plot
plot_data <- data.frame(Expression = log2(gene_data), conditions = marselli_conditions)

# Create a box plot with overlaid data points
box_plot <- ggplot(plot_data, aes(x = conditions, y = Expression, color = conditions)) +
  geom_boxplot(outlier.shape = NA) +  # Hide outliers
  geom_jitter(width = 0.2, size = 3, alpha = 0.7) +
  labs(x = NA, y = "Normalized Expression",
       title = gene_name,
  ) +
  theme_minimal() +
  theme(axis.text = element_text(size = 14, face = "bold", color = "black"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.title.x = element_markdown(),
        axis.title.y = element_markdown(),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 14),
        axis.line = element_line(color = "black", size = 1),
        panel.grid = element_blank(), 
        legend.position = "none"
  ) +
  scale_color_manual(values = c("ND" = "blue", "T2D" = "red"))

print(box_plot)

ggsave(paste0("output/gene_plots/",gene_name,"_merged_bulk_log2cpm_boxplot.tiff"), box_plot, width = 3, height = 3, device = "tiff", dpi = 300)

# Calculate summary statistics for error bars
error_bar_data <- plot_data %>%
  group_by(conditions) %>%
  summarize(mean = mean(Expression), sd = sd(Expression))

summary_data <- plot_data %>%
  group_by(conditions) %>%
  summarize(mean = mean(Expression), sd = sd(Expression),
            q25 = quantile(Expression, 0.25), q75 = quantile(Expression, 0.75))

# Create a vertical violin plot with overlaid data points, error bars, crossbars, colors, and quartiles
violin_plot <- ggplot() +
  geom_violin(data = plot_data, aes(x = conditions, y = Expression, fill = conditions, alpha = 0.7), trim = TRUE) +
  geom_point(data = plot_data, aes(x = conditions, y = Expression), position = position_jitter(width = 0.08), size = 2, alpha = 0.7, color = "black") +
  geom_errorbar(data = summary_data, aes(x = conditions, y = mean, ymin = mean - sd, ymax = mean + sd),
                width = 0.2, position = position_dodge(0.5), color = "black") +
  geom_crossbar(data = summary_data, aes(x = conditions, y = mean, ymin = q25, ymax = q75),
                width = 0.2, position = position_dodge(0.5), color = "black") +
  labs(x = "Group", y = "Normalized Expression \n (log2CPM)",
       title = paste("Customized Violin Plot with Data Points, Error Bars, and Quartiles for", unique(plot_data$gene_name)),
       subtitle = "Comparison of expression between groups") +
  theme_minimal() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 14),
        legend.title = element_blank(),
        legend.position = "none",
        axis.line = element_line(color = "black", size = 1),
        panel.grid = element_blank()) + # Remove gridlines)  # Add black solid axis lines)
  
  scale_fill_manual(values = c("ND" = "blue", "T2D" = "red"))

print(violin_plot)

# 3. Marselli hits vs Asplund genes ------------------------------------
#'[' ==================================================================================================================
#'[' Compare Marselli edgeR DEGs with T2D genes from Asplund study             
#'[' ==================================================================================================================

# Marselli data objects already loaded
# de_genes_lrt_marselli_cpm
# de_genes_lrt_marselli_counts
# de_genes_lrt_marselli

de_genes_lrt_marselli_filtered <- de_genes_lrt_marselli %>% dplyr::select(gene_symbol = genes, log2fc = logFC, pval = PValue, padj = FDR)

# Annotate according to differential expression
de_genes_lrt_marselli_filtered <- de_genes_lrt_marselli_filtered %>% mutate(diffexpressed = case_when(
  log2fc > 0 & padj < 0.05 ~ 'UP',
  log2fc < 0 & padj < 0.05 ~ 'DOWN',
  padj > 0.05 ~ 'NO'
))

#load Asplund supplemental table
asplund <- read_xlsx("data/asplund/inline-supplementary-material-2.xlsx")

# merge marselli and asplund log2FC data
marselli_asplund_merge <- inner_join(de_genes_lrt_marselli_filtered, asplund, by = c('gene_symbol' = "HGNC symbol"))

marselli_asplund_merge$point_color <- ifelse(marselli_asplund_merge$log2fc >= 0.58, "red", "purple")

marselli_v_asplund_logFC_plot <- ggplot(marselli_asplund_merge, aes(x = log2fc, y = logFC)) +
  geom_point(data = subset(marselli_asplund_merge, abs(log2fc) >= 0.58),  # Filter points outside the range
             aes(color = point_color), 
             alpha = 0.7,
             size = 3) +
  geom_text_repel(data = subset(marselli_asplund_merge, abs(log2fc) >= 0.58), 
                  aes(label = gene_symbol), color = "black",
                  size = 5,
                  max.overlaps = 25,  # Increase if needed
                  box.padding = 0.6,  # Increase padding around each label
                  point.padding = 0.5,  # Increase padding around each point
                  nudge_x = 0.2,  # Nudge labels slightly in the x direction
                  nudge_y = 0.2)  + # Nudge labels slightly in the y direction+
  
  geom_point(data = subset(marselli_asplund_merge, padj < 0.05), # Highlight points with padj < 0.05 in Marselli data - all asplund points are under FDR cutoff already
             alpha = 1,
             size = 3, 
             shape = 21, # Use a filled circle shape
             stroke = 1,
             color = "black") + 
  
  scale_color_manual(values = c("red" = "red", "purple" = "purple"), name = "Point Color") +  # Manually specify colors
  theme_minimal() +
  xlab("log<sub>2</sub> (T2D/ND)<br>Marselli") + 
  ylab("log<sub>2</sub> (T2D/ND)<br>Asplund") + 
  scale_x_continuous(breaks = seq(-5, 5, by = 1)) +  # Adjust the breaks as needed
  theme(axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.title.x = element_markdown(),
        axis.title.y = element_markdown(),
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 14),
        axis.line = element_line(color = "black", size = 1),
        panel.grid = element_blank(),
  )

# Display the plot
print(marselli_v_asplund_logFC_plot)

ggsave("analysis/marselli_asplund_log2fc_comparison.svg", marselli_v_asplund_logFC_plot, width = 8, height = 6, device = "svg")
ggsave("analysis/marselli_asplund_log2fc_comparison.tiff", marselli_v_asplund_logFC_plot,width = 8.5, height = 6, device = "tiff", dpi = "retina")

marselli_asplund_merge_fixlabels <- dplyr::select(marselli_asplund_merge,
                "gene_symbol" = "gene_symbol", 
                "ensemblid" = "genes",
                "marselli_log2fc" = "log2fc",
                "marselli_pval" = "pval",
                "marselli_padj" = "padj",
                "asplund_log2fc" = "logFC",
                "asplund_pval" = "P.Value",
                "asplund_fdr" = "FDR"
                )

write.table(marselli_asplund_merge_fixlabels, "analysis/marselli_asplund_deg_merge.txt", row.names = F, quote = F, sep = '\t')


# 4. Marselli DEGs vs T2DKP list of T2D effectors ------------------------------------
#'[' ==================================================================================================================
#'[' Load the table from T2DKP and compare the intersection of consistent hit genes from Marselli (or marselli+asplund)   
#'[' ==================================================================================================================

t2dkp <- read_csv("data/t2dkp/t2d_effector_summary_251_filtered.csv")

marselli_filtered <- de_genes_lrt_marselli %>%
  dplyr::select(gene_symbol = genes, log2fc = logFC, logCPM = logCPM, pval = PValue, padj = FDR) %>%
  filter(padj <= 0.05)

t2dkp_v_marselli <- inner_join(marselli_filtered, t2dkp, by = c('gene_symbol' = 'Gene'))

write.table(t2dkp_v_marselli, "analysis/t2dkp_vs_marselli.txt", row.names = F, quote = F, sep = '\t')
