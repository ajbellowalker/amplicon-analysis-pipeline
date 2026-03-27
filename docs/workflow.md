### Bash workflow
pwd

conda env create -n qiime2-amplicon-2024.10 --file https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.10-py310-osx-conda.yml
conda activate qiime2-amplicon-2024.10
qiime --version
qiime --help
conda deactivate

python qiime_manifest.py
head -5 qiime_manifest_v2.tsv
cat -vet qiime_manifest_v2.tsv | head

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path qiime_manifest_v2.tsv \
  --output-path paired-end-demux.qza \
  --input-format PairedEndFastqManifestPhred33V2

qiime demux summarize \
  --i-data paired-end-demux.qza \
  --o-visualization demux.qzv

qiime tools view demux.qzv

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs paired-end-demux.qza \
  --p-trunc-len-f 240 \
  --p-trunc-len-r 220 \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-n-threads 0 \
  --p-n-reads-learn 250000 \
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats stats.qza

qiime metadata tabulate \
  --m-input-file stats.qza \
  --o-visualization stats.qzv

qiime tools view stats.qzv

qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv

qiime tools view table.qzv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv

qiime tools view rep-seqs.qzv

qiime feature-table filter-samples \
  --i-table table.qza \
  --p-min-frequency 1000 \
  --o-filtered-table table-filtered.qza

qiime feature-table summarize \
  --i-table table-filtered.qza \
  --o-visualization table-filtered.qzv

qiime tools view table-filtered.qzv

qiime tools export \
  --input-path table-filtered.qza \
  --output-path exported-table

biom convert \
  -i exported-table/feature-table.biom \
  -o exported-table/feature-table.tsv \
  --to-tsv

head -1 exported-table/feature-table.tsv

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza

qiime metadata tabulate \
  --m-input-file metadata_clean.tsv \
  --o-visualization metadata_check.qzv

qiime tools view metadata_check.qzv

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table table.qza \
  --p-sampling-depth 18000 \
  --m-metadata-file metadata_clean.tsv \
  --output-dir core-metrics-results

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/shannon_vector.qza \
  --m-metadata-file metadata_clean.tsv \
  --o-visualization shannon_groups.qzv

qiime tools view shannon_groups.qzv

qiime tools export \
  --input-path core-metrics-results/bray_curtis_distance_matrix.qza \
  --output-path exported-bray

conda env create \
  -n qiime2-2023.9 \
  --file https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2023.9-py38-osx-conda.yml
conda activate qiime2-2023.9

qiime --version

wget https://data.qiime2.org/2023.9/common/silva-138-99-515-806-nb-classifier.qza

qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-99-515-806-nb-classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza \
  --p-n-jobs 2

qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv

qiime tools view taxonomy.qzv

qiime taxa barplot \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file metadata_clean.tsv \
  --o-visualization taxa-barplot.qzv

qiime tools view taxa-barplot.qzv

qiime tools export \
  --input-path taxonomy.qza \
  --output-path taxonomy_export

conda deactivate
conda activate qiime2-amplicon-2024.10

qiime tools export \
  --input-path table.qza \
  --output-path table_export

biom convert \
  -i table_export/feature-table.biom \
  -o table_export/feature-table.tsv \
  --to-tsv

### R workflow
install.packages("vegan")  
install.packages("readr")
install.packages("BiocManager")
BiocManager::install("DESeq2", ask = FALSE, update = FALSE)
BiocManager::install("apeglm", ask = FALSE, update = FALSE)

library(phyloseq)
library(apeglm)
library(biomformat)
library(tidyverse)
library(DESeq2)
library(vegan)
library(ggplot2)
library(readr)

bray <- read.table("exported-bray/distance-matrix.tsv",
                   header = TRUE,
                   row.names = 1,
                   check.names = FALSE)

bray_dist <- as.dist(bray)

metadata <- read_tsv("metadata_clean.tsv")
metadata$Cow       <- factor(metadata$Cow)
metadata$Period    <- factor(metadata$Period)
metadata$Silage    <- factor(metadata$Silage)
metadata$Starch    <- factor(metadata$Starch)
metadata$Timepoint <- factor(metadata$Timepoint)
metadata$SampleType <- factor(metadata$SampleType)
metadata$Group     <- factor(metadata$Group)
metadata$Interaction <- factor(metadata$Interaction)
str(metadata)

metadata <- metadata[match(labels(bray_dist), metadata$`sample-id`), ]
all(labels(bray_dist) == metadata$`sample-id`)

metadata_clean <- metadata[metadata$SampleType == "EPM", ]
bray_mat <- as.matrix(bray_dist)
bray_clean <- bray_mat[
  metadata_clean$`sample-id`,
  metadata_clean$`sample-id`
]
bray_clean <- as.dist(bray_clean)
head(rownames(bray_mat))
head(metadata_clean$`sample-id`)

adonis2(bray_clean ~ Starch * Silage + Timepoint,
        data = metadata_clean,
        permutations = 999,
        strata = metadata_clean$Cow)

adonis2(bray_clean ~ Starch * Silage + Timepoint,
        data = metadata_clean,
        permutations = 999,
        strata = metadata_clean$Cow,
        by = "terms")

bd_starch <- betadisper(bray_clean, metadata_clean$Starch)
anova(bd_starch)

bd_silage <- betadisper(bray_clean, metadata_clean$Silage)
anova(bd_silage)

set.seed(123)

nmds <- metaMDS(bray_clean, k = 2, trymax = 100)
nmds <- metaMDS(bray_clean, k = 3)
nmds$stress

nmds_scores <- as.data.frame(scores(nmds))
nmds_scores$Starch <- metadata_clean$Starch
nmds_scores$Silage <- metadata_clean$Silage
nmds_scores$Timepoint <- metadata_clean$Timepoint
nmds_scores$Cow <- metadata_clean$Cow

ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = Starch)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(aes(group = Starch), linetype = 2) +
  theme_classic() +
  labs(title = "NMDS – Bray-Curtis (Epimural Microbiome)",
       subtitle = "Colored by Starch Level")

# Read taxonomy
tax_raw <- read_tsv("taxonomy_export/taxonomy.tsv", show_col_types = FALSE)

# Split taxonomy string into ranks
tax_split <- tax_raw %>%
  separate(Taxon,
           into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),
           sep = ";",
           fill = "right")

# Remove prefix (e.g. D_0__, k__, etc.)
tax_split <- tax_split %>%
  mutate(across(Kingdom:Species, ~ gsub("^[A-Za-z]_[0-9]+__", "", .))) %>%
  mutate(across(Kingdom:Species, ~ gsub("^[a-z]__", "", .)))

# Convert to matrix
tax_mat <- as.matrix(tax_split[,2:8])
rownames(tax_mat) <- tax_split$`Feature ID`
tax_table_ps <- tax_table(tax_mat)

# Read feature table
otu_raw <- read_tsv(
  "table_export/feature-table.tsv",
  skip = 1   # skip the "# Constructed from biom file" line
)

# Convert to matrix
otu_mat <- as.matrix(otu_raw[,-1])
rownames(otu_mat) <- otu_raw[[1]]

# Make sure values are numeric
otu_mat <- apply(otu_mat, 2, as.numeric)
rownames(otu_mat) <- otu_raw[[1]]
otu_table_ps <- otu_table(otu_mat, taxa_are_rows = TRUE)

# Make sample IDs rownames
metadata_df <- as.data.frame(metadata)
rownames(metadata_df) <- metadata_df$`sample-id`

# Remove the sample-id column (phyloseq doesn't need it twice)
metadata_df$`sample-id` <- NULL

# Convert to phyloseq sample_data object
sampledata_ps <- sample_data(metadata_df)
sampledata_ps

ps <- phyloseq(otu_table_ps, tax_table_ps, sampledata_ps)
ps

ps_clean <- subset_samples(ps, SampleType == "EPM")
ps_clean <- prune_taxa(taxa_sums(ps_clean) > 0, ps_clean)

dds <- phyloseq_to_deseq2(ps_clean, ~ Cow + Starch + Silage + Timepoint)
dds <- DESeq(dds)

res_starch <- results(dds, contrast = c("Starch","High","Low"))
res_starch <- lfcShrink(dds,
                        coef = "Starch_Low_vs_High",
                        type = "apeglm")

sig_asvs <- res_starch %>%
  as.data.frame() %>%
  rownames_to_column("ASV") %>%
  filter(padj < 0.05)

sig_starch <- res_starch[which(res_starch$padj < 0.05), ]
sig_starch <- as.data.frame(sig_starch)
sig_starch <- cbind(sig_starch,
                    as(tax_table(ps_clean)[rownames(sig_starch), ], "matrix"))

head(sig_starch)

sum(res_starch$padj < 0.05, na.rm = TRUE)

res_df <- as.data.frame(res_starch)
res_df$ASV <- rownames(res_df)

tax_df <- as.data.frame(tax_table(ps_clean))
tax_df$ASV <- rownames(tax_df)

res_annotated <- merge(res_df, tax_df, by = "ASV")

# Remove leading spaces and __ prefix
res_annotated$Genus <- gsub(" g__", "", res_annotated$Genus)
res_annotated$Genus <- gsub("^g__", "", res_annotated$Genus)
res_annotated$Genus <- trimws(res_annotated$Genus)

res_annotated$Family <- gsub(" f__", "", res_annotated$Family)
res_annotated$Family <- gsub("^f__", "", res_annotated$Family)
res_annotated$Family <- trimws(res_annotated$Family)

res_annotated$Order <- gsub(" o__", "", res_annotated$Order)
res_annotated$Order <- gsub("^o__", "", res_annotated$Order)
res_annotated$Order <- trimws(res_annotated$Order)

res_annotated$Class <- gsub(" c__", "", res_annotated$Class)
res_annotated$Class <- gsub("^c__", "", res_annotated$Class)
res_annotated$Class <- trimws(res_annotated$Class)

res_annotated$Phylum <- gsub(" p__", "", res_annotated$Phylum)
res_annotated$Phylum <- gsub("^p__", "", res_annotated$Phylum)
res_annotated$Phylum <- trimws(res_annotated$Phylum)

res_annotated$Species <- gsub(" s__", "", res_annotated$Species)
res_annotated$Species <- gsub("^s__", "", res_annotated$Species)
res_annotated$Species <- trimws(res_annotated$Species)

# Extract Strach Significant Group
sig_starch <- subset(res_annotated, padj < 0.05)
head(sig_starch[order(sig_starch$padj), ])
table(sig_starch$log2FoldChange > 0)

sig_starch %>%
  group_by(Genus) %>%
  summarise(n_ASVs = n()) %>%
  arrange(desc(n_ASVs)) %>%
  print(n = 20) # Print top 20 genera

sig_starch %>%
  group_by(Phylum) %>%
  summarise(n = n(),
            mean_LFC = mean(log2FoldChange)) %>%
  arrange(desc(n))
  
res_df$Significant <- ifelse(res_df$padj < 0.05, "Yes", "No")

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = Significant), alpha = 0.6) +
  theme_classic() +
  labs(title = "Differential Abundance: Low vs High Starch",
       x = "Log2 Fold Change (Low vs High)",
       y = "-log10 adjusted p-value")

prevotella_hits <- subset(res_annotated, Genus == "Prevotella")
prevotella_hits[order(prevotella_hits$padj), ]
table(prevotella_hits$log2FoldChange > 0)

# Extract Urea Significant Group
urea_candidates <- c("Succinivibrio",
                     "Selenomonas",
                     "Campylobacter",
                     "Neisseriaceae",
                     "Ruminobacter",
                     "Treponema",
                     "Desulfovibrio")
urea_candidates_hits <- subset(res_annotated, Genus %in% urea_candidates)
urea_candidates_hits[order(urea_candidates_hits$padj), ]
table(urea_candidates_hits$log2FoldChange > 0)

urea_hits <- subset(res_annotated,
                    Genus %in% c("Ruminobacter","Selenomonas","Succinimonas","Treponema"))

urea_hits[, c("Genus","log2FoldChange","padj")]
table(urea_hits$log2FoldChange > 0)

prev_sig <- subset(res_annotated,
                   Genus == "Prevotella" & padj < 0.05)

table(prev_sig$log2FoldChange > 0)

urea_sig <- subset(res_annotated,
                   Genus %in% c("Ruminobacter","Selenomonas",
                                "Succinimonas","Treponema") &
                   padj < 0.05)

table(urea_sig$log2FoldChange > 0)

median(prev_sig$log2FoldChange)
median(urea_sig$log2FoldChange)

proteo_sig <- subset(sig_starch, Phylum == "Proteobacteria")
proteo_sig[, c("Genus","log2FoldChange","padj")]
table(proteo_sig$log2FoldChange > 0)

mean(abs(proteo_sig$log2FoldChange))
mean(abs(prev_sig$log2FoldChange))

# Aggregate phyloseq to Genus
ps_genus <- tax_glom(ps_clean, taxrank = "Genus")

# Convert to DESeq2
dds_genus <- phyloseq_to_deseq2(ps_genus, ~ Cow + Starch + Silage + Timepoint)

dds_genus <- DESeq(dds_genus)

res_genus <- results(dds_genus, contrast = c("Starch","High","Low"))

res_genus <- lfcShrink(
  dds_genus,
  coef = "Starch_Low_vs_High",
  type = "apeglm"
)

res_genus_df <- as.data.frame(res_genus)
res_genus_df$Genus <- tax_genus$Genus
res_genus_hits <- subset(res_genus_df, padj < 0.05)

res_genus_hits %>%
  mutate(direction = ifelse(log2FoldChange > 0, "Low starch", "High starch")) %>%
  group_by(direction) %>%
  summarise(mean_abs_LFC = mean(abs(log2FoldChange)))

res_genus_hits %>%
  arrange(log2FoldChange) %>%
  head(10)

# Attach taxonomy
tax_genus <- as.data.frame(tax_table(ps_genus))
res_genus_df$Genus <- tax_genus$Genus

# Remove leading spaces and __ prefix
res_genus_df$Genus <- gsub(" g__", "", res_genus_df$Genus)
res_genus_df$Genus <- gsub("^g__", "", res_genus_df$Genus)
res_genus_df$Genus <- trimws(res_genus_df$Genus)

head(res_genus_df[order(res_genus_df$padj), ])
colnames(res_genus_df)
sum(res_genus_df$padj < 0.05, na.rm = TRUE)
res_genus_df_hits <- subset(res_genus_df, padj < 0.05)[order(res_genus_df$padj), ]
table(res_genus_df_hits$log2FoldChange > 0)

res_genus_df_hits %>%
  arrange(log2FoldChange) %>%
  select(Genus, log2FoldChange, padj)

res_genus_df_hits %>%
  mutate(direction = ifelse(log2FoldChange > 0, "Low starch", "High starch")) %>%
  group_by(direction) %>%
  summarise(mean_abs_LFC = mean(abs(log2FoldChange)),
            n = n()) %>%
  arrange(desc(n))

tax <- tax_table(ps_clean)

# Replace NA in Phylum
tax[,"Phylum"][is.na(tax[,"Phylum"])] <- "Unassigned"

tax_table(ps_clean) <- tax

# Transform to relative abundance
ps_rel <- transform_sample_counts(ps_clean, function(x) x / sum(x))

tax_df <- as.data.frame(tax_table(ps_rel))
tax_df$Phylum <- gsub(" p__", "", tax_df$Phylum)
tax_df$Phylum <- gsub("^p__", "", tax_df$Phylum)
tax_df$Phylum <- trimws(tax_df$Phylum)

proteo_taxa <- rownames(tax_df)[tax_df$Phylum == "Proteobacteria"]
length(proteo_taxa)
ps_proteo <- prune_taxa(proteo_taxa, ps_rel)

# Sum Proteobacteria abundance per sample
proteo_abund <- sample_sums(ps_proteo)

# Convert to dataframe
proteo_df <- data.frame(
  sample = names(proteo_abund),
  Proteobacteria_RA = as.numeric(proteo_abund)
)

# assuming phyloseq object is called ps
alpha_df <- estimate_richness(ps_clean, measures = c("Shannon", "Observed"))
alpha_df$sample_id <- rownames(alpha_df)

alpha_df <- merge(alpha_df,
                  metadata_clean,
                  by.x = "sample_id",
                  by.y = "sample-id")
colnames(alpha_df)

model_alpha <- lmer(Shannon ~ Starch + Silage + Timepoint + (1|Cow),
                    data = alpha_df)
summary(model_alpha)

aov_shannon <- aov(Shannon ~ Starch * Silage + Timepoint, data = alpha_df)
summary(aov_shannon)
TukeyHSD(aov_shannon)

