#LOAD library ----
# Core packages
install.packages("vegan")
install.packages("phyloseq")
install.packages("ggfortify")
install.packages("caret")
install.packages("adegenet")
install.packages("igraph")
install.packages("randomForest")
library(tidyverse)
library(vegan)
library(phyloseq)
library(adegenet)
library(igraph)
library(caret)
library(randomForest)
library(ggplot2)
library(ggfortify)

# For biomarker (LEfSe-style)
library(microbiomeMarker)

# For BiocManager
install.packages("BiocManager")
BiocManager::install("phyloseq")

library(phyloseq)
library(MASS)
library(dplyr)

## Load Data ----
# Read cleaned data table (tab separated). Replace path if needed.
data <- read_xslx("C:/Users/ADMIN/Desktop/GENUS.xlsx")

# Identify microbial columns (adjust if your column names differ) ----
microbe_start <- which(colnames(GENUS) == "Alistipes")
microbe_end <- which(colnames(GENUS) == "Victivallaceae")
microbe_cols <- microbe_start:microbe_end

# Microbe abundance table ----
microbe_abund <- GENUS %>% dplyr::select(all_of(microbe_cols))
rownames(microbe_abund) <- GENUS$Name

# Relative abundance ----
rel_abund <- microbe_abund / rowSums(microbe_abund)

# Metadata ----
meta <- GENUS %>% dplyr::select(Name, Breed, Sex, BW, WL, SL, BL, TL, KL, Chao1, Shannon, Simpson, Pielou)
rownames(meta) <- meta$Name

# Quick check ----
cat("Loaded", nrow(GENUS), "samples and", ncol(microbe_abund), "taxa.\n")

## Microbiome Composition ----
# Top 10 genera across all samples
top10 <- rel_abund %>% dplyr::summarise_all(mean, na.rm=TRUE) %>%
  pivot_longer(cols = dplyr::everything(), names_to = "GENUS1", values_to = "MeanRA") %>%
  dplyr::arrange(desc(MeanRA)) %>% dplyr::slice(1:10)

p_top10 <- ggplot(top10, aes(x=reorder(GENUS1, -MeanRA), y=MeanRA)) +
  geom_col(fill = "steelblue") +
  coord_flip() + theme_minimal() +
  labs(title = "Top 10 Genera in FUNAAB Alpha Chickens", y="Mean Relative Abundance", x="GENUS1")

p_top10

## Diversity Indices by Breed ----
div <- GENUS %>% dplyr::select(Breed, Chao1, Shannon, Simpson, Pielou)
div_long <- div %>% pivot_longer(-Breed, names_to = "Index", values_to = "Value")
p_div <- ggplot(div_long, aes(Breed, Value, fill = Breed)) +
  geom_boxplot() +
  facet_wrap(~Index, scales = "free_y") +
  theme_minimal() +
  labs(title = "Diversity Indices by Breed")
p_div

# Statistical tests: Kruskal-Wallis ----
alpha_metrics <- c("Chao1", "Shannon", "Simpson", "Pielou")
kw_res <- map(alpha_metrics, ~kruskal.test(as.formula(paste(.x, "~ Breed")), data=div))
names(kw_res) <- alpha_metrics
kw_res

## PCA of Microbiome ----
# PCA on centered log-ratio transformed data (add pseudocount)
clr_abund <- apply(rel_abund, 2, function(x) x)

# For PCA we use Hellinger or CLR. We'll use CLR after adding small pseudocount to raw counts.
# If microbe_abund are counts, prefer CLR on counts+1. If proportions, use Hellinger.
if(all(rel_abund <= 1)) {
  hell <- decostand(as.matrix(microbe_abund), method="hellinger")
  pca <-prcomp(hell, center=TRUE, scale.=TRUE)
} else {
  clr_mat <- clr(as.matrix(microbe_abund + 1))
  pca <- prcomp(clr_mat, center=TRUE, scale.=TRUE)
}

pca_plot <- autoplot(pca, data = GENUS, colour = "Breed", shape = "Sex") +
  labs(title = "PCA of Microbiome Composition")

pca_plot

## DAPC (Breed as Group) ----
library(compositions)

# Factor for grouping
grp <- as.factor(GENUS$Breed)

# Prevalence filter
otu_samp <- as.data.frame(microbe_abund)
prev.prop <- apply(otu_samp, 2, function(x) mean(x > 0, na.rm = TRUE))
keep <- which(prev.prop >= 0.2)

# CLR transform
otu_keep <- otu_samp[, keep, drop=FALSE]
otu_clr <- as.data.frame(clr(as.matrix(otu_keep + 1)))

# Run DAPC
set.seed(123)
dapc_res <- dapc(otu_clr, grp = grp, n.pca = 50, n.da = length(unique(grp)) -1)

# DAPC scatterplot
scatter(
  dapc_res, 
  scree.da = FALSE,
  posi.da = "bottomleft", 
  bg = "white",
  pch = 19,
  col = rainbow(length(unique(grp)))
)
compoplot(
  dapc_res, 
  col = rainbow(length(unique(grp))),
  show.lab = FALSE
)

# Screeplot for PCA eigenvalues 
barplot(
  dapc_res$eig,
  main = "Discriminant analysis eigenvalues",
  col = "skyblue"
)

## Microbiome–Growth Associations ----
# Use relative abundances for correlations
traits <- GENUS %>% dplyr::select(BW, WL, SL, BL, TL, KL)
abund_df <- as.data.frame(rel_abund)


# Spearman correlations of taxa vs BW (example) ----
cor_res <- apply(abund_df, 2, function(x) cor.test(x, traits$BW, method = "spearman", exact = FALSE))
cor_p <- sapply(cor_res, function(x) x$p.value)
cor_rho <- sapply(cor_res, function(x) x$estimate)
cor_df <- data.frame(taxa = names(cor_rho), rho = unlist(cor_rho), pval = cor_p, padj = p.adjust(cor_p, method = "fdr"))
cor_sig <- cor_df %>% dplyr::arrange(padj) %>% filter(padj < 0.05)
cor_sig
cors <- cor(abund_df, traits, method = "spearman")

# Heatmap of correlations (taxa x traits) ----
install.packages("pheatmap")
library(pheatmap)
pheatmap::pheatmap(cors, main = "Spearman correlations (taxa vs traits)")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("microbiomeMarker")
install.packages("microbiomeMarker")
library(microbiomeMarker)
library(phyloseq)

## Biomarker Identification (LEfSe-style with microbiomeMarker) ----
# Ensure microbe_abund is non-empty
OTU <- otu_table(as.matrix(microbe_abund), taxa_are_rows=FALSE)

meta_df <- as.data.frame(meta)
rownames(meta_df) <-meta_df$Name
meta_df$Name <- NULL
SAM <- sample_data(meta_df)

# Taxonomy table ----
tax_df <- dplyr::data_frame(GENUS = colnames(microbe_abund))
rownames(tax_df) <- tax_df$GENUS
TAX <- tax_table(as.matrix(tax_df))

# Build phyloseq object
physeq <- phyloseq(OTU, TAX, SAM)

# Run LEfSe
args(run_lefse)
tax_table(physeq)
dim(tax_table(physeq))
colnames(tax_table(physeq))
library(tidyr)
library(dplyr)
tax_df <- as.data.frame(tax_table(physeq))
tax_df_split <- tax_df %>%
  tidyr::separate(GENUS, 
                  into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                  sep = ",",
                  fill = "right", 
                  remove = TRUE)
tax_table(physeq) <- as.matrix(tax_df_split)
colnames(tax_table(physeq))
lefse_res <- run_lefse(
  physeq, 
  group = "Breed", 
  lda_cutoff = 2, 
  norm = "CPM"
)

# Extract results
res_df <- marker_table(lefse_res)

# Plot
library(ggplot2)
ggplot(res_df, aes(x = reorder(feature, ef_lda), y = ef_lda, fill = enrich_group)) +
  geom_col() +
  coord_flip() +
  labs(x = "Taxa", y = "LDA score(log10)", title = "LEFSe Discriminant Taxa") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top")

# Save results
write.csv(res_df, "lefse_markers_table.csv")

## Predict Growth Traits from Microbiome
set.seed(123)
BiocManager::install("microbiome")
library(microbiome)
library(compositions)
library(caret)

# Prepare features: CLR-transformed taxa (use otu_keep from earlier) ----
X <- as.data.frame(compositions::clr(as.matrix(otu_keep +1)))
Y <- GENUS$BW
df_model <-data.frame(BW = Y, X)
train_ctrl <- trainControl(method = 'repeatedcv',number = 5, repeats = 3)

# Random Forest
rf_fit <- train(BW ~ .,data = df_model, method ='rf', trControl = train_ctrl, importance =TRUE)
print(rf_fit)

# Predictions (resubstitution here; for proper test do train/test split)
preds <- predict(rf_fit, df_model)
resid <- df_model$BW - preds
rmse <- sqrt(mean(resid^2, na.rm=TRUE))
cat('RF RMSE (resubstitution):', rmse, '\n')
library(tibble)

# Variable importance
imp <- varImp(rf_fit)$importance
imp_top <- imp %>% rownames_to_column(var = 'taxa') %>% arrange(desc(Overall)) %>% slice(1:20)
p_imp <- ggplot(imp_top,aes(x=reorder(taxa, Overall), y=Overall)) + geom_col() + coord_flip() + theme_minimal() + labs(title='Top 20 Important Taxa for BW')
print(p_imp)

## Microbe–Microbe Network ----
# Compute spearman correlations
cors_micro <- cor(as.matrix(rel_abund), method = 'spearman')

# keep strong edges
thr <- 0.6
adj <- cors_micro
adj[abs(adj) < thr] <- 0
diag(adj) <- 0
library(igraph)

# Build igraph
g <- graph_from_adjacency_matrix(adj, mode='undirected', weighted=TRUE, diag=FALSE)
V(g)$name <- colnames(rel_abund)
V(g)$abund <- colMeans(rel_abund)

# Save network plot
png('microbe_network.png', width=1000, height=800)
plot(g, vertex.size = (V(g))$abund * 30) + 3, vertex.label.cex = 0.8, edge.width = abs(E(g)$weight)*3, main='Microbe Co-occurrence Network')
cat('Network plot saved to microbe_network.png\n')

#Linear Regression Models: Growth Traits vs. Genera ----
# Loop for summary of multiple dependent variables (BW, WL, etc.)
dep_vars <- c("BW", "WL", "SL", "BL", "TL", "KL")  
regmodels <- lapply(dep_vars, function (y) { 
  formula <- as.formula(paste(y, "~GENUS1"))
  lm(formula, data = GENUS1)
})
lapply(regmodels, summary)
dep_vars_1 <- lm(BW ~ GENUS, data = GENUS)
str(GENUS)
library(dplyr)
GENUS2 <- GENUS[, 10:46]
GENUS2 <- GENUS %>% dplyr::select(10:46)
names(GENUS2)
GENUS1 <- GENUS[, 4:9]
GENUS1 <- GENUS %>% dplyr::select(4, 9)
sub_GENUS <- data.frame(GENUS1 = GENUS$GENUS1, GENUS2 = GENUS$GENUS2)
model_GENUS <- lm(GENUS1 ~ GENUS2, data = sub_GENUS)
subdata1 <- GENUS[, c("BW", "Alistipes":"Victivallaceae"),
                  na.action = na.omit()]
mod_Gen <- lm(BW ~ GENUS[10:46], data = GENUS)
reg_GENUS <- lm(GENUS1 ~ GENUS2, data = GENUS)
head(GENUS)
str(GENUS)
model <- lm(BW ~ ., data = GENUS[c("BW", Alistipes + Anaerostipes + Bacteroides + Barnesiella + Blautia)])
model1 <- lm(BW ~ Alistipes + Anaerostipes + Bacteroides + Barnesiella + Blautia + Campylobacter + CHKCI001 + Christensenellaceae_R-7_group + Clostridia_UCG-014 + Clostridia_vadinBB60_group + Colidextribacter + Collinsella + Enterococcus + Eubacterium_hallii_group + Faecalibacterium + Frisingicoccus + Fusicatenibacter + Gastranaerophilales + Helicobacter + Lachnoclostridium + Lactobacillus + NK4A214_group + Olsenella + Parabacteroides + Phascolarctobacterium + RF39 + Rikenella + Rothia + Ruminococcus_torques_group + Slackia + Streptococcus + Subdoligranulum + Sutterella + Synergistes + UCG-005 + Uncultured + Victivallaceae, data = GENUS)

#Detailed Linear Models for Specific Traits ----
# Models 1-6 focusing on individual growth traits
model1 <- lm(GENUS$BW ~ ., data = GENUS[, 10:46])
summary(model1)
model2 <- lm(GENUS$WL ~ ., data = GENUS[, 10:46])
summary(model2)
model3 <- lm(GENUS$SL ~ ., data = GENUS[, 10:46])
summary(model3)
model4 <- lm(GENUS$BL ~ ., data = GENUS[, 10:46])
summary(model4)
model5 <- lm(GENUS$TL ~ ., data = GENUS[, 10:46])
summary(model5)
model6 <- lm(GENUS$KL ~ ., data = GENUS[, 10:46])
summary(model6)

#Microbiome Composition & PCA Analysis ----
# Data cleaning (removing zero variance) and PCA calculation
hell_clean <- hell[, apply(hell, 2, var) != 0]
pca <- prcomp(hell_clean, center = TRUE, scale. = TRUE)
K <- min(nrow(x), ncol(x))

#Visualization ----
# PCA plot colored by Breed and Sex
pca_plot <- autoplot(pca, data=data, colour="Breed",shape="Sex") +
  labs(title="PCA of Microbiome composition")

#Relative Abundance Calculations ----
# (Your sweep and abundance table code)
str(GENUS)
abund_table <- GENUS$Abundance
abund_table <- as.data.frame(GENUS$Abundance)
rel_abund <- sweep(as.matrix(abund_table), 2, colSums(abund_table), FUN = "/")
GENUS$relative_abundance <- rel_abund


pca_plot

traits <- GENUS %>% dplyr::select(BW, WL, SL, BL, TL, KL)
abund_df <- as.data.frame