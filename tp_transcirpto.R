# -------------------------------
# 1. Librairies et options
# -------------------------------
library(dplyr)
library(Seurat)
library(future)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)

# Options pour le calcul parallèle
plan(sequential)
options(future.globals.maxSize=10*1024**3)

# -------------------------------
# 2. Chargement des données
# -------------------------------
p5 <- read.csv("CAFs/GSM4805570_CountsMatrix_20G00953M_TN.txt.gz", sep="\t")
p4a <- read.csv("CAFs/GSM4805566_CountsMatrix_19G02977A_TN.txt.gz", sep="\t")
p4b <- read.csv("CAFs/GSM4805568_CountsMatrix_19G02977B_TN.txt.gz", sep="\t")

caf.data <- data.matrix(cbind(p5, p4a, p4b))

ens <- read.csv("ensembl-38-108-genes.txt", sep="\t")
ens2symb <- setNames(ens$Gene.name, ens$Gene.stable.ID)
ens2type <- setNames(ens$Gene.type, ens$Gene.stable.ID)

symbols <- ens2symb[rownames(caf.data)]
types <- ens2type[rownames(caf.data)]

# Filtrage des gènes codants uniques
good <- types=="protein_coding" & !is.na(symbols) & !duplicated(symbols)
caf.data <- caf.data[good,]
rownames(caf.data) <- symbols[good]

# -------------------------------
# 3. Création de l'objet Seurat
# -------------------------------
caf <- CreateSeuratObject(counts=caf.data, project="cafs",
                          min.cells=0.01*ncol(caf.data), min.features=1000)

# Pourcentage de gènes mitochondriaux
caf[["percent.mt"]] <- PercentageFeatureSet(caf, pattern="^MT-")

# Vérification de la qualité des cellules
VlnPlot(caf, features=c("nFeature_RNA","nCount_RNA","percent.mt"))

# Filtrage des cellules
caf <- subset(caf, subset = nFeature_RNA > 1000 &
                          nCount_RNA < 50000 &
                          percent.mt < 50)

# -------------------------------
# 4. Normalisation et PCA
# -------------------------------
caf <- NormalizeData(caf)
caf <- FindVariableFeatures(caf)
caf <- ScaleData(caf)
caf <- RunPCA(caf)

# Elbow plot classique
ElbowPlot(caf)

# Elbow plot personnalisé
var_explained <- caf[["pca"]]@stdev^2
percent_var <- var_explained / sum(var_explained) * 100

df_elbow <- data.frame(PC = 1:30, Variance = percent_var[1:30])

ggplot(df_elbow, aes(x = PC, y = Variance)) +
  geom_point(size = 3) +
  geom_line(size = 1) +
  theme_minimal() +
  labs(title = "Elbow Plot (avec ligne)", x = "Composante Principale", y = "% variance expliquée")

# -------------------------------
# 5. Clustering (objectif : 4 clusters)
# -------------------------------
caf <- FindNeighbors(caf, dims=1:20)
caf <- FindClusters(caf, resolution = 0.1)
table(Idents(caf))

# -------------------------------
# 6. t-SNE
# -------------------------------
caf <- RunTSNE(caf, dims=1:20)
DimPlot(caf, reduction="tsne", label=TRUE)

# -------------------------------
# 7. Marqueurs
# -------------------------------
markers <- FindAllMarkers(caf, only.pos = TRUE)
head(markers)

# -------------------------------
# 8. Enrichissement GO pour tous les clusters
# -------------------------------
ego_list <- list()
clusters <- 0:3  # 4 clusters

for (clust in clusters) {
  genes <- markers %>% filter(cluster == clust) %>% pull(gene)
  
  ego <- enrichGO(
    gene = genes,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    readable = TRUE
  )
  
  ego_list[[paste0("cluster", clust)]] <- ego
}

# Affichage des dotplots pour tous les clusters
if (!requireNamespace("gridExtra", quietly = TRUE)) install.packages("gridExtra")
library(gridExtra)

plots <- list()
for (clust in 0:3) {
  p <- dotplot(ego_list[[paste0("cluster", clust)]]) +
       ggplot2::ggtitle(paste("Enrichissement GO - Cluster", clust))
  plots[[clust + 1]] <- p
}

grid.arrange(grobs = plots, ncol = 2, nrow = 2)

# -------------------------------
# 9. Signatures d’origine des CAF
# -------------------------------
vsmc  <- c("PLN","SORBS2","PHLDA2","SNCG","MT1M","MYH11")
hsc   <- c("FABP5","HIGD1B","AGT","RGS5","CPE","SSTR2")
sames <- c("LUM","PTGDS","DCN","COL1A1","FBLN1","LTBP2")

FeaturePlot(caf, features=vsmc)
FeaturePlot(caf, features=hsc)
FeaturePlot(caf, features=sames)

# -------------------------------
# 10. Tableau résumé simple
# -------------------------------
summary_simple <- data.frame(
  cluster = 0:3,
  source_probable = c("VSMC", "HSC", "CAF commun", "VSMC"),  # à remplir selon FeaturePlot / AverageExpression
  fonction_principale = c(
    "Matrice extracellulaire / adhésion", 
    "Métabolisme lipidique", 
    "Remodelage ECM / angiogenèse", 
    "Contraction musculaire"
  )
)

# Afficher le tableau
summary_simple
