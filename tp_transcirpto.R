library(dplyr)
library(Seurat)
# parallel computation and memory allocation for Seurat
library(future)
plan(sequential)
options(future.globals.maxSize=10*1024**3)
p5 <- read.csv("CAFs/GSM4805570_CountsMatrix_20G00953M_TN.txt.gz",
               sep="\t")
p4a <- read.csv("CAFs/GSM4805566_CountsMatrix_19G02977A_TN.txt.gz",
                sep="\t")
p4b <- read.csv("CAFs/GSM4805568_CountsMatrix_19G02977B_TN.txt.gz",
                sep="\t")
caf.data <- data.matrix(cbind(p5,p4a,p4b))
ens <- read.csv("ensembl-38-108-genes.txt", sep="\t")
ens2symb <- setNames(ens$Gene.name, ens$Gene.stable.ID)
ens2type <- setNames(ens$Gene.type, ens$Gene.stable.ID)
symbols <- ens2symb[rownames(caf.data)]
types <- ens2type[rownames(caf.data)]
good <- types=="protein_coding" & !is.na(symbols) & !duplicated(symbols)
sum(good)
symb <- symbols[good]
caf.data <- caf.data[good,]
rownames(caf.data) <- symb

#objet
caf <- CreateSeuratObject(counts=caf.data, project="cafs",
                          min.cells=0.01*ncol(caf.data), min.features=1000)
caf



caf[["percent.mt"]] <- PercentageFeatureSet(caf, pattern="^MT-")

VlnPlot(caf, features=c("nFeature_RNA","nCount_RNA","percent.mt"))

caf <- subset(
  caf,
  subset = nFeature_RNA > 1000 &
    nCount_RNA < 50000 &
    percent.mt < 50
)

### --- 3. Normalisation & PCA ----
### --- Normalisation, HVG, scaling, PCA ----
caf <- NormalizeData(caf)
caf <- FindVariableFeatures(caf)
caf <- ScaleData(caf)
caf <- RunPCA(caf)

### --- ElbowPlot classique (Seurat)
ElbowPlot(caf)

### --- ElbowPlot personnalisé AVEC LIGNE ----
# Extraire la variance expliquée
var_explained <- caf[["pca"]]@stdev^2
percent_var <- var_explained / sum(var_explained) * 100

# Prendre les 30 premiers PCs
df_elbow <- data.frame(
  PC = 1:30,
  Variance = percent_var[1:30]
)

# Graphe avec points + ligne
library(ggplot2)

ggplot(df_elbow, aes(x = PC, y = Variance)) +
  geom_point(size = 3) +
  geom_line(size = 1) +
  theme_minimal() +
  labs(
    title = "Elbow Plot (avec ligne)",
    x = "Composante Principale",
    y = "% variance expliquée"
  )


### --- 4. Clustering (objectif : 4 clusters) ----
caf <- FindNeighbors(caf, dims=1:20)
caf <- FindClusters(caf, resolution = 0.1)
table(Idents(caf))



### --- 6. t-SNE ----
caf <- RunTSNE(caf, dims=1:20)
DimPlot(caf, reduction="tsne", label=TRUE)

### --- 7. Marqueurs ----
markers <- FindAllMarkers(caf, only.pos = TRUE)
head(markers)

### --- 8. GO Enrichissement (ex : cluster 0)
genes_c0 <- markers %>% filter(cluster==0) %>% pull(gene)

ego <- enrichGO(
  gene = genes_c0,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  readable = TRUE
)

head(ego)
dotplot(ego)

### --- 7. Signatures d’origine des CAF ----
vsmc  <- c("PLN","SORBS2","PHLDA2","SNCG","MT1M","MYH11")
hsc   <- c("FABP5","HIGD1B","AGT","RGS5","CPE","SSTR2")
sames <- c("LUM","PTGDS","DCN","COL1A1","FBLN1","LTBP2")

FeaturePlot(caf, features=vsmc)
FeaturePlot(caf, features=hsc)
FeaturePlot(caf, features=sames)
