options(repos = c(CRAN = "https://cloud.r-project.org"))

# Installation des packages nécessaires
install.packages("dplyr")
install.packages("Seurat")
install.packages("future")

# Chargement des bibliothèques
library(dplyr)
library(Seurat)
library(future)

# Configuration pour Seurat
plan(sequential)
options(future.globals.maxSize = 10*1024^3)

# Chargement des données
cat("Chargement des données...\n")
p5 <- read.csv("CAFs/GSM4805570_CountsMatrix_20G00953M_TN.txt.gz", sep = "\t")
p4a <- read.csv("CAFs/GSM4805566_CountsMatrix_19G02977A_TN.txt.gz", sep = "\t")
p4b <- read.csv("CAFs/GSM4805568_CountsMatrix_19G02977B_TN.txt.gz", sep = "\t")

caf.data <- data.matrix(cbind(p5, p4a, p4b))
cat("Données chargées: ", dim(caf.data), "\n")

# Préparation des données
cat("Préparation des données...\n")
ens <- read.csv("ensembl-38-108-genes.txt", sep = "\t")

ens2symb <- setNames(ens$Gene.name, ens$Gene.stable.ID)
ens2type <- setNames(ens$Gene.type, ens$Gene.stable.ID)

symbols <- ens2symb[rownames(caf.data)]
types   <- ens2type[rownames(caf.data)]

good <- types == "protein_coding" & !is.na(symbols) & !duplicated(symbols)
cat("Gènes retenus:", sum(good), "\n")

symb <- symbols[good]
caf.data <- caf.data[good,]
rownames(caf.data) <- symb

# Création de l'objet Seurat
cat("Création de l'objet Seurat...\n")
caf <- CreateSeuratObject(
  counts = caf.data,
  project = "cafs",
  min.cells = 0.01 * ncol(caf.data),
  min.features = 1000
)

print(caf)
cat("Analyse terminée avec succès!\n")
