setwd("home_dir")

library(Seurat)
library(gtools)
library(dplyr)
library(patchwork)
library(scran)
library(RColorBrewer)

library(viridis)

library(scCustomize)
library(scDblFinder)
library(stringr)

library(Matrix)
library(writexl)
library(readxl)

library(scater)
library(ggplot2)
library(ggthemes)
library(tidyverse)
library(SingleCellExperiment)
library(ggpubr)

library(glmGamPoi)

library(readxl)
library(gtools)

library(SingleR)
library(pals)
library(MAST)
library(harmony)

options(future.globals.maxSize = 60000 * 1024^2)


######################
######################
df.metadata.Xenium<-read_excel("xenium_CO_metadata.xlsx",sheet=2)
######################
######################


spatial_10x_objects<-list()
seurat_10x_objects<-list()

###pre-processing for loop
for (sampleno in 1:nrow(df.metadata.Xenium))   {
print(sampleno)

xenium.obj <- NULL
xenium.obj <- LoadXenium(paste0(df.metadata.Xenium$xenium.maindir[sampleno],"/",df.metadata.Xenium$xenium.subdir[sampleno]), fov = "fov")

xenium.obj@meta.data[["orig.ident.new"]]<-df.metadata.Xenium$orig.ident.new[sampleno]

spatial_10x_objects[[sampleno]]<-xenium.obj
seurat_10x_objects[[sampleno]] <- subset(xenium.obj,nCount_Xenium > 10)

}


for (sampleno in 1:nrow(df.metadata.Xenium))   {

	print(sampleno)

seurat_10x_objects[[sampleno]]@assays[["RNA"]]<-seurat_10x_objects[[sampleno]]@assays[["Xenium"]]

seurat_10x_objects[[sampleno]] <- SCTransform(seurat_10x_objects[[sampleno]], method = "glmGamPoi", verbose = TRUE)

}



##########
##########
##########
##########
integ_features <- SelectIntegrationFeatures(object.list = seurat_10x_objects, nfeatures = 2000)

# Merge normalized samples
merged_seurat <- merge(x = seurat_10x_objects[[1]], y = seurat_10x_objects[2:length(seurat_10x_objects)], merge.data = TRUE)
DefaultAssay(merged_seurat) <- "SCT"

# Manually set variable features of merged Seurat object
VariableFeatures(merged_seurat) <- integ_features

# Calculate PCs using manually set variable features
merged_seurat <- RunPCA(merged_seurat, assay = "SCT", npcs = 20)
##########
##########

harmonized_seurat <- RunHarmony(merged_seurat,group.by.vars = c("orig.ident.new"),reduction = "pca", assay.use = "SCT", reduction.save = "harmony",verbose = TRUE)
##########
##########


############
############
# UMAP and Clustering
harmonized_seurat <- RunUMAP(harmonized_seurat, assay = "SCT", reduction = "harmony", dims = 1:20)
harmonized_seurat <- FindNeighbors(harmonized_seurat,, assay = "SCT", reduction = "harmony", dims = 1:20)
harmonized_seurat <- FindClusters(harmonized_seurat, resolution = 0.3)

#dim(harmonized_seurat@meta.data)
#[1] 101268     10

###########
###########
harmonized_seurat@assays$RNA<-JoinLayers(harmonized_seurat@assays$RNA)
harmonized_seurat<-NormalizeData(harmonized_seurat,assay="RNA")
###########
###########
DefaultAssay(harmonized_seurat)<-"RNA"
###########
###########
Idents(harmonized_seurat)<-"seurat_clusters"
DE_all_celltype.predicted<-NULL
#Identifying the cluster specific markers
DE_all_celltype.predicted<-FindAllMarkers(object = harmonized_seurat, min.pct = 0.10,min.cells.group=50,min.diff.pct=0.1,only.pos =TRUE) #
###########
###########


