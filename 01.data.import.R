.libPaths('/data/00.software/R/4.3.3/library/')
getwd()
# setwd('/data/02.project/00.other.lab/01.hualun/03.mus.lp_hp.adipo/03.analysis/')

options(scipen = 3)
options(future.globals.maxSize = 100000 * 1024^2)
library(Seurat)
library(scDblFinder)
library(SoupX)
library(rio)
library(tidyverse)
library(future)
library(ggplot2)
library(plotly)
library(ggsci)
library(stringr)
library(data.table)
library(org.Mm.eg.db)
nbrOfWorkers()
plan(multisession, workers=50)
use.cols <- c('#DA87B6','#A0CC58','#8F4C9A','#E47C7B','#658DC9','#D6251F','#1E2D5B','#EF7D1C','#4EA74A','#F5D23A','#A15528','#BA1D7C','#0D7C7C')
npg.cols <- pal_npg(palette = 'nrc')(10)
d3.cols <- pal_d3("category10")(10)

# 01. data import & clean -----
dir='../../02.align'
samples <- c(paste0('Young_',1:3),
             paste0('Aged9_',1:3),
             paste0('Aged13_',1:3))

soupx_dbl <- function(sp){
  samples <- get(sp)
  
  sceList_clean <- lapply(samples,function(pro){
    print(pro)
    sce_raw <- Read10X(file.path(dir,pro,paste0(pro,'-1_S1_L001_outs'),'raw_cell_gene_matrix'))
    
    sce_filtered <- Read10X(file.path(dir,pro,paste0(pro,'-1_S1_L001_outs'),'filtered_cell_gene_matrix'))
    
    print(paste0(pro,' working on seurat'))
    sce_seu <- CreateSeuratObject(counts =  sce_filtered,
                                  project =  pro) %>%
      SCTransform(verbose = F) %>%
      RunPCA(verbose = F) %>%
      FindNeighbors(dims = 1:30, verbose = F) %>%
      FindClusters( resolution = 0.5, verbose = F)
    
    print(paste0(pro,' working on soupx'))
    sce_souxp <- SoupChannel(sce_raw,sce_filtered, calcSoupProfile = FALSE) %>%
      estimateSoup() %>%
      setClusters(clusters = sce_seu$seurat_clusters) %>%
      autoEstCont() %>%
      adjustCounts() %>%
      CreateSeuratObject() %>%
      SCTransform(verbose = F) %>%
      RunPCA(verbose = F) %>%
      FindNeighbors(dims = 1:30, verbose = F) %>%
      FindClusters( resolution = 0.5, verbose = F)
    
    print(paste0(pro,' working on dbf'))
    sce_soupx_dbl <- as.SingleCellExperiment(sce_souxp) %>%
      scDblFinder(clusters = 'seurat_clusters')
    sce_soupx_dbl$doublet_logic <- ifelse(sce_soupx_dbl$scDblFinder.class == "doublet", TRUE, FALSE)
    sce_soupx_dbl <- as.Seurat(sce_soupx_dbl, data = NULL) %>%
      subset(doublet_logic == 'FALSE')
    sce_soupx_dbl <- CreateSeuratObject( GetAssayData(sce_soupx_dbl,assay = 'RNA',slot = 'counts'),
                                         meta.data = sce_soupx_dbl@meta.data,
                                         project = pro)
    
    return(sce_soupx_dbl)
  })
}

for (tissue in 'samples') {
  sceList <- soupx_dbl(tissue)
  samples <- get(tissue)
  names(sceList) = samples
  for (samp in samples) {
    sceList[[samp]]$sample <- samp
  }
  assign('Clean_qc',sceList)
}
rm(sceList,tissue,samp)

# 02.delete low-quality cells -----

for (samples in c('Clean_qc')) {
  tmp <-  paste0(samples,'.merge')
  print(paste0('start merge ', samples))
  tmp.mer <-  merge(get(samples)[[1]], y = get(samples)[-1],add.cell.ids = names(get(samples)))
  tmp.mer$mitoRatio <- (PercentageFeatureSet(object = tmp.mer, pattern = '^mt-')) / 100
  tmp.mer$rpRatio <- (PercentageFeatureSet(object = tmp.mer, pattern = "^Rps|Rpl")) / 100
  
  tmp.mer$cells <- colnames(tmp.mer)
  tmp.mer$condition <- NA
  tmp.mer$condition[which(str_detect(tmp.mer$cells, "^Young"))] <- "Young"
  tmp.mer$condition[which(str_detect(tmp.mer$cells, "^Aged9"))] <- "Aged_9"
  tmp.mer$condition[which(str_detect(tmp.mer$cells, "^Aged13"))] <- "Aged_13"
  assign(tmp,tmp.mer)
  print(paste0('Merge samples & ratio finished'))
}
 rm(tmp,tmp.mer,sceList, tmp.mer)

## filter low-quality cells and low-expression genes -----
dim(Clean_qc.merge)
mito.mad <- mad(Clean_qc.merge$mitoRatio, constant=1)
mito.mad.upper <- median(Clean_qc.merge$mitoRatio) + 5*mito.mad
dim(Clean_qc.merge@meta.data[Clean_qc.merge$mitoRatio > mito.mad.upper,])
rp.mad <- mad(Clean_qc.merge$rpRatio, constant=1)
rp.mad.upper <- median(Clean_qc.merge$rpRatio) + 5*rp.mad
dim(Clean_qc.merge@meta.data[Clean_qc.merge$rpRatio > rp.mad.upper,])
Clean_qc.merge$log10GenesPerUMI <- log10(Clean_qc.merge$nFeature_RNA) / log10(Clean_qc.merge$nCount_RNA)
dim(Clean_qc.merge@meta.data[Clean_qc.merge$log10GenesPerUMI < 0.8,])

## filter cells
p1 <- Clean_qc.merge@meta.data %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

p2 <- Clean_qc.merge@meta.data %>% 
  ggplot(aes(color=sample, x=nFeature_RNA, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = c(350,1000))
p1 + p2

Clean_qc.merge@meta.data %>% 
  ggplot(aes(color=sample, x=nCount_RNA, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = c(650) )

Clean_qc.merge@meta.data %>% 
  ggplot(aes(x=sample, y=log10(nFeature_RNA), fill=sample)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")

p5 <- Clean_qc.merge@meta.data %>% 
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "lightgray", high = "red") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)

p6 <- Clean_qc.merge@meta.data %>% 
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=rpRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "lightgray", high = "red") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)
p5 + p6

p7 <- Clean_qc.merge@meta.data %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.09)

p8 <- Clean_qc.merge@meta.data %>% 
  ggplot(aes(x=rpRatio, color = sample, fill=sample)) +
  geom_density(alpha = 0.5) +
  theme_classic() +
  geom_vline(xintercept = rp.mad.upper)

p9 <- Clean_qc.merge@meta.data %>% 
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.5) +
  theme_classic() +
  geom_vline(xintercept = 0.8)

Clean_qc.merge.filtered <- subset(Clean_qc.merge, 
                                  cells = unique(c(rownames(Clean_qc.merge@meta.data[Clean_qc.merge$condition == 'Aged_13' & 
                                                                                       Clean_qc.merge$nCount_RNA > 350 & 
                                                                                       Clean_qc.merge$nFeature_RNA > 600 &
                                                                                       Clean_qc.merge$mitoRatio < 0.5,] ),
                                                   rownames(Clean_qc.merge@meta.data[Clean_qc.merge$condition == 'Aged_9' &
                                                                                       Clean_qc.merge$nCount_RNA > 350 &
                                                                                       Clean_qc.merge$nFeature_RNA > 600 &
                                                                                       Clean_qc.merge$mitoRatio < mito.mad.upper,] ),
                                                   rownames(Clean_qc.merge@meta.data[Clean_qc.merge$condition == 'Young' &
                                                                                       Clean_qc.merge$nFeature_RNA > 350 &
                                                                                       Clean_qc.merge$nCount_RNA > 600 &
                                                                                       Clean_qc.merge$mitoRatio < mito.mad.upper,]))
                                  ))

plan(multisession, workers=1)

cc_file <- read.csv('/data/02.project/00.other.lab/01.hualun/02.mus.age.placenta/03.analysis/cc_file',header = T)
g2m_genes <- select(org.Mm.eg.db, 
                    keys=as.character(cc_file[cc_file$phase == 'G2/M',]$geneID),
                    columns="SYMBOL", keytype="ENSEMBL")
s_genes <- select(org.Mm.eg.db, 
                  keys=as.character(cc_file[cc_file$phase == 'S',]$geneID),
                  columns="SYMBOL", keytype="ENSEMBL")
Clean_qc.merge.filtered <- JoinLayers(Clean_qc.merge.filtered) %>% NormalizeData() %>% 
  CellCycleScoring(s.features = s_genes$SYMBOL,
                   g2m.features = g2m_genes$SYMBOL)
dim(Clean_qc.merge.filtered)

## filter genes
Clean_qc.merge.filtered <- CreateSeuratObject(GetAssayData(object = Clean_qc.merge.filtered, layer = 'counts'), 
                                              meta.data = Clean_qc.merge.filtered@meta.data,
                                              min.cells = 20,
                                              min.features = 300)
metadata <- Clean_qc.merge.filtered@meta.data
Clean_qc.merge.filtered[['RNA']] <- split(Clean_qc.merge.filtered[['RNA']], f = Clean_qc.merge.filtered$sample)
Clean_qc.merge.filtered@meta.data <- metadata

rm(list = ls(pattern = 'p'))
save.image(compress = F)
