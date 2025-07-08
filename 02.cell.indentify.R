# intergrate data ------
plan(multisession, workers=30)
Clean_sct.inte <- Clean_qc.merge.filtered %>%
  SCTransform(vars.to.regress = c('mitoRatio','G2M.Score','S.Score')) %>% 
  RunPCA() %>% 
  IntegrateLayers(method = CCAIntegration,
                  k.anchor = 10,
                  normalization.method = "SCT") %>%
  FindNeighbors( reduction = "integrated.dr", 
                 dims = 1:30) %>% 
  FindClusters(resolution = .5) %>% 
  RunUMAP( reduction = "integrated.dr", 
           dims = 1:30)
Clean_sct.inte <- FindClusters(Clean_sct.inte, resolution = c(.2,.4,.6,1))
DimPlot(Clean_sct.inte, reduction = "umap", label = T,group.by = 'SCT_snn_res.1',split.by = 'condition')
VlnPlot(Clean_sct.inte, group.by = 'SCT_snn_res.1', pt.size = 0,features = c('nFeature_RNA','nCount_RNA')) # 10 14 24

Clean_sct.inte <- RunUMAP(Clean_sct.inte, reduction = "integrated.dr", 
                          dims = 1:30)
DimPlot(Clean_sct.inte, reduction = "umap", label = T,group.by = 'SCT_snn_res.0.2',split.by = 'condition')

# cell types identify----
Idents(Clean_sct.inte) <- Clean_sct.inte$SCT_snn_res.0.2
DefaultAssay(Clean_sct.inte) <- 'SCT'
Clean_sct.inte <- PrepSCTFindMarkers(Clean_sct.inte)
deg.res.0.2 <- FindAllMarkers(Clean_sct.inte, only.pos = T,logfc.threshold = 0.25,min.pct = .2) %>% dplyr::filter(p_val_adj < 0.05)
deg.res.0.2.top10 <- deg.res.0.2 %>% group_by(cluster) %>% top_n(n = 10,wt = avg_log2FC)

# Clean_sct.inte <- FindSubCluster(Clean_sct.inte, cluster = '10',graph.name = 'SCT_snn',resolution = .1)
# DimPlot(Clean_sct.inte, reduction = "umap", label = T,group.by = 'sub.cluster',split.by = 'condition')
# Idents(Clean_sct.inte) <- Clean_sct.inte$sub.cluster
# deg.res.0.2 <- FindAllMarkers(Clean_sct.inte, only.pos = T,logfc.threshold = 0.25,min.pct = .2) %>% dplyr::filter(p_val_adj < 0.05)
# deg.res.0.2.top10 <- deg.res.0.2 %>% group_by(cluster) %>% top_n(n = 10,wt = avg_log2FC)

## major cell type

# epi 9

FeaturePlot(Clean_sct.inte,split.by = 'condition', 
            features = c('Epcam','Krt19','Car3'),
            order = T)
# endo 1 3 8 11 LEC-12
FeaturePlot(Clean_sct.inte, 
            features = c('Pecam1','Kdr','Tek','Vwf','Cldn5','Reln'),min.cutoff = 'q1',
            order = T)
# pericyte na
FeaturePlot(Clean_sct.inte, 
            features = c('Rgs5','Notch3','Abcc9','Myh3','Mybpc1'),
            order = T)
# MSC 0 2 5 7 10 myo na
FeaturePlot(Clean_sct.inte, 
            features = c('Pdgfrb','Col3a1','Pdgfra','Acta2','Myh11'),
            order = T)
# immune-myeloid 1 3 8 11 --- 4 6
FeaturePlot(Clean_sct.inte, 
            features = c('Ptprc','Itgam','Lyz2','Cd3e'),
            order = T)
# Blood NA
FeaturePlot(Clean_sct.inte, 
            features = c('Hbb-bt','Hba-a1'),
            order = T)

# cycling 10
FeaturePlot(Clean_sct.inte,
            features = c('Mki67','Top2a','Ccnb1'), cols = c('gray90','red3'), split.by = 'condition',
            order = T)
# test special gene
FeaturePlot(Clean_sct.inte,
            features = c('Clec3b','Spp1','Msln','Igha'), cols = c('gray90','red3'), split.by = 'condition',
            order = T)

Clean_sct.inte <- RenameIdents(Clean_sct.inte,
                               '0' = 'Stromal cells',
                               '2' = 'Stromal cells',
                               '5' = 'Stromal cells',
                               '7' = 'Stromal cells',
                               '10' = 'Stromal cells',
                               '12' = 'Endothelial cells',
                               '9' = 'Epithelial cells',
                               '1' = 'Immune cells',
                               '3' = 'Immune cells',
                               '8' = 'Immune cells',
                               '11' = 'Immune cells',
                               '13' = 'Immune cells',
                               # '12' = 'Endothelial cells'
                               '4' = 'Immune cells',
                               '6' = 'Immune cells'
)
Clean_sct.inte$cell.type.major <- Idents(Clean_sct.inte)
DimPlot(Clean_sct.inte, group.by = 'cell.type.major',label = T)

## ICs subcluster ----

### integrate 
ICs.sub <- subset(Clean_sct.inte, cells = rownames(Clean_sct.inte@meta.data[Clean_sct.inte$cell.type.major == 'Immune cells',]))

ICs.inte <- ICs.sub %>%
  FindNeighbors( reduction = "integrated.dr", 
                 dims = 1:30) %>% 
  FindClusters(resolution = c(.2,.5,1)) %>% 
  RunUMAP( reduction = "integrated.dr", 
           dims = 1:30)

DimPlot(ICs.inte, reduction = "umap", label = T,group.by = 'SCT_snn_res.0.5')
DimPlot(ICs.inte, reduction = "umap", label = T,group.by = 'SCT_snn_res.1')

ICs.inte <- FindSubCluster(ICs.inte, cluster = c('15'),graph.name = 'SCT_snn',resolution = .1)
DimPlot(ICs.inte, reduction = "umap", label = T,group.by = 'sub.cluster')

### HVGs
Idents(ICs.inte) <- ICs.inte$SCT_snn_res.1
ICs.deg.res.1 <- FindAllMarkers(ICs.inte, only.pos = T,logfc.threshold = 0.25,min.pct = .1,recorrect_umi = F) %>% dplyr::filter(p_val_adj < 0.05)
ICs.deg.res.1.top10 <- ICs.deg.res.1 %>% group_by(cluster) %>% top_n(n = 10,wt = avg_log2FC)

Idents(ICs.inte) <- ICs.inte$SCT_snn_res.0.5
DefaultAssay(ICs.inte) <- 'SCT'
ICs.deg.res.0.5 <- FindAllMarkers(ICs.inte, only.pos = T,logfc.threshold = 0.25,min.pct = .1,recorrect_umi = F) %>% 
  dplyr::filter(p_val_adj < 0.05)
ICs.deg.res.0.5.top10 <- ICs.deg.res.0.5 %>% group_by(cluster) %>% top_n(n = 10,wt = avg_log2FC)

## major markers
FeaturePlot(ICs.inte, 
            features = c('Nkg7','Cd14','Cd8a','Cd3e'),
            order = T)

#NKT NA
FeaturePlot(ICs.inte,
            features = c('Gzmb','Klrb1c','Tbx21','Ncr1','Foxp3','Il2ra','Cd4','Cd8a','Cd3e'),
            order = T)
#NK 0 1 4 9 
FeaturePlot(ICs.inte, 
            features = c('Klrb1c','Itga2','Klrk1'),
            order = T)
# macro 7 11
FeaturePlot(ICs.inte, 
            features = c('Adgre1','C1qc','Apoe','C1qa'))

# mono 12 sub 15_0
FeaturePlot(ICs.inte, 
            features = c('Ly6c2','Plac8',"Cybb","Pou2f2"), order = T)

# B 13-Early B sub15_1-Naive B
FeaturePlot(ICs.inte, 
            features = c('Cd79a','Cd79b','Ms4a1','Mki67'),order = T)

# T __DP Memory T-5 cd4_17 cd8_14 γδ T_2 6 
FeaturePlot(ICs.inte,
            features = c('Ccl5','Cd3g','Cd8b1','Cd3e','Foxp3'),order = T)
FeaturePlot(ICs.inte,
            features = c('Il23r','Il17re','Il17a','Npnt'),order = T)

# neut 3 10
FeaturePlot(ICs.inte,
            features = c('S100a9','S100a8','Asprv1','G0s2','Lcn2','Rsad2'),order = T)

FeaturePlot(ICs.inte,
            features = c('Ccl3','Ccl4'),order = T)
# mast NA
FeaturePlot(ICs.inte,
            features = c('Ifitm1','Ccr3','Cd200r3','Cpa3'),order = T)
FeaturePlot(ICs.inte,
            features = c('Kit','Fcer1g','Fcer2','Kit'),order = F)

FeaturePlot(ICs.inte,
            features = c('Ptprc','Itgam'),order = T)

# GMP na
FeaturePlot(ICs.inte,
            features = c('Stfa1','Ngp','Mpo'),order = T)

# Immature dendritic cell 8
FeaturePlot(ICs.inte,
            features = c('Cd14','Cd83','Cd74','H2−Ab1'),order = T,cols = c('gray90','red3'))

VlnPlot(ICs.inte, features = c('nCount_RNA','nFeature_RNA'),pt.size = 0)
DimPlot(Clean_sct.inte, cells.highlight = rownames(ICs.inte@meta.data[ICs.inte$SCT_snn_res.0.5 == '16',]))

Idents(ICs.inte) <- ICs.inte$sub.cluster
ICs.inte <- RenameIdents(ICs.inte,
                         '3' = 'Neutrophils',
                         '10' = 'Neutrophils',
                         '7' = 'Macrophages',
                         '11' = 'Macrophages',
                         '12' = 'Monocytes',
                         '15_0' = 'Monocytes',
                         '8' = 'DCs',
                         '0' = 'NK',
                         '1' = 'NK',
                         '4' = 'NK',
                         '9' = 'NK',
                         '13' = 'Early B',
                         '15_1' = 'Naive B',
                         '5' = 'DPT',
                         '17' = 'Cd4+ T',
                         '14' = 'Cd8+ T',
                         '2' = 'γδT',
                         '6' = 'γδT',
                         '16' = 'Stromal cells')
ICs.inte$cell.type <- Idents(ICs.inte)
DimPlot(ICs.inte, group.by = 'cell.type',label = T,repel = T) + scale_color_simpsons()
Idents(ICs.inte) <- ICs.inte$cell.type

## redefine CD4
FeaturePlot(ICs.inte,features = c('Cd4','Cd8a'),order = T)
Idents(ICs.inte) <- ICs.inte$SCT_snn_res.0.5
ICs.inte <- FindSubCluster(ICs.inte, cluster = c('5'),graph.name = 'SCT_snn',resolution = .1)
DimPlot(ICs.inte,label = T,repel = T,group.by = 'sub.cluster')
cd4.cells <- c(rownames(ICs.inte@meta.data[ICs.inte$sub.cluster == '5_1',]))
cd8.cells <- c(rownames(ICs.inte@meta.data[ICs.inte$sub.cluster == '5_0',]),
               rownames(ICs.inte@meta.data[ICs.inte$sub.cluster == '5_2',]),
               rownames(ICs.inte@meta.data[ICs.inte$sub.cluster == '14',])
               )

Nuocytes.cells <- c(rownames(ICs.inte@meta.data[ICs.inte$sub.cluster == '17',]))

ICs.inte$cell.type <- as.character(ICs.inte$cell.type)
ICs.inte$cell.type[cd4.cells] <- 'Cd4+ T'
ICs.inte$cell.type[cd8.cells] <- 'Cd8+ T'
ICs.inte$cell.type[Nuocytes.cells] <- 'Nuocytes'

## SCs subcluster ----

### integrate 
SCs.sub <- subset(Clean_sct.inte, cells = c(rownames(Clean_sct.inte@meta.data[Clean_sct.inte$cell.type.major == 'Stromal cells',]), 
                                            rownames(ICs.inte@meta.data[ICs.inte$cell.type == 'Stromal cells',])))

SCs.inte <- SCs.sub %>%
  FindNeighbors( reduction = "integrated.dr", 
                 dims = 1:30) %>% 
  FindClusters(resolution = c(.2,.5,1)) %>% 
  RunUMAP( reduction = "integrated.dr", 
           dims = 1:30)
SCs.inte <- FindClusters(SCs.inte,resolution = .1)

DimPlot(SCs.inte, reduction = "umap", label = T,group.by = 'SCT_snn_res.0.1')
DimPlot(SCs.inte, reduction = "umap", label = T,group.by = 'SCT_snn_res.0.2')
DimPlot(SCs.inte, reduction = "umap", label = T,group.by = 'SCT_snn_res.0.5')
DimPlot(SCs.inte, reduction = "umap", label = T,group.by = 'SCT_snn_res.1')

SCs.inte <- FindSubCluster(SCs.inte, cluster = c('15'),graph.name = 'SCT_snn',resolution = .1)
DimPlot(SCs.inte, reduction = "umap", label = T,group.by = 'sub.cluster')

### HVGs
Idents(SCs.inte) <- SCs.inte$SCT_snn_res.1
SCs.deg.res.1 <- FindAllMarkers(SCs.inte, only.pos = T,logfc.threshold = 0.25,min.pct = .1,recorrect_umi = F) %>% dplyr::filter(p_val_adj < 0.05)
SCs.deg.res.1.top10 <- SCs.deg.res.1 %>% group_by(cluster) %>% top_n(n = 10,wt = avg_log2FC)

Idents(SCs.inte) <- SCs.inte$SCT_snn_res.0.5
DefaultAssay(SCs.inte) <- 'SCT'
SCs.deg.res.0.5 <- FindAllMarkers(SCs.inte, only.pos = T,logfc.threshold = 0.25,min.pct = .1,recorrect_umi = F) %>% 
  dplyr::filter(p_val_adj < 0.05)
SCs.deg.res.0.5.top10 <- SCs.deg.res.0.5 %>% group_by(cluster) %>% top_n(n = 10,wt = avg_log2FC)

Idents(SCs.inte) <- SCs.inte$SCT_snn_res.0.2
DefaultAssay(SCs.inte) <- 'SCT'
SCs.deg.res.0.2 <- FindAllMarkers(SCs.inte, only.pos = T,logfc.threshold = 0.25,min.pct = .1,recorrect_umi = F) %>% dplyr::filter(p_val_adj < 0.05)
SCs.deg.res.0.2.top10 <- SCs.deg.res.0.2 %>% group_by(cluster) %>% top_n(n = 10,wt = avg_log2FC)

Idents(SCs.inte) <- SCs.inte$SCT_snn_res.0.1
DefaultAssay(SCs.inte) <- 'SCT'
SCs.deg.res.0.1 <- FindAllMarkers(SCs.inte, only.pos = T,logfc.threshold = 0.15,min.pct = .1,recorrect_umi = F) %>% dplyr::filter(p_val_adj < 0.05)
SCs.deg.res.0.1.top10 <- SCs.deg.res.0.1 %>% group_by(cluster) %>% top_n(n = 10,wt = avg_log2FC)

## major markers
FeaturePlot(SCs.inte, 
            features = c('Pdgfra','Fn1'),
            order = T)

# col1a1 col13a1 col14a1 all
FeaturePlot(SCs.inte,
            features = c('Col13a1','Col14a1','Col1a1','Col3a1','Dcn'),min.cutoff = 'q1',
            order = T)
# Myofibro all
FeaturePlot(SCs.inte, 
            features = c('Tgfbi','Agt','Acta2'),
            order = T)
# SMC na
FeaturePlot(SCs.inte, 
            features = c('Acta2','Actc1','Des','Myh11'),order = T)

# pericyte na
FeaturePlot(SCs.inte, 
            features = c('Rgs5','Notch3','Abcc9'), order = T)
# ds es
FeaturePlot(SCs.inte,
            features = c('Igf1','Pcolce','Mmp11','Sfrp1'),order = T)

VlnPlot(SCs.inte, features = c('nCount_RNA','nFeature_RNA','mitoRatio'),pt.size = 0)
DimPlot(Clean_sct.inte, cells.highlight = rownames(SCs.inte@meta.data[SCs.inte$SCT_snn_res.0.1 == '3',]))
DimPlot(SCs.inte, cells.highlight = rownames(ICs.inte@meta.data[ICs.inte$cell.type == 'Stromal cells',]))

SCs.inte$res_0.2_sub <- SCs.sub$SCT_snn_res.0.2
DimPlot(SCs.inte, reduction = "umap", label = T,group.by = 'res_0.2_sub')

Idents(SCs.inte) <- SCs.inte$SCT_snn_res.0.1
SCs.inte <- RenameIdents(SCs.inte,
                         '0' = 'Nalf1+ Fibroblast',
                         '1' = 'Clec3b+ Fibroblast',
                         '2' = 'Ccl2+ Fibroblast',
                         '3' = 'Polluted',
                         '4' = 'Top2a+ Fibroblast')
SCs.inte$cell.type <- Idents(SCs.inte)
DimPlot(SCs.inte, group.by = 'cell.type',label = T,repel = T)

summary(SCs.inte$nFeature_RNA)
summary(SCs.inte$nCount_RNA)
summary(SCs.inte@meta.data[SCs.inte$SCT_snn_res.0.1 == '3',]$nFeature_RNA)
summary(SCs.inte@meta.data[SCs.inte$SCT_snn_res.0.1 == '3',]$nCount_RNA)

summary(Clean_sct.inte@meta.data[Clean_sct.inte$condition == 'Aged_13',]$nFeature_RNA)
summary(Clean_sct.inte@meta.data[Clean_sct.inte$condition == 'Aged_9',]$nFeature_RNA)
summary(Clean_sct.inte@meta.data[Clean_sct.inte$condition == 'Young',]$nFeature_RNA)

summary(SCs.inte@meta.data[SCs.inte$SCT_snn_res.0.1 == '0',]$nFeature_RNA)
summary(ICs.inte@meta.data[ICs.inte$cell.type == 'Neutrophils',]$nFeature_RNA)
summary(ICs.inte@meta.data[ICs.inte$SCT_snn_res.0.5 == '3',]$nFeature_RNA)
summary(ICs.inte@meta.data[ICs.inte$SCT_snn_res.0.5 == '10',]$nFeature_RNA)

## ECs subcluster ----

### integrate 
ECs.sub <- subset(Clean_sct.inte, cells = c(rownames(Clean_sct.inte@meta.data[Clean_sct.inte$cell.type.major == 'Endothelial cells',])))

ECs.inte <- ECs.sub %>%
  FindNeighbors( reduction = "integrated.dr", 
                 dims = 1:30) %>% 
  FindClusters(resolution = c(.2,.5,1)) %>% 
  RunUMAP( reduction = "integrated.dr", 
           dims = 1:30)
ECs.inte <- FindClusters(ECs.inte,resolution = .1)

DimPlot(ECs.inte, reduction = "umap", label = T,group.by = 'SCT_snn_res.0.1')
DimPlot(ECs.inte, reduction = "umap", label = T,group.by = 'SCT_snn_res.0.2')
DimPlot(ECs.inte, reduction = "umap", label = T,group.by = 'SCT_snn_res.0.5')
DimPlot(ECs.inte, reduction = "umap", label = T,group.by = 'SCT_snn_res.1')

### HVG
Idents(ECs.inte) <- ECs.inte$SCT_snn_res.0.1
DefaultAssay(ECs.inte) <- 'SCT'
ECs.deg.res.0.1 <- FindAllMarkers(ECs.inte, only.pos = T,logfc.threshold = 0.15,min.pct = .1,recorrect_umi = F) %>% dplyr::filter(p_val_adj < 0.05)
ECs.deg.res.0.1.top10 <- ECs.deg.res.0.1 %>% group_by(cluster) %>% top_n(n = 10,wt = avg_log2FC)

## major markers
FeaturePlot(ECs.inte, 
            features = c('Ackr1','Cxcl12','Gja5'),
            order = T)

## Epi subcluster ----

### integrate 
Epi.sub <- subset(Clean_sct.inte, cells = rownames(Clean_sct.inte@meta.data[Clean_sct.inte$cell.type.major == 'Epithelial cells',]))

Epi.inte <- Epi.sub %>%
  FindNeighbors( reduction = "integrated.dr", 
                 dims = 1:30) %>% 
  FindClusters(resolution = c(.2,.5,1,.1)) %>% 
  RunUMAP( reduction = "integrated.dr", 
           dims = 1:30)

DimPlot(Epi.inte, reduction = "umap", label = T,group.by = 'SCT_snn_res.0.2')
DimPlot(Epi.inte, reduction = "umap", label = T,group.by = 'SCT_snn_res.0.5')
DimPlot(Epi.inte, reduction = "umap", label = T,group.by = 'SCT_snn_res.0.1')
DimPlot(Epi.inte, reduction = "umap", label = T,group.by = 'SCT_snn_res.1')

### HVGs
Idents(Epi.inte) <- Epi.inte$SCT_snn_res.0.1
Epi.deg.res.0.1 <- FindAllMarkers(Epi.inte, only.pos = T,logfc.threshold = 0.25,min.pct = .1,recorrect_umi = F) %>% dplyr::filter(p_val_adj < 0.05)
Epi.deg.res.0.1.top10 <- Epi.deg.res.0.1 %>% group_by(cluster) %>% top_n(n = 10,wt = avg_log2FC)

Idents(Epi.inte) <- Epi.inte$SCT_snn_res.0.5
DefaultAssay(Epi.inte) <- 'SCT'
Epi.deg.res.0.5 <- FindAllMarkers(Epi.inte, only.pos = T,logfc.threshold = 0.25,min.pct = .1,recorrect_umi = F) %>% 
  dplyr::filter(p_val_adj < 0.05)
Epi.deg.res.0.5.top10 <- Epi.deg.res.0.5 %>% group_by(cluster) %>% top_n(n = 10,wt = avg_log2FC)

## major markers
FeaturePlot(Epi.inte, 
            features = c('Sox9','Mmp7','Lgr5','Ptgs1','Foxf1','Tppp3','Scgb3a1'),
            order = T)

Idents(Epi.inte) <- Epi.inte$SCT_snn_res.0.1
Epi.inte <- RenameIdents(Epi.inte,
                         '0' = 'Lumenal',
                         '1' = 'Sox9+',
                         '2' = 'Sox9+')
Epi.inte$cell.type <- Idents(Epi.inte)
DimPlot(Epi.inte, group.by = 'cell.type',label = T,repel = T)

DimPlot(Clean_sct.inte, cells.highlight = rownames(Epi.inte@meta.data[Epi.inte$SCT_snn_res.0.1 == '2',]))

## return cell types ----
Clean_sct.inte$cell.type.minor <- as.character( Clean_sct.inte$cell.type.major) 
Clean_sct.inte@meta.data[Clean_sct.inte$cell.type.minor == 'Immune cells',]$cell.type.minor <- as.character(ICs.inte$cell.type)
Clean_sct.inte@meta.data[Clean_sct.inte$cell.type.minor == 'Stromal cells',]$cell.type.minor <- as.character(SCs.inte$cell.type)
Clean_sct.inte@meta.data[Clean_sct.inte$cell.type.minor == 'Epithelial cells',]$cell.type.minor <- as.character(Epi.inte$cell.type)

DimPlot(Clean_sct.inte, group.by = 'cell.type.minor',label = T,cols = c(npg.cols, use.cols), repel = T) + NoAxes() + ggtitle("") 
Idents(Clean_sct.inte) <- Clean_sct.inte$cell.type.minor
table(Clean_sct.inte$cell.type.minor)

# remove polluted cells -----
Clean_sct.inte.rm.pt <- subset(Clean_sct.inte, cells = rownames(Clean_sct.inte@meta.data[Clean_sct.inte$cell.type.minor == 'Polluted',]), invert = T)
Clean_sct.inte.rm.pt <- RunUMAP(Clean_sct.inte.rm.pt,dims = 1:20)
DimPlot(Clean_sct.inte.rm.pt,group.by = 'cell.type.minor',label = T, cols = c(npg.cols, use.cols))
