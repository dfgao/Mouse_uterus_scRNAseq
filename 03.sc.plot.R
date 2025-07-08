
############################################################################## scRNA-seq section ########################################################################################
source('function.R')

# 01. umap plot ----
Clean_sct.inte.rm.pt$cell.type.percise.new <- Clean_sct.inte.rm.pt$cell.type.minor
Clean_sct.inte.rm.pt$cell.type.percise.new <- as.character(Clean_sct.inte.rm.pt$cell.type.percise.new)
Clean_sct.inte.rm.pt$cell.type.percise.new[cd4.cells] <- 'Cd4+ T'
Clean_sct.inte.rm.pt$cell.type.percise.new[cd8.cells] <- 'Cd8+ T'
Clean_sct.inte.rm.pt$cell.type.percise.new[Nuocytes.cells] <- 'Nuocytes'
Clean_sct.inte.rm.pt$condition <- factor(Clean_sct.inte.rm.pt$condition, levels = c('Young','Aged_9','Aged_13'))

p1 <- DimPlot(Clean_sct.inte.rm.pt, group.by = 'cell.type.percise.new', cols = c( use.cols,npg.cols), label = T) + ggtitle("") + NoAxes()
p2 <- DimPlot(Clean_sct.inte.rm.pt, group.by = 'cell.type.percise.new', cols = c(use.cols,npg.cols), split.by = 'condition',label = T,pt.size = .1) + ggtitle("") + NoAxes() + NoLegend()
p1 + p2

Clean_sct.inte.rm.pt$condition[which(str_detect(Clean_sct.inte.rm.pt$cells, "^Young"))] <- "Young"
Clean_sct.inte.rm.pt$condition[which(str_detect(Clean_sct.inte.rm.pt$cells, "^Aged9"))] <- "Aged_9"
Clean_sct.inte.rm.pt$condition[which(str_detect(Clean_sct.inte.rm.pt$cells, "^Aged13"))] <- "Aged_13"

Clean_sct.inte.rm.pt$cell.type.percise.new <- factor(Clean_sct.inte.rm.pt$cell.type.percise.new, levels = c('Nalf1+ Fibroblast',
                                                                                                              'Clec3b+ Fibroblast',
                                                                                                              'Ccl2+ Fibroblast',
                                                                                                              'Top2a+ Fibroblast',
                                                                                                              'Sox9+',
                                                                                                              'Lumenal',
                                                                                                              'Endothelial cells',
                                                                                                              'Monocytes',
                                                                                                            'Macrophages',
                                                                                                              'DCs',
                                                                                                              'Neutrophils',
                                                                                                              'Early B',
                                                                                                              'Naive B',
                                                                                                            'Cd4+ T',
                                                                                                            'Cd8+ T',
                                                                                                            'γδT',
                                                                                                            'Nuocytes',
                                                                                                            'NK'))
## major cell types
DimPlot(Clean_sct.inte.rm.pt, group.by = 'cell.type.major', cols = c( use.cols,npg.cols), label = T) + ggtitle("") + NoAxes()
Clean_sct.inte.rm.pt$cell.type.minor.fibro <- Clean_sct.inte.rm.pt$cell.type.percise.new 
Clean_sct.inte.rm.pt$cell.type.minor.fibro <-
  sub(".*Fibroblast.*", "Fibroblast", Clean_sct.inte.rm.pt$cell.type.percise.new)
Clean_sct.inte.rm.pt$cell.type.minor.fibro <- factor(Clean_sct.inte.rm.pt$cell.type.minor.fibro, levels = c('Fibroblast',
                                                                                                            'Sox9+',
                                                                                                            'Lumenal',
                                                                                                            'Endothelial cells',
                                                                                                            'Monocytes',
                                                                                                            'Macrophages',
                                                                                                            'DCs',
                                                                                                            'Neutrophils',
                                                                                                            'Early B',
                                                                                                            'Naive B',
                                                                                                            'Cd4+ T',
                                                                                                            'Cd8+ T',
                                                                                                            'γδT',
                                                                                                            'Nuocytes',
                                                                                                            'NK'))
## cirlize umap
library(plot1cell)

### cell.type.percise.new
Idents(Clean_sct.inte.rm.pt) <- Clean_sct.inte.rm.pt$cell.type.percise.new
circ_data <- prepare_circlize_data_test(Clean_sct.inte.rm.pt, scale = 0.7 )
circ_data$Cluster <- factor(circ_data$Cluster, levels = c(levels(Clean_sct.inte.rm.pt$cell.type.percise.new)))
set.seed(1234)
circ_data$x_polar2 <- circ_data$x_polar*100
cluster_colors<- c(use.cols,npg.cols)[1:18]
con_colors <- c('#0B996F','#fb8500','#f47068') 
circlize::circos.clear()
plot.cir.test(circ_data,do.label = T, pt.size = 0.1, col.use = cluster_colors ,bg.color = 'white', kde2d.n = 1000, repel = T, label.cex = 1)
add_track(circ_data, group = "condition", colors = con_colors, track_num = 2) 

### cell.type.minor.fibro
Idents(Clean_sct.inte.rm.pt) <- Clean_sct.inte.rm.pt$cell.type.minor.fibro
circ_data <- prepare_circlize_data_test(Clean_sct.inte.rm.pt, scale = 0.7 )
circ_data$Cluster <- factor(circ_data$Cluster, levels = c(levels(Clean_sct.inte.rm.pt$cell.type.minor.fibro)))
set.seed(1234)
circ_data$x_polar2 <- circ_data$x_polar*100
cluster_colors<- c(use.cols,npg.cols)[c(1,5:18)]
con_colors <- c('#0B996F','#fb8500','#f47068') 
circlize::circos.clear()
plot.cir.test(circ_data,do.label = T, pt.size = 0.1, col.use = cluster_colors ,bg.color = 'white', kde2d.n = 1000, repel = T, label.cex = 1)
add_track(circ_data, group = "condition", colors = con_colors, track_num = 2) 

DimPlot(Clean_sct.inte.rm.pt, group.by = 'cell.type.minor.fibro', cols = cluster_colors, split.by = 'condition',label = T,pt.size = .1) + ggtitle("") + NoAxes() + NoLegend()

# 03.ltsr composition of cell types ----

library(magrittr)
library(lme4)
library(numDeriv)

source('/data/02.project/00.other.lab/01.hualun/02.mus.age.placenta/03.analysis/LTSR.raw.code.R')

run_celltype_composition <- function(
    seu,
    sample_col        = "sample",
    celltype_col      = "cell.type.percise.new",
    treatment_levels  = c("Young","Aged9","Aged13"),
    reps              = c("one","two","three"),
    ltsr_vars         = c("Treatment","Rep"),
    FC                = 2,
    var_num1_name     = "Var_Num1"
) {

  levels(seu[[sample_col]]) <- paste0(
    rep(treatment_levels, each = length(reps)), "_", reps
  )

  cell_number <- FetchData(seu, vars = c(sample_col, celltype_col)) %>%
    as_tibble() %>%
    dplyr::count(!!sym(sample_col), !!sym(celltype_col)) %>%
    pivot_wider(
      names_from  = !!sym(sample_col),
      values_from = n,
      values_fill = 0
    )

  sample_ids         <- colnames(cell_number)[-1]
  cell_types         <- cell_number[[celltype_col]]
  n_cells_per_sample <- colSums(cell_number[,-1])

  sample_meta <- tibble(
    Sample_ID = sample_ids,
    Treatment = rep(treatment_levels, each = length(reps)),
    Rep       = rep(reps, times = length(treatment_levels)),
    cell.num  = n_cells_per_sample
  )

  obs_tbl <- data.frame(
    Sample_ID = rep(sample_ids, times = n_cells_per_sample),
    Treatment = rep(sample_meta$Treatment, times = n_cells_per_sample),
    Rep       = rep(sample_meta$Rep, times = n_cells_per_sample),
    Var_Num1  = rep(1, sum(n_cells_per_sample))
  )

  obs_tbl$Cell_type <- unlist(
    lapply(sample_ids, function(sid) {
      rep(cell_number[[celltype_col]], times = cell_number[[sid]])
    })
  )
  
  results <- CellTypeCompositionAnalysis(
    obs_tbl,
    colSample   = "Sample_ID",
    colCelltype = "Cell_type",
    colVarCats  = c(unlist(ltsr_vars)),
    colVarNums  = var_num1_name
  )
  ranef_tbl <- results$ranef
  sdse_tbl  <- results$sdse
  
  ranef_plot <- plot_ranef.new(
    ranef_tbl,
    vars           = list(Treatment = treatment_levels),
    celltypes      = cell_types,
    celltype_order = rev(cell_types),
    maxFC          = FC,
    LTSR2p         = FALSE
  ) +
    xlab("Condition") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  sdse_plot <- plot_sdse(
    sdse_tbl,
    colSample = "Sample_ID",
    ci                = 0.95,
    xlim              = c(0,1)
  )

  list(
    sample_meta = sample_meta,
    obs_tbl     = obs_tbl,
    ranef_tbl   = ranef_tbl,
    sdse_tbl    = sdse_tbl,
    ranef_plot  = ranef_plot,
    sdse_plot   = sdse_plot,
    combined    = ranef_plot + sdse_plot
  )
}

res <- run_celltype_composition(
  seu                   = Clean_sct.inte.rm.pt,
  sample_col            = "sample",
  celltype_col          = "cell.type.percise.new",
  treatment_levels      = c("Young","Aged9","Aged13"),
  reps                  = c("1","2","3"),
  ltsr_vars             = c("Treatment","Rep")
)

print(res$combined)


### add fibro
obs_tbl.fibro <- obs_tbl 
obs_tbl.fibro$Cell_type <- sub(".*Fibroblast.*", "Fibroblast", obs_tbl$Cell_type)
results <- CellTypeCompositionAnalysis(obs_tbl.fibro, "Sample_ID", "Cell_type", c("Treatment",'Rep'), "Var_Num1")
ranef_tbl <- results$ranef
sdse_tbl <- results$sdse
vars1 <- list(Treatment = c('Young','Aged9','Aged13'))

cell_types.fibro <- levels(Clean_sct.inte.rm.pt$cell.type.minor.fibro)
ranef_plot <- plot_ranef.new(ranef_tbl, vars = vars1, celltypes = cell_types.fibro, celltype_order = rev(cell_types.fibro),
                             maxFC = 2, LTSR2p = FALSE) + xlab('Condition') + theme(axis.text.x = element_text(angle = 45,hjust = 1))
sdse_plot <- plot_sdse(sdse_tbl, "Sample_ID", ci = 0.95, xlim = c(0, 1))
ranef_plot + sdse_plot

### cell ratio
library(ggalluvial)

plot.data.num <- FetchData(Clean_sct.inte.rm.pt, 
                           vars = c("condition", "cell.type.minor.fibro")) %>%
  dplyr::count(condition, cell.type.minor.fibro) %>% 
  tidyr::spread(condition, n) 

plot.data.num[is.na(plot.data.num)] <- 0

plot.data <- FetchData(Clean_sct.inte.rm.pt, 
                       vars = c("condition", "cell.type.minor.fibro")) %>%
  dplyr::count(condition, cell.type.minor.fibro, name = "n") %>% 
  dplyr::group_by(condition) %>%
  dplyr::mutate(Prop = n / sum(n)) %>% 
  dplyr::ungroup()

plot.data$levels <- c(rep(1:15,3))
plot.data$Prop <- plot.data$Prop*100
ggplot(plot.data,
       aes(x = condition, stratum = cell.type.minor.fibro, alluvium = levels,
           y = Prop,
           fill = cell.type.minor.fibro, label = cell.type.minor.fibro)) +
  scale_fill_manual(values = cluster_colors) +  
  geom_flow(stat = "alluvium", 
            lode.guidance = "frontback",
            width = 0.3,
            color = "darkgray") +
  geom_stratum(alpha = .8) +
  guides(fill = guide_legend(title = 'Cell type')) +
  ylab('Proportion (%)')+
  theme_bw() +
  theme(legend.position = "right",
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15))

# 04.Augur:prioritize the cell types most responsive to biological perturbations -----
library(Augur)
table(Clean_sct.inte.rm.pt$condition)
young_age9 <- subset(Clean_sct.inte.rm.pt, condition == 'Aged_13',invert = T)
dim(young_age9)
young_age9_augur <- calculate_auc(young_age9, cell_type_col = 'cell.type.percise.new', label_col = 'condition', n_threads = 80)
young_age9_augur.fibro <- calculate_auc(input = young_age9, cell_type_col = "cell.type.minor.fibro",     
  label_col  = "condition", n_subsamples = 15,  folds = 2,
  # min_cells  = 50, 
  var_quantile  = 0.3,
  feature_perc = 0.3, classifier = "rf",n_threads  = 110, show_progress= T)

age9_13 <- subset(Clean_sct.inte.rm.pt, condition == 'Young',invert = T)
dim(age9_13)
# age9_13_augur <- calculate_auc(age9_13, cell_type_col = 'cell.type.percise.new', label_col = 'condition', n_threads = 80)
age9_13_augur.fibro <- calculate_auc(input = age9_13_augur, cell_type_col = "cell.type.minor.fibro",     
                                        label_col  = "condition", n_subsamples = 15,  subsample_size= 50, folds = 2,
                                        min_cells  = 50, var_quantile  = 0.3,
                                        feature_perc = 0.3, classifier = "lr",n_threads  = 110, show_progress= T)


young_age13 <- subset(Clean_sct.inte.rm.pt, condition == 'Aged_9',invert = T)
dim(young_age13)
DimPlot(young_age13,group.by = 'cell.type.percise.new',label = T)
young_age13_augur <- calculate_auc(young_age13, cell_type_col = 'cell.type.percise.new', label_col = 'condition', n_threads = 90)
young_age13_augur.fibro <- calculate_auc(input = young_age13, cell_type_col = "cell.type.minor.fibro",     
                                           label_col  = "condition", n_subsamples = 15,  subsample_size= 50, folds = 2,
                                           min_cells  = 50, var_quantile  = 0.3,
                                           feature_perc = 0.3, classifier = "lr",n_threads  = 110, show_progress= T)
plot_loli <- function (augur) {
  aucs = augur$AUC
  size_sm = 10
  size_lg = 10
  range = range(aucs$auc)
  expand = abs(diff(range)) * 0.1
  p = aucs %>% ggplot(aes(x = reorder(cell_type, auc), y = auc)) + 
    geom_hline(aes(yintercept = 0.5), linetype = "dotted", 
               size = 0.8) +
    geom_point(size = 2) + 
    geom_text(aes(label = format(auc, digits = 3), 
                  y = ifelse(auc < 0.5, 0.5, auc)), 
              size = 4, nudge_y = expand, hjust = 0.5) + 
    geom_segment(aes(xend = cell_type, yend = 0.5)) + 
    scale_y_continuous("AUC", limits = c(min(range[1] - expand, 0.5), range[2] + expand * 1.5)) + 
    coord_flip() + 
    theme_bw() + 
    theme(axis.text.x = element_text(size = size_sm), 
          axis.text.y = element_text(size = size_sm + 2),
          axis.title.x = element_text(size = size_lg),
          axis.title.y = element_blank(), panel.grid = element_blank(), 
          strip.text = element_text(size = size_lg), strip.background = element_blank(),
          axis.line.y = element_blank(), axis.line.x = element_blank(), 
          legend.position = "top", legend.text = element_text(size = size_sm), 
          legend.title = element_text(size = size_sm), 
          legend.key.size = unit(0.6, "lines"),
          legend.margin = margin(rep(0, 4)),
          legend.background = element_blank(), 
          plot.title = element_text(size = size_lg, hjust = 0.5))
  p
}

plot_loli(young_age9_augur)
plot_loli(young_age13_augur)
plot_loli(age9_13_augur)

setdiff(age9_13_augur$AUC$cell_type,young_age13_augur$AUC$cell_type)
setdiff(age9_13_augur$AUC$cell_type,young_age9_augur$AUC$cell_type)

augur_auc <- young_age13_augur$AUC %>% left_join(y = young_age9_augur$AUC, by = 'cell_type') %>% left_join(y = age9_13_augur$AUC, by = 'cell_type')
augur_auc <- augur_auc[,c(1,3,4,2)]
colnames(augur_auc) <- c('Cell_type','MM/YM','OM/MM','OM/YM')

df_long <- augur_auc %>%
  pivot_longer(cols = -Cell_type,
               names_to  = "Comparison",
               values_to = "Value")

my_colors <- c(use.cols,npg.cols)[1:18]
names(my_colors) <- levels(Clean_sct.inte.rm.pt$cell.type.percise.new)
ggplot(df_long, aes(x = Comparison, y = Value,
                    group = Cell_type, color = Cell_type)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  scale_color_manual(values = my_colors) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  ) +
  labs(x = NULL, y = "Augur score")

### add fibro
augur_auc.fibro <- augur_auc
augur_auc.fibro$Cell_type <- as.character(augur_auc.fibro$Cell_type)
augur_auc.fibro$Cell_type[grep("Fibroblast", augur_auc.fibro$Cell_type, ignore.case = TRUE)] <- "Fibroblast"
augur_auc.fibro <- aggregate(
  . ~ Cell_type,
  data = augur_auc.fibro,
  FUN  = function(x) mean(x, na.rm = TRUE)
)
augur_auc.fibro$Cell_type <- factor(augur_auc.fibro$Cell_type, levels = levels(Clean_sct.inte.rm.pt$cell.type.minor.fibro))

df_long <- augur_auc.fibro %>%
  pivot_longer(cols = -Cell_type,
               names_to  = "Comparison",
               values_to = "Value")

my_colors <- cluster_colors
names(my_colors) <- levels(Clean_sct.inte.rm.pt$cell.type.minor.fibro)
ggplot(df_long, aes(x = Comparison, y = Value,
                    group = Cell_type, color = Cell_type)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  scale_color_manual(values = my_colors) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  ) +
  labs(x = NULL, y = "Augur score")

# 05.DEGs in cell types ----
library(ClusterGVis)
Clean_sct.inte.rm.pt <- PrepSCTFindMarkers(Clean_sct.inte.rm.pt)
rm.pt.ct.deg <- FindAllMarkers(Clean_sct.inte.rm.pt, only.pos = T,min.pct = 0.1,logfc.threshold = 0.25,recorrect_umi = F) %>% 
  dplyr::filter(p_val_adj < 0.05)
rm.pt.ct.deg.top10 <- rm.pt.ct.deg %>% dplyr::group_by(cluster) %>% top_n(n = 10,wt = avg_log2FC)
ht.data <- prepareDataFromscRNA(Clean_sct.inte.rm.pt,
                                diffData = rm.pt.ct.deg,
                                showAverage = T,
                                assays = 'SCT',slot = 'data',
                                group.by = 'cell.type.percise.new',keep.uniqGene = F,
                                scale.data = T)

enrich.go <- enrichCluster(object = ht.data,
                           OrgDb = org.Mm.eg.db,
                           type = "BP",
                           organism = "mmu",
                           pvalueCutoff = 0.01,
                           topn = 5,
                           seed = 1234)
enrich.go$ratio <- -log(enrich.go$pvalue)
rio::export(rm.pt.ct.deg, file = './uterus.sc/table1.scRNAseq cell types DEGs.xlsx')

pdf('./uterus.sc/cell.type.degs.heatmap.pdf',height = 12,width = 14,onefile = F)

visCluster(object = ht.data,
           ht.col.list = list(col_range = c(-4, 0, 4)),
           plot.type = "both",
           column_names_rot = 45,
           show_row_dend = F,
           # markGenes = markGenes,
           markGenes.side = "left",
           annoTerm.data = enrich.go,
           line.side = "left",
           cluster.order = c(1:18),
           go.col = rep(c(use.cols,npg.cols)[1:18],each = 5),
           sample.col = c(use.cols,npg.cols)[1:18],
           ctAnno.col = c(use.cols,npg.cols)[1:18],
           go.size = 8,
           add.bar = T)
dev.off()

rm.pt.ct.deg_pos_neg <- FindAllMarkers(Clean_sct.inte.rm.pt, only.pos = F,min.pct = 0.1,logfc.threshold = 0.25) %>% 
  dplyr::filter(p_val_adj < 0.05)
rio::export(rm.pt.ct.deg_pos_neg, file = './uterus.sc/rm.pt.ct.deg_pos_neg.xlsx')


### add fibro
Idents(Clean_sct.inte.rm.pt) <- Clean_sct.inte.rm.pt$cell.type.minor.fibro
rm.pt.ct.deg.add.fibro <- FindAllMarkers(Clean_sct.inte.rm.pt, only.pos = T,min.pct = 0.1,logfc.threshold = 0.25,recorrect_umi = F) %>% dplyr::filter(p_val_adj < 0.05)
rm.pt.ct.deg.add.fibro.top10 <- rm.pt.ct.deg.add.fibro %>% dplyr::group_by(cluster) %>% top_n(n = 10,wt = avg_log2FC)
ht.data <- prepareDataFromscRNA(Clean_sct.inte.rm.pt,
                                diffData = rm.pt.ct.deg.add.fibro,
                                showAverage = T,
                                assays = 'SCT',slot = 'data',
                                group.by = 'cell.type.minor.fibro',keep.uniqGene = F,
                                scale.data = T)

enrich.go <- enrichCluster(object = ht.data,
                           OrgDb = org.Mm.eg.db,
                           type = "BP",
                           organism = "mmu",
                           pvalueCutoff = 0.01,
                           topn = 5,
                           seed = 1234)
enrich.go$ratio <- -log(enrich.go$pvalue)
rio::export(rm.pt.ct.deg, file = './2025.7.add/table1.scRNAseq cell types DEGs.xlsx')

pdf('./2025.7.add/cell.type.degs.heatmap.pdf',height = 12,width = 14,onefile = F)

visCluster(object = ht.data,
           ht.col.list = list(col_range = c(-4, 0, 4)),
           plot.type = "both",
           column_names_rot = 45,
           show_row_dend = F,
           # markGenes = markGenes,
           markGenes.side = "left",
           annoTerm.data = enrich.go,
           line.side = "left",
           cluster.order = c(1:15),
           go.col = rep(cluster_colors,each = 5),
           sample.col = cluster_colors,
           ctAnno.col = cluster_colors,
           go.size = 8,
           add.bar = T)
dev.off()

# 06.dotplot of key genes -----
library(plot1cell)

Idents(Clean_sct.inte.rm.pt) <- Clean_sct.inte.rm.pt$cell.type.percise.new
com.dot.new(Clean_sct.inte.rm.pt, feature = c('Nalf1','Adarb2','Aox3',
                                               'Clec3b','Ces2g','Aspn',
                                               'Ccl2','Ccl7','Cxcl10',
                                               'Top2a','Mki67','Anln',
                                               'Sox9','Fxyd4','Sprr2f',
                                               'Upk3b','Lrrn4','Msln',
                                               'Cldn5','Mmrn1','Reln',
                                              'Msr1','Adgre4',"Cybb",
                                               'C1qa','C1qc','C1qb',
                                               'Cd74','Ccr7','Cd209a',
                                               'S100a9','S100a8','Il36g',
                                               'Jchain','Igha','Mzb1',
                                               'Bank1','Cd19','Ms4a1',
                                              'Cd4','Trat1','Cd40lg',
                                              'Cd8a','Cd6','Cd8b1',
                                              'Scart1','Trgc1','Il23r',
                                              'Il1rl1','Hs3st1','Il5',
                                              'Klra4','Klrb1c','Klra8')
            ,groups = "condition",strip.color = c(use.cols,npg.cols)[1:18])

nad.enzyme <- c('Nampt','Nmnat1','Nmnat2','Nmnat3','Ido1','Cd38','Bst1','Sirt1','Sirt3','Sirt5','Parp2')
com.dot.new(Clean_sct.inte.rm.pt, feature = nad.enzyme,
            groups = "condition",strip.color = c(use.cols,npg.cols)[1:18])

### add fibro
Idents(Clean_sct.inte.rm.pt) <- Clean_sct.inte.rm.pt$cell.type.minor.fibro
com.dot.new(Clean_sct.inte.rm.pt, feature = c('Col1a1','Col3a1','Col6a4',
                                              'Sox9','Fxyd4','Epcam',
                                              'Upk3b','Lrrn4','Msln',
                                              'Cldn5','Mmrn1','Reln',
                                              'Msr1','Adgre4',"Cybb",
                                              'C1qa','C1qc','C1qb',
                                              'Cd74','Ccr7','Cd209a',
                                              'S100a9','S100a8','Il36g',
                                              'Jchain','Igha','Mzb1',
                                              'Bank1','Cd19','Ms4a1',
                                              'Cd4','Trat1','Cd40lg',
                                              'Cd8a','Cd6','Cd8b1',
                                              'Scart1','Trgc1','Il23r',
                                              'Il1rl1','Hs3st1','Il5',
                                              'Klra4','Klrb1c','Klra8')
            ,groups = "condition",strip.color = cluster_colors)

com.dot.new(Clean_sct.inte.rm.pt, feature = nad.enzyme,
            groups = "condition",strip.color = cluster_colors)

# 07.trans noise ----

## cell levels
library(MASS)
library(ggpubr)
celltypes <- unique(Clean_sct.inte.rm.pt@meta.data$cell.type.percise.new)
celltypes <- celltypes[which(!is.na(celltypes))]
Idents(Clean_sct.inte.rm.pt) <- Clean_sct.inte.rm.pt$cell.type.percise.new
set.seed(12345)

res <- lapply(celltypes, function(x) getEuclideanDistance(x, 
                                                          obj = Clean_sct.inte.rm.pt,
                                                          assay = 'SCT',
                                                          slot = 'counts',
                                                          group.by = 'condition',
                                                          ident1 = 'Young',
                                                          ident2 = 'Aged_9',
                                                          ident3 = 'Aged_13',
                                                          lowcv = T))
names(res) <- celltypes
res.df <- data.frame(TN.value = unlist(do.call(c, res))) %>% 
  rownames_to_column(var = 'info') %>% 
  separate(col = info,
           into = c('Celltype','Condition','Sample_cells'),
           sep = '\\.',
           remove = T,
           extra = "merge")
res.df$Condition <- factor(ifelse(res.df$Condition == 'ident1','Young',
                                  ifelse(res.df$Condition == 'ident2','Aged_9','Aged_13')), levels = c('Young','Aged_9','Aged_13'))
res.df$Celltype <- factor(res.df$Celltype, levels = levels(Clean_sct.inte.rm.pt$cell.type.percise.new))

my_comparisons <- list(c('Young','Aged_9'),c('Aged_9','Aged_13'))
ggplot(res.df, aes(x=Condition, y=TN.value, fill=Condition)) + 
  geom_boxplot(alpha = .7) + 
  theme_bw() +
  labs(x = '', y = 'Transcriptional \nheterogeneity') +
  theme(legend.position = "none") +
  scale_fill_manual(values = con_colors) +
  facet_wrap(~Celltype, scale="free",nrow = 3) +
  ggpubr::stat_compare_means(comparisons = my_comparisons,paired = F,
                             method = "wilcox.test")

### add fibro
celltypes <- unique(Clean_sct.inte.rm.pt@meta.data$cell.type.minor.fibro)
celltypes <- celltypes[which(!is.na(celltypes))]
Idents(Clean_sct.inte.rm.pt) <- Clean_sct.inte.rm.pt$cell.type.minor.fibro
set.seed(12345)

res.fibro <- lapply(celltypes, function(x) getEuclideanDistance(x, 
                                                          obj = Clean_sct.inte.rm.pt,
                                                          assay = 'SCT',
                                                          slot = 'counts',
                                                          group.by = 'condition',
                                                          ident1 = 'Young',
                                                          ident2 = 'Aged_9',
                                                          ident3 = 'Aged_13',
                                                          lowcv = T))
names(res.fibro) <- celltypes
res.fibro.df <- data.frame(TN.value = unlist(do.call(c, res.fibro))) %>% 
  rownames_to_column(var = 'info') %>% 
  separate(col = info,
           into = c('Celltype','Condition','Sample_cells'),
           sep = '\\.',
           remove = T,
           extra = "merge")
res.fibro.df$Condition <- factor(ifelse(res.fibro.df$Condition == 'ident1','Young',
                                  ifelse(res.fibro.df$Condition == 'ident2','Aged_9','Aged_13')), levels = c('Young','Aged_9','Aged_13'))
res.fibro.df$Celltype <- factor(res.fibro.df$Celltype, levels = levels(Clean_sct.inte.rm.pt$cell.type.minor.fibro))

my_comparisons <- list(c('Young','Aged_9'),c('Aged_9','Aged_13'))
ggplot(res.fibro.df, aes(x=Condition, y=TN.value, fill=Condition)) + 
  geom_boxplot(alpha = .7,width = .5) + 
  theme_bw() +
  labs(x = '', y = 'Transcriptional heterogeneity') +
  theme(legend.position = "none",axis.text.x = element_text(angle = 45,hjust = 1)) +
  scale_fill_manual(values = con_colors) +
  facet_wrap(~Celltype, scale="free",nrow = 3) +
  ggpubr::stat_compare_means(comparisons = my_comparisons,paired = F,
                             method = "wilcox.test")

## sample levels
sum_df <- res.fibro.df %>%
  dplyr::group_by(Celltype, Condition) %>%
  dplyr::summarise(total_TN = sum(TN.value, na.rm = TRUE), .groups = "drop")

wide_df <- sum_df %>%
  pivot_wider(
    names_from  = Condition,
    values_from = total_TN
  )
ratio_df <- wide_df %>%
  mutate(
    log2_Aged9_vs_Young      = log2(Aged_9    / Young),
    log2_Aged13_vs_Aged9     = log2(Aged_13   / Aged_9),
    log2_Aged13_vs_Young     = log2(Aged_13   / Young)
  )
plot_df <- ratio_df %>%
  select(Celltype, starts_with("log2_")) %>%
  pivot_longer(
    cols      = -Celltype,
    names_to  = "comparison",
    values_to = "log2ratio"
  )

plot_df$Celltype <- factor(plot_df$Celltype, levels = rev(levels(plot_df$Celltype)))
ggplot(plot_df, aes(x = log2ratio, y = Celltype, color = comparison)) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0) + 
  scale_color_brewer(palette = "Set1") +
  labs(
    x = "log2 ratio",
    y = "Cell type",
    color = "Comparison"
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8),
    legend.position = "right"
  )

p1 <- plot_df[plot_df$comparison == 'log2_Aged9_vs_Young',] %>% 
ggplot( aes(x = log2ratio, y = Celltype, color = Celltype)) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0) + 
  labs(
    x = "log2 ratio",
    y = "Cell type",
    color = "cluster_colors"
  ) +
  scale_color_manual(values = rev(cluster_colors)) + 
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8),
    legend.position = ""
  ) + ggtitle('Aged9_vs_Young')

p2 <- plot_df[plot_df$comparison == 'log2_Aged13_vs_Aged9',] %>% 
  ggplot( aes(x = log2ratio, y = Celltype, color = Celltype)) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0) + 
  labs(
    x = "log2 ratio",
    y = "Cell type",
    color = "Celltype"
  ) +
  scale_color_manual(values = rev(cluster_colors)) + 
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8),
    legend.position = ""
  ) + ggtitle('Aged13_vs_Aged9')

p3 <- plot_df[plot_df$comparison == 'log2_Aged13_vs_Young',] %>% 
  ggplot( aes(x = log2ratio, y = Celltype, color = Celltype)) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0) + 
  labs(
    x = "log2 ratio",
    y = "Cell type",
    color = "Celltype"
  ) +
  scale_color_manual(values = rev(cluster_colors)) + 
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8),
    legend.position = ""
  ) + ggtitle('Aged13_vs_Young')
p1 + p2 + p3
# 08.umap of key genes ----
p1 <- FeaturePlot(Clean_sct.inte.rm.pt, features = c('Cd38'),split.by = 'condition',order = T ,combine = F, cols = c('gray90','red3'))
for(i in 1:length(p1)) {
  p1[[i]] <- p1[[i]] + NoLegend() + NoAxes() + theme(panel.background=element_rect(fill='transparent', color='black'), title = element_text(size = 8))
}
patchwork::wrap_plots(c(p1),nrow = 1)

# 08.1 fibroblast subtypes ----

### umap
fibro.seu <- subset(Clean_sct.inte.rm.pt, subset = (cell.type.minor.fibro == 'Fibroblast'))
DimPlot(fibro.seu, group.by = 'cell.type.percise.new',split.by = 'condition')
fibro.seu <- fibro.seu %>% RunUMAP( reduction = "integrated.dr", dims = 1:30)
DimPlot(fibro.seu, group.by = 'cell.type.percise.new', label = T) + ggtitle("") + NoAxes() + scale_color_d3()
DimPlot(fibro.seu, group.by = 'cell.type.percise.new', split.by = 'condition',label = T,pt.size = .1) + ggtitle("") + NoAxes() + NoLegend() + scale_color_d3()

FeaturePlot(fibro.seu, features = c('Top2a','Mki67'), order = T,cols = c('gray90','red3'),split.by = 'condition')
### dotplot
Idents(fibro.seu) <- fibro.seu$cell.type.percise.new
DefaultAssay(fibro.seu) <- 'SCT'
fibro.deg<- FindAllMarkers(fibro.seu, only.pos = T,logfc.threshold = 0.25,min.pct = .1,recorrect_umi = F) %>% dplyr::filter(p_val_adj < 0.05)
fibro.deg.top10 <- fibro.deg %>% group_by(cluster) %>% top_n(n = 10,wt = avg_log2FC)
com.dot.new(fibro.seu, feature = c('Nalf1','Adarb2','Aox3',
                                              'Clec3b','Ces2g','Aspn',
                                              'Ccl2','Ccl7','Cxcl10',
                                              'Top2a','Mki67','Anln')
            ,groups = "condition",strip.color = d3.cols[1:4])
rio::export(fibro.deg, file = '2025.7.add/fibro.cell.type.degs.xlsx')

### heatmap
ht.data <- prepareDataFromscRNA(fibro.seu,
                                diffData = fibro.deg,
                                showAverage = T,
                                assays = 'SCT',slot = 'data',
                                group.by = 'cell.type.percise.new',keep.uniqGene = F,
                                scale.data = T)
enrich.go <- enrichCluster(object = ht.data,
                           OrgDb = org.Mm.eg.db,
                           type = "BP",
                           organism = "mmu",
                           pvalueCutoff = 0.01,
                           topn = 10,
                           seed = 1234)
enrich.go$ratio <- -log(enrich.go$pvalue)

pdf('./2025.7.add/fibro.cell.type.degs.heatmap.pdf',height = 12,width = 14,onefile = F)
visCluster(object = ht.data,
           ht.col.list = list(col_range = c(-4, 0, 4)),
           plot.type = "both",
           column_names_rot = 45,
           show_row_dend = F,
           # markGenes = markGenes,
           markGenes.side = "left",
           annoTerm.data = enrich.go,
           line.side = "left",
           cluster.order = c(1:4),
           go.col = rep(d3.cols[1:4],each = 10),
           sample.col = d3.cols[1:4],
           ctAnno.col = d3.cols[1:4],
           go.size = 8,bar.width = 5,
           add.bar = T)
dev.off()

enrich.go <- enrichCluster(object = ht.data,
                           OrgDb = org.Mm.eg.db,
                           type = "BP",
                           organism = "mmu",
                           pvalueCutoff = 0.01,
                           topn = 50,
                           seed = 1234)
enrich.go[enrich.go$group == 'C1',]$group <- 'Nalf1'
enrich.go[enrich.go$group == 'C2',]$group <- 'Celc3b'
enrich.go[enrich.go$group == 'C3',]$group <- 'Ccl2'
enrich.go[enrich.go$group == 'C4',]$group <- 'Top2a'
rio::export(enrich.go, file = '2025.7.add/fibro.cell.type.degs.GOBP.top50.xlsx')

### LTSR
res <- run_celltype_composition(
  seu                   = fibro.seu,
  sample_col            = "sample",
  celltype_col          = "cell.type.percise.new",
  treatment_levels      = c("Young","Aged9","Aged13"),
  reps                  = c("1","2","3"),
  ltsr_vars             = c("Treatment","Rep"),
  FC = 1.2
)
print(res$combined)

### ca_degs using findmarkers
comlist <- t(combn(c("Young","Aged_9","Aged_13"), 3))
deg_corss_condition_age_up <- list()
Idents(fibro.seu) <- 'cell.type.percise.new'

#### ca_DEG-up
for (ct in unique(fibro.seu$cell.type.percise.new)) {
  cat("==> Processing cell type:", ct, "\n")
  
  sub.seu.1 <- subset(fibro.seu, idents = ct)
  Idents(sub.seu.1) <- 'condition'
  for (gp in nrow(comlist)) {
    Young <- comlist[gp,1]
    Aged_9 <- comlist[gp,2]
    Aged_13 <- comlist[gp,3]
    DEGs_young_age9 <- FindMarkers(sub.seu.1,
                                   ident.1 = Aged_9,
                                   ident.2 = Young,
                                   logfc.threshold = 0.1,
                                   group.by = 'condition',
                                   recorrect_umi = F,only.pos = T) 
    
    DEGs_age9_age13 <- FindMarkers(sub.seu.1,
                                   ident.1 = Aged_13,
                                   ident.2 = Aged_9,
                                   logfc.threshold = 0.2,
                                   group.by = 'condition',
                                   recorrect_umi = F,only.pos = T) %>% dplyr::filter(p_val_adj < 0.05)
    
    ca_DEGs_age_up <- intersect(rownames(DEGs_young_age9), rownames(DEGs_age9_age13))
    deg_corss_condition_age_up[[ct]] <- ca_DEGs_age_up
  }
}
DEGs_young_age9['Nme1',]
DEGs_age9_age13['Nme1',]
FeaturePlot(fibro.seu, features = 'Nme1',split.by = 'condition',order = T,min.cutoff = 'q50')

library(purrr)
max_len <- max(map_int(deg_corss_condition_age_up, length))
deg_corss_condition_age_up.df <- map_dfr(deg_corss_condition_age_up, ~ { length(.x) <- max_len; .x })
deg_corss_condition_age_up.df <- list2DF(lapply(deg_corss_condition_age_up, `length<-`, max(lengths(deg_corss_condition_age_up))))
colnames(deg_corss_condition_age_up.df)
rio::export(deg_corss_condition_age_up.df, file = '2025.7.add/fibro.ca_DGEs_age_up_by_findmarkers.xlsx')

## down
deg_corss_condition_age_down <- list()
for (ct in unique(fibro.seu$cell.type.percise.new)) {
  cat("==> Processing cell type:", ct, "\n")
  
  sub.seu.1 <- subset(fibro.seu, idents = ct)
  Idents(sub.seu.1) <- 'condition'
  for (gp in nrow(comlist)) {
    Young <- comlist[gp,1]
    Aged_9 <- comlist[gp,2]
    Aged_13 <- comlist[gp,3]
    DEGs_young_age9 <- FindMarkers(sub.seu.1,
                                   ident.1 = Young,
                                   ident.2 = Aged_9,
                                   logfc.threshold = 0.1,
                                   group.by = 'condition',
                                   recorrect_umi = F,only.pos = T) 
    
    DEGs_age9_age13 <- FindMarkers(sub.seu.1,
                                   ident.1 = Aged_9,
                                   ident.2 = Aged_13,
                                   logfc.threshold = 0.2,
                                   group.by = 'condition',
                                   recorrect_umi = F,only.pos = T) %>% dplyr::filter(p_val_adj < 0.05)
    
    ca_DEGs_age_down <- intersect(rownames(DEGs_young_age9), rownames(DEGs_age9_age13))
    deg_corss_condition_age_down[[ct]] <- ca_DEGs_age_down
  }
}

max_len <- max(map_int(deg_corss_condition_age_down, length))
deg_corss_condition_age_down.df <- map_dfr(deg_corss_condition_age_down, ~ { length(.x) <- max_len; .x })
deg_corss_condition_age_down.df <- list2DF(lapply(deg_corss_condition_age_down, `length<-`, max(lengths(deg_corss_condition_age_down))))
colnames(deg_corss_condition_age_down.df)
rio::export(deg_corss_condition_age_down.df, file = '2025.7.add/fibro.ca_DGEs_age_down_by_findmarkers.xlsx')

### ca-DEGs from pseudobulk
library(stringr)

agg <- AggregateExpression(
  object     = fibro.seu,
  assay      = "RNA",
  slot       = "counts", 
  group.by   = c("sample", "cell.type.percise.new", "condition")
)
pb_counts <- agg$RNA %>% as.data.frame() %>% round(0)

coldata <- data.frame(
  sample_celltype = colnames(pb_counts),
  stringsAsFactors = FALSE
) %>%
  tidyr::separate(
    col = sample_celltype,
    into = c("sample", "cell.type.percise.new", "condition"),
    sep  = "_",
    remove = FALSE
  )
pb_list <- split(seq_len(ncol(pb_counts)), coldata$cell.type)
pseudo_bulk_by_ct <- lapply(pb_list, function(idxs){
  list(
    counts  = pb_counts[, idxs, drop = FALSE],
    coldata = coldata[idxs, ]
  )
})

library(DESeq2)
results_cross <- map(pseudo_bulk_by_ct, function(x){
  mat <- x$counts
  cd  <- x$coldata
  
  keep <- cd$condition %in% c("Young", "Aged-9", "Aged-13")
  mat  <- mat[, keep]
  cd   <- cd[keep, ]
  
  dds <- DESeqDataSetFromMatrix(
    countData = mat,
    colData   = cd %>% mutate(condition = factor(condition, 
                                                 levels = c("Young","Aged-9","Aged-13"))),
    design    = ~ condition
  )
  
  keep_gene <- rowSums(counts(dds) >= 10) >= 2
  dds <- dds[keep_gene, ]
  
  dds <- DESeq(dds, quiet = TRUE)
  # Aged_9 vs Young
  res1 <- results(dds,
                  contrast = c("condition","Aged-9","Young"),
                  alpha    = 0.05) %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    as_tibble() %>% 
    na.omit()
  
  # Aged_13 vs Aged_9
  res2 <- results(dds,
                  contrast = c("condition","Aged-13","Aged-9"),
                  alpha    = 0.05) %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    as_tibble() %>% 
    na.omit()
  
  up1   <- dplyr::filter(res1, log2FoldChange > 0.5) %>% pull(gene)
  down1 <- dplyr::filter(res1, log2FoldChange < -0.5) %>% pull(gene)
  up2   <- dplyr::filter(res2, pvalue < 0.05, log2FoldChange > 0.5) %>% pull(gene)
  down2 <- dplyr::filter(res2, pvalue < 0.05, log2FoldChange < -0.5) %>% pull(gene)
  
  list(
    up_intersect   = intersect(up1,   up2),
    down_intersect = intersect(down1, down2)
  )
})

summary_df <- imap_dfr(results_cross, ~ tibble(
  cell.type      = .y,
  direction      = c("up","down"),
  n.genes        = c(length(.x$up_intersect),
                     length(.x$down_intersect)),
  genes          = c(paste(.x$up_intersect, collapse = ","),
                     paste(.x$down_intersect, collapse = ","))
))
rio::export(summary_df, file = '2025.7.add//fibro.ca_DGEs_by_pseudobulk.xlsx')
FeaturePlot(fibro.seu,features = 'Blnk',order = T,split.by = 'condition')

### ca_DEG young or aged_9 vs aged13
results_cross <- map(pseudo_bulk_by_ct, function(x){
  mat <- x$counts
  cd  <- x$coldata
  
  keep <- cd$condition %in% c("Young", "Aged-9", "Aged-13")
  mat  <- mat[, keep]
  cd   <- cd[keep, ]
  
  dds <- DESeqDataSetFromMatrix(
    countData = mat,
    colData   = cd %>% mutate(condition = factor(condition, 
                                                 levels = c("Young","Aged-9","Aged-13"))),
    design    = ~ condition
  )
  
  keep_gene <- rowSums(counts(dds) >= 10) >= 2
  dds <- dds[keep_gene, ]
  
  dds <- DESeq(dds, quiet = TRUE)
  # Aged_9 vs Young
  res1 <- results(dds,
                  contrast = c("condition","Aged-9","Young"),
                  alpha    = 0.05) %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    as_tibble() %>% 
    na.omit()
  
  # Aged_13 vs Aged_9
  res2 <- results(dds,
                  contrast = c("condition","Aged-13","Young"),
                  alpha    = 0.05) %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    as_tibble() %>% 
    na.omit()
  
  up1   <- dplyr::filter(res1, pvalue < 0.05, log2FoldChange > 0.5) %>% pull(gene)
  down1 <- dplyr::filter(res1, pvalue < 0.05, log2FoldChange < -0.5) %>% pull(gene)
  up2   <- dplyr::filter(res2, pvalue < 0.05, log2FoldChange > 0.5) %>% pull(gene)
  down2 <- dplyr::filter(res2, pvalue < 0.05, log2FoldChange < -0.5) %>% pull(gene)
  
  list(
    up_intersect   = intersect(up1,   up2),
    down_intersect = intersect(down1, down2)
  )
})

summary_df <- imap_dfr(results_cross, ~ tibble(
  cell.type      = .y,
  direction      = c("up","down"),
  n.genes        = c(length(.x$up_intersect),
                     length(.x$down_intersect)),
  genes          = c(paste(.x$up_intersect, collapse = ","),
                     paste(.x$down_intersect, collapse = ","))
))
rio::export(summary_df, file = '2025.7.add//fibro.ca_DGEs_by_pseudobulk_age13-young & age13-age9 jiaoji.xlsx')
FeaturePlot(fibro.seu,features = 'Blnk',order = T,split.by = 'condition')


### cytotrace
library(CytoTRACE2)
fibro.top2a.seu <- subset(Clean_sct.inte.rm.pt, subset = (cell.type.percise.new == 'Top2a+ Fibroblast'))
dim(fibro.top2a.seu)

fibro.top2a.cyto <- cytotrace2(fibro.top2a.seu,is_seurat = T,ncores = 50,species = 'mouse')
ggboxplot(fibro.top2a.cyto@meta.data, x="condition", y="CytoTRACE2_Score", width = 0.6, 
          color = "black",
          fill="condition",
          xlab = F,
          bxp.errorbar=T,
          bxp.errorbar.width=0.5, 
          size=.1,
          outlier.shape=NA,
          legend = "right",
          alpha = 0.8) + 
  ylab('Potency score')  + ggtitle('Top2a+ Fibroblast') +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)) + 
  scale_fill_manual(values = con_colors)

# 09.SASP ----
sasp.genenset <- rio::import('/data/02.project/00.other.lab/01.hualun/02.mus.age.placenta/03.analysis/SASP.geneset.xlsx')
sasp.geneuse <- sasp.genenset$`Gene ID` %>% tolower() %>% stringr::str_to_title() 
sasp.geneuse <- list(c(sasp.geneuse))
Clean_sct.inte.rm.pt <- AddModuleScore(Clean_sct.inte.rm.pt, features = sasp.geneuse,name = 'SASP')

p1 <- FeaturePlot(Clean_sct.inte.rm.pt, features = 'SASP1', split.by = 'condition', order = T, min.cutoff = 'q50', cols = c('grey90', "red3"), combine = F)
for(i in 1:length(p1)) {
  p1[[i]] <- p1[[i]] + NoLegend() + NoAxes() + theme(panel.background=element_rect(fill='transparent', color='black'), title = element_text(size = 8))
}
patchwork::wrap_plots(c(p1),nrow = 1)

sasp.score <- Clean_sct.inte.rm.pt@meta.data[,c('cell.type.percise.new','SASP1','condition')] %>% na.omit()
sasp.score$condition <- factor(sasp.score$condition, levels = c('Young','Aged_9','Aged_13'))

ggplot(sasp.score, aes(x=condition, y=SASP1, fill=condition)) + 
  geom_boxplot(alpha = .7) + 
  theme_bw()+
  labs(x = '', y = 'SASP score') + 
  theme(legend.position = "none") +
  scale_fill_manual(values = con_colors) +
  facet_wrap(~cell.type.percise.new, scale="free",nrow = 3) +
  ggpubr::stat_compare_means(comparisons = my_comparisons,paired = F,
                             method = "wilcox.test")

ggdensity(sasp.score, 
          x = "SASP1",
          add = "median", rug = F,
          color = "condition", fill = "condition",
          palette = c('#0B996F',"#D6570D")) +
  ylab("Density") + 
  xlab('SASP score') + 
  ggtitle('SASP') + 
  theme_bw() + 
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 14)) 
wilcox.test(sasp.score[sasp.score$condition == 'LP',]$SASP1, sasp.score[sasp.score$condition == 'HP',]$SASP1)
median(sasp.score[sasp.score$condition == 'LP' & sasp.score$cell.type.percise.new == 'DSCs',]$SASP1)
median(sasp.score[sasp.score$condition == 'HP' & sasp.score$cell.type.percise.new == 'DSCs',]$SASP1)


# 10.inflammatory reponse----
infm.genenset <- rio::import('/data/02.project/00.other.lab/01.hualun/02.mus.age.placenta/03.analysis/inflammatory_response.xlsx')
infm.geneuse <- infm.genenset$inflammatory.response
infm.geneuse <- list(c(infm.geneuse))
Clean_sct.inte.rm.pt <- AddModuleScore(Clean_sct.inte.rm.pt, features = infm.geneuse,name = 'infm')

p1 <- FeaturePlot(Clean_sct.inte.rm.pt, features = 'infm1', split.by = 'condition', order = T, min.cutoff = 'q50', cols = c('grey90', "red3"), combine = F)
for(i in 1:length(p1)) {
  p1[[i]] <- p1[[i]] + NoLegend() + NoAxes() + theme(panel.background=element_rect(fill='transparent', color='black'), title = element_text(size = 8))
}
patchwork::wrap_plots(c(p1),nrow = 1)

infm.score <- Clean_sct.inte.rm.pt@meta.data[,c('cell.type.percise.new','infm1','condition')] %>% na.omit()
infm.score$condition <- factor(infm.score$condition, levels = c('Young','Aged_9','Aged_13'))

ggplot(infm.score, aes(x=condition, y=infm1, fill=condition)) + 
  geom_boxplot(alpha = .7) + 
  theme_bw()+
  labs(x = '', y = 'infm score') + 
  theme(legend.position = "none") +
  scale_fill_manual(values = con_colors) +
  facet_wrap(~cell.type.percise.new, scale="free",nrow = 3) +
  ggpubr::stat_compare_means(comparisons = my_comparisons,paired = F,
                             method = "wilcox.test")

# 11.ca-DEGs cross condition -----

## up
comlist <- t(combn(c("Young","Aged_9","Aged_13"), 3))
deg_corss_condition_age_up <- list()
Idents(Clean_sct.inte.rm.pt) <- 'cell.type.percise.new'

for (ct in unique(Clean_sct.inte.rm.pt$cell.type.percise.new)) {
  cat("==> Processing cell type:", ct, "\n")

  sub.seu.1 <- subset(Clean_sct.inte.rm.pt, idents = ct)
  Idents(sub.seu.1) <- 'condition'
  for (gp in nrow(comlist)) {
    Young <- comlist[gp,1]
    Aged_9 <- comlist[gp,2]
    Aged_13 <- comlist[gp,3]
    DEGs_young_age9 <- FindMarkers(sub.seu.1,
                        ident.1 = Aged_9,
                        ident.2 = Young,
                        logfc.threshold = 0.1,
                        group.by = 'condition',
                        recorrect_umi = F,only.pos = T) %>% dplyr::filter(p_val_adj < 0.05)
    
    DEGs_age9_age13 <- FindMarkers(sub.seu.1,
                                   ident.1 = Aged_13,
                                   ident.2 = Aged_9,
                                   logfc.threshold = 0.2,
                                   group.by = 'condition',
                                   recorrect_umi = F,only.pos = T) %>% dplyr::filter(p_val_adj < 0.05)
    
    ca_DEGs_age_up <- intersect(rownames(DEGs_young_age9), rownames(DEGs_age9_age13))
    deg_corss_condition_age_up[[ct]] <- ca_DEGs_age_up
  }
}

library(purrr)
max_len <- max(map_int(deg_corss_condition_age_up, length))
deg_corss_condition_age_up.df <- map_dfr(deg_corss_condition_age_up, ~ { length(.x) <- max_len; .x })

deg_corss_condition_age_up.df <- list2DF(lapply(deg_corss_condition_age_up, `length<-`, max(lengths(deg_corss_condition_age_up))))
colnames(deg_corss_condition_age_up.df)
rio::export(deg_corss_condition_age_up.df, file = 'uterus.sc/ca_DGEs_age_up_by_findmarkers.xlsx')

## down
deg_corss_condition_age_down <- list()
for (ct in unique(Clean_sct.inte.rm.pt$cell.type.percise.new)) {
  cat("==> Processing cell type:", ct, "\n")
  
  sub.seu.1 <- subset(Clean_sct.inte.rm.pt, idents = ct)
  Idents(sub.seu.1) <- 'condition'
  for (gp in nrow(comlist)) {
    Young <- comlist[gp,1]
    Aged_9 <- comlist[gp,2]
    Aged_13 <- comlist[gp,3]
    DEGs_young_age9 <- FindMarkers(sub.seu.1,
                                   ident.1 = Young,
                                   ident.2 = Aged_9,
                                   logfc.threshold = 0.1,
                                   group.by = 'condition',
                                   recorrect_umi = F,only.pos = T) %>% dplyr::filter(p_val_adj < 0.05)
    
    DEGs_age9_age13 <- FindMarkers(sub.seu.1,
                                   ident.1 = Aged_9,
                                   ident.2 = Aged_13,
                                   logfc.threshold = 0.2,
                                   group.by = 'condition',
                                   recorrect_umi = F,only.pos = T) %>% dplyr::filter(p_val_adj < 0.05)
    
    ca_DEGs_age_down <- intersect(rownames(DEGs_young_age9), rownames(DEGs_age9_age13))
    deg_corss_condition_age_down[[ct]] <- ca_DEGs_age_down
  }
}

deg_corss_condition_age_down <- lapply(deg_corss_condition_age_down, function(vec) {
  data.frame(gene = vec, stringsAsFactors = FALSE)
})
rio::export(deg_corss_condition_age_down, file = 'uterus.sc/ca_DGEs_age_down_by_findmarkers.xlsx')


### ca-DEGs from pseudobulk
library(stringr)

agg <- AggregateExpression(
  object     = Clean_sct.inte.rm.pt,
  assay      = "RNA",
  slot       = "counts", 
  group.by   = c("sample", "cell.type.percise.new", "condition")
)
pb_counts <- agg$RNA %>% as.data.frame() %>% round(0)

coldata <- data.frame(
  sample_celltype = colnames(pb_counts),
  stringsAsFactors = FALSE
) %>%
  tidyr::separate(
    col = sample_celltype,
    into = c("sample", "cell.type.percise.new", "condition"),
    sep  = "_",
    remove = FALSE
  )
pb_list <- split(seq_len(ncol(pb_counts)), coldata$cell.type)
pseudo_bulk_by_ct <- lapply(pb_list, function(idxs){
  list(
    counts  = pb_counts[, idxs, drop = FALSE],
    coldata = coldata[idxs, ]
  )
})

library(DESeq2)
results_cross <- map(pseudo_bulk_by_ct, function(x){
  mat <- x$counts
  cd  <- x$coldata
  
  keep <- cd$condition %in% c("Young", "Aged-9", "Aged-13")
  mat  <- mat[, keep]
  cd   <- cd[keep, ]
  
  dds <- DESeqDataSetFromMatrix(
    countData = mat,
    colData   = cd %>% mutate(condition = factor(condition, 
                                                 levels = c("Young","Aged-9","Aged-13"))),
    design    = ~ condition
  )
  
  keep_gene <- rowSums(counts(dds) >= 10) >= 2
  dds <- dds[keep_gene, ]
  
  dds <- DESeq(dds, quiet = TRUE)
  # Aged_9 vs Young
  res1 <- results(dds,
                  contrast = c("condition","Aged-9","Young"),
                  alpha    = 0.05) %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    as_tibble() %>% 
    na.omit()
  
  # Aged_13 vs Aged_9
  res2 <- results(dds,
                  contrast = c("condition","Aged-13","Aged-9"),
                  alpha    = 0.05) %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    as_tibble() %>% 
    na.omit()
  
  up1   <- dplyr::filter(res1, pvalue < 0.05, log2FoldChange > 0.5) %>% pull(gene)
  down1 <- dplyr::filter(res1, pvalue < 0.05, log2FoldChange < -0.5) %>% pull(gene)
  up2   <- dplyr::filter(res2, pvalue < 0.05, log2FoldChange > 0.5) %>% pull(gene)
  down2 <- dplyr::filter(res2, pvalue < 0.05, log2FoldChange < -0.5) %>% pull(gene)
  
  list(
    up_intersect   = intersect(up1,   up2),
    down_intersect = intersect(down1, down2)
  )
})

results_cross$DCs$up_intersect
results_cross$DCs$down_intersect

summary_df <- imap_dfr(results_cross, ~ tibble(
  cell.type      = .y,
  direction      = c("up","down"),
  n.genes        = c(length(.x$up_intersect),
                     length(.x$down_intersect)),
  genes          = c(paste(.x$up_intersect, collapse = ","),
                     paste(.x$down_intersect, collapse = ","))
))
rio::export(summary_df, file = 'uterus.sc/ca-DEGs.pseudobulk.xlsx')

## pathway
hub<-AnnotationHub::AnnotationHub()
annolist <- AnnotationHub::query(hub, "Mus musculus")
ensdb110 <- hub[["AH113713"]]

plot.list <- list()
library(clusterProfiler)

for (ct in unique(summary_df$cell.type)) {
  
  genes <- results_cross[[ct]]$up_intersect
  gene2entrzid <- bitr(genes, fromType = 'GENENAME', toType = "ENTREZID", OrgDb = ensdb110)
  erich.go.BP <- enrichGO(gene=gene2entrzid$ENTREZID,
                          'org.Mm.eg.db',
                          pvalueCutoff = 0.01,
                          qvalueCutoff = 0.05,
                          keyType = 'ENTREZID',
                          readable = T,
                          ont = "BP")
  
  df2 <- erich.go.BP@result %>%
    filter(!is.na(GeneRatio), !is.na(p.adjust), !is.na(Count)) %>%
    slice_head(n = 10) %>%
    separate(GeneRatio, into = c("num","denom"), sep = "/", convert = TRUE) %>%
    mutate(GeneRatioNum = num / denom) %>%
    arrange(-GeneRatioNum) %>%
    mutate(Description = factor(Description, levels = rev(Description)))
  
  er.plot <- ggplot(df2, aes(x = GeneRatioNum, y = Description)) +
    geom_point(aes(size = Count, color = p.adjust)) +
    scale_size_continuous(name = "Count") +
    scale_color_continuous(name = "adj.P") +
    labs(x = "Gene Ratio", y = NULL, title = "Top 10 enriched BP terms") +
    theme_minimal()
  erich.name <- paste0(ct,'_plot')
  assign(erich.name, er.plot)
  plot.list[[erich.name]] <- er.plot
}
patchwork::wrap_plots(c(plot.list),nrow = 3)

# 12.pseodutime -----
library(monocle3)
library(ggpmisc)

## run mono
Idents(Clean_sct.inte.rm.pt)
getMonocds <- function(obj, assay, slot) {
  tab <- table(Clean_sct.inte.rm.pt$cell.type.percise.new)
  cell_types <- names(which(tab > 500))
  
  cds_list <- lapply(cell_types, function(celltype) {
    print(paste("Working on", celltype))
    
    tmp <- subset(obj, cells = WhichCells(obj, idents = celltype))
    mono <- GetAssayData(tmp, assay = assay, slot = slot)
    mono.cell_meta <- tmp@meta.data
    mono.gene_annotation <- data.frame(gene_short_name = rownames(mono))
    rownames(mono.gene_annotation) <- rownames(mono)
    
    cds <- new_cell_data_set(mono,
                             cell_metadata = mono.cell_meta,
                             gene_metadata = mono.gene_annotation) %>% 
      preprocess_cds( num_dim = 30) %>%
      reduce_dimension(preprocess_method = "PCA", cores = 90) %>% 
      cluster_cells(reduction_method = 'UMAP',
                    resolution = 0.001) %>% 
      learn_graph(use_partition = F)
    
    return(cds)
  })
  names(cds_list) <- cell_types
  return(cds_list)
}
sc.cds_list <- getMonocds(obj = Clean_sct.inte.rm.pt,assay = 'RNA',slot = 'counts')
save.image(compress = F)

## then analysis at each cell types
plot.ct <- 'Top2a+ Fibroblast'
cds <- sc.cds_list[[plot.ct]]
p1 <- plot_cells(cds,
                 reduction_method="UMAP", 
                 color_cells_by="condition",
                 label_cell_groups = F,
                 show_trajectory_graph = F,
                 cell_stroke = .2,
                 group_label_size = 5,
                 cell_size = .6) +
  ggtitle(plot.ct) +
  # scale_color_manual(values = scales::alpha(colour = con_colors,alpha = 0.7)) + 
  NoAxes()
p1
cds <- order_cells(cds)
p2 <- plot_cells(cds, color_cells_by = "pseudotime", 
                 label_cell_groups = F, 
                 label_leaves = FALSE,  
                 label_branch_points = F,
                 group_label_size = 1,
                 cell_size = .3,
                 label_roots = F,
                 cell_stroke = .2,
                 trajectory_graph_segment_size = 1
) + NoAxes() 

p1 + p2

## bin plot 
cdsmeta <- do.call(cbind, lapply(cds@colData@listData, as.data.frame))
colnames(cdsmeta) <- names(cds@colData@listData)
pseudotime.df <- data.frame(pseudotime = pseudotime(cds)) %>% 
  rownames_to_column(var = 'cells') %>% 
  left_join(y = cdsmeta[,c(21,22)], by = 'cells')

pseudotime.df$bin <- cut(pseudotime.df$pseudotime, 
                         breaks = seq(0, ceiling(max(pseudotime.df$pseudotime)), by = 1), 
                         include.lowest = TRUE, right = FALSE)
result <- pseudotime.df %>%
  dplyr::group_by(bin, condition) %>%
  dplyr::summarise(count = n()) %>%
  spread(key = condition, value = count, fill = 0) %>%
  dplyr::mutate(treat_ratio = Aged_13 / (Aged_13 + Aged_9 + Young)) 

result$bin_numeric <- as.numeric(result$bin)

model <- lm(treat_ratio ~ bin_numeric, data = result)
summary_model <- summary(model)
r_squared <- summary_model$r.squared
R <- cor(result$bin_numeric,result$treat_ratio, method = 'spearman')
p_value <- summary_model$coefficients[2, 4] 

p3 <- ggplot(result, aes(x = bin_numeric, y = treat_ratio)) +
  geom_point(size = 2.5,color = con_colors[3]) + 
  geom_smooth(method = "lm", se = TRUE, color = "gray70",alpha = .2) +  
  annotate("text", x = 5, y = 0.8, 
           label = paste("R = ", round(R, 3), "\np = ", format(p_value, digits = 3)),
           color = "black", size = 5) +  
  labs(x = "Pseudotime", y = "Aged 13 cell ratio") + 
  ggtitle(plot.ct) +
  theme_minimal() + 
  theme(axis.text = element_text(size = 15),axis.title = element_text(size = 18))


result <- pseudotime.df %>%
  dplyr::group_by(bin, condition) %>%
  dplyr::summarise(count = n()) %>%
  spread(key = condition, value = count, fill = 0) %>%
  dplyr::mutate(treat_ratio = Aged_9 / (Aged_13 + Aged_9 + Young)) 

result$bin_numeric <- as.numeric(result$bin)

model <- lm(treat_ratio ~ bin_numeric, data = result)
summary_model <- summary(model)
r_squared <- summary_model$r.squared
R <- cor(result$bin_numeric,result$treat_ratio, method = 'spearman')
p_value <- summary_model$coefficients[2, 4] 

p2 <- ggplot(result, aes(x = bin_numeric, y = treat_ratio)) +
  geom_point(size = 2.5,color = con_colors[2]) + 
  geom_smooth(method = "lm", se = TRUE, color = "gray70",alpha = .2) +  
  annotate("text", x = 5, y = 0.8, 
           label = paste("R = ", round(R, 3), "\np = ", format(p_value, digits = 3)),
           color = "black", size = 5) +  
  labs(x = "Pseudotime", y = "Aged 9 cell ratio") + 
  ggtitle(plot.ct) +
  theme_minimal() + 
  theme(axis.text = element_text(size = 15),axis.title = element_text(size = 18))


result <- pseudotime.df %>%
  dplyr::group_by(bin, condition) %>%
  dplyr::summarise(count = n()) %>%
  spread(key = condition, value = count, fill = 0) %>%
  dplyr::mutate(treat_ratio = Young / (Aged_13 + Aged_9 + Young)) 

result$bin_numeric <- as.numeric(result$bin)

model <- lm(treat_ratio ~ bin_numeric, data = result)
summary_model <- summary(model)
r_squared <- summary_model$r.squared
R <- cor(result$bin_numeric,result$treat_ratio, method = 'spearman')
p_value <- summary_model$coefficients[2, 4] 

p1 <- ggplot(result, aes(x = bin_numeric, y = treat_ratio)) +
  geom_point(size = 2.5,color = con_colors[1]) + 
  geom_smooth(method = "lm", se = TRUE, color = "gray70",alpha = .2) +  
  annotate("text", x = 5, y = 0.8, 
           label = paste("R = ", round(R, 3), "\np = ", format(p_value, digits = 3)),
           color = "black", size = 5) +  
  labs(x = "Pseudotime", y = "Young cell ratio") + 
  ggtitle(plot.ct) +
  theme_minimal() + 
  theme(axis.text = element_text(size = 15),axis.title = element_text(size = 18))
p1 + p2 + p3

## heatmap
library(Mfuzz)
modulated_genes <- graph_test(cds, neighbor_graph = "principal_graph", cores = 40)
genes <- row.names(subset(modulated_genes, q_value <= 0.0001 & morans_I > 0.1))
mat <- pre_pseudotime_matrix(cds_obj = cds,
                             gene_list = genes)

ck <- clusterData_new(data = mat,cluster.num = 5,method = 'kmeans')
# names(ck)[4] <- 'type'
sasp.geneuse[[1]][match(sasp.geneuse[[1]], rownames(mat))]
pdf(paste0('./25.6.add/monocle3.heatmap.',plot.ct,'.pdf'),height = 10,width = 8,onefile = F)
visCluster(object = ck,
           plot.type = "both",
           add.sampleanno = F,
           markGenes = c('Inha','Cxcl8','Fas','Cxcl12','Cxcl1','Il12b','Il6r','Igf1','Cd38'))
dev.off()
cluster.list <- lapply(ck$cluster.list, function(vec) {
  data.frame(gene = vec, stringsAsFactors = FALSE)
})
names(cluster.list) <- paste0("Cluster_", names(cluster.list))
rio::export(cluster.list, file = paste0('25.6.add/',plot.ct,'.heatmap.cluster.genes.xlsx'))

# 14.NAD pathway -----
library(KEGGREST)
nad.pathway <- keggGet('mmu00760')
nad.pathway <- sapply(seq(2,86,2), function(x){
  gns <- unlist(strsplit(nad.pathway[[1]][['GENE']][x],";"))[1]
})

nad.geneuse <- list(c(nad.pathway))
Clean_sct.inte.rm.pt <- AddModuleScore(Clean_sct.inte.rm.pt, features = nad.geneuse,name = 'nad')

p1 <- FeaturePlot(Clean_sct.inte.rm.pt, features = 'nad1', split.by = 'condition', order = T, min.cutoff = 'q10', cols = c('grey90', "red3"), combine = F)
for(i in 1:length(p1)) {
  p1[[i]] <- p1[[i]] + NoLegend() + NoAxes() + theme(panel.background=element_rect(fill='transparent', color='black'), title = element_text(size = 8))
}
patchwork::wrap_plots(c(p1),nrow = 1)

nad.score <- Clean_sct.inte.rm.pt@meta.data[,c('cell.type.percise.new','nad1','condition')] %>% na.omit()
nad.score$condition <- factor(nad.score$condition, levels = c('Young','Aged_9','Aged_13'))

ggplot(nad.score, aes(x=condition, y=nad1, fill=condition)) + 
  geom_boxplot(alpha = .7,outlier.size = .5) + 
  theme_bw()+
  labs(x = '', y = 'nad score') + 
  theme(legend.position = "none") +
  scale_fill_manual(values = con_colors) +
  facet_wrap(~cell.type.percise.new, scale="free",nrow = 3) +
  ggpubr::stat_compare_means(comparisons = my_comparisons,paired = F,
                             method = "wilcox.test")

# 16.cellchat -----
library(CellChat)
cc.cols <- c(use.cols,npg.cols)[1:18]
names(cc.cols) <- levels(Clean_sct.inte.rm.pt$cell.type.percise.new)
names(cluster_colors) <- levels(Clean_sct.inte.rm.pt$cell.type.minor.fibro)

# Young
split.seu <- SplitObject(Clean_sct.inte.rm.pt,split.by = 'condition')
CellChatDB <- CellChatDB.mouse
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") 

Young.cc <- createCellChat(object = split.seu$Young, meta = split.seu$Young@meta.data, group.by = "cell.type.percise.new")
Young.cc@DB <- CellChatDB.use
Young.cc <- subsetData(Young.cc) %>% 
  identifyOverExpressedGenes() %>% 
  identifyOverExpressedInteractions() %>% 
  computeCommunProb(type = "triMean") %>% 
  filterCommunication( min.cells = 10) %>% 
  computeCommunProbPathway() %>% 
  aggregateNet()

groupSize <- as.numeric(table(Young.cc@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(Young.cc@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(Young.cc@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

## young cc add fibro
Young.cc.add.fibro <- createCellChat(object = split.seu$Young, meta = split.seu$Young@meta.data, group.by = "cell.type.minor.fibro")
Young.cc.add.fibro@DB <- CellChatDB.use
Young.cc.add.fibro <- subsetData(Young.cc.add.fibro) %>% 
  identifyOverExpressedGenes() %>% 
  identifyOverExpressedInteractions() %>% 
  computeCommunProb(type = "triMean") %>% 
  filterCommunication( min.cells = 10) %>% 
  computeCommunProbPathway() %>% 
  aggregateNet()

### aged9
Aged_9.cc <- createCellChat(object = split.seu$Aged_9, meta = split.seu$Aged_9@meta.data, group.by = "cell.type.percise.new")
Aged_9.cc@DB <- CellChatDB.use
Aged_9.cc <- subsetData(Aged_9.cc) %>% 
  identifyOverExpressedGenes() %>% 
  identifyOverExpressedInteractions() %>% 
  computeCommunProb(type = "triMean") %>% 
  filterCommunication( min.cells = 10) %>% 
  computeCommunProbPathway() %>% 
  aggregateNet()

## Aged_9 cc add fibro
Aged_9.cc.add.fibro <- createCellChat(object = split.seu$Aged_9, meta = split.seu$Aged_9@meta.data, group.by = "cell.type.minor.fibro")
Aged_9.cc.add.fibro@DB <- CellChatDB.use
Aged_9.cc.add.fibro <- subsetData(Aged_9.cc.add.fibro) %>% 
  identifyOverExpressedGenes() %>% 
  identifyOverExpressedInteractions() %>% 
  computeCommunProb(type = "triMean") %>% 
  filterCommunication( min.cells = 10) %>% 
  computeCommunProbPathway() %>% 
  aggregateNet()

### aged13
Aged_13.cc <- createCellChat(object = split.seu$Aged_13, meta = split.seu$Aged_13@meta.data, group.by = "cell.type.percise.new")
Aged_13.cc@DB <- CellChatDB.use
Aged_13.cc <- subsetData(Aged_13.cc) %>% 
  identifyOverExpressedGenes() %>% 
  identifyOverExpressedInteractions() %>% 
  computeCommunProb(type = "triMean") %>% 
  filterCommunication( min.cells = 10) %>% 
  computeCommunProbPathway() %>% 
  aggregateNet()

## Aged_13 cc add fibro
Aged_13.cc.add.fibro <- createCellChat(object = split.seu$Aged_13, meta = split.seu$Aged_13@meta.data, group.by = "cell.type.minor.fibro")
Aged_13.cc.add.fibro@DB <- CellChatDB.use
Aged_13.cc.add.fibro <- subsetData(Aged_13.cc.add.fibro) %>% 
  identifyOverExpressedGenes() %>% 
  identifyOverExpressedInteractions() %>% 
  computeCommunProb(type = "triMean") %>% 
  filterCommunication( min.cells = 10) %>% 
  computeCommunProbPathway() %>% 
  aggregateNet()

### merge--young_age9
cc.list <- list(Young = Young.cc, Aged_9 = Aged_9.cc)
young_aged9.cc <- mergeCellChat(cc.list, add.names = names(cc.list))

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(young_aged9.cc, weight.scale = T, measure = "weight",color.use = cc.cols,sources.use = 'Macrophages')
netVisual_diffInteraction(young_aged9.cc, weight.scale = T, measure = "weight",color.use = cc.cols,sources.use = 'Macrophages')
gg1 <- netVisual_heatmap(young_aged9.cc,color.use = cc.cols)
gg2 <- netVisual_heatmap(young_aged9.cc, measure = "weight",color.use = cc.cols)
gg1 + gg2

### merge--young_age9 add fibro
cc.list <- list(Young = Young.cc.add.fibro, Aged_9 = Aged_9.cc.add.fibro)
young_aged9.cc.add.fibro <- mergeCellChat(cc.list, add.names = names(cc.list))

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(young_aged9.cc.add.fibro, weight.scale = T, measure = "weight",color.use = cluster_colors,sources.use = 'Macrophages')
netVisual_diffInteraction(young_aged9.cc.add.fibro, weight.scale = T, measure = "weight",color.use = cluster_colors,sources.use = 'Fibroblast')
gg1 <- netVisual_heatmap(young_aged9.cc.add.fibro,color.use = cluster_colors)
gg2 <- netVisual_heatmap(young_aged9.cc.add.fibro, measure = "weight",color.use = cluster_colors)
gg1 + gg2

### merge cell type
group.cellType <- c(rep("FIB", 1), rep("Sox9", 1), rep("Lumenal", 1), rep('EC',1),rep('Meyoid',4),rep('Lymph',7))
group.cellType <- factor(group.cellType, levels = c("FIB", "Sox9", "Lumenal",'EC','Meyoid','Lymph'))
group.cellType <- c(rep("FIB", 1), rep("Sox9", 1), rep("Lumenal", 1), rep('EC',1),rep('Immune',11))
group.cellType <- factor(group.cellType, levels = c("FIB", "Sox9", "Lumenal",'EC','Immune'))

cc.list <- list(Young = Young.cc.add.fibro, Aged_9 = Aged_9.cc.add.fibro)
cc.list <- lapply(cc.list, function(x) {mergeInteractions(x, group.cellType)})
young_aged9.cc.add.fibro.merge <- mergeCellChat(cc.list, add.names = names(cc.list))
netVisual_diffInteraction(young_aged9.cc.add.fibro.merge, weight.scale = T, measure = "count.merged", label.edge = F)
netVisual_diffInteraction(young_aged9.cc.add.fibro.merge, weight.scale = T, measure = "weight.merged", label.edge = F)

### merge--young_aged13
cc.list <- list(Young = Young.cc, Aged_13 = Aged_13.cc)
young_aged13.cc <- mergeCellChat(cc.list, add.names = names(cc.list))
gg1 <- compareInteractions(young_aged13.cc, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(young_aged13.cc, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(young_aged13.cc, weight.scale = T, measure = "weight",color.use = cc.cols,sources.use = 'Macrophages')

gg1 <- netVisual_heatmap(young_aged13.cc,color.use = cc.cols,)
gg2 <- netVisual_heatmap(young_aged13.cc, measure = "weight",color.use = cc.cols)
gg1 + gg2

### merge--young_age13 add fibro
cc.list <- list(Young = Young.cc.add.fibro, Aged_13 = Aged_13.cc.add.fibro)
young_aged13.cc.add.fibro <- mergeCellChat(cc.list, add.names = names(cc.list))

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(young_aged13.cc.add.fibro, weight.scale = T, measure = "weight",color.use = cluster_colors,sources.use = 'Macrophages')
netVisual_diffInteraction(young_aged13.cc.add.fibro, weight.scale = T, measure = "weight",color.use = cluster_colors,sources.use = 'Fibroblast')
gg1 <- netVisual_heatmap(young_aged13.cc.add.fibro,color.use = cluster_colors)
gg2 <- netVisual_heatmap(young_aged13.cc.add.fibro, measure = "weight",color.use = cluster_colors)
gg1 + gg2

### merge cell type
cc.list <- list(Young = Young.cc.add.fibro, Aged_13 = Aged_13.cc.add.fibro)
cc.list <- lapply(cc.list, function(x) {mergeInteractions(x, group.cellType)})
young_aged13.cc.add.fibro.merge <- mergeCellChat(cc.list, add.names = names(cc.list))
netVisual_diffInteraction(young_aged13.cc.add.fibro.merge, weight.scale = T, measure = "weight.merged", label.edge = F)

### merge--aged9_aged13
cc.list <- list( Aged_9 = Aged_9.cc, Aged_13 = Aged_13.cc)
aged9_aged13.cc <- mergeCellChat(cc.list, add.names = names(cc.list))

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(aged9_aged13.cc, weight.scale = T, measure = "weight",color.use = cc.cols,sources.use = 'Macrophages')
gg1 <- netVisual_heatmap(aged9_aged13.cc,color.use = cc.cols,)
gg2 <- netVisual_heatmap(aged9_aged13.cc, measure = "weight",color.use = cc.cols)
gg1 + gg2

### merge--aged9_age13 add fibro
cc.list <- list(Aged_9 = Aged_9.cc.add.fibro, Aged_13 = Aged_13.cc.add.fibro)
age9_aged13.cc.add.fibro <- mergeCellChat(cc.list, add.names = names(cc.list))

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(age9_aged13.cc.add.fibro, weight.scale = T, measure = "weight",color.use = cluster_colors,sources.use = 'Macrophages')
netVisual_diffInteraction(age9_aged13.cc.add.fibro, weight.scale = T, measure = "weight",color.use = cluster_colors,sources.use = 'Fibroblast')
gg1 <- netVisual_heatmap(age9_aged13.cc.add.fibro,color.use = cluster_colors)
gg2 <- netVisual_heatmap(age9_aged13.cc.add.fibro, measure = "weight",color.use = cluster_colors)
gg1 + gg2

### merge cell type
cc.list <- list(Aged_9 = Aged_9.cc.add.fibro, Aged_13 = Aged_13.cc.add.fibro)
cc.list <- lapply(cc.list, function(x) {mergeInteractions(x, group.cellType)})
age9_aged13.cc.add.fibro.merge <- mergeCellChat(cc.list, add.names = names(cc.list))
netVisual_diffInteraction(age9_aged13.cc.add.fibro.merge, weight.scale = T, measure = "weight.merged", label.edge = F)

# 19. young cell type score ----
young.seu <- subset(Clean_sct.inte.rm.pt, condition == 'Young')
table(young.seu$cell.type.percise.new)
DimPlot(young.seu,group.by = 'cell.type.percise.new',label = T)
young.ct.hvgs <- FindAllMarkers(Clean_sct.inte.rm.pt, only.pos = T,min.pct = 0.1,logfc.threshold = 0.75,recorrect_umi = F) %>% 
  dplyr::filter(p_val < 0.05)
young.ct.hvgs <- young.ct.hvgs[young.ct.hvgs$avg_log2FC > 1,]

for (ct in as.character(levels(young.ct.hvgs$cluster))) {
  cat("working on", ct, "\n")
  Clean_sct.inte.rm.pt <- AddModuleScore(Clean_sct.inte.rm.pt, features = list(c(young.ct.hvgs[young.ct.hvgs$cluster == ct,]$gene)),name = paste0('Young_',ct,'_score'))
}

young.score.ct <- data.frame(condition = 0,cell.type = 0, score = 0)
for (ct in as.character(levels(young.ct.hvgs$cluster))) {
  cat("working on", ct, "\n")
  tmp <- Clean_sct.inte.rm.pt@meta.data[Clean_sct.inte.rm.pt$cell.type.percise.new == ct, c('condition','cell.type.percise.new',paste0('Young_',ct,'_score1'))]
  colnames(tmp) <- c('condition','cell.type','score')
  young.score.ct <- rbind(young.score.ct,tmp)
}
young.score.ct <- young.score.ct[-1,]

ggplot(young.score.ct, aes(x=condition, y=score, fill=condition)) + 
  geom_boxplot(alpha = .7,outlier.size = .5) + 
  theme_bw()+
  labs(x = '', y = 'Young score') + 
  theme(legend.position = "none") +
  scale_fill_manual(values = con_colors) +
  facet_wrap(~cell.type, scale="free",nrow = 3) +
  ggpubr::stat_compare_means(comparisons = my_comparisons,paired = F,
                             method = "wilcox.test")

# 20.NC metabolism pathway ----
kegg.pathways <- read.gmt('../../KEGG_metabolism_nc.gmt')
kegg.pathways$gene <- tolower(kegg.pathways$gene ) %>% stringr::str_to_title() 
# kegg.pathways <- kegg.pathways[kegg.pathways$term == 'Pentose phosphate pathway',]

for (pathway in unique(kegg.pathways$term)) {
  cat("working on", pathway, "\n")
  Clean_sct.inte.rm.pt <- AddModuleScore(Clean_sct.inte.rm.pt, features = list(c(kegg.pathways[kegg.pathways$term == pathway,]$gene)),name = pathway)
}

meta <- Clean_sct.inte.rm.pt@meta.data %>%
  as.data.frame() 

pathways <- paste0(unique(kegg.pathways$term),'1')

df_long <- meta %>%
  pivot_longer(
    cols      = all_of(pathways),
    names_to  = "pathway",
    values_to = "score"
  )

df_stat <- df_long %>%
  group_by(condition, pathway) %>%
  summarise(
    avg_score    = mean(score, na.rm = TRUE),        
    prop_nonzero = mean(score > 0, na.rm = TRUE),  
    .groups      = "drop"
  )

ggplot(df_stat, aes(x = condition, y = pathway)) +
  geom_point(aes(
    color = avg_score,
    size  = prop_nonzero*0.7
  )) +
  scale_color_gradient2(
    low     = "blue",
    mid     = "white",
    high    = "red",
    midpoint= 0,         
    name    = "Avg score"
  ) +
  scale_size_continuous(
    range = c(0.5, 4),
    name  = "Prop non-zero"
  ) +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    axis.title   = element_blank(),
    panel.grid   = element_blank()
  ) 

## add cell type 
library(tidyverse)
library(ggplot2)

df_stat <- df_long %>%
  group_by(cell.type.percise.new, condition, pathway) %>%
  summarise(
    avg_score    = mean(score, na.rm = TRUE),
    prop_nonzero = mean(score > 0, na.rm = TRUE),
    .groups      = "drop"
  ) %>%
  
  mutate(
    cell.type.percise.new = factor(cell.type.percise.new,
                                   levels = unique(cell.type.percise.new)),
    condition = factor(condition,
                       levels = unique(condition)),
    pathway   = factor(pathway,
                       levels = unique(pathway))
  )


ggplot(df_stat, aes(x = cell.type.percise.new, y = condition)) +
  geom_point(aes(
    color = avg_score,
    size  = prop_nonzero
  ), alpha = 0.8) +
  scale_color_gradient2(
    low      = "blue",
    mid      = "white",
    high     = "red",
    midpoint = median(df_stat$avg_score, na.rm = TRUE),
    name     = "Avg score"
  ) +
  scale_size_continuous(
    range = c(0.3, 4),
    name  = "Prop non-zero"
  ) +
  facet_wrap(~ pathway, ncol = 4) +   
  theme_minimal() +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    axis.title   = element_blank(),
    strip.text   = element_text(size = 8),
    panel.grid   = element_blank()
  )


## scale 
df_scaled <- meta %>%
  pivot_longer(all_of(pathways),
               names_to  = "pathway",
               values_to = "raw_score") %>%
  group_by(pathway) %>%
  mutate(scaled = (raw_score - mean(raw_score, na.rm=TRUE)) / sd(raw_score, na.rm=TRUE)) %>%
  ungroup()

df_stat2 <- df_scaled %>%
  group_by(cell.type.percise.new, condition, pathway) %>%
  summarise(
    avg_scaled   = mean(scaled, na.rm = TRUE),
    prop_nonzero = mean(raw_score > 0, na.rm = TRUE),
    .groups      = "drop"
  )

ggplot(df_stat2, aes(x = cell.type.percise.new, y = condition)) +
  geom_point(aes(
    color = avg_scaled,
    size  = prop_nonzero
  ), alpha = 0.8) +
  scale_color_gradient2(
    low      = "blue",
    mid      = "white",
    high     = "red",
    midpoint = 0,
    name     = "Avg scaled\nscore"
  ) +
  scale_size_continuous(
    range = c(0.3, 4),
    name  = "Prop non-zero"
  ) +
  facet_wrap(~ pathway, ncol = 4) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title  = element_blank(),
    strip.text  = element_text(size = 8),
    panel.grid  = element_blank()
  )
