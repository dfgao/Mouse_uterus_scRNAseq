
## LTSR ----

plot_ranef.new <- function(ranef_tbl, vars, celltypes = NULL, celltype_order = "hclust", references = NULL,
                           maxFC = 3, LTSR2p = F, highlightLtsr = 0.0, filterLtsr = 0.0, swap_axes = F) {
  ranef_tbl <- .getCondValLtsr(ranef_tbl, vars, celltypes = celltypes, references = references)
  save(ranef_tbl, file = "ranef_tbl.RData")
  condval_mat <- ranef_tbl %>%
    dplyr::select(
      Celltype, grpval, condval
    ) %>%
    spread(
      "grpval", "condval"
    ) %>%
    column_to_rownames(
      var = "Celltype"
    ) %>%
    as.matrix()
  if (length(celltype_order) == 1 && celltype_order == "hclust") {
    dendy <- hclust(dist(condval_mat))
    ordered_celltype <- rownames(condval_mat)[dendy$ord]
  } else if (!is.null(celltype_order) && length(celltype_order) == dim(condval_mat)[1]) {
    ordered_celltype <- celltype_order
  }
  
  ranef_tbl <- ranef_tbl %>% mutate(
    Celltype = factor(Celltype, levels = ordered_celltype),
    condval = condval %>% pmin(log(maxFC)) %>% pmax(log(1 / maxFC)),
    ltsr = ltsr %>% pmin(0.9999) %>% pmax(0.5)
  )
  
  if (swap_axes) {
    ranef_tbl$Celltype <- factor(ranef_tbl$Celltype, levels = rev(levels(ranef_tbl$Celltype)))
    ranef_tbl$grpval <- factor(ranef_tbl$grpval, levels = rev(levels(ranef_tbl$grpval)))
  }
  
  if (filterLtsr > 0) {
    filtered_celltypes <- ranef_tbl %>%
      group_by(Celltype) %>%
      summarise(maxLtsr = max(ltsr)) %>%
      dplyr::filter(maxLtsr >= filterLtsr) %>%
      dplyr::select(Celltype) %>%
      unlist(use.names = F)
    ranef_tbl <- ranef_tbl %>% dplyr::filter(Celltype %in% filtered_celltypes)
  }
  
  geom_dots <- geom_point(
    aes(
      fill = log2(exp(condval)),
      size = -log10(1 - ltsr)
    ),
    color = "white",
    shape = 21
  )
  
  if (swap_axes) {
    p <- (
      ggplot(ranef_tbl, aes(y = grpval, x = Celltype)) +
        facet_grid(grpvar ~ ., scales = "free_y", space = "free_y", switch = "x") +
        geom_dots
    )
  } else {
    p <- (
      ggplot(ranef_tbl, aes(x = grpval, y = Celltype)) +
        facet_grid(. ~ grpvar, scales = "free_x", space = "free_x", switch = "x") +
        geom_dots
    )
  }
  
  p <- (
    p + scale_fill_distiller(
      palette = "RdBu",
      limits = log2(c(1 / maxFC, maxFC)),
      breaks = log2(c(1 / maxFC, maxFC)),
      labels = c(paste0("1/", maxFC), maxFC),
      oob = squish,
      guide = guide_colorbar(
        title = "Fold change", title.position = "top", direction = "horizontal",
        barwidth = 5, barheight = 0.75, raster = F, order = 1)
    )
    + scale_size(
      limits = -log10(1 - c(0.5, 0.9999)),
      breaks = -log10(1 - c(0.5, 0.7, 0.9, 0.99)),
      range = c(0.5, 9),
      labels = ifelse(
        rep(LTSR2p, 4),
        c("0.5", "0.3", "0.1", "<0.01"),
        c("0.5", "0.7", "0.9", ">0.99")
      ),
      guide = guide_legend(
        title = ifelse(LTSR2p, "p", "LTSR"), reverse = T, order = 2,
        override.aes = list(fill = "black", color = "white")
      )
    )
    + theme_bw()
    + theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      strip.placement = "outside",
      strip.background = element_blank(),
      legend.spacing.y = unit(0.5, "line")
    )
  )
  
  if (highlightLtsr > 0) {
    p <- (
      p + geom_point(
        aes(
          color = ltsr > highlightLtsr,
          alpha = ltsr > highlightLtsr
        ),
        shape = 21, size = 9
      )
      + scale_color_manual(
        label = c(
          "", ifelse(LTSR2p, paste0("p < ", 1 - highlightLtsr), paste0("LTSR > ", highlightLtsr))
        ),
        values = c("white", "red"),
        guide = guide_legend(title = NULL, override.aes = list(size = 9), reverse = T, order = 3)
      )
      + scale_alpha_manual(values = c(0, 1), guide = F)
    )
  }
  
  p
}


## circle plot ----

plot.cir.test <- function (data_plot, do.label = T, contour.levels = c(0.2, 0.4, 0.6), 
                           pt.size = 0.5, kde2d.n = 1000, contour.nlevels = 100, bg.color = "#F9F2E4", 
                           col.use = NULL, label.cex = 0.5, repel = FALSE) 
{
  centers <- data_plot %>% dplyr::group_by(Cluster) %>% summarise(x = median(x = x), 
                                                                  y = median(x = y))
  z <- MASS::kde2d(data_plot$x, data_plot$y, n = kde2d.n)
  celltypes <- names(table(data_plot$Cluster))
  cell_colors <- (scales::hue_pal())(length(celltypes))
  if (!is.null(col.use)) {
    cell_colors = col.use
    col_df <- data.frame(Cluster = celltypes, color2 = col.use)
    cells_order <- rownames(data_plot)
    data_plot <- merge(data_plot, col_df, by = "Cluster")
    rownames(data_plot) <- data_plot$cells
    data_plot <- data_plot[cells_order, ]
    data_plot$Colors <- data_plot$color2
  }
  circos.clear()
  par(bg = bg.color)
  circos.par(cell.padding = c(0, 0, 0, 0), track.margin = c(0.01, 0), track.height = 0.01, gap.degree = c(rep(2, (length(celltypes) -  1)), 12), points.overflow.warning = FALSE)
  # circos.clear()
  # circos.par(
  #   cell.padding = c(0, 0, 0, 0),
  #   track.margin = c(0.01, 0),
  #   canvas.xlim = c(-1.2, 1.2),
  #   canvas.ylim = c(-1.2, 1.2)
  # )
  circos.initialize(sectors = data_plot$Cluster, x = data_plot$x_polar2)
  circos.track(data_plot$Cluster, data_plot$x_polar2, y = data_plot$dim2, 
               bg.border = NA, panel.fun = function(x, y) {
                 circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + 
                               mm_y(4), CELL_META$sector.index, cex = 0.5, 
                             col = "black", facing = "bending.inside", niceFacing = T)
                 breaks = seq(0, 100, by = 50)
                 circos.axis(labels.cex = 0.3, col = "black", labels.col = "black",major.at = breaks,
                             labels = paste0(breaks, "%"))
               })
  for (i in 1:length(celltypes)) {
    dd <- data_plot[data_plot$Cluster == celltypes[i], ]
    circos.segments(x0 = min(dd$x_polar2), y0 = 0, x1 = max(dd$x_polar2), 
                    y1 = 0, col = cell_colors[i], lwd = 3, sector.index = celltypes[i])
  }
  text(x = 1, y = 0.1, labels = "Cluster", cex = 0.4, col = "black", 
       srt = -90)
  points(data_plot$x, data_plot$y, pch = 19, col = alpha(data_plot$Colors, 
                                                         0.2), cex = pt.size)
  contour(z, drawlabels = F, nlevels = 100, levels = contour.levels, 
          col = "#ae9c76", add = TRUE)
  if (do.label) {
    if (repel) {
      textplot(x = centers$x, y = centers$y, words = centers$Cluster, 
               cex = label.cex, new = F, show.lines = F)
    }
    else {
      text(centers$x, centers$y, labels = centers$Cluster, 
           cex = label.cex, col = "black")
    }
  }
}

cell_order_test <- function (dat) 
{
  celltypes <- names(table(dat$Cluster))
  new_dat <- list()
  for (i in 1:length(celltypes)) {
    dat$Cluster <- as.character(dat$Cluster)
    dat1 <- dat[dat$Cluster == celltypes[i], ]
    dat1$x_polar <- (1:nrow(dat1))/nrow(dat1)
    new_dat[[i]] <- dat1
  }
  new_dat <- do.call("rbind", new_dat)
  new_dat
}

prepare_circlize_data_test <- function (seu_obj, scale = 1) 
{
  celltypes <- levels(seu_obj)
  cell_colors <- (scales::hue_pal())(length(celltypes))
  data_plot <- get_metadata(seu_obj, color = cell_colors, 
                            coord_scale = scale)
  data_plot <- cell_order_test(data_plot)
  data_plot$x_polar2 <- log10(data_plot$x_polar)
  data_plot
}


## dotplot ----

complex_dotplot_single <- function (seu_obj, feature, celltypes = NULL, groups, splitby = NULL, 
                                    color.palette = NULL, font.size = 12, strip.color = NULL, 
                                    do.scale = T, scale.by = "radius") 
{
  if (is.null(color.palette)) {
    color.palette <- colorRampPalette(c("grey80", "lemonchiffon1", 
                                        "indianred1", "darkred"))(255)
  }
  scale.func <- switch(EXPR = scale.by, size = scale_size, 
                       radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
  if (is.null(celltypes)) {
    celltypes <- levels(seu_obj)
  }
  if (length(groups) == 1) {
    groups_level <- levels(seu_obj@meta.data[, groups])
    if (is.null(groups_level)) {
      seu_obj@meta.data[, groups] <- factor(seu_obj@meta.data[, 
                                                              groups], levels = names(table(seu_obj@meta.data[, 
                                                                                                              groups])))
      groups_level <- levels(seu_obj@meta.data[, groups])
    }
    if (!is.null(splitby)) {
      if (is.null(levels(seu_obj@meta.data[, splitby]))) {
        seu_obj@meta.data[, splitby] <- factor(seu_obj@meta.data[, 
                                                                 splitby], levels = names(table(seu_obj@meta.data[, 
                                                                                                                  splitby])))
      }
      splitby_level <- levels(seu_obj@meta.data[, splitby])
      count_df <- extract_gene_count(seu_obj, features = feature, 
                                     cell.types = celltypes, meta.groups = c(groups, 
                                                                             splitby))
      count_df$new_group <- paste(count_df[, groups], 
                                  count_df[, "celltype"], count_df[, splitby], 
                                  sep = "___")
      exp_df <- aggregate(. ~ new_group, data = count_df[, 
                                                         c("new_group", feature)], FUN = function(x) {
                                                           mean(expm1(x))
                                                         })
      pct_df <- aggregate(. ~ new_group, data = count_df[, 
                                                         c("new_group", feature)], FUN = function(x) {
                                                           length(x[x > 0])/length(x)
                                                         })
      colnames(exp_df)[2] <- "avg.exp"
      colnames(pct_df)[2] <- "pct.exp"
      data_plot <- merge(exp_df, pct_df, by = "new_group")
      data_plot$groups <- as.character(lapply(X = strsplit(data_plot$new_group, 
                                                           split = "___"), FUN = function(x) {
                                                             x[[1]]
                                                           }))
      data_plot$celltype <- as.character(lapply(X = strsplit(data_plot$new_group, 
                                                             split = "___"), FUN = function(x) {
                                                               x[[2]]
                                                             }))
      data_plot$splitby <- as.character(lapply(X = strsplit(data_plot$new_group, 
                                                            split = "___"), FUN = function(x) {
                                                              x[[3]]
                                                            }))
      data_plot$groups <- factor(data_plot$groups, levels = groups_level)
      data_plot$splitby <- factor(data_plot$splitby, levels = splitby_level)
      data_plot$celltype <- factor(data_plot$celltype, 
                                   levels = rev(celltypes))
    }
    else {
      count_df <- extract_gene_count(seu_obj, features = feature, 
                                     cell.types = celltypes, meta.groups = groups)
      count_df$new_group <- paste(count_df[, groups], 
                                  count_df[, "celltype"], sep = "___")
      exp_df <- aggregate(. ~ new_group, data = count_df[, 
                                                         c("new_group", feature)], FUN = function(x) {
                                                           mean(expm1(x))
                                                         })
      pct_df <- aggregate(. ~ new_group, data = count_df[, 
                                                         c("new_group", feature)], FUN = function(x) {
                                                           length(x[x > 0])/length(x)
                                                         })
      colnames(exp_df)[2] <- "avg.exp"
      colnames(pct_df)[2] <- "pct.exp"
      data_plot <- merge(exp_df, pct_df, by = "new_group")
      data_plot$groups <- as.character(lapply(X = strsplit(data_plot$new_group, 
                                                           split = "___"), FUN = function(x) {
                                                             x[[1]]
                                                           }))
      data_plot$celltype <- as.character(lapply(X = strsplit(data_plot$new_group, 
                                                             split = "___"), FUN = function(x) {
                                                               x[[2]]
                                                             }))
      data_plot$groups <- factor(data_plot$groups, levels = groups_level)
      data_plot$celltype <- factor(data_plot$celltype, 
                                   levels = rev(celltypes))
    }
    data_plot$pct.exp <- round(100 * data_plot$pct.exp, 
                               2)
    data_plot$avg.exp <- scale(data_plot$avg.exp)
    p <- ggplot(data_plot, aes(y = celltype, x = groups)) + 
      geom_tile(fill = "white", color = "white") + geom_point(aes(colour = avg.exp, 
                                                                  size = pct.exp)) + scale_color_gradientn(colours = color.palette) + 
      theme(panel.background = element_rect(fill = "white", 
                                            colour = "black"), axis.text.x = element_text(angle = 45, 
                                                                                          hjust = 1, size = font.size), plot.title = element_text(size = (font.size + 
                                                                                                                                                            2), hjust = 0.5, face = "bold"), axis.text = element_text(size = font.size), 
            legend.text = element_text(size = (font.size - 
                                                 2)), legend.title = element_text(size = (font.size)), 
            strip.text = element_text(size = font.size), 
            legend.position = "right") + ylab("") + xlab("") + 
      ggtitle(feature)
    if (do.scale) {
      p = p + scale_size(range = c(0, 10))
    }
    else {
      if (max(data_plot$pct.exp) >= 20) {
        p = p + scale_size(range = c(0, 10))
      }
      else {
        p = p + scale.func(range = c(0, 10), limits = c(0, 
                                                        20))
      }
    }
    if (!is.null(splitby)) {
      p <- p + facet_wrap(~splitby, scales = "free_x")
      g <- change_strip_background(p, type = "top", strip.color = strip.color)
      print(grid.draw(g))
    }
    else {
      p
    }
  }
  else {
    gene_count <- extract_gene_count(seu_obj = seu_obj, 
                                     features = feature, cell.types = celltypes, meta.groups = c(groups, 
                                                                                                 splitby))
    allgroups <- c(groups, splitby)
    for (i in 1:length(allgroups)) {
      if (is.null(levels(seu_obj@meta.data[, allgroups[i]]))) {
        seu_obj@meta.data[, allgroups[i]] <- factor(seu_obj@meta.data[, 
                                                                      allgroups[i]], levels = names(table(seu_obj@meta.data[, 
                                                                                                                            allgroups[i]])))
      }
      group_level <- levels(seu_obj@meta.data[, allgroups[i]])
      gene_count[, allgroups[i]] <- factor(gene_count[, 
                                                      allgroups[i]], levels = group_level)
    }
    gene_count$celltype <- factor(gene_count$celltype, levels = celltypes)
    all_levels <- list()
    for (i in 1:length(groups)) {
      if (is.null(levels(seu_obj@meta.data[, groups[i]]))) {
        seu_obj@meta.data[, groups[i]] <- factor(seu_obj@meta.data[, 
                                                                   groups[i]], levels = names(table(seu_obj@meta.data[, 
                                                                                                                      groups[i]])))
      }
      group_level <- levels(seu_obj@meta.data[, groups[i]])
      all_levels[[i]] <- group_level
    }
    all_levels <- as.character(unlist(all_levels))
    data_plot <- list()
    for (i in 1:length(groups)) {
      count_df <- gene_count
      count_df$new_group <- paste(gene_count[, groups[i]], 
                                  gene_count[, "celltype"], sep = "___")
      exp_df <- aggregate(. ~ new_group, data = count_df[, 
                                                         c("new_group", feature)], FUN = function(x) {
                                                           mean(expm1(x))
                                                         })
      pct_df <- aggregate(. ~ new_group, data = count_df[, 
                                                         c("new_group", feature)], FUN = function(x) {
                                                           length(x[x > 0])/length(x)
                                                         })
      colnames(exp_df)[2] <- "avg.exp"
      colnames(pct_df)[2] <- "pct.exp"
      df1 <- merge(exp_df, pct_df, by = "new_group")
      df1$groupID <- groups[i]
      data_plot[[i]] <- df1
    }
    data_plot <- do.call("rbind", data_plot)
    data_plot$groups <- as.character(lapply(X = strsplit(data_plot$new_group, 
                                                         split = "___"), FUN = function(x) {
                                                           x[[1]]
                                                         }))
    data_plot$celltype <- as.character(lapply(X = strsplit(data_plot$new_group, 
                                                           split = "___"), FUN = function(x) {
                                                             x[[2]]
                                                           }))
    data_plot$groups <- factor(data_plot$groups, levels = all_levels)
    data_plot$celltype <- factor(data_plot$celltype, levels = rev(celltypes))
    data_plot$groupID <- factor(data_plot$groupID, levels = groups)
    data_plot$pct.exp <- round(100 * data_plot$pct.exp, 
                               2)
    data_plot$avg.exp <- scale(data_plot$avg.exp)
    if (is.null(splitby)) {
      p <- ggplot(data_plot, aes(y = celltype, x = groups)) + 
        geom_tile(fill = "white", color = "white") + 
        geom_point(aes(colour = avg.exp, size = pct.exp)) + 
        scale_color_gradientn(colours = color.palette) + 
        theme(panel.background = element_rect(fill = "white", 
                                              colour = "black"), axis.text.x = element_text(angle = 45, 
                                                                                            hjust = 1, size = font.size), plot.title = element_text(size = (font.size + 
                                                                                                                                                              2), hjust = 0.5, face = "bold"), axis.text = element_text(size = font.size), 
              legend.text = element_text(size = (font.size - 
                                                   2)), legend.title = element_text(size = (font.size)), 
              strip.text = element_text(size = font.size), 
              legend.position = "right") + ylab("") + xlab("") + 
        ggtitle(feature) + facet_wrap(~groupID, scales = "free_x")
      if (do.scale) {
        p = p + scale_size(range = c(0, 10))
      }
      else {
        if (max(data_plot$pct.exp) >= 20) {
          p = p + scale_size(range = c(0, 10))
        }
        else {
          p = p + scale.func(range = c(0, 10), limits = c(0, 
                                                          20))
        }
      }
      g <- change_strip_background(p, type = "top", strip.color = strip.color)
      print(grid::grid.draw(g))
    }
    else {
      df2 <- reshape2::melt(gene_count[, c(groups, splitby)], 
                            measure.vars = groups)
      df2 <- df2[!duplicated(df2$value), ]
      colnames(df2)[colnames(df2) == "value"] <- "groups"
      data_plot2 <- list()
      for (i in 1:length(groups)) {
        df3 <- data_plot[data_plot$groupID == groups[i], 
        ]
        df4 <- df2[df2$variable == groups[i], c("groups", 
                                                splitby[i])]
        colnames(df4)[2] <- "split"
        df5 <- merge(df3, df4, by = "groups")
        data_plot2[[i]] <- df5
      }
      data_plot2 <- do.call("rbind", data_plot2)
      fill_x1 <- grDevices::rainbow(length(groups), alpha = 0.5)
      fill_x2 <- list()
      for (i in 1:length(splitby)) {
        n_col <- unique(gene_count[, splitby[i]])
        fill_x2[[i]] <- (scales::hue_pal(l = 90))(length(n_col))
      }
      fill_x2 <- as.character(unlist(fill_x2))
      fill_x <- c(fill_x1, fill_x2)
      p <- ggplot(data_plot2, aes(y = celltype, x = groups)) + 
        geom_tile(fill = "white", color = "white") + 
        geom_point(aes(colour = avg.exp, size = pct.exp)) + 
        scale_color_gradientn(colours = color.palette) + 
        theme(panel.background = element_rect(fill = "white", 
                                              colour = "black"), axis.text.x = element_text(angle = 45, 
                                                                                            hjust = 1, size = font.size), plot.title = element_text(size = (font.size + 
                                                                                                                                                              2), hjust = 0.5, face = "bold"), axis.text = element_text(size = font.size), 
              legend.text = element_text(size = (font.size - 
                                                   2)), legend.title = element_text(size = (font.size)), 
              strip.text = element_text(size = font.size), 
              legend.position = "right") + ylab("") + xlab("") + 
        ggtitle(feature) + facet_nested(~groupID + split, 
                                        scales = "free_x", strip = strip_nested(background_x = elem_list_rect(fill = fill_x)))
      if (do.scale) {
        p = p + scale_size(range = c(0, 10))
      }
      else {
        if (max(data_plot$pct.exp) >= 20) {
          p = p + scale_size(range = c(0, 10))
        }
        else {
          p = p + scale.func(range = c(0, 10), limits = c(0, 
                                                          20))
        }
      }
      p
    }
  }
}
com.dot.new <- com.dot.new <- function (seu_obj, features, celltypes = NULL, groups, color.palette = NULL, 
                                        strip.color = NULL) 
{
  pb <- progress_bar$new(format = "  Ploting [:bar] :percent eta: :eta", 
                         clear = FALSE, total = length(features), width = 100)
  plot_list <- list()
  for (i in 1:length(features)) {
    pp <- invisible(complex_dotplot_single(seu_obj = seu_obj, 
                                           feature = features[i], groups = groups, celltypes = celltypes))
    pp <- pp$data
    pp$gene <- features[i]
    plot_list[[i]] <- pp
    pb$tick()
    Sys.sleep(1/length(features))
  }
  all_data <- do.call("rbind", plot_list)
  all_data$gene <- factor(all_data$gene, levels = rev(features))
  all_data$celltype <- factor(all_data$celltype, levels = levels(seu_obj))
  if (is.null(color.palette)) {
    color.palette <- colorRampPalette(c("grey80", "lemonchiffon1",
                                        "indianred1", "darkred"))(255)
  }
  p <- invisible(ggplot(all_data, aes(x = groups, y = gene)) + 
                   geom_tile(fill = "white", color = "white") + 
                   geom_point(aes(colour = avg.exp, size = pct.exp), alpha = 0.9) + 
                   scale_color_gradientn(colours = color.palette) + 
                   scale_size(range = c(0, 5)) +
                   theme(
                     panel.background = element_rect(fill = "white", colour = "black"), 
                     axis.text.x = element_text(angle = 45,hjust = 1),
                     axis.text.y = element_text(face = 'italic'),
                     plot.title = element_text(size = 10, hjust = 0.5,face = "bold"), 
                     axis.text = element_text(size = 12), 
                     axis.title = element_text(size = 8), 
                     legend.text = element_text(size = 8), 
                     legend.title = element_text(size = 12),
                     legend.position = "right", 
                     strip.text = element_text(size = 8, colour = "black",face = "bold")) + 
                   ylab("") + xlab("") + ggtitle("") + 
                   facet_wrap(~celltype, ncol = length(levels(seu_obj))))
  g <- change_strip_background(p, type = "top", strip.color = strip.color)
  print(grid.draw(g))
}
extract_gene_count <- function (seu_obj, features, cell.types = NULL, data.type = "data", 
                                meta.groups = NULL) 
{
  if (is.null(cell.types)) {
    cell.types = levels(seu_obj)
  }
  seu_obj@meta.data$celltype <- as.character(seu_obj@active.ident)
  if (is.null(meta.groups)) {
    meta.groups = colnames(seu_obj@meta.data)
  }
  if (!is.null(cell.types)) {
    new_seu <- subset(seu_obj, idents = cell.types)
  }
  feature_count <- Seurat::FetchData(new_seu, slot = data.type, 
                                     vars = c(features, meta.groups, "celltype"))
  umap_data <- data.frame(new_seu[["umap"]]@cell.embeddings)
  feature_count$UMAP1 <- umap_data$UMAP_1
  feature_count$UMAP2 <- umap_data$UMAP_2
  feature_count
}
change_strip_background <- function (ggplt_obj, type = "top", strip.color = NULL) 
{
  g <- ggplot_gtable(ggplot_build(ggplt_obj))
  if (type == "top") {
    strip_both <- which(grepl("strip-t", g$layout$name))
    fills <- strip.color
    if (is.null(fills)) {
      fills <- (scales::hue_pal(l = 90))(length(strip_both))
    }
  }
  else if (type == "right") {
    strip_both <- which(grepl("strip-r", g$layout$name))
    fills <- strip.color
    if (is.null(fills)) {
      fills <- (scales::hue_pal(l = 90))(length(strip_both))
    }
  }
  else {
    strip_t <- which(grepl("strip-t", g$layout$name))
    strip_r <- which(grepl("strip-r", g$layout$name))
    strip_both <- c(strip_t, strip_r)
    fills <- strip.color
    if (is.null(fills)) {
      fills <- c((scales::hue_pal(l = 90))(length(strip_t)), 
                 (scales::hue_pal(l = 90))(length(strip_r)))
    }
  }
  k <- 1
  for (i in strip_both) {
    j <- which(grepl("rect", g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k + 1
  }
  g
}


## trans noise ----
getEuclideanDistance <- function(celltype, obj, assay, slot, ident1, ident2, ident3, group.by, lowcv = T){
  print(paste("Working on", celltype))
  library(hopach)
  tmp <- subset(obj, cells = WhichCells(obj, idents = celltype))
  
  counts <- GetAssayData(object = tmp, slot = slot, assay = assay)
  nonzero <- counts > 0
  keep_genes <- Matrix::rowSums(nonzero) > 0
  expr <- counts[keep_genes, ]
  
  ifelse(min(table(tmp@meta.data[[group.by]])) > 300,
         expr <- expr[,c(rownames(tmp@meta.data[tmp@meta.data[[group.by]] == ident1,])[sample(1:nrow(tmp@meta.data[tmp@meta.data[[group.by]] == ident1,]),300)],
                         rownames(tmp@meta.data[tmp@meta.data[[group.by]] == ident2,])[sample(1:nrow(tmp@meta.data[tmp@meta.data[[group.by]] == ident2,]),300)],
                         rownames(tmp@meta.data[tmp@meta.data[[group.by]] == ident3,])[sample(1:nrow(tmp@meta.data[tmp@meta.data[[group.by]] == ident3,]),300)] )
         ],
         expr <- expr)
  tmp <- subset(tmp,cells = colnames(expr))
  
  Down_Sample_Matrix <-function (expr_mat) {
    min_lib_size <- min(colSums(expr_mat))
    down_sample <- function(x) {
      prob <- min_lib_size/sum(x)
      return(unlist(lapply(x, function(y) {
        rbinom(1, y, prob)
      })))
    }
    down_sampled_mat <- apply(expr_mat, 2, down_sample)
    return(down_sampled_mat)
  }
  ds_expr <- Down_Sample_Matrix(expr)
  
  nsample <- min(table(tmp@meta.data[[group.by]])[c(ident1,ident2,ident3)])
  
  if(nsample < 10){
    print("Not enough cells")
    return(NULL)
  } 
  print(nsample)
  ident3_r <- sample(rownames(tmp@meta.data)[which(tmp@meta.data[[group.by]] == ident3)], nsample)
  ident2_r <- sample(rownames(tmp@meta.data)[which(tmp@meta.data[[group.by]] == ident2)], nsample)
  ident1_r <- sample(rownames(tmp@meta.data)[which(tmp@meta.data[[group.by]] == ident1)], nsample)
  ds_expr_r <- ds_expr[, c(ident1_r, ident2_r, ident3_r)]
  
  if(lowcv){
    getLowCVgenes <- function(matr){
      means <- Matrix::rowMeans(matr)
      bins <- quantile(means, c(seq(from = 0, to = 1, length = 11)))
      mean_bin <- unlist(lapply(means, function(x) min(which(bins >= x))))
      asplit <- split(names(means), mean_bin)
      genes <- unique(unlist(lapply(asplit[setdiff(names(asplit), c("1", "11"))], function(x){
        coef_var <- apply(matr, 1, function(x) sd(x)/mean(x))
        bottom10percent <- names(head(sort(coef_var), round(10*length(coef_var))))
      })))
      genes
    }
    genes <- getLowCVgenes(ds_expr_r)
  }
  else{
    genes <- rownames(ds_expr_r)
  }
  
  calcEuclDist <- function(matr, ident1, ident2, ident3){
    tmp <- data.matrix(sqrt(matr[genes, ident1]))
    mean <- rowMeans(sqrt(matr[genes, ident1]))
    d_ident1 <- distancevector(t(tmp), mean , d="euclid")
    names(d_ident1) <- ident1
    
    tmp <- data.matrix(sqrt(matr[genes, ident2]))
    mean <- rowMeans(sqrt(matr[genes, ident2]))
    d_ident2 <- distancevector(t(tmp), mean , d="euclid")
    names(d_ident2) <- ident2
    
    tmp <- data.matrix(sqrt(matr[genes, ident3]))
    mean <- rowMeans(sqrt(matr[genes, ident3]))
    d_ident3 <- distancevector(t(tmp), mean , d="euclid")
    names(d_ident3) <- ident3
    
    list(ident1 = d_ident1, ident2 = d_ident2, ident3 = d_ident3)
  }
  ds <- calcEuclDist(matr = ds_expr_r, ident2 = ident2_r, ident1 = ident1_r, ident3 = ident3_r)
  ds
}

