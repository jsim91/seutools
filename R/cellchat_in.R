cellchat_netAnalysis_signalingRole_heatmap <- function(object, signaling = NULL, pattern = c("outgoing", "incoming", "all"),
                                                       slot.name = "netP", color.use = NULL, color.heatmap = "BuGn", title = NULL,
                                                       font.size.expansion = 1, cluster.rows = FALSE, cluster.cols = FALSE) {
  pattern <- match.arg(pattern)
  if (length(slot(object, slot.name)$centr) == 0) {
    stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores! ")
  }
  centr <- slot(object, slot.name)$centr
  outgoing <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
  incoming <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
  dimnames(outgoing) <- list(levels(object@idents), names(centr))
  dimnames(incoming) <- dimnames(outgoing)
  for (i in 1:length(centr)) {
    outgoing[, i] <- centr[[i]]$outdeg
    incoming[, i] <- centr[[i]]$indeg
  }
  if (pattern == "outgoing") {
    mat <- t(outgoing)
    legend.name <- "Outgoing"
  }
  else if (pattern == "incoming") {
    mat <- t(incoming)
    legend.name <- "Incoming"
  } else if (pattern == "all") {
    mat <- t(outgoing + incoming)
    legend.name <- "Overall"
  }
  if (is.null(title)) {
    title <- paste0(legend.name, " signaling patterns")
  } else {
    title <- paste0(paste0(legend.name, " signaling patterns"),
                    " - ", title)
  }
  if (!is.null(signaling)) {
    mat1 <- mat[rownames(mat) %in% signaling, , drop = FALSE]
    mat <- matrix(0, nrow = length(signaling), ncol = ncol(mat))
    idx <- match(rownames(mat1), signaling)
    mat[idx[!is.na(idx)], ] <- mat1
    dimnames(mat) <- list(signaling, colnames(mat1))
  }
  mat.ori <- mat
  mat <- sweep(mat, 1L, apply(mat, 1, max), "/", check.margin = FALSE)
  mat[mat == 0] <- NA
  if (is.null(color.use)) {
    color.use <- scPalette(length(colnames(mat)))
  }
  color.heatmap.use = (grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, name = color.heatmap))))(length(unique(unlist(lapply(X = mat, FUN = function(x) x)))) + 1)
  df <- data.frame(group = colnames(mat))
  rownames(df) <- colnames(mat)
  names(color.use) <- colnames(mat)
  col_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use),
                                      which = "column", show_legend = FALSE, show_annotation_name = FALSE,
                                      simple_anno_size = grid::unit(0.2, "cm"))
  ha2 = HeatmapAnnotation(Strength = anno_barplot(colSums(mat.ori), border = FALSE, gp = gpar(fill = color.use, col = color.use)),
                          show_annotation_name = FALSE)
  pSum <- rowSums(mat.ori)
  pSum.original <- pSum
  pSum <- -1/log(pSum)
  pSum[is.na(pSum)] <- 0
  idx1 <- which(is.infinite(pSum) | pSum < 0)
  if (length(idx1) > 0) {
    values.assign <- seq(max(pSum) * 1.1, max(pSum) * 1.5, length.out = length(idx1))
    position <- sort(pSum.original[idx1], index.return = TRUE)$ix
    pSum[idx1] <- values.assign[match(1:length(idx1), position)]
  }
  ha1 = rowAnnotation(Strength = anno_barplot(pSum, border = FALSE), show_annotation_name = FALSE)
  if (min(mat, na.rm = T) == max(mat, na.rm = T)) {
    legend.break <- max(mat, na.rm = T)
  } else {
    legend.break <- c(round(min(mat, na.rm = T), digits = 1), round(max(mat, na.rm = T), digits = 1))
  }
  ht1 = Heatmap(mat, col = color.heatmap.use, na_col = "white",
                name = "Relative strength", bottom_annotation = col_annotation,
                top_annotation = ha2, right_annotation = ha1, cluster_rows = cluster.rows,
                cluster_columns = cluster.rows, row_names_side = "left",
                row_names_rot = 0, row_names_gp = gpar(fontsize = 8*font.size.expansion),
                column_names_gp = gpar(fontsize = 8*font.size.expansion), column_title = title,
                column_title_gp = gpar(fontsize = 10*font.size.expansion),
                rect_gp = gpar(lwd = 0.2, col = "black"), border = "black",
                column_names_rot = 90, heatmap_legend_param = list(title_gp = gpar(fontsize = 8*font.size.expansion, fontface = "plain"),
                                                                   title_position = "leftcenter-rot", border = "black",
                                                                   at = legend.break, legend_height = unit(20, "mm"), grid_width = unit(2.5*font.size.expansion, "mm"),
                                                                   legend_height = unit(24*font.size.expansion, "mm"), labels_gp = gpar(fontsize = 10*font.size.expansion),
                                                                   labels_gp = gpar(fontsize = 8*font.size.expansion)))
  return(ht1)
}


cellchat_netAnalysis_signalingRole_network <- function(object, signaling, slot.name = "netP", color.use = NULL,
                                                       measure = c("outdeg", "indeg", "flowbet", "info"),
                                                       measure.name = c("Sender", "Receiver", "Mediator", "Influencer"),
                                                       color.heatmap = "BuGn", font.size.expansion = 1,
                                                       cluster.rows = FALSE, cluster.cols = FALSE, heatmap_title = NULL)
  # width = 6.5, height = 1.4, font.size.title = 10,
{
  # testing
  # object = cc_list[[i]]
  # signaling = names(path_maps)[j]
  # slot.name = "netP"
  # color.use = NULL
  # measure = c("outdeg", "indeg", "flowbet", "info")
  # measure.name = c("Sender", "Receiver", "Mediator", "Influencer")
  # color.heatmap = "BuGn"
  # font.size.expansion = 1.25
  # cluster.rows = FALSE
  # cluster.cols = FALSE
  # heatmap_title = paste0(names(cc_list)[i]," ",names(path_maps)[j]," signaling network")

  if (length(slot(object, slot.name)$centr) == 0) {
    stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores! ")
  }
  centr <- slot(object, slot.name)$centr[signaling]
  plot_list <- vector("list", length(centr))
  for (i in 1:length(centr)) {
    centr0 <- centr[[i]]
    mat <- matrix(unlist(centr0), ncol = length(centr0),
                  byrow = FALSE)
    mat <- t(mat)
    rownames(mat) <- names(centr0)
    colnames(mat) <- names(centr0$outdeg)
    if (!is.null(measure)) {
      mat <- mat[measure, ]
      if (!is.null(measure.name)) {
        rownames(mat) <- measure.name
      }
    }
    mat <- sweep(mat, 1L, apply(mat, 1, max), "/", check.margin = FALSE)
    if (is.null(color.use)) {
      color.use <- scPalette(length(colnames(mat)))
    }
    color.heatmap.use = (grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, name = color.heatmap))))(100)
    df <- data.frame(group = colnames(mat))
    rownames(df) <- colnames(mat)
    cell.cols.assigned <- setNames(color.use, unique(as.character(df$group)))
    # col_annotation <- HeatmapAnnotation(df = df, col = list(group = cell.cols.assigned),
    #                                     which = "column", show_legend = FALSE, show_annotation_name = FALSE,
    #                                     simple_anno_size = grid::unit(0.2, "cm"))
    row_annotation <- HeatmapAnnotation(df = df, col = list(group = cell.cols.assigned),
                                        which = "row", show_legend = FALSE, show_annotation_name = FALSE,
                                        simple_anno_size = grid::unit(0.3, "cm"))
    hm_mat <- t(mat)
    ht1 = Heatmap(matrix = hm_mat, col = color.heatmap.use, na_col = "white",
                  name = "Importance", left_annotation = row_annotation, column_names_centered = TRUE,
                  # bottom_annotation = col_annotation,
                  cluster_rows = cluster.rows, cluster_columns = cluster.cols,
                  row_names_side = "left", row_names_rot = 0, row_names_gp = gpar(fontsize = 10*font.size.expansion, fontface = "bold"),
                  column_names_gp = gpar(fontsize = 12*font.size.expansion),
                  # width = unit(width, "cm"), height = unit(height, "cm"),
                  column_title = ifelse(!is.null(heatmap_title), heatmap_title, paste0(names(centr[i]), " signaling network")),
                  column_title_gp = gpar(fontsize = 12*font.size.expansion),
                  rect_gp = gpar(lwd = 0.5, col = "black"), border = "black",
                  column_names_rot = 0, heatmap_legend_param = list(title = "Importance",
                                                                    title_gp = gpar(fontsize = 10*font.size.expansion, fontface = "plain"),
                                                                    title_position = "leftcenter-rot", border = "black",
                                                                    at = c(round(min(mat, na.rm = T), digits = 1),
                                                                           round(max(mat, na.rm = T), digits = 1)),
                                                                    legend_height = unit(24*font.size.expansion, "mm"), labels_gp = gpar(fontsize = 10*font.size.expansion),
                                                                    grid_width = unit(2.5*font.size.expansion, "mm")))
    tmp_draw <- ComplexHeatmap::draw(ht1)
    # plot_list[[i]] <- tmp_draw
    tmp_grob <- grid::grid.grabExpr(draw(tmp_draw))
    plot_list[[i]] <- tmp_grob
    # plot_list[[i]] <- ht1
    # plot_list[[i]] <- hm_mat
  }

  return(plot_list)
  # return(ggpubr::ggarrange(plotlist = plot_list, nrow = 1))
}


cellchat_compare_network_genes <- function()
{
  # this function should take (1) pathway(s) or (2) detect pathways that are different between two groups and plot representation(s) of these genes such as
  # violin plots, correlation plots of L vs R, that shows differences in expression/expression profiles for the genes identified by cellchat
  return()
}
