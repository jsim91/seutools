cellchat_netAnalysis_signalingRole_heatmap <- function(object, signaling = NULL, pattern = c("outgoing", "incoming", "all"),
                                                       slot.name = "netP", color.use = NULL, color.heatmap = "BuGn", title = NULL,
                                                       font.size.expansion = 1, cluster.rows = FALSE, cluster.cols = FALSE) {
  # testing
  # object = readRDS(file = "J:/10x/TB_sc/scbp2/final_annotations_outs/cellchat/GRV_p_cluster_less4_more3/cellchat_p_GRV.rds")
  # signaling = NULL
  # pattern = c("outgoing", "incoming", "all")[1]
  # slot.name = c("net","netP")[2]
  # color.use = NULL
  # color.heatmap = "BuGn"
  # title = NULL
  # font.size.expansion = 1
  # cluster.rows = FALSE
  # cluster.cols = FALSE

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
                                                       consider_cell_types = "all", row_plotting_threshold = 0.2,
                                                       measure = c("outdeg", "indeg", "flowbet", "info"),
                                                       measure.name = c("Sender", "Receiver", "Mediator", "Influencer"),
                                                       color.heatmap = "BuGn", font.size.expansion = 1,
                                                       cluster.rows = FALSE, cluster.cols = FALSE, heatmap_title = NULL)
{
  require(grid)
  require(gridExtra)
  # testing
  # object = object.list[[1]]
  # signaling = netp[1:2]
  # slot.name = "netP"
  # row_plotting_threshold = 0.2
  # heatmap_title = names(object.list)[1]
  # color.use = NULL
  # measure = c("outdeg", "indeg", "flowbet", "info")
  # measure.name = c("Sender", "Receiver", "Mediator", "Influencer")
  # color.heatmap = "BuGn"
  # font.size.expansion = 1
  # cluster.rows = FALSE
  # cluster.cols = FALSE
  # consider_cell_types <- "all"

  if(all(consider_cell_types[1]=="all",length(consider_cell_types)==1)) {
    consider_cell_types <- unique(names(object.list[[1]]@netP[["centr"]][[1]]$outdeg_unweighted))
  } else {
    consider_cell_types <- unique(consider_cell_types)
  }
  if (is.null(color.use)) {
    color.use <- scPalette(length(consider_cell_types))
  }
  cell.cols.assigned <- setNames(color.use, consider_cell_types)
  if (length(slot(object, slot.name)$centr) == 0) {
    stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores! ")
  }
  centr <- slot(object, slot.name)$centr[signaling]
  plot_list <- vector("list", length(centr))
  for (i in 1:length(centr)) {
    if(is.null(centr[[i]])) {
      plot_list[[i]] <- ggplot(data = cars, mapping = aes(x = speed, y = dist)) +
        geom_point(alpha = 0) +
        theme_void() +
        ggtitle(paste0(names(centr[i]), " signaling network")) +
        theme(plot.title = element_text(size = 12*font.size.expansion, hjust = 0.5))
      next
    } else {
      centr0 <- centr[[i]]
    }
    mat <- matrix(unlist(centr0), ncol = length(centr0), byrow = FALSE)
    mat <- t(mat)
    rownames(mat) <- names(centr0)
    colnames(mat) <- names(centr0$outdeg)
    mat <- mat[,which(colnames(mat) %in% consider_cell_types)]
    if (!is.null(measure)) {
      mat <- mat[measure, ]
      if (!is.null(measure.name)) {
        rownames(mat) <- measure.name
      }
    }
    mat_dim2_before <- ncol(mat)
    mat <- sweep(mat, 1L, apply(mat, 1, max), "/", check.margin = FALSE)
    mat[is.na(mat)] <- 0
    # if (is.null(color.use)) {
    #   color.use <- scPalette(length(colnames(mat)))
    # }
    if (!is.null(row_plotting_threshold)) {
      which_rm <- which(apply(mat, 2, max, na.rm = TRUE)<row_plotting_threshold)
      if (length(which_rm)!=0) {
        mat <- mat[,-which_rm]
      }
    }
    mat_dim2_after <- ncol(mat)
    if(is.null(mat_dim2_after)) {
      plot_list[[i]] <- ggplot() + theme_void()
    } else {
      num_lost_rows <- mat_dim2_before - mat_dim2_after
      ratio_lost <- num_lost_rows/mat_dim2_before

      color.heatmap.use = (grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, name = color.heatmap))))(100)
      # force color.heatmap.use limits c(0,1); even if min value in heatmap is >0 (because of thresholding)
      df <- data.frame(group = colnames(mat))
      rownames(df) <- colnames(mat)
      # cell.cols.assigned <- setNames(color.use, unique(as.character(df$group)))
      cell.cols.assigned <- cell.cols.assigned[which(names(cell.cols.assigned) %in% colnames(mat))]
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
                    column_title = ifelse(!is.null(heatmap_title), paste0(names(centr[i]), " signaling network (",heatmap_title,")"),
                                          paste0(names(centr[i]), " signaling network")),
                    column_title_gp = gpar(fontsize = 12*font.size.expansion),
                    rect_gp = gpar(lwd = 0.5, col = "black"), border = "black",
                    column_names_rot = 0, heatmap_legend_param = list(title = "Importance",
                                                                      title_gp = gpar(fontsize = 10*font.size.expansion, fontface = "plain"),
                                                                      title_position = "leftcenter-rot", border = "black",
                                                                      # at = c(round(min(mat, na.rm = T), digits = 1),
                                                                      #        round(max(mat, na.rm = T), digits = 1)),
                                                                      at = c(0, 1),
                                                                      legend_height = unit(24*font.size.expansion, "mm"),
                                                                      labels_gp = gpar(fontsize = 10*font.size.expansion),
                                                                      grid_width = unit(2.5*font.size.expansion, "mm")))
      tmp_draw <- ComplexHeatmap::draw(ht1)
      tmp_grob <- grid::grid.grabExpr(draw(tmp_draw))
      tmp_plt <- ggpubr::ggarrange(plotlist = tmp_grob)
      void_plt <- ggplot() + theme_void()
      blank_grob <- ggplotGrob(void_plt)
      tmp_out <- gridExtra::grid.arrange(tmp_grob, blank_grob,
                                         heights = c(0.13 + (1-(num_lost_rows/mat_dim2_before)), (num_lost_rows/mat_dim2_before)))
      plot_list[[i]] <- ggpubr::ggarrange(tmp_out)
    }
  }
  return(plot_list)
}


cellchat_compare_network_genes <- function()
{
  # this function should take (1) pathway(s) or (2) detect pathways that are different between two groups and plot representation(s) of these genes such as
  # violin plots, correlation plots of L vs R, that shows differences in expression/expression profiles for the genes identified by cellchat
  return()
}


cellchat_netAnalysis_signalingRole_merged_heatmap <- function(cellchat_object,
                                                              slot.name = "netP",
                                                              color.use = NULL,
                                                              name_vector_key = c("progr_more_3" = "progressor >12mo", "progr_less_4" = "progressor <12mo"),
                                                              font.size.expansion = 1) {
  require(CellChat)

  # testing
  # cellchat_object <- readRDS(file = "J:/10x/TB_sc/scbp2/final_annotations_outs/cellchat/Media_cl_cluster_controller11_less5/cellchat_cl_Media.rds")
  # slot.name = "net"
  # color.use = NULL
  # font.size.expansion = 0.5
  # name_vector_key = c("controller_11" = "controller", "progr_less_5" = "progressor <12mo")

  if(is.null(color.use)) {
    grad2 <- c("#ea2e2e","#f537f4","#4542ff")
    grad2[2] <- "#ffffff"
  } else {
    grad2 <- color.use
  }

  cc_1 <- slot(cellchat_object,slot.name)[[1]]$centr
  cc_2 <- slot(cellchat_object,slot.name)[[2]]$centr

  name1 <- names(cellchat_object@var.features)[1]
  name2 <- names(cellchat_object@var.features)[2]
  if(!is.null(name_vector_key)) {
    name1_adj <- as.character(name_vector_key[name1])
    name2_adj <- as.character(name_vector_key[name2])
    if(sum(c(is.na(name1_adj),is.na(name2_adj)))!=0) {
      stop("error in name_vector_key, no names found; check formatting")
    }
  } else {
    name1_adj <- name1
    name2_adj <- name2
  }

  if(slot.name=="netP") {
    path1 <- cellchat_object@netP[[1]]$pathways
    path2 <- cellchat_object@netP[[2]]$pathways
    pathunion <- union(path1,path2)
  } else if(slot.name=="net") {
    lr1 <- cellchat_object@net[[1]]$LRs
    lr2 <- cellchat_object@net[[2]]$LRs
    pathunion <- union(lr1,lr2)
  } else {
    stop("use either 'netP' or 'net' for 'slot.name'")
  }
  signaling <- pathunion

  get_outdeg <- function(obj) {
    return(obj$outdeg)
  }
  get_indeg <- function(obj) {
    return(obj$indeg)
  }

  outdeg_id1 <- t(sapply(X = cc_1, FUN = get_outdeg))
  indeg_id1 <- t(sapply(X = cc_1, FUN = get_indeg))
  capture_cname <- colnames(indeg_id1)
  outdeg_id1 <- outdeg_id1[which(row.names(outdeg_id1) %in% pathunion),]
  indeg_id1 <- indeg_id1[which(row.names(indeg_id1) %in% pathunion),]
  if(mean(colnames(outdeg_id1)==colnames(indeg_id1))!=1) {
    stop("column (celltype) order mismatch")
  }
  if(mean(row.names(outdeg_id1)==row.names(indeg_id1))!=1) {
    stop("column (celltype) order mismatch")
  }

  outdeg_id2 <- t(sapply(X = cc_2, FUN = get_outdeg))
  indeg_id2 <- t(sapply(X = cc_2, FUN = get_indeg))
  outdeg_id2 <- outdeg_id2[which(row.names(outdeg_id2) %in% pathunion),]
  indeg_id2 <- indeg_id2[which(row.names(indeg_id2) %in% pathunion),]
  if(mean(colnames(outdeg_id2)==colnames(indeg_id2))!=1) {
    stop("column (celltype) order mismatch")
  }
  if(mean(row.names(outdeg_id2)==row.names(indeg_id2))!=1) {
    stop("column (celltype) order mismatch")
  }

  unique_paths <- unique(c(row.names(outdeg_id1),row.names(outdeg_id2)))

  pathlist <- vector("list", length = length(unique_paths)); names(pathlist) <- unique_paths;
  pathlist <- lapply(pathlist, function(arg1, reptimes) return(rep(0,reptimes)), reptimes = ncol(outdeg_id1))

  outdeg_id1 <- split(x = outdeg_id1, f = row.names(outdeg_id1));
  outdeg_id1 <- append(outdeg_id1, pathlist[which(!names(pathlist) %in% names(outdeg_id1))])
  indeg_id1 <- split(x = indeg_id1, f = row.names(indeg_id1));
  indeg_id1 <- append(indeg_id1, pathlist[which(!names(pathlist) %in% names(indeg_id1))])

  outdeg_id2 <- split(x = outdeg_id2, f = row.names(outdeg_id2));
  outdeg_id2 <- append(outdeg_id2, pathlist[which(!names(pathlist) %in% names(outdeg_id2))])
  indeg_id2 <- split(x = indeg_id2, f = row.names(indeg_id2));
  indeg_id2 <- append(indeg_id2, pathlist[which(!names(pathlist) %in% names(indeg_id2))])

  outdeg_id1 <- outdeg_id1[unique_paths]; indeg_id1 <- indeg_id1[unique_paths]
  outdeg_id2 <- outdeg_id2[unique_paths]; indeg_id2 <- indeg_id2[unique_paths]

  outdeg_id1 <- lapply(X = outdeg_id1, function(x) return(data.frame('v1'=x)))
  indeg_id1 <- lapply(X = indeg_id1, function(x) return(data.frame('v1'=x)))

  outdeg_id2 <- lapply(X = outdeg_id2, function(x) return(data.frame('v1'=x)))
  indeg_id2 <- lapply(X = indeg_id2, function(x) return(data.frame('v1'=x)))

  scale_paths <- function(arg1) {
    if(sum(arg1)!=0) {
      arg2 <- scales::rescale(x = arg1[,1], to = c(0,1))
      df <- data.frame('v1' = arg2)
    } else {
      arg2 <- arg1
      df <- data.frame('v1' = arg2)
    }
    return(df)
  }

  outdeg_id1 <- lapply(X = outdeg_id1, FUN = scale_paths); outdeg_cdf <- do.call(cbind, outdeg_id1);
  colnames(outdeg_cdf) <- names(outdeg_id1); row.names(outdeg_cdf) <- colnames(outdeg_id1)
  indeg_id1 <- lapply(X = indeg_id1, FUN = scale_paths); indeg_cdf <- do.call(cbind, indeg_id1);
  colnames(indeg_cdf) <- names(indeg_id1); row.names(indeg_cdf) <- colnames(indeg_id1)

  outdeg_id2 <- lapply(X = outdeg_id2, FUN = scale_paths); outdeg_pdf <- do.call(cbind, outdeg_id2);
  colnames(outdeg_pdf) <- names(outdeg_id1); row.names(outdeg_pdf) <- colnames(outdeg_id1)
  indeg_id2 <- lapply(X = indeg_id2, FUN = scale_paths); indeg_pdf <- do.call(cbind, indeg_id2);
  colnames(indeg_pdf) <- names(indeg_id1); row.names(indeg_pdf) <- colnames(indeg_id1)

  if(mean(row.names(outdeg_pdf)==row.names(outdeg_cdf))!=1) {
    stop("mismatch outgoing rows")
  }
  if(mean(colnames(outdeg_pdf)==colnames(outdeg_cdf))!=1) {
    stop("mismatch outgoing columns")
  }

  if(mean(row.names(indeg_pdf)==row.names(indeg_cdf))!=1) {
    stop("mismatch incoming rows")
  }
  if(mean(colnames(indeg_pdf)==colnames(indeg_cdf))!=1) {
    stop("mismatch incoming columns")
  }

  outdeg_diff <- outdeg_pdf - outdeg_cdf
  indeg_diff <- indeg_pdf - indeg_cdf

  min_out <- min(outdeg_diff); max_out <- max(outdeg_diff)
  colmap_out <- colorRamp2(breaks = c(min_out, 0, max_out), colors = rev(grad2), space = "LAB")
  outdeg_diff_hm <- t(outdeg_diff)
  colnames(outdeg_diff_hm) <- capture_cname

  min_in <- min(indeg_diff); max_in <- max(indeg_diff)
  colmap_in <- colorRamp2(breaks = c(min_in, 0, max_in), colors = rev(grad2), space = "LAB")
  indeg_diff_hm <- t(indeg_diff)
  colnames(indeg_diff_hm) <- capture_cname

  sum_mat <- outdeg_diff_hm + indeg_diff_hm
  mean_rm0 <- function(arg1){
    drop0 <- arg1[union(which(arg1!=0),which(is.na(arg1)))]
    if(length(drop0)==0) {
      return(0)
    } else {
      return(mean(drop0))
    }
  }
  # sum_rmean <- rowMeans(x = sum_mat, na.rm = TRUE)
  sum_rmean <- apply(X = sum_mat, MARGIN = 1, FUN = mean_rm0)
  reorder <- order(sum_rmean)

  hm_diff_out <- Heatmap(matrix = outdeg_diff_hm[reorder,], col = colmap_out, name = 'hm_out', na_col = "black", cluster_columns = FALSE,
                         cluster_rows = FALSE, row_names_side = "left",
                         # heatmap_legend_param=list(legend_height=unit(3,"cm"), border = "black", labels = c("controller", "equal", "progressor"), at = c(min_out, 0, max_out),
                         #                           grid_width=unit(0.6,"cm"),title_position="topleft", legend_direction = "vertical",
                         #                           labels_gp=gpar(fontsize=11),title_gp=gpar(alpha = 0, fontsize = 0.1)),
                         row_names_gp=gpar(fontsize=13*font.size.expansion,fontface="bold"), column_title = "Outgoing Signal",
                         column_title_gp = gpar(fontface = "bold", cex = 2*font.size.expansion),
                         column_names_gp=gpar(fontsize=14*font.size.expansion,fontface="bold"), show_heatmap_legend = FALSE,
                         row_gap=unit(1,"mm"),column_gap=unit(1,"mm"),row_dend_gp=gpar(lwd=1.2),row_dend_width=unit(1,"cm"),
                         column_dend_gp = gpar(lwd=1.2), column_dend_height = unit(1,"cm"), rect_gp = gpar(lwd = 0.5, col = "black"), border = "black")

  hm_diff_in <- Heatmap(matrix = indeg_diff_hm[reorder,], col = colmap_in, name = 'hm_in', na_col = "black", cluster_columns = FALSE,
                        cluster_rows = FALSE, row_names_side = "left",
                        # heatmap_legend_param=list(legend_height=unit(3,"cm"), border = "black", labels = c("controller", "equal", "progressor"), at = c(min_in, 0, max_in),
                        #                           grid_width=unit(0.6,"cm"),title_position="topleft", legend_direction = "vertical",
                        #                           labels_gp=gpar(fontsize=11),title_gp=gpar(alpha = 0, fontsize = 0.1)),
                        row_names_gp=gpar(fontsize=13*font.size.expansion,fontface="bold"), column_title = "Incoming Signal",
                        column_title_gp = gpar(fontface = "bold", cex = 2*font.size.expansion),
                        column_names_gp=gpar(fontsize=14*font.size.expansion,fontface="bold"), show_heatmap_legend = FALSE,
                        row_gap=unit(1,"mm"),column_gap=unit(1,"mm"),row_dend_gp=gpar(lwd=1.2),row_dend_width=unit(1,"cm"),
                        column_dend_gp = gpar(lwd=1.2), column_dend_height = unit(1,"cm"), rect_gp = gpar(lwd = 0.5, col = "black"), border = "black")

  min_out_diff <- min(outdeg_diff); max_out_diff <- max(outdeg_diff)
  colmap_out_diff <- colorRamp2(breaks = c(min_out_diff, 0, max_out_diff), colors = rev(grad2), space = "LAB")
  lgd_out_diff <- Legend(col_fun = colmap_out_diff,
                         # title = paste0("Relative Interaction Strength\n",name2_adj," - ",name1_adj),
                         title = "Relative Interaction Strength",
                         border = "black", direction = "horizontal",
                         at = c(min_out_diff, 0, max_out_diff), legend_width = unit(8, "cm"),
                         legend_height = unit(1.2, "cm"), title_gp = gpar(fontface = "bold", cex = 1.7*font.size.expansion),
                         title_position = "topcenter", labels = c(name1_adj,"equal",name2_adj),
                         labels_gp = gpar(fontface = "bold", cex = 1.7*font.size.expansion))

  grob_in_diff_hm = grid::grid.grabExpr(draw(hm_diff_in))
  grob_out_diff_hm = grid::grid.grabExpr(draw(hm_diff_out))
  arr_diff_hm <- ggpubr::ggarrange(plotlist = list(grob_out_diff_hm, grob_in_diff_hm), nrow = 1)
  grob_diff_leg = grid::grid.grabExpr(draw(lgd_out_diff))
  arr_diff_leg <- cowplot::plot_grid(grob_diff_leg, ncol = 1, nrow = 1)
  arr_diff_out <- ggpubr::ggarrange(plotlist = list(arr_diff_leg, arr_diff_hm), nrow = 2, ncol = 1, heights = c(0.1,0.9))
  return(arr_diff_out)
  # pdf(file = paste0("/media/MPEdge16/TB_sc/sc/py/scbp2/cellchat/",stim_condition,"_pc_cluster_variableN/cellchat_pc_",
  #                   stim_condition,"_diff_hm",ifelse(use_equal_sampling,"_equal_sampling",""),".pdf"), width = 18, height = 12)
  # arr_diff_out
  # dev.off()
}
