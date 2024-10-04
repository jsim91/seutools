heatmap_calculate <- function(seurat_obj, gene_set, set_name, clusters)
{
  require(ComplexHeatmap)
  require(FCSimple)
  require(ggplot2)

  # testing
  # seurat_obj = seu
  # gene_set = c("")
  # set_name =
  # clusters =


  ct = seurat_obj@assays[["RNA"]]@data

  gene_row <- which(row.names(ct) %in% gsub("_.+$","",gene_set))
  pared_data <- as.matrix(ct[gene_row,])
  cluster_nums <- as.character(sobj@meta.data$seurat_clusters)
  cluster_nums[which(!cluster_nums %in% clusters)] <- "other"
  uclus <- unique(cluster_nums[-which(cluster_nums=="other")]); uclus <- uclus[order(as.numeric(uclus))]
  uclus <- c(uclus,"other")

  fcsobj <- list(data = t(pared_data), source = 1:ncol(pared_data),
                 cluster_nums = list(clusters = as.numeric(as.character(sobj@meta.data$seurat_clusters))))
  fcsobj <- fcs_cluster_heatmap(fcs_join_obj = fcsobj, algorithm = "cluster_nums", include_parameters = "all", cluster_heatmap = FALSE)
  ggplot2::ggsave(filename = "test.pdf",
                  plot = grid::grid.grabExpr(draw(fcsobj[["cluster_nums_heatmap"]][["heatmap"]])),
                  device = "pdf", width = 12, height = 12, units = "in", dpi = 900)

  # grid::grid.grabExpr(draw(fcsobj[["cluster_nums_heatmap"]][["heatmap"]]))

  in_mat <- matrix(data = NA, nrow = nrow(pared_data), ncol = length(uclus))
  colnames(in_mat) <- uclus; row.names(in_mat) <- row.names(pared_data)

  for(i in 1:nrow(in_mat)) {
    gene_ct <- pared_data[which(row.names(pared_data)==row.names(in_mat)[i]),]
    names(gene_ct) <- cluster_nums
    for(j in 1:ncol(in_mat)) {
      in_mat[i,j] <- mean(gene_ct[which(names(gene_ct)==colnames(in_mat)[j])])
    }
  }

  hm_height <- (2.25) + nrow(in_mat)/7
  hm_width <- (3.25) + ncol(in_mat)/4

  out_hm <- Heatmap(matrix = in_mat)

  # ggsave(filename = "cytokine_signaling.pdf",
  #        plot = grid::grid.grabExpr(draw(out_hm)),
  #        device = "pdf", width = hm_width,
  #        height = hm_height,
  #        units = "in", dpi = 900)
  # return(grid::grid.grabExpr(draw(out_hm)))
}

seurat_tile_reduction <- function(seurat_object, condition_column, cluster_column, reduction = "umap",
                                  color_clusters = "all", label_clusters = "all",
                                  pt_alpha = 0.05, text_expansion = 1, pt_size = 1, color_seed = 123,
                                  # outline_method = c("nudge","fontsize"),
                                  postfix_title_string = NA, force_colors = FALSE,
                                  force_xlim = FALSE, force_ylim = FALSE, return_as_list = FALSE,
                                  plot_order = c(2,1,3), annotation_method = "repel", # c("repel","text","shadowtext","none")
                                  override_color_aes = NA, frameon = FALSE)
{
  require(ggplot2)
  require(ggpubr)
  require(shadowtext)
  require(ggrastr)
  require(ggrepel)

  # testing
  # seurat_object = seu_small
  # condition_column = "condition"
  # cluster_column = "cell_type"
  # reduction = "umap"
  # color_clusters = "all"
  # label_clusters = "all"
  # pt_alpha = 0.05
  # text_expansion = 1
  # pt_size = 1
  # color_seed = 123
  # outline_method = c("nudge","fontsize")
  # postfix_title_string = NA
  # force_xlim = FALSE
  # force_ylim = FALSE
  # return_as_list = FALSE
  # plot_order = c(2,1,3)
  # annotation_method = "repel" # c("repel","text","shadowtext","none")
  # override_color_aes = NA
  # frameon = FALSE


  coords <- seurat_object@reductions[[tolower(reduction)]]@cell.embeddings

  plot_data <- data.frame(redx = coords[,1], redy = coords[,2],
                          cluster = as.character(seurat_object@meta.data[,cluster_column]),
                          condition = seurat_object@meta.data[,condition_column])
  unique_clus <- unique(plot_data$cluster)
  set.seed(color_seed)
  plot_data$cluster <- factor(x = plot_data$cluster, levels = sample(unique_clus,length(unique_clus),replace=F))
  xrange <- range(plot_data$redx); yrange = range(plot_data$redy)

  uclus <- unique(plot_data$cluster); uclus <- uclus[order(uclus,decreasing=F)]
  clusx <- rep(NA,length(uclus)); names(clusx) <- uclus; clusy <- clusx
  for(i in 1:length(clusx)) {
    clusx[i] <- median(plot_data$redx[which(plot_data$cluster==names(clusx)[i])])
    clusy[i] <- median(plot_data$redy[which(plot_data$cluster==names(clusx)[i])])
  }
  if(label_clusters[1]!="all") {
    lab_ind <- which(names(clusx) %in% as.character(label_clusters))
    if(length(lab_ind)!=0) {
      clusx <- clusx[lab_ind]
      clusy <- clusy[lab_ind]
    } else {
      clusx <- c(); clusy <- c()
    }
  }

  spl_data <- split(x = plot_data, f = plot_data$condition)
  ds_to <- min(sapply(X = spl_data, FUN = function(x) return(nrow(x))))
  for(i in 1:length(spl_data)) {
    set.seed(123)
    spl_data[[i]] <- spl_data[[i]][sample(x = 1:nrow(spl_data[[i]]), size = ds_to, replace = FALSE),]
  }

  plot_red <- function(input, color_clus = color_clusters, xanno = clusx, yanno = clusy,
                       palpha = pt_alpha, texp = text_expansion, psize = pt_size,
                       plimx = xrange, plimy = yrange, amethod = annotation_method,
                       cseed = color_seed,#text_omethod=outline_method,
                       oca = override_color_aes, fcol = force_colors,
                       flimx = force_xlim, flimy = force_ylim, pts = postfix_title_string,
                       fo = frameon)
  {
    # testing
    # input = spl_data[[1]]
    # color_clus = color_clusters
    # xanno = clusx
    # yanno = clusy
    # palpha = pt_alpha
    # texp = text_expansion
    # psize = pt_size
    # plimx = xrange
    # plimy = yrange
    # amethod = annotation_method
    # cseed = color_seed
    # oca = override_color_aes
    # flimx = force_xlim
    # flimy = force_ylim
    # pts = postfix_title_string
    # fo = frameon

    gg_color_hue <- function(n) {
      hues = seq(15, 375, length = n + 1)
      hcl(h = hues, l = 65, c = 100)[1:n]
    }

    if(length(xanno)!=0) {
      text_add <- data.frame(UMAP1 = xanno, UMAP2 = yanno, cluster = names(xanno))
    } else {
      text_add <- as.data.frame(matrix(data=NA,nrow=0,ncol=2))
    }
    # if(length(color_clus)!=0) {
    if(color_clus[1]!="none") {
      if(color_clus[1]=="all") {
        subs_rows <- 1:nrow(text_add)
      } else {
        subs_rows <- which(text_add$cluster %in% as.character(color_clus))
      }
      if(length(subs_rows)!=0) {
        color_text_add <- text_add[subs_rows,]
        text_add <- text_add[subs_rows,]

        color_text_add$cluster <- factor(color_text_add$cluster)
        text_add$cluster <- factor(text_add$cluster)
        # set.seed(cseed)
        # color_text_add <- color_text_add[sample(1:nrow(color_text_add),nrow(color_text_add),replace=F),]
      } else {
        color_text_add <- matrix(data=NA,nrow=0,ncol=2)
      }
    } else {
      if(nrow(text_add)!=0) {
        text_add$cluster <- factor(text_add$cluster)
      }
      color_text_add <- matrix(data=NA,nrow=0,ncol=2)
    }

    input$cluster <- factor(input$cluster)
    if(color_clus[1]!="none") {
      if(color_clus[1]!="all") {
        foreground <- input[which(input$cluster %in% color_clus),]
        input <- input[which(!input$cluster %in% color_clus),]
      } else {
        foreground <- input
      }
    } else {
      plt <- ggplot(data = input, mapping = aes(x = redx, y = redy)) +
        geom_point(pch = 19, alpha = palpha, size = psize) + theme_void() +
        xlim(plimx) + ylim(plimy)
      if(!isFALSE(flimx[1])) {
        plt <- plt + xlim(flimx)
      }
      if(!isFALSE(flimy[1])) {
        plt <- plt + ylim(flimy)
      }
    }
    if(color_clus[1]!="none") {
      if(color_clus[1]!="all") {
        plt <- ggplot(data = input, mapping = aes(x = redx, y = redy)) +
          geom_point(pch = 19, alpha = palpha, size = psize) + theme_void() +
          xlim(plimx) + ylim(plimy)
        plt <- plt + geom_point(data = foreground,
                                mapping = aes(x = redx, y = redy, color = cluster),
                                alpha = palpha, size = psize) +
          guides(color = guide_legend(override.aes = list(size = 6, alpha = 1))) +
          theme(legend.position = "none",
                legend.title = element_blank(),
                legend.text = element_text(size = 16*texp))
        if(!isFALSE(flimx[1])) {
          plt <- plt + xlim(flimx)
        }
        if(!isFALSE(flimy[1])) {
          plt <- plt + ylim(flimy)
        }
      } else {
        # plt <- ggplot(data = foreground, mapping = aes(x = UMAP1, y = UMAP2)) +
        #   geom_point(pch = 19, alpha = palpha, size = psize) + theme_void() +
        #   xlim(plimx) + ylim(plimy)
        plt <- ggplot() + theme_void() + xlim(plimx) + ylim(plimy) +
          geom_point(data = foreground,
                     mapping = aes(x = redx, y = redy, color = cluster),
                     alpha = palpha, size = psize) +
          guides(color = guide_legend(override.aes = list(size = 6, alpha = 1))) +
          theme(legend.position = "none",
                legend.title = element_blank(),
                legend.text = element_text(size = 16*texp))
        if(!isFALSE(flimx[1])) {
          plt <- plt + xlim(flimx)
        }
        if(!isFALSE(flimy[1])) {
          plt <- plt + ylim(flimy)
        }
      }
      if(!is.na(oca[1])) {
        plt <- plt + scale_color_manual(values = oca)
      }
    }
    if(amethod[1]!="none") {
      if(nrow(text_add)!=0) {
        if(amethod[1]=="shadowtext") {
          plt <- plt + annotate("shadowtext", x = color_text_add$UMAP1, y = color_text_add$UMAP2,
                                label = color_text_add$cluster, size = 3*texp)
        } else if(amethod[1]=="repel") {
          plt <- plt + ggrepel::geom_text_repel(data = color_text_add, mapping = aes(x = UMAP1, y = UMAP2, color = cluster, label = cluster),
                                                size = 3*texp, bg.color = "black", bg.r = 0.05, seed = 123)
        } else if(amethod[1]=="text") {
          plt <- plt + geom_text(data = color_text_add, mapping = aes(x = UMAP1, y = UMAP2, color = cluster, label = cluster),
                                 size = 3*texp)
        }
      } else if(nrow(color_text_add)!=0) {
        if(amethod[1]=="shadowtext") {
          plt <- plt + annotate("shadowtext", x = color_text_add$UMAP1, y = color_text_add$UMAP2,
                                label = color_text_add$cluster, size = 3*texp)
        } else if(amethod[1]=="repel") {
          plt <- plt + ggrepel::geom_text_repel(data = color_text_add, mapping = aes(x = UMAP1, y = UMAP2, color = cluster, label = cluster),
                                                size = 3*texp, bg.color = "black", bg.r = 0.05, seed = 123)
        } #else if(amethod[1]=="text") {
        #   if(text_omethod[1]=="nudge") {
        #     plt <- plt + geom_text(data = color_text_add, mapping = aes(x = UMAP1, y = UMAP2, label = cluster),
        #                            size = anno_ts, color = "black", nudge_x = 0.005, nudge_y = 0.005) +
        #       plt <- plt + geom_text(data = color_text_add, mapping = aes(x = UMAP1, y = UMAP2, color = cluster, label = cluster),
        #                              size = anno_ts)
        #   } else if(text_omethod[1]=="fontsize") {
        #     plt <- plt + geom_text(data = color_text_add, mapping = aes(x = UMAP1, y = UMAP2, label = cluster),
        #                            size = anno_ts*1.05, color = "black") +
        #       plt <- plt + geom_text(data = color_text_add, mapping = aes(x = UMAP1, y = UMAP2, color = cluster, label = cluster),
        #                              size = anno_ts)
        #   }
        # }
      }
    } else {
      plt <- plt + guides(color = guide_legend(override.aes = list(alpha = 1, stroke = 0.1, size = 6))) +
        theme(legend.position = "bottom", legend.text = element_text(size = 16*texp), legend.title = element_blank())
    }
    plt <- plt + ggtitle(ifelse(!is.na(pts), paste0(input$condition[1]," - ",pts), input$condition[1])) +
      theme(axis.text = element_blank(),
            axis.title = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 26*texp))
    if(fo) {
      plt <- plt + theme_bw() + theme(axis.ticks = element_blank(),
                                      axis.title = element_blank(),
                                      axis.text = element_blank(),
                                      legend.title = element_blank(),
                                      panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank())
    }
    if(!isFALSE(fcol)) {
      plt <- plt +
        scale_color_manual(values = fcol)
    }
    return(plt)
  }

  out_plots <- lapply(X = spl_data, FUN = plot_red)

  if(return_as_list) {
    return(out_plots[plot_order])
  }

  if(!is.na(plot_order[1])) {
    arr_plot <- ggpubr::ggarrange(plotlist = out_plots[plot_order], nrow = 1)
  } else {
    arr_plot <- ggpubr::ggarrange(plotlist = out_plots, nrow = 1)
  }
  return(arr_plot)
}


seurat_footprint <- function(seurat_object, cluster_column, pid_column, condition_column,
                             media_condition, subtract_media, color_by_column, which_clusters = "all",
                             reduction = "umap", scale.factor = 1000, pca_fraction_variance = 0.95,
                             cluster_data = FALSE, umap_n_neighbors = 5, leiden_resolution = 0.2,
                             umap_min_dist = 0.3, report_values_as = "fraction",
                             feature_reduction_method = "pca", return_as_list = FALSE,
                             return_data = FALSE, manual_colors = NULL, umap_spread = 1,
                             label_points = TRUE) # prepare feature input with either 'pca' or 'boruta' as feature reduction method
{
  suppressPackageStartupMessages({
    require(ggplot2)
    require(ggpubr)
    require(Boruta)
    require(uwot)
    require(gtools)
    require(ggrepel)
    require(FCSimple)
    require(combinat)
  })

  # testing
  # seurat_object = seu
  # cluster_column = "cell_type"
  # pid_column = "pid"
  # condition_column = "condition"
  # media_condition = "media"
  # subtract_media = TRUE
  # color_by_column = "age_group"
  # scale.factor = 1000
  # pca_fraction_variance = 1
  # umap_n_neighbors = 5
  # cluster_data = FALSE
  # leiden_resolution = 0.5
  # umap_min_dist = 0.25
  # report_values_as = "frequency"
  # feature_reduction_method = "none"
  # which_clusters = selected_clusters
  # return_data = FALSE
  # return_as_list = FALSE
  # manual_colors = NULL
  # # #
  # reduction = "umap"

  metadata <- seurat_object@meta.data
  # if(which_clusters[1]!="all") {
  #   metadata <- metadata[which(metadata[,cluster_column] %in% which_clusters),]
  #   if(nrow(metadata)==0) {
  #     stop("no data left after filtering for 'which_clusters'")
  #   }
  # }
  metadata[,cluster_column] <- as.character(metadata[,cluster_column])
  metadata[,pid_column] <- as.character(metadata[,pid_column])
  small_meta <- metadata[!duplicated(metadata[,pid_column]),][,c(pid_column, color_by_column)]
  colnames(small_meta)[2] <- "meta_col"
  unique_condition <- unique(metadata[,condition_column])

  if(subtract_media) {
    upid_list <- vector("list", length = length(unique_condition))
    for(i in 1:length(unique_condition)) {
      upid_list[[i]] <- unique(metadata[which(metadata[,condition_column]==unique_condition[i]),][,pid_column])
    }
    upid_all <- unlist(upid_list)
    upid_table <- table(upid_all)
    keep_id <- names(upid_table)[which(upid_table==max(upid_table))]
    if(length(keep_id)<length(unique(metadata[,pid_column]))) {
      metadata <- metadata[which(metadata[,pid_column] %in% keep_id),]
    }
  }

  clus_list <- vector("list", length = length(unique_condition)); names(clus_list) <- unique_condition

  for(i in 1:length(unique_condition)) {
    tmp_meta <- metadata[metadata[,condition_column]==unique_condition[i],]
    if(report_values_as=="normalized counts") {
      freq_table <- table(tmp_meta[,pid_column], tmp_meta[,cluster_column])

      freq_matrix <- matrix(freq_table, nrow = length(unique(tmp_meta[,pid_column])), ncol = length(unique(tmp_meta[,cluster_column])))

      rownames(freq_matrix) <- sort(unique(tmp_meta[,pid_column]))
      colnames(freq_matrix) <- sort(unique(tmp_meta[,cluster_column]))
      freq_df <- as.data.frame.matrix(freq_matrix)
      row_sums <- rowSums(freq_matrix)
      normalized_matrix <- sweep(freq_matrix, 1, row_sums, FUN="/") * scale.factor

      # clus_list[[i]] <- normalized_matrix
    } else {
      upid <- unique(tmp_meta[,pid_column]); upid <- upid[order(upid)]; uclus <- unique(tmp_meta[,cluster_column]); uclus <- uclus[order(uclus)]
      freq_matrix <- matrix(data = NA, nrow = length(upid), ncol = length(uclus))
      row.names(freq_matrix) <- upid; colnames(freq_matrix) <- uclus
      if(report_values_as=="fraction") {
        for(k in 1:nrow(freq_matrix)) {
          tmpnum <- tmp_meta[,cluster_column][tmp_meta[,pid_column]==row.names(freq_matrix)[k]]
          for(j in 1:ncol(freq_matrix)) {
            freq_matrix[k,j] <- mean(tmpnum==colnames(freq_matrix)[j])
          }
        }
      } else if(report_values_as=="frequency") {
        for(k in 1:nrow(freq_matrix)) {
          tmpnum <- tmp_meta[,cluster_column][tmp_meta[,pid_column]==row.names(freq_matrix)[k]]
          for(j in 1:ncol(freq_matrix)) {
            freq_matrix[k,j] <- mean(tmpnum==colnames(freq_matrix)[j]) * 100
          }
        }
      }
    }
    clus_list[[i]] <- freq_matrix
  }

  rnames <- lapply(X = clus_list, FUN = row.names); rname_table <- table(unlist(rnames)); keepid <- names(which(rname_table==max(rname_table)))
  cnames <- lapply(X = clus_list, FUN = colnames); cname_table <- table(unlist(cnames)); keepcl <- names(which(cname_table==max(cname_table)))
  for(i in 1:length(clus_list)) {
    clus_list[[i]] <- clus_list[[i]][which(row.names(clus_list[[i]]) %in% keepid),]
    clus_list[[i]] <- clus_list[[i]][,which(colnames(clus_list[[i]]) %in% keepcl)]
  }

  for(i in 1:length(clus_list)) {
    if(any(i==1, length(clus_list)==1)) {
      next
    } else {
      if(mean(row.names(clus_list[[i]])==row.names(clus_list[[i-1]]))!=1) {
        stop("row mismatch")
      }
    }
  }

  if(subtract_media) {
    which_media <- which(names(clus_list)==media_condition)
    perms <- permutations(n = length(clus_list), r = 2, v = 1:length(clus_list), repeats.allowed = TRUE)
    drop_perms <- intersect(which(perms[,1]==which_media), which(perms[,2]!=which_media))
    if(length(drop_perms)!=0) {
      perms <- perms[-drop_perms,]
    }
    subtr_list <- vector("list", length = nrow(perms))
    for(i in 1:length(subtr_list)) {
      if(perms[i,1]==perms[i,2]) {
        subtr_list[[i]] <- clus_list[[perms[i,1]]]
        names(subtr_list)[i] <- names(clus_list)[perms[i,1]]
      } else {
        subtr_list[[i]] <- clus_list[[perms[i,1]]] - clus_list[[perms[i,2]]]
        names(subtr_list)[i] <- paste0(names(clus_list)[perms[i,1]], " - ", names(clus_list)[perms[i,2]])
      }
    }
    umap_inputs <- subtr_list[order(nchar(names(subtr_list)))]
  } else {
    umap_inputs <- clus_list
  }

  if(which_clusters[1]!="all") {
    if(mean(colnames(umap_inputs[[1]])==colnames(umap_inputs[[2]]))!=1) {
      stop("data matrix column mismatch")
    }
    if(mean(row.names(umap_inputs[[1]])==row.names(umap_inputs[[2]]))!=1) {
      stop("data matrix row mismatch")
    }
    for(i in 1:length(umap_inputs)) {
      umap_inputs[[i]] <- umap_inputs[[i]][,which(colnames(umap_inputs[[i]]) %in% which_clusters)]
    }
  }

  if(return_data) {
    return(umap_inputs)
  }

  small_meta_vector <- small_meta[,2]; names(small_meta_vector) <- small_meta[,1]
  bor_class <- small_meta_vector[row.names(umap_inputs[[1]])]

  bor_fun <- function(arg1, split_group = factor(as.character(bor_class))) {
    # testing
    # arg1 = umap_inputs[[2]]
    # split_group = factor(as.character(bor_class))

    bor_out <- Boruta::Boruta(x = arg1, y = split_group, maxRuns = 300, doTrace = 0, holdHistory = TRUE)

    imp <- bor_out$ImpHistory; colnames(imp) <- gsub("^(X|x)","",colnames(imp))
    imp <- imp[,order(colMeans(imp, na.rm = TRUE))]
    imp_score <- colMeans(imp)
    imp_df <- data.frame(cluster = names(imp_score), importance = as.numeric(imp_score))

    dec <- data.frame(cluster = names(bor_out$finalDecision), decision = as.character(bor_out$finalDecision))
    dec$cluster <- gsub("^(X|x)","",dec$cluster)

    bor_df <- merge(x = imp_df, y = dec, by = "cluster", sort = FALSE)
    bor_df <- bor_df[order(bor_df$importance, decreasing = TRUE),]

    bor_att <- Boruta::attStats(bor_out); row.names(bor_att) <- gsub("^(X|x)","",row.names(bor_att))
    bor_imp <- as.data.frame(bor_out$ImpHistory); colnames(bor_imp) <- gsub("^(X|x)","",colnames(bor_imp))
    colmed1 <- unlist(sapply(X = bor_imp, FUN = median))
    bor_imp <- bor_imp[,gsub("^c","",dec$cluster[which(dec$decision=="Confirmed")])]
    if(class(bor_imp)!="data.frame") {
      bor_imp <- data.frame(val1 = bor_imp); colnames(bor_imp) <- dec$cluster[which(dec$decision=="Confirmed")]
    }
    colmed <- unlist(sapply(X = bor_imp, FUN = median))
    bor_imp <- reshape2::melt(data = bor_imp)
    bor_imp$variable <- factor(x = bor_imp$variable, levels = names(colmed)[order(colmed)])
    bor_imp <- bor_imp[!grepl(pattern = "^shadowM", x = bor_imp$variable),]

    borplot_expansion <- 1

    borplt <- ggplot(data = bor_imp) +
      geom_boxplot(mapping = aes(x = variable, y = value), fill = "grey") +
      ylab("Importance") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, size = 14*borplot_expansion, face = "bold", vjust = 0.8),
            axis.text.y = element_text(size = 12*borplot_expansion),
            axis.title.y = ggtext::element_markdown(size = 16*borplot_expansion, margin = margin(r = 10)),
            axis.title.x = element_blank(),
            legend.title = element_blank())

    borscore <- data.frame(cluster = row.names(bor_att), score = bor_att$medianImp, decision = bor_att$decision)

    return(list(boruta_box = borplt, boruta_scores = borscore))
  }

  bor_list <- lapply(X = umap_inputs, FUN = bor_fun, factor(as.character(bor_class)))

  bor_plots <- lapply(X = bor_list, FUN = function(x) return(x[[1]]))
  bor_plot_arr <- ggpubr::ggarrange(plotlist = bor_plots, nrow = 1)

  if(nrow(metadata)==nrow(seurat_object@reductions[[tolower(reduction)]]@cell.embeddings)) {
    check1 <- TRUE
    seu_embeds <- as.data.frame(seurat_object@reductions[[tolower(reduction)]]@cell.embeddings); colnames(seu_embeds) <- c("redx","redy")
    seu_embeds$cluster <- metadata[,cluster_column]

    umap_by_boruta <- function(arg1, embeds = seu_embeds, byred = reduction) {
      # testing
      # arg1 <- bor_list[[1]]
      # embeds <- seu_embeds
      # byred = reduction

      arg1 <- arg1[[2]]
      df <- merge(x = embeds, y = arg1, by.x = 'cluster', by.y = 'cluster', all.x = TRUE, sort = FALSE)
      df <- df[!is.na(df$decision),]
      df$score <- ifelse(df$decision=="Confirmed", df$score, NA)

      text_expansion_factor <- 1

      mapplt <- ggplot(data = df, mapping = aes(x = redx, y = redy, color = score)) +
        geom_point(size = 0.5, alpha = 0.2) +
        scale_colour_gradientn(colours = c("blue", "red"),
                               na.value = "grey",
                               limits = c(min(df$score, na.rm = TRUE), max(df$score, na.rm = TRUE))) +
        xlab(toupper(paste0(byred,"1"))) + ylab(toupper(paste0(byred,"2"))) +
        theme_minimal() +
        guides(color = guide_colourbar(title.position = "top", frame.colour = "black", ticks.colour = "black",
                                       draw.ulim = T, draw.llim = T, ticks.linewidth = 0.5)) +
        theme(axis.title = element_text(size = 20*text_expansion_factor, face = "bold"),
              plot.title = element_text(size = 22*text_expansion_factor, face = "bold", hjust = 0.5),
              axis.text = element_blank(),
              legend.title = element_text(size = 14*text_expansion_factor, hjust = 0.5),
              legend.key.height = unit(5, "mm"),
              legend.key.width = unit(15, "mm"),
              legend.direction = "horizontal",
              legend.position = "bottom",
              legend.text = element_text(size = 12*text_expansion_factor))
      return(mapplt)
    }
    map_plots <- lapply(X = bor_list, FUN = umap_by_boruta, embeds = seu_embeds, byred = reduction)
    map_plot_arr <- ggpubr::ggarrange(plotlist = map_plots, nrow = 1)
  } else {
    check1 <- FALSE
  }

  if(feature_reduction_method[1]=="pca") {
    pca_list <- vector("list", length = length(umap_inputs)); names(pca_list) <- names(umap_inputs)
    for(i in 1:length(umap_inputs)) {
      pca_result <- prcomp(umap_inputs[[i]], scale. = TRUE)
      var_explained <- (pca_result$sdev)^2
      total_variance <- sum(var_explained)
      cumulative_var_explained <- cumsum(var_explained)
      fraction_cumulative_var_explained <- cumulative_var_explained / total_variance
      pca_list[[i]] <- pca_result$x[,1:which(fraction_cumulative_var_explained>pca_fraction_variance)[1]]
    }
  } else {
    pca_list <- umap_inputs
  }

  analysis_list <- vector("list", length = length(pca_list)); names(analysis_list) <- names(pca_list)
  for(i in 1:length(analysis_list)) {
    umap_out <- uwot::umap(X = pca_list[[i]], n_neighbors = umap_n_neighbors, init = "spca",
                           min_dist = umap_min_dist, batch = TRUE, seed = 123, spread = umap_spread)
    tmp_obj <- list(data = pca_list[[i]],
                    raw = "",
                    source = row.names(pca_list[[i]]))
    tmp_obj <- FCSimple::fcs_update(fcs_join_obj = tmp_obj)
    tmp_obj <- FCSimple::fcs_cluster(fcs_join_obj = tmp_obj, language = "r", algorithm = "leiden",
                                     leiden_louvain_resolution = leiden_resolution,
                                     adjacency_knn = umap_n_neighbors, search_only = FALSE,
                                     search_method = "RANN", num_cores = 1)
    colnames(umap_out) <- c("UMAP1","UMAP2"); umap_out <- as.data.frame(umap_out)
    umap_out$pid <- row.names(pca_list[[i]])
    umap_out <- merge(x = umap_out, y = small_meta, by.x = 'pid', by.y = pid_column, sort = FALSE, all.x = TRUE)
    umap_out$group <- as.character(tmp_obj$leiden$clusters)
    umap_out$type <- names(analysis_list)[i]
    analysis_list[[i]] <- umap_out
  }

  plot_umap_embed <- function(arg1, text_expansion_factor = 1, labp = label_points,
                              point_expansion_factor = 1, consider_clustered_clusters = cluster_data,
                              mancol = manual_colors) {
    # testing
    # arg1 = analysis_list[[1]]
    # text_expansion_factor = 1
    # labp = TRUE
    # point_expansion_factor = 1

    plot_title <- arg1$type[1]

    plt <- ggplot() +
      geom_point(data = arg1, mapping = aes(x = UMAP1, y = UMAP2, fill = meta_col),
                 size = 6*point_expansion_factor, pch = 21, alpha = 0.7) +
      theme_minimal() +
      guides(fill = guide_legend(override.aes = list(size = ifelse((6*point_expansion_factor)>5,(6*point_expansion_factor),5), alpha = 1))) +
      ggtitle(plot_title)
    if(labp) {
      tmpdf <- arg1; tmpdf$ptlab <- gsub("_.+$","",arg1$pid)
      if(all(consider_clustered_clusters, "group" %in% colnames(arg1))) {
        plt <- plt + ggrepel::geom_label_repel(data = tmpdf, mapping = aes(x = UMAP1, y = UMAP2, label = ptlab, color = group),
                                               seed = 123, max.overlaps = Inf, size = 3.5*text_expansion_factor, verbose = FALSE,
                                               fontface = "bold")
      } else {
        plt <- plt + ggrepel::geom_label_repel(data = tmpdf, mapping = aes(x = UMAP1, y = UMAP2, label = ptlab), color = "black",
                                               seed = 123, max.overlaps = Inf, size = 3.5*text_expansion_factor, verbose = FALSE,
                                               fontface = "bold")
      }
    }
    if(!is.null(mancol)) {
      plt <- plt + scale_fill_manual(values = mancol)
    }
    plt <- plt +
      theme(axis.title = element_text(size = 20*text_expansion_factor, face = "bold"),
            plot.title = element_text(size = 22*text_expansion_factor, face = "bold", hjust = 0.5),
            axis.text = element_blank(),
            legend.position = "bottom",
            legend.title = element_blank(),
            legend.text = element_text(size = 20*text_expansion_factor))
    return(plt)
  }
  umaps <- lapply(X = analysis_list, FUN = plot_umap_embed)
  umaps_arr <- ggpubr::ggarrange(plotlist = umaps, nrow = 1, common.legend = TRUE, legend = "bottom")

  if(return_as_list) {
    if(check1) {
      return(list(umaps_arr, map_plot_arr, bor_plot_arr))
    } else {
      return(list(umaps_arr, bor_plot_arr))
    }
  } else {
    if(check1) {
      all_arr <- ggpubr::ggarrange(plotlist = list(umaps_arr, map_plot_arr, bor_plot_arr), nrow = 3, heights = c(0.5, 0.5, 0.2))
    } else {
      all_arr <- ggpubr::ggarrange(plotlist = list(umaps_arr, bor_plot_arr), nrow = 2, heights = c(0.6, 0.4))
    }
    return(all_arr)
  }
}


seurat_mean_count_hm <- function(seurat_object,
                                 assay = 'RNA',
                                 cluster_column = 'cell_type',
                                 plot_clusters = c('CM_CD4','EM_CD8'),
                                 pid_column = 'pid',
                                 gene_set = c('CD8A','TNF','IL7R','IL32','SELL','LEF1'),
                                 low_mid_high_cols = c("#DA29D9","black","#fff176"),
                                 scale_per_gene = TRUE,
                                 cluster_rows = FALSE,
                                 cluster_columns = FALSE,
                                 auto_order_genes = TRUE,
                                 get_legend = FALSE,
                                 split_by_pid = TRUE,
                                 text_expansion_factor = 1,
                                 cluster_annotation_color = NULL,
                                 cluster_annotation_ref = "none",
                                 pid_order = NULL,
                                 gene_order = NULL)
{
  require(ComplexHeatmap)
  require(grid)
  require(circlize)
  # need to update so that bar annotation across top is working correctly

  # testing
  # seurat_object = seu
  # assay = 'RNA'
  # cluster_column = 'cell_type'
  # plot_clusters = c('CM_CD4','Mono','EM_CD8','MAIT/gd')
  # pid_column = 'pid'
  # gene_set = c('CD8A','TNF','IL7R','IL32','SELL','LEF1','CD14','CLEC12A','TRAV1-2')
  # low_mid_high_cols = c("#1976D2","black","#FF1D23")
  # scale_per_gene = TRUE
  # cluster_rows = FALSE
  # auto_order_genes = TRUE
  # get_legend = FALSE
  # split_by_pid = TRUE
  # text_expansion_factor = 0.5
  # cluster_annotation_color = annocol # or NULL
  # cluster_annotation_ref = ref

  cluster_numbers <- seurat_object@meta.data[,cluster_column]
  exprs_matrix <- seurat_object[[assay]]$data
  pid <- seurat_object@meta.data[,pid_column]

  include_clusters <- plot_clusters[order(plot_clusters)]
  matrix_keep_row <- which(row.names(exprs_matrix) %in% gene_set)
  gene_set_in <- gene_set[which(gene_set %in% row.names(exprs_matrix))]
  if(length(gene_set_in)==0) {
    stop("no genes from gene_set found in matrix")
  } else {
    print(paste0(length(unique(gene_set_in))," of ",length(unique(gene_set))," features from gene_set found in exprs_matrix"))
  }
  condense_matrix1 <- as.matrix(exprs_matrix[matrix_keep_row,])
  condense_matrix <- t(t(condense_matrix1)[,gene_set_in])
  uclus <- unique(cluster_numbers); uclus <- uclus[order(uclus)]
  if(include_clusters[1]=="all") {
    unclus <- c()
  } else {
    unclus <- uclus[which(!uclus %in% include_clusters)]
  }
  if(length(unclus)!=0) {
    if(split_by_pid) {
      upid <- unique(pid); upid <- upid[order(upid)]
      tmplate_mat <- matrix(data = NA, ncol = length(upid), nrow = nrow(condense_matrix))
      colnames(tmplate_mat) <- upid; row.names(tmplate_mat) <- row.names(condense_matrix)
      split_clus_mats <- vector("list", length = length(include_clusters)+1); names(split_clus_mats) <- c(include_clusters,"other")
      for(i in 1:length(split_clus_mats)) {
        tmp_mat <- tmplate_mat
        if(names(split_clus_mats)[i]=="other") {
          which_clus_match <- which(cluster_numbers %in% unclus)
        } else {
          which_clus_match <- which(cluster_numbers==names(split_clus_mats)[i])
        }
        pid_match <- pid[which_clus_match]
        condense_matrix_match <- condense_matrix[,which_clus_match]
        for(j in 1:nrow(tmp_mat)) {
          gene_values <- condense_matrix_match[which(row.names(condense_matrix_match)==row.names(tmp_mat)[j]),]
          for(k in 1:ncol(tmp_mat)) {
            tmp_mat[j,k] <- mean(gene_values[which(pid_match==colnames(tmp_mat)[k])])
          }
        }
        split_clus_mats[[i]] <- tmp_mat
      }
      large_mat <- do.call(cbind,split_clus_mats)
      # large_mat[is.na(large_mat)] <- 0 # sub in 0 for NaN values
      # large_mat[is.na(large_mat)] <- mean(x = large_mat, na.rm = TRUE)
      for(i in 1:nrow(large_mat)) {
        #   # large_mat[i,is.na(large_mat[i,])] <- mean(x = range(large_mat[i,], na.rm = TRUE))
        large_mat[i,is.na(large_mat[i,])] <- min(x = range(large_mat[i,], na.rm = TRUE))
      }
      # large_mat[is.na(large_mat)] <- mean(x = range(large_mat, na.rm = TRUE))
      # large_mat[is.na(large_mat)] <- 0
      for(i in 1:length(include_clusters)) {
        if(i==1) {
          hm_groups <- rep(x = include_clusters[i], times = length(upid))
        } else {
          hm_groups <- append(hm_groups, rep(x = include_clusters[i], times = length(upid)))
        }
      }
      hm_groups <- c(hm_groups, rep(x = "other", times=length(upid)))
      hm_groups <- factor(hm_groups,levels=c(include_clusters,"other"))
      # if(!show_other) {
      #   other_ind <- which(hm_groups=="other")
      #   if(length(other_ind)!=0) {
      #     large_mat <- large_mat[,-other_ind]
      #     hm_groups <- hm_groups[-other_ind]
      #   }
      # }
      if(auto_order_genes) {
        unique_col <- as.character(unique(hm_groups))
        mean_expr_matrix <- matrix(data = NA, nrow = nrow(large_mat), ncol = length(unique_col))
        row.names(mean_expr_matrix) <- row.names(large_mat); colnames(mean_expr_matrix) <- unique_col
        for(i in 1:nrow(mean_expr_matrix)) {
          for(j in 1:ncol(mean_expr_matrix)) {
            mean_expr_matrix[i,j] <- mean(large_mat[which(row.names(large_mat)==row.names(mean_expr_matrix)[i]),
                                                    which(as.character(hm_groups)==colnames(mean_expr_matrix)[j])])
          }
        }
        row_order_list <- vector("list", length = length(unique_col)); names(row_order_list) <- unique_col
        for(i in 1:nrow(mean_expr_matrix)) {
          which_top <- which.max(mean_expr_matrix[i,])
          row_order_list[[colnames(mean_expr_matrix)[which_top]]] <- append(row_order_list[[colnames(mean_expr_matrix)[which_top]]],
                                                                            row.names(mean_expr_matrix)[i])
        }
        for(i in 1:length(row_order_list)) {
          subs_mem <- t(t(mean_expr_matrix)[,row_order_list[[i]]])
          if(nrow(subs_mem)==0) {
            next
          } else {
            row_order_list[[i]] <- row_order_list[[i]][order(subs_mem[,i], decreasing = TRUE)]
          }
        }
        genes_ordered <- unlist(row_order_list)
        large_mat <- t(t(large_mat)[,genes_ordered])
      }

      # hm_groups <- paste0("cluster ",hm_groups)
      # hm_groups <- c(hm_groups, rep(x = "other", times=length(upid)))
      # hm_groups <- factor(hm_groups,levels=c(paste0("cluster ",include_clusters),"other"))
      # if(!show_other) {
      #   other_ind <- which(hm_groups=="other")
      #   if(length(other_ind)!=0) {
      #     large_mat <- large_mat[,-other_ind]
      #     hm_groups <- hm_groups[-other_ind]
      #   }
      # }
      # hm_groups <- rep(x = paste0("cluster ",include_clusters), times = length(upid))
      # hm_groups <- hm_groups[order(rep(include_clusters, times = length(upid)))]
      # hm_groups <- c(hm_groups, rep(x = "other", times=length(upid)))
    } else {
      unclus_ind <- which(cluster_numbers %in% unclus)
      hm_matrix <- matrix(data = NA,
                          ncol = length(include_clusters)+1,
                          nrow = nrow(condense_matrix))
      colnames(hm_matrix) <- c(include_clusters,"other"); row.names(hm_matrix) <- row.names(condense_matrix)
      for(i in 1:nrow(hm_matrix)) {
        tmp_exprs_count <- condense_matrix[which(row.names(condense_matrix)==row.names(hm_matrix)[i]),]
        for(j in 1:ncol(hm_matrix)) {
          if(colnames(hm_matrix)[j]!="other") {
            hm_matrix[i,j] <- mean(tmp_exprs_count[which(cluster_numbers==colnames(hm_matrix)[j])])
          } else {
            hm_matrix[i,j] <- mean(tmp_exprs_count[unclus_ind])
          }
        }
      }
    }
  } else {
    if(split_by_pid) {
      upid <- unique(pid); upid <- upid[order(upid)]
      tmplate_mat <- matrix(data = NA, ncol = length(upid), nrow = nrow(condense_matrix))
      colnames(tmplate_mat) <- upid; row.names(tmplate_mat) <- row.names(condense_matrix)
      split_clus_mats <- vector("list", length = length(include_clusters)); names(split_clus_mats) <- include_clusters
      for(i in 1:length(split_clus_mats)) {
        tmp_mat <- tmplate_mat
        which_clus_match <- which(cluster_numbers==names(split_clus_mats)[i])
        pid_match <- pid[which_clus_match]
        condense_matrix_match <- condense_matrix[,which_clus_match]
        for(j in 1:nrow(tmp_mat)) {
          gene_values <- condense_matrix_match[which(row.names(condense_matrix_match)==row.names(tmp_mat)[j]),]
          for(k in 1:ncol(tmp_mat)) {
            tmp_mat[j,k] <- mean(gene_values[which(pid_match==colnames(tmp_mat)[k])])
          }
        }
        if(!is.null(pid_order)) {
          pid_ord <- pid_order[which(pid_order %in% colnames(tmp_mat))]
          tmp_mat <- tmp_mat[,pid_ord]
        }
        if(!is.null(gene_order)) {
          gene_ord <- gene_order[which(gene_order %in% row.names(tmp_mat))]
          tmp_mat <- tmp_mat[gene_ord,]
        }
        split_clus_mats[[i]] <- tmp_mat
      }
      large_mat <- do.call(cbind,split_clus_mats)
      # large_mat[is.na(large_mat)] <- 0 # sub in 0 for NaN values
      # large_mat[is.na(large_mat)] <- mean(x = large_mat, na.rm = TRUE)
      for(i in 1:nrow(large_mat)) {
        #   # large_mat[i,is.na(large_mat[i,])] <- mean(x = range(large_mat[i,], na.rm = TRUE))
        large_mat[i,is.na(large_mat[i,])] <- min(x = range(large_mat[i,], na.rm = TRUE))
      }
      # large_mat[is.na(large_mat)] <- mean(x = range(large_mat, na.rm = TRUE))
      # large_mat[is.na(large_mat)] <- 0
      for(i in 1:length(include_clusters)) {
        if(i==1) {
          hm_groups <- rep(x = include_clusters[i], times = length(upid))
        } else {
          hm_groups <- append(hm_groups, rep(x = include_clusters[i], times = length(upid)))
        }
      }
      # if(prefix_cluster) {
      #   hm_groups <- paste0("cluster ",hm_groups)
      #   hm_groups <- factor(hm_groups,levels=paste0("cluster ",include_clusters))
      # } else {
        hm_groups <- factor(hm_groups,levels=include_clusters)
      # }
      # hm_groups <- rep(x = paste0("cluster ",include_clusters), times = length(upid))
      # hm_groups <- hm_groups[order(rep(include_clusters, times = length(upid)))]
      if(auto_order_genes) {
        unique_col <- as.character(unique(hm_groups))
        mean_expr_matrix <- matrix(data = NA, nrow = nrow(large_mat), ncol = length(unique_col))
        row.names(mean_expr_matrix) <- row.names(large_mat); colnames(mean_expr_matrix) <- unique_col
        for(i in 1:nrow(mean_expr_matrix)) {
          for(j in 1:ncol(mean_expr_matrix)) {
            mean_expr_matrix[i,j] <- mean(large_mat[which(row.names(large_mat)==row.names(mean_expr_matrix)[i]),
                                                    which(as.character(hm_groups)==colnames(mean_expr_matrix)[j])])
          }
        }
        row_order_list <- vector("list", length = length(unique_col)); names(row_order_list) <- unique_col
        for(i in 1:nrow(mean_expr_matrix)) {
          which_top <- which.max(mean_expr_matrix[i,])
          row_order_list[[colnames(mean_expr_matrix)[which_top]]] <- append(row_order_list[[colnames(mean_expr_matrix)[which_top]]],
                                                                            row.names(mean_expr_matrix)[i])
        }
        for(i in 1:length(row_order_list)) {
          subs_mem <- t(t(mean_expr_matrix)[,row_order_list[[i]]])
          row_order_list[[i]] <- row_order_list[[i]][order(subs_mem[,i], decreasing = TRUE)]
        }
        genes_ordered <- unlist(row_order_list)
        large_mat <- t(t(large_mat)[,genes_ordered])
      }
    } else {
      hm_matrix <- matrix(data = NA,
                          ncol = length(uclus),
                          nrow = nrow(condense_matrix))
      colnames(hm_matrix) <- uclus; row.names(hm_matrix) <- row.names(condense_matrix)
      for(i in 1:nrow(hm_matrix)) {
        tmp_exprs_count <- condense_matrix[which(row.names(condense_matrix)==row.names(hm_matrix)[i]),]
        for(j in 1:ncol(hm_matrix)) {
          hm_matrix[i,j] <- mean(tmp_exprs_count[which(cluster_numbers==colnames(hm_matrix)[j])])
        }
      }
    }
  }
  if(scale_per_gene) {
    if(split_by_pid) {
      for(i in 1:nrow(large_mat)) {
        large_mat[i,] <- scales::rescale(x = large_mat[i,], to = c(0,1), from = range(large_mat[i,]))
      }
    } else {
      for(i in 1:nrow(hm_matrix)) {
        hm_matrix[i,] <- scales::rescale(x = hm_matrix[i,], to = c(0,1), from = range(hm_matrix[i,]))
      }
    }
  }
  if(is.null(cluster_annotation_color[1])) {
    gg_color_hue <- function(n) {
      hues = seq(15, 375, length = n + 1)
      hcl(h = hues, l = 65, c = 100)[1:n]
    }
    cluster_annotation_color <- gg_color_hue(length(split_clus_mats))
    for(i in 1:length(cluster_annotation_color)) {
      if(i>length(ref)) {
        names(cluster_annotation_color)[i] <- "other"
      } else {
        names(cluster_annotation_color)[i] <- names(ref)[i]
      }
    }
  } #else {
  #   for(i in 1:length(cluster_annotation_color)) {
  #     if(i>length(ref)) {
  #       names(cluster_annotation_color)[i] <- "other"
  #     } else {
  #       names(cluster_annotation_color)[i] <- names(ref)[i]
  #     }
  #   }
  # }
  if(split_by_pid) {
    col_fun = colorRamp2(breaks = c(min(range(large_mat)), mean(range(large_mat)), max(range(large_mat))),
                         colors = low_mid_high_cols)
    if(unlist(cluster_annotation_ref)[1]!="none") {
      type_anno <- rep(NA,length(hm_groups))
      for(i in 1:length(cluster_annotation_ref)) {
        anno_pos_match <- which(hm_groups %in% cluster_annotation_ref[[i]])
        type_anno[anno_pos_match] <- names(cluster_annotation_ref)[i]
      }
      na_pos <- which(is.na(type_anno))
      if(length(na_pos)!=0) {
        type_anno[na_pos] <- "other"
      }

      # ha <- HeatmapAnnotation(bar = factor(type_anno, levels = names(cluster_annotation_color)),
      ha <- HeatmapAnnotation(bar = factor(type_anno),
                              col = list(bar = cluster_annotation_color), show_legend = FALSE,
                              show_annotation_name = FALSE, which = "column", simple_anno_size = unit(3, "mm"),
                              annotation_legend_param = list(legend_height=unit(3,"cm"),
                                                             grid_width=unit(0.6,"cm"),title_position="topleft",
                                                             labels_gp=gpar(fontsize=16*text_expansion_factor),
                                                             title_gp=gpar(fontsize=0*text_expansion_factor)))
      outhm <- ComplexHeatmap::Heatmap(matrix = large_mat, name = "GEX", column_split = hm_groups, # ComplexHeatmap doc 2.7.2 split single heatmap
                                       na_col = "black", col = col_fun,
                                       heatmap_legend_param=list(at=round(c(min(range(large_mat)),max(range(large_mat))),2),
                                                                 labels = c("low", "high"),
                                                                 legend_height=unit(3,"cm"),
                                                                 grid_width=unit(0.6,"cm"),title_position="topleft",
                                                                 labels_gp=gpar(fontsize=14*text_expansion_factor),
                                                                 title_gp=gpar(fontsize=0*text_expansion_factor,fontface="bold")),
                                       show_row_dend = FALSE, show_column_dend = FALSE, cluster_rows = cluster_rows,
                                       cluster_columns = cluster_columns, top_annotation = ha, show_heatmap_legend = FALSE,
                                       row_names_gp = gpar(fontsize=12*text_expansion_factor),
                                       column_title_gp = gpar(fontsize=14*text_expansion_factor, vjust = 20), column_title_side = "top",
                                       show_column_names = FALSE, row_names_side = "left", cluster_column_slices = FALSE)
    } else {
      outhm <- ComplexHeatmap::Heatmap(matrix = large_mat, name = "GEX", column_split = hm_groups, # ComplexHeatmap doc 2.7.2 split single heatmap
                                       na_col = "black", col = col_fun,
                                       heatmap_legend_param=list(at=round(c(min(range(large_mat)),max(range(large_mat))),2),
                                                                 labels = c("low", "high"),
                                                                 legend_height=unit(3,"cm"),
                                                                 grid_width=unit(0.6,"cm"),title_position="topleft",
                                                                 labels_gp=gpar(fontsize=14*text_expansion_factor),
                                                                 title_gp=gpar(fontsize=0*text_expansion_factor,fontface="bold")),
                                       show_row_dend = FALSE, show_column_dend = FALSE, cluster_rows = cluster_rows,
                                       cluster_columns = cluster_columns, show_heatmap_legend = FALSE,
                                       row_names_gp = gpar(fontsize=12*text_expansion_factor),
                                       column_title_gp = gpar(fontsize=14*text_expansion_factor, vjust = 20), column_title_side = "top",
                                       show_column_names = FALSE, row_names_side = "left", cluster_column_slices = FALSE)
    }
    if(get_legend) {
      scaled_leg <- Legend(at=round(c(min(range(large_mat)),max(range(large_mat))),2),
                           labels = c("low", "high"), col_fun = col_fun,
                           legend_height=unit(3,"cm"),
                           grid_width=unit(0.6,"cm"),title_position="topleft",
                           labels_gp=gpar(fontsize=14*text_expansion_factor),
                           title_gp=gpar(fontsize=0*text_expansion_factor,fontface="bold"))
      discrete_leg <- Legend(at = names(cluster_annotation_color), title = "none",
                             title_position="topleft", labels_gp=gpar(fontsize=14*text_expansion_factor),
                             # title_gp=gpar(fontsize=14,fontface="bold"),
                             title_gp = gpar(fontsize=0*text_expansion_factor, alpha = 0),
                             legend_gp = gpar(fill = cluster_annotation_color),
                             grid_height = unit(7, "mm"), grid_width = unit(0.6,"cm"))
      return(list(discrete_leg,scaled_leg))
    } else {
      return(list(out_figure = grid::grid.grabExpr(ComplexHeatmap::draw(outhm)),
                  hm = outhm,
                  tile_data = large_mat))
    }
  } else {
    col_fun = colorRamp2(breaks = c(min(range(hm_matrix)), mean(range(hm_matrix)), max(range(hm_matrix))),
                         colors = low_mid_high_cols)
    outhm <- ComplexHeatmap::Heatmap(name = "feature\nmean\ncount", matrix = hm_matrix, cluster_rows = cluster_rows,
                                     col = col_fun,
                                     heatmap_legend_param=list(at=round(c(min(range(hm_matrix)),max(range(hm_matrix))),2),
                                                               legend_height=unit(3,"cm"),
                                                               grid_width=unit(0.6,"cm"),title_position="topleft",
                                                               labels_gp=gpar(fontsize=14*text_expansion_factor),
                                                               title_gp=gpar(fontsize=15*text_expansion_factor,fontface="bold")),
                                     show_row_dend = FALSE, row_names_gp=gpar(fontsize=13*text_expansion_factor,fontface="bold"),
                                     row_names_gp = gpar(fontsize=12*text_expansion_factor),
                                     column_names_gp = gpar(fontsize=14,fontface="bold"), column_names_side = "top")
    if(get_legend) {
      scaled_leg <- Legend(at=round(c(min(range(large_mat)),max(range(large_mat))),2),
                           labels = c("low", "high"), col_fun = col_fun,
                           legend_height=unit(3,"cm"),
                           grid_width=unit(0.6,"cm"),title_position="topleft",
                           labels_gp=gpar(fontsize=14*text_expansion_factor),
                           title_gp=gpar(fontsize=0*text_expansion_factor,fontface="bold"))
      return(scaled_leg)
    } else {
      return(list(out_figure = grid::grid.grabExpr(ComplexHeatmap::draw(outhm)),
                  hm = outhm,
                  tile_data = hm_matrix))
    }
  }
}

seurat_test_clusters <- function(seurat_object, test_by_column = "condition", pid_column = "pid",
                                 cluster_column = "cell_type", cell_barcode_column = "barcode",
                                 order_x = c("MTB300","Media","GRV"),
                                 bin_colors = c("coral","azure4","seagreen"),
                                 subtract_background = TRUE, comparison_list = NULL,
                                 backround_condition = "Media", y_axis_subset = "T/NK",
                                 connect_points = TRUE, data_paired = TRUE,
                                 return_plots = TRUE, return_plot_data = FALSE,
                                 coord_stretch_factor = 0.11, text_size_factor = 1,
                                 shape_key = NULL, compare_step_distance = 0,
                                 y_as_log = FALSE, coordinate_title_color = FALSE)
{
  require(ggplot2)
  require(ggpubr)

  # testing
  # seurat_object <- seu
  # test_by_column <- "condition"
  # pid_column <- "pid"
  # cluster_column <- "cell_type"
  # cell_barcode_column <- "barcode"
  # order_x <- c("media","stim")
  # bin_colors <- c("coral","azure4","seagreen")
  # subtract_background = FALSE
  # tile_plots = TRUE
  # backround_condition = "media"
  # y_axis_subset = "PBMC"
  # connect_points = TRUE
  # data_paired = TRUE
  # return_plots = TRUE
  # return_plot_data = FALSE
  # coord_stretch_factor = 0.11
  # text_size_factor = 1
  # shape_key = NULL

  names(bin_colors) <- order_x

  test_data <- seurat_object@meta.data[,c(pid_column,test_by_column,cluster_column,cell_barcode_column)]
  upid <- unique(test_data[,pid_column])
  ucondition <- order_x
  uclus <- unique(test_data[,cluster_column]); uclus <- uclus[order(uclus)]

  freq_matrix <- matrix(data = NA, nrow = length(upid), ncol = length(uclus))
  row.names(freq_matrix) <- upid; colnames(freq_matrix) <- uclus
  freq_list <- vector("list", length = length(ucondition)); names(freq_list) <- ucondition
  for(i in 1:length(freq_list)) {
    freq_matrix_copy <- freq_matrix
    tmp_data <- test_data[which(test_data[,test_by_column]==names(freq_list)[i]),]
    for(j in 1:nrow(freq_matrix_copy)) {
      tmp_clus <- tmp_data[,cluster_column][which(tmp_data[,pid_column]==row.names(freq_matrix_copy)[j])]
      for(k in 1:ncol(freq_matrix_copy)) {
        freq_matrix_copy[j,k] <- mean(tmp_clus==colnames(freq_matrix_copy)[k]) * 100
      }
    }
    freq_list[[i]] <- freq_matrix_copy
  }
  if(!return_plots) {
    return(freq_list)
  }

  if(subtract_background) {
    background_index <- which(names(freq_list)==backround_condition)
    if(length(background_index)!=1) {
      stop("error (1)")
    }
    background_data <- freq_list[[background_index]]
    for(i in 1:length(freq_list)) {
      if(any(mean(row.names(freq_list[[i]])==row.names(background_data))!=1,mean(colnames(freq_list[[i]])==colnames(background_data))!=1)) {
        stop("error (2)")
      }
      freq_list[[i]] <- freq_list[[i]] - background_data
    }
    freq_list <- freq_list[-background_index]
  }

  cluster_list <- vector("list", length = ncol(freq_list[[1]])); names(cluster_list) <- colnames(freq_list[[1]])
  for(i in 1:length(cluster_list)) {
    for(j in 1:length(freq_list)) {
      if(j==1) {
        tmp_df <- data.frame(var1 = freq_list[[j]][,i])
      } else {
        tmp_df <- cbind(tmp_df,data.frame(tmp_var = freq_list[[j]][,i]))
        colnames(tmp_df)[ncol(tmp_df)] <- paste0('var',j)
      }
    }
    cluster_list[[i]] <- cbind(tmp_df, data.frame(cluster = rep(names(cluster_list)[i],nrow(freq_list[[1]])),
                                                  pid = row.names(freq_list[[1]])))
    colnames(cluster_list[[i]])[1:length(freq_list)] <- names(freq_list)
  }
  if(!is.null(shape_key[1])) {
    for(i in 1:length(cluster_list)) {
      cluster_list[[i]] <- merge(x = cluster_list[[i]], y = shape_key, by = pid_column, sort = FALSE)
    }
  }
  if(return_plot_data) {
    return(cluster_list)
  }

  iterate_boxplot <- function(arg1, b_colors = bin_colors,
                              tbc = test_by_column,
                              pidc = pid_column,
                              clc = cluster_column,
                              cellbc = cell_barcode_column,
                              should_connect = connect_points,
                              anchor_stim = backround_condition,
                              comps = comparison_list,
                              # tp = tile_plots,
                              dp = data_paired,
                              yax_lab = y_axis_subset,
                              csf = coord_stretch_factor,
                              tsf = text_size_factor,
                              shp = shape_key,
                              csd = compare_step_distance,
                              yal = y_as_log,
                              ctcol = coordinate_title_color)
  {
    require(reshape2)
    require(combinat)
    require(ggplot2)
    require(ggpubr)
    require(scales)

    # testing
    # arg1 = cluster_list[[1]]
    # b_colors = bin_colors
    # tbc = test_by_column
    # pidc = pid_column
    # clc = cluster_column
    # cellbc = cell_barcode_column
    # should_connect = connect_points
    # anchor_stim = backround_condition
    # comps = comparison_list
    # dp = data_paired
    # yax_lab = y_axis_subset
    # csf = coord_stretch_factor
    # tsf = text_size_factor
    # shp = shape_key
    # csd = compare_step_distance

    if(!is.null(shp[1])) {
      shp <- colnames(shape_key)[which(colnames(shape_key)!=pidc)]
    } else {
      shp <- NULL
    }

    arg1$group <- as.character(1:nrow(arg1))
    arg_melt <- reshape2::melt(data = arg1, value.name = clc)
    # colnames(arg_melt) <- c("cluster","pid","group","condition","frequency") # condition = variable
    colnames(arg_melt)[which(colnames(arg_melt)==clc)] <- "frequency"
    colnames(arg_melt)[ncol(arg_melt)] <- "frequency"

    if(!dp) {
      arg_melt <- arg_melt[!is.na(arg_melt$frequency),]
    }

    if(is.null(comps)) {
      compare_these <- combinat::permn(unique(arg_melt$variable))
      for(i in 1:length(compare_these)) {
        compare_these[[i]] <- as.character(compare_these[[i]][1:2])
      }
      find_dupls <- rep(NA,length=length(compare_these))
      for(i in 1:length(find_dupls)) {
        find_dupls[i] <- paste0(compare_these[[i]][order(compare_these[[i]])],collapse="")
      }
      rm_elements <- which(duplicated(find_dupls))
      if(length(rm_elements)>0) {
        compare_these <- compare_these[-rm_elements]
      }
    } else {
      compare_these <- comps
    }
    if(!is.null(shp)) {
      # arg_melt <- merge(x = arg_melt, y = shape_key, by = "pid", all.x = TRUE, all.y = FALSE)
      colnames(arg_melt)[which(colnames(arg_melt)==shp)] <- "shape_group"
    }
    pl <- ggplot(data = arg_melt, mapping = aes(x = variable, y = frequency, color = variable)) +
      geom_boxplot(fill = "#bfbfbf", lwd = 0.5, alpha = 0.4, width = 0.45)
    if(should_connect) {
      pl <- pl + geom_line(mapping = aes(group = group), color = "black")
      if("shape_group" %in% colnames(arg_melt)) {
        pl <- pl + geom_point(mapping = aes(fill = variable, shape = shape_group), size = 4*tsf, color = "black") +
          scale_shape_manual(values = 21:25)
      } else {
        pl <- pl + geom_point(mapping = aes(fill = variable), pch = 21, size = 3*tsf, color = "black")
      }
    } else {
      pl <- pl + geom_dotplot(data = arg_melt, aes(fill = variable), color = "black", stroke = 1.1,
                              binaxis = "y", stackdir = "center", position = "dodge", binpositions="all")
    }
    if(dp) {
      pl <- pl +
        stat_compare_means(paired = TRUE, method = "wilcox", comparisons = compare_these, size = 7.6*tsf, step.increase = csd)
    } else {
      pl <- pl +
        stat_compare_means(paired = FALSE, method = "wilcox", comparisons = compare_these, size = 7.6*tsf, step.increase = csd)
    }
    if(!isFALSE(ctcol[1])) {
      get_title_col <- ctcol[arg_melt[,"cluster"][1]]
    } else {
      get_title_col = "black"
    }
    pl <- pl + coord_cartesian(ylim = c(min(arg_melt$frequency, na.rm = TRUE), max(arg_melt$frequency, na.rm = TRUE)*(1 + length(compare_these)*csf)))
    pl <- pl +
      guides(fill = "none", shape = guide_legend(position = "bottom", override.aes = list(size = 7*tsf),
                                                 theme = theme(legend.text = element_text(size = 18*tsf, margin = margin(t = -2, r = 15, b = 0, l = 1)),
                                                               legend.title = element_blank())),
             color = "none") +
      scale_fill_manual(values = b_colors) +
      scale_color_manual(values = b_colors) +
      ylab(paste0("% of ",yax_lab)) +
      # ggtitle(ifelse(pl_type=="numeric", paste0("cluster ",arg_melt[,"cluster"][1]), arg_melt[,"cluster"][1])) +
      ggtitle(arg_melt[,"cluster"][1]) +
      theme_minimal() +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_text(size = 23*tsf, color = "black"),
            axis.text.x = element_text(size = 22*tsf, face = "bold", color = "black"),
            axis.text.y = element_text(size = 20*tsf),
            plot.title = element_text(size = 25*tsf, hjust = 0.5, face = "bold", color = get_title_col))
    if(yal) {
      pl <- pl + scale_y_continuous(trans=scales::pseudo_log_trans(base = 10))
    }
    if(is.null(shp[1])) {
      pl <- pl + theme(legend.position = "none")
    }
    return(pl)
  }
  test_plots <- lapply(X = cluster_list, FUN = iterate_boxplot, b_colors = bin_colors, dp = data_paired)
  return(test_plots)
}

seurat_size_bar <- function(seurat_object, group_column = "pid", condition_column = "condition", cluster_column = "cell_type", rm_legend = TRUE, bar_outline_lwd = 0.1)
{
  require(ggplot2)
  require(ggpubr)
  require(data.table)

  # testing
  # patient_id <- meta_clus$pid
  # condition <- meta_clus$condition
  # cluster <- meta_clus$cluster

  # testing
  # seurat_object = seu
  # pid_column = "pid"
  # condition_column = "condition"
  # cluster_column = "cell_type"

  plot_data <- data.frame(pid = seurat_object@meta.data[,group_column],
                          condition = seurat_object@meta.data[,condition_column],
                          cluster = seurat_object@meta.data[,cluster_column])
  plot_data$pid_condition <- paste0(plot_data$pid,"_",plot_data$condition)
  plot_data$group_var <- factor(plot_data$cluster)

  cluster_table <- as.matrix(table(plot_data$cluster,plot_data$pid_condition))
  # pre_melt <- cbind(as.matrix(data.frame(cluster = row.names(cluster_table))),cluster_table)
  # data.m <- reshape2::melt(pre_melt, id.vars="cluster")
  data.m <- as.data.frame(cluster_table)
  colnames(data.m) <- c("cluster","PID","size")

  rm_row <- which(data.m$PID=="cluster")
  if(length(rm_row)!=0) {
    data.m <- data.m[-rm_row,]
  }
  data.m$cluster <- factor(data.m$cluster)

  stack_bars <- ggplot(data = data.m, mapping = aes(x = cluster, y = size)) +
    geom_bar(aes(fill = PID), stat = "identity", color = "black", linewidth = bar_outline_lwd) +
    theme_minimal() +
    ylab("number of cells") +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = 12),
          axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5),
          # axis.title = element_text(size = 15, face = "bold"))
          axis.title = element_blank())
  if(rm_legend) {
    stack_bars <- stack_bars + theme(legend.position = "none")
  }
  return(stack_bars)
}

seurat_feature_overlay <- function(seurat_object,
                                   reduction = "umap",
                                   assay = "ADT",
                                   layer = "data", # or "counts"
                                   color_meta = "cell_type",
                                   color_by = "cluster", # c("cluster","cell")
                                   feature_reference = NULL, # allow for csv reference where features in col1 are renamed to features in col2, row-wise for titling plots for readability
                                   plot_features = NULL,
                                   pt_expansion = 0.7,
                                   pt_alpha = 0.25,
                                   text_expansion = 1,
                                   text_expansion_annotate = 1,
                                   label_clusters = "none", # "all", "none", or a set of seurat_object@meta.data[,color_meta] to label
                                   add_title = TRUE,
                                   draw_legend = TRUE,
                                   downsample_size = NA){
  require(ggplot2)
  require(shadowtext)
  require(viridis)
  library(ggrastr)
  require(ggpubr)

  # testing
  # seurat_object <- seu_adt
  # reduction <- "umap"
  # assay <- "ADT"
  # layer <- "data" # or "counts"
  # color_meta = "cell_type"
  # color_by = "cluster" # c("cluster","cell")
  # feature_reference = NULL # allow for csv reference where features in col1 are renamed to features in col2, row-wise for titling plots for readability
  # plot_features = NULL
  # pt_expansion = 0.7
  # pt_alpha = 0.25
  # text_expansion = 1
  # text_expansion_annotate = 1
  # label_clusters = "all" # "all", "none", or a set of seurat_object@meta.data[,color_meta] to label
  # add_title = TRUE
  # downsample_size = NA

  exprs <- seurat_object
  DefaultAssay(object = seurat_object) <- assay
  available_features <- row.names(seurat_object)
  if(is.null(plot_features)) {
    plot_features <- row.names(seurat_object)
  } else {
    plot_features <- plot_features[which(plot_features %in% row.names(seurat_object))]
    if(length(plot_features)==0) {
      stop("None of plot_features found in seurat_object.")
    }
  }
  if(!is.na(downsample_size[1])) {
    sample_loc <- sample(x = 1:dim(seurat_object)[2], size=downsample_size, replace=FALSE)
    seurat_object <- subset(x = seurat_object, cells = sample_loc)
    metadata <- seurat_object@meta.data
    exprs <- seurat_object@assays[[assay]]@layers[["data"]] # assumes normed data
  } else {
    metadata <- seurat_object@meta.data
    exprs <- seurat_object@assays[[assay]]@layers[["data"]]
  }
  row.names(exprs) <- row.names(seurat_object)
  if(assay=="ADT") {
    cd_features <- plot_features[grep("CD[0-9]*",plot_features)]; other_features <- plot_features[-grep("CD[0-9]*",plot_features)]
    cdnum <- as.numeric(gsub("^.+CD","",stringr::str_extract(cd_features,"^.+CD[0-9]*")))
    cd_features <- cd_features[order(cdnum)]; other_features <- other_features[order(other_features)]
    plot_features <- c(cd_features,other_features)
  }
  plt_list <- vector("list", length = length(plot_features)); names(plt_list) <- plot_features
  for(i in 1:length(plt_list)){
    get_feature <- exprs[grep(pattern = paste0("^",names(plt_list)[i],"(-|$)"), x = row.names(seurat_object)),]
    intens <- rep(NA,times=nrow(metadata))
    uclus <- unique(metadata[,color_meta])
    if(color_by=="cluster") {
      for(j in 1:length(uclus)){
        clus_val <- get_feature[which(metadata[,color_meta]==uclus[j])]
        cluster_val <- mean(clus_val)
        intens[which(metadata[,color_meta]==uclus[j])] <- cluster_val
      }
      plt_list[[i]] <- metadata
      plt_list[[i]]$value <- intens
      colnames(plt_list[[i]])[ncol(plt_list[[i]])] <- names(plt_list)[i]
    } else if(color_by=="cell"){
      plt_list[[i]] <- metadata
      plt_list[[i]]$value <- get_feature
      plt_list[[i]] <- plt_list[[i]][order(plt_list[[i]]$value, decreasing = FALSE),] # draw highest values last
      colnames(plt_list[[i]])[ncol(plt_list[[i]])] <- names(plt_list)[i]
      plt_list[[i]]$value <- get_feature
    }
    plt_list[[i]]$feature <- names(plt_list)[i]
    plt_list[[i]]$redx <- seurat_object@reductions[[tolower(reduction)]]@cell.embeddings[,1]
    plt_list[[i]]$redy <- seurat_object@reductions[[tolower(reduction)]]@cell.embeddings[,2]
  }

  adt_on_umap <- function(arg1,
                          lab_pl = add_title,
                          pex = pt_expansion,
                          pal = pt_alpha,
                          metac = color_meta,
                          incl_leg = draw_legend,
                          tex = text_expansion,
                          tex_annotate = text_expansion_annotate,
                          lab_clus = label_clusters){
    # testing
    # arg1 = plt_list[[1]]
    # lab_pl = add_title
    # pex = pt_expansion
    # pal = pt_alpha
    # metac = color_meta
    # incl_leg = FALSE
    # tex = text_expansion
    # tex_annotate = text_expansion_annotate
    # lab_clus = label_clusters

    pl_title <- arg1$feature[1]
    if(!is.factor(arg1[,metac])){
      arg1[,metac] <- factor(arg1[,metac])
    }
    uclus1 <- as.character(unique(arg1[,metac]))
    xval <- rep(NA,times=length(uclus1)); names(xval) <- uclus1[order(uclus1)]; yval <- xval
    for(i in 1:length(xval)){
      xval[i] <- median(arg1$redx[which(arg1[,metac]==names(xval)[i])])
      yval[i] <- median(arg1$redy[which(arg1[,metac]==names(xval)[i])])
    }
    cluslab <- data.frame(xval = xval, yval = yval, labl = names(xval))
    capture_feature <- arg1[,'feature'][1]
    colnames(arg1)[which(colnames(arg1)==metac)] <- "pop"
    colnames(arg1)[which(colnames(arg1)==capture_feature)] <- "cby"
    pl <- ggplot(data = arg1, mapping = aes(x=redx,y=redy)) +
      geom_point_rast(aes(color=cby),pch=19, alpha = pal, cex = pex) +
      scale_color_viridis(option = "D", name = capture_feature) +
      guides(color = guide_colourbar(title.position = "top", frame.colour = "black", ticks.colour = "black",
                                     draw.ulim = T, draw.llim = T, ticks.linewidth = 0.5)) +
      theme_void() +
      theme(legend.key.height = unit(5, "mm"),
            legend.key.width = unit(15, "mm"),
            legend.direction = "horizontal",
            legend.position = "bottom",
            legend.text = element_text(size = 12*tex),
            legend.title = element_blank())
    if(!incl_leg) {
      pl <- pl + theme(legend.position = "none")
    }
    if(lab_clus[1]=="all") {
      pl <- pl + annotate("shadowtext", x = xval, y = yval, label = names(xval), size = 8*tex_annotate)
    } else if(lab_clus[1]!="none") {
      which_anno <- which(names(xval) %in% lab_clus)
      annox <- xval[which_anno]
      annoy <- yval[which_anno]
      lab_df <- data.frame(xval = annox, yval = annoy, id = names(annox))
      pl <- pl + ggrepel::geom_text_repel(data = lab_df, mapping = aes(x = xval, y = yval, label = id),
                                          size = 4*tex_annotate, bg.color = "black", bg.r = 0.075, color = "white", seed = 123)
    }
    if(lab_pl) {
      pl <- pl + ggtitle(capture_feature) +
        theme(plot.title = element_text(hjust = 0.5, size = 24*tex, face = "bold", vjust = -1))
    }
    return(pl)
  }

  out_plots <- lapply(plt_list, adt_on_umap, pex = pt_expansion, pal = pt_alpha,
                      metac = color_meta, incl_leg = draw_legend, tex = text_expansion,
                      tex_annotate = text_expansion_annotate, lab_clus = label_clusters)
  return(out_plots)
}


seurat_reduction_by_value <- function(seurat_object, reduction, values, feature_name, plot_title = NULL,
                                      text_expansion = 1, point_alpha = 0.25, point_expansion = 0.7, return_legend = FALSE,
                                      plot_w_legend = TRUE, legx = 0.5, legy = 0.9, legend_orientation = "horizontal",
                                      leg_height = unit(0.6,"cm"), leg_width = unit(3,"cm"))
{
  require(ggplot2)
  require(ggpubr)
  require(shadowtext)
  require(ggrastr)
  require(viridis)

  get_dr <- as.data.frame(seurat_object@reductions[[tolower(reduction)]]@cell.embeddings)

  if(ncol(get_dr)!=2) {
    stop("reduction cell.embeddings must be a matrix or data frame with two columns")
  }

  capture_col1 <- colnames(get_dr)[1]
  capture_col2 <- colnames(get_dr)[2]
  colnames(get_dr) <- c("redx", "redy")

  indata <- cbind(get_dr, data.frame(value = values))

  plt <- ggplot(data = indata, mapping = aes(x=redx,y=redy)) +
    geom_point_rast(aes(color=value),pch=19, alpha = point_alpha, cex = point_expansion) +
    scale_color_viridis(option = "D", name = feature_name) +
    # guides(color = guide_colorbar(title.position="top", title.hjust = 0.5, title.vjust = 0.5,
    #                               title = feature_name, barwidth = grid::unit(x = 0.44, units = "npc"),
    #                               barheight = grid::unit(x = 6, units = "mm"), direction = "horizontal")) +
    guides(color = guide_colourbar(title.position = "top", frame.colour = "black", ticks.colour = "black",
                                   draw.ulim = T, draw.llim = T, ticks.linewidth = 0.5, direction = "horizontal")) +
    xlab(capture_col1) + ylab(capture_col2) +
    theme_void() +
    theme(plot.title = element_blank(),
          legend.key.height = unit(5, "mm"),
          legend.key.width = unit(15, "mm"),
          legend.direction = "horizontal",
          legend.position = "bottom",
          legend.text = element_text(size = 12*text_expansion),
          legend.title = element_blank())
  if(any(is.null(legx), is.null(legy))) {
    plt <- plt + theme(legend.position = "bottom")
  } else {
    plt <- plt + theme(legend.position = c(legx, legy))
  }
  if(!is.null(plot_title)) {
    plt <- plt + ggtitle(plot_title) +
      theme(plot.title = element_text(hjust = 0.5, size = 24*text_expansion, face = "bold", vjust = -1))
  }
  if(!plot_w_legend) {
    plt <- plt + theme(legend.position = "none")
  }
  if(return_legend) {
    require(ComplexHeatmap)
    leg_lab <- round(seq(from=floor(min(range(indata$value))),to=ceiling(max(range(indata$value))),length.out=5),0)
    col_fun <- colorRamp2(breaks = leg_lab,
                          colors = viridis(n = length(leg_lab)))
    if(legend_orientation=="horizontal") {
      leg <- Legend(at = leg_lab,
                    labels = leg_lab,
                    col_fun = col_fun, direction = "horizontal",
                    legend_height = leg_height, title = feature_name,
                    legend_width = leg_width, title_position = "topcenter",
                    labels_gp = gpar(fontsize=14*text_expansion),
                    title_gp = gpar(fontsize = 15*text_expansion, fontface = "bold"))
    } else if(legend_orientation=="vertical") {
      leg <- Legend(at = leg_lab,
                    labels = leg_lab, col_fun = col_fun,
                    legend_height = leg_height, title = feature_name,
                    grid_width = leg_width, title_position = "topleft",
                    labels_gp = gpar(fontsize=14*text_expansion),
                    title_gp = gpar(fontsize=15*text_expansion, fontface="bold"))
    }
    return(leg)
  } else {
    return(plt)
  }
}


seurat_feature_violin <- function(seurat_object, plot_features, categorical_column = "leiden2",
                                  plot_categorical_types = "all", assay = "RNA", text_expansion = 1,
                                  nudge_nonzero = 0.35, y_limit_expansion_factor = 0.5,
                                  condition = "media", condition_cat = "condition") {
  # testing
  # seurat_object = seu_adt
  # plot_features = feature_prot
  # categorical_column = "cell_type"
  # plot_categorical_types = "all"
  # assay = "ADT"
  # text_expansion = 1
  # nudge_nonzero = 0.35
  # y_limit_expansion_factor = 0.5
  # condition = "media"
  # condition_cat = "stim"

  require(ggplot2)

  # feature_counts <- GetAssayData(object = seurat_object, assay = assay, layer = "counts")
  feature_counts <- seurat_object@assays[[assay]]@layers$data
  row.names(feature_counts) <- row.names(seurat_object)
  colnames(feature_counts) <- colnames(seurat_object)
  meta <- seurat_object@meta.data

  plot_features <- plot_features[which(plot_features %in% row.names(feature_counts))]
  if(length(plot_features)==0) {
    stop("no 'plot_features' found in counts matrix")
  } else{
    rm_features <- plot_features[which(!plot_features %in% row.names(feature_counts))]
    if(length(rm_features)!=0) {
      warning(paste0("one or more 'plot_features' not found in counts matrix: ",paste0(rm_features, collapse = ", ")))
    }
  }

  spl_mat2 <- vector("list", length = length(plot_features)); names(spl_mat2) <- plot_features
  for(i in 1:length(plot_features)) {
    ct_index <- which(row.names(feature_counts)==plot_features[i])
    spl_mat2[[i]] <- feature_counts[ct_index,]
  }
  if(plot_categorical_types[1]=="all") {
    plot_categorical_types <- unique(meta[,categorical_column])
  } else {
    plot_categorical_types <- plot_categorical_types[which(plot_categorical_types %in% meta[,categorical_column])]
  }
  spl_mat <- vector("list",length(spl_mat2)); names(spl_mat) <- names(spl_mat2)
  for(i in 1:length(spl_mat2)) {
    spl_mat[[i]] <- data.frame(ct = spl_mat2[[i]], cat = meta[,categorical_column],
                               condition = meta[,condition_cat])
    if(!is.na(condition[1])) {
      spl_mat[[i]] <- spl_mat[[i]][which(spl_mat[[i]]$condition %in% condition),]
    }
    spl_mat[[i]] <- spl_mat[[i]][which(spl_mat[[i]]$cat %in% plot_categorical_types),]
    spl_mat[[i]]$gene <- names(spl_mat)[i]
  }

  violin_internal <- function(indata, texp = text_expansion,
                              nudge_nz = nudge_nonzero, yle = y_limit_expansion_factor) {
    # testing
    # indata <- spl_mat[[1]]
    # texp = text_expansion
    # nudge_nz = nudge_nonzero
    # yle = y_limit_expansion_factor

    require(ggrastr)

    freqs <- data.frame(cat = unique(indata$cat), freq = rep(NA,length(unique(indata$cat))))
    for(i in 1:nrow(freqs)) {
      freqs$freq[i] <- round(mean(indata$ct[which(indata$cat==freqs$cat[i])]!=0)*100,3)
    }
    indata_jitter <- indata[which(indata$ct!=0),]

    expand_out <- max(indata$ct)*(yle/max(indata$ct) + 1)

    plt <- ggplot(data = indata, aes(x = cat, y = ct, fill = cat)) +
      geom_violin(scale = "width", trim = TRUE, alpha = 0.7) +
      # geom_jitter(data = indata_jitter, width = 0.2, height = 0, size = 1, alpha = 0.5) +
      ggrastr::geom_jitter_rast(data = indata_jitter, width = 0.2, height = 0, size = 1, alpha = 0.5) +
      scale_fill_viridis_d() +
      coord_cartesian(ylim = c(0, max(indata$ct))) +
      coord_flip() +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 0),
            legend.position = "none") +
      geom_text(data = freqs, aes(x = cat, y = max(indata$ct), label = freq), size = 5*texp,
                nudge_y = ifelse(max(freqs$freq>=10), nudge_nz*1, nudge_nz)) +
      expand_limits(y = c(0, ifelse(max(freqs$freq>=10), expand_out*1.025, expand_out))) +
      ylab("normalized count") +
      ggtitle(indata$gene[1]) +
      theme(axis.text.y = element_text(size = 15*texp, face = "bold"),
            axis.title.y = element_blank(),
            axis.text.x = element_text(size = 14*texp),
            axis.title.x = element_text(size = 15*texp, face = "bold"),
            plot.title = element_text(size = 18*texp, face = "bold", hjust = 0.5))
    return(plt)
  }
  vplots <- lapply(X = spl_mat, FUN = violin_internal, texp = text_expansion,
                   nudge_nz = nudge_nonzero, yle = y_limit_expansion_factor)
  return(vplots)
}

seurat_feature_violin_test <- function(seurat_object,
                                       plot_features,
                                       categorical_column,
                                       plot_categorical_types,
                                       assay = "RNA",
                                       text_expansion = 1,
                                       condition = NULL,
                                       condition_cat = "condition",
                                       test_cat = "time_to_group",
                                       test_method = "wilcox",
                                       rotate_x = TRUE,
                                       comparison_list = NULL,
                                       apply_p.adjust = TRUE,
                                       add_pvalue_step = FALSE,
                                       show_quantiles = TRUE)
{
  # testing
  # assay = "RNA"
  # text_expansion = 1
  # condition = NULL
  # condition_cat = "condition"
  # test_cat = "time_to_group"
  # test_method = "wilcox"
  # rotate_x = TRUE
  # comparison_list = NULL
  # apply_p.adjust = TRUE
  # add_pvalue_step = FALSE
  # show_quantiles = TRUE
  #
  # seurat_object = seu_rna
  # plot_features = fg3
  # categorical_column = "annotation"
  # plot_categorical_types = vc3
  # condition = "MTB300"
  # rotate_x = F
  # show_quantiles = T
  # comparison_list = list(c(paste0(vc3[1],"\ncontroller"),paste0(vc3[1],"\nprogressor")),
  #                        c(paste0(vc3[2],"\ncontroller"),paste0(vc3[2],"\nprogressor")),
  #                        c(paste0(vc3[3],"\ncontroller"),paste0(vc3[3],"\nprogressor")))
  # test_cat = "progr_contr"
  # add_pvalue_step = FALSE
  # text_expansion = 0.75
  # assay = "RNA"


  require(ggplot2)
  require(rstatix)
  require(ggtext)
  require(ggpubr)
  require(ggrastr)

  if(!is.null(condition)) {
    seurat_object <- subset(x = seurat_object, subset = condition == condition)
  }

  feature_counts <- seurat_object@assays[[assay]]@layers$data
  row.names(feature_counts) <- row.names(seurat_object)
  colnames(feature_counts) <- colnames(seurat_object)
  meta <- seurat_object@meta.data

  plot_features <- plot_features[which(plot_features %in% row.names(feature_counts))]
  if(length(plot_features)==0) {
    stop("no 'plot_features' found in counts matrix")
  } else{
    rm_features <- plot_features[which(!plot_features %in% row.names(feature_counts))]
    if(length(rm_features)!=0) {
      warning(paste0("one or more 'plot_features' not found in counts matrix: ",paste0(rm_features, collapse = ", ")))
    }
  }

  spl_mat2 <- vector("list", length = length(plot_features)); names(spl_mat2) <- plot_features
  for(i in 1:length(plot_features)) {
    ct_index <- which(row.names(feature_counts)==plot_features[i])
    spl_mat2[[i]] <- feature_counts[ct_index,]
  }
  if(plot_categorical_types[1]=="all") {
    plot_categorical_types <- unique(meta[,categorical_column])
  } else {
    plot_categorical_types <- plot_categorical_types[which(plot_categorical_types %in% meta[,categorical_column])]
  }
  spl_mat <- vector("list",length(spl_mat2)); names(spl_mat) <- names(spl_mat2)
  for(i in 1:length(spl_mat2)) {
    spl_mat[[i]] <- data.frame(ct = spl_mat2[[i]], cat = meta[,categorical_column],
                               condition = meta[,condition_cat],
                               test_cat = meta[,test_cat])
    if(!is.null(condition[1])) {
      spl_mat[[i]] <- spl_mat[[i]][which(spl_mat[[i]]$condition %in% condition),]
    }
    spl_mat[[i]] <- spl_mat[[i]][which(spl_mat[[i]]$cat %in% plot_categorical_types),]
    spl_mat[[i]]$gene <- names(spl_mat)[i]
  }

  violin_internal <- function(indata,
                              texp = text_expansion,
                              tm = test_method,
                              rotx = rotate_x,
                              comps = comparison_list,
                              cond = condition,
                              apa = apply_p.adjust,
                              astep = add_pvalue_step,
                              shq = show_quantiles) {
    # testing
    # indata = spl_mat[[1]]
    # texp = text_expansion
    # tm = test_method
    # rotx = rotate_x
    # comps = comparison_list
    # apa = apply_p.adjust
    # astep = add_pvalue_step
    # shq = show_quantiles
    # cond = condition

    indata$add_col <- paste0(indata$cat, "\n", indata$test_cat)
    concat_table_names <- names(table(indata$add_col))

    # comps <- split(concat_table_names, f = gsub(pattern = "\n.+$", replacement = "", x = concat_table_names))
    if(is.null(comps)) {
      compare_these <- combinat::permn(concat_table_names)
      for(i in 1:length(compare_these)) {
        compare_these[[i]] <- as.character(compare_these[[i]][1:2])
      }
      find_dupls <- rep(NA,length=length(compare_these))
      for(i in 1:length(find_dupls)) {
        find_dupls[i] <- paste0(compare_these[[i]][order(compare_these[[i]])],collapse="")
      }
      rm_elements <- which(duplicated(find_dupls))
      if(length(rm_elements)>0) {
        compare_these <- compare_these[-rm_elements]
      }
    } else {
      compare_these <- comps
    }

    # freqs <- data.frame(cat = unique(indata$cat), freq = rep(NA,length(unique(indata$cat))))
    freqs <- data.frame(add_col = unique(indata$add_col), freq = rep(NA,length(unique(indata$add_col))))
    for(i in 1:nrow(freqs)) {
      # freqs$freq[i] <- round(mean(indata$ct[which(indata$cat==freqs$cat[i])]!=0)*100,3)
      freqs$freq[i] <- round(mean(indata$ct[which(indata$add_col==freqs$add_col[i])]!=0)*100,3)
    }
    indata_jitter <- indata[which(indata$ct!=0),]

    spl_indata <- split(x = indata, f = indata$add_col)
    check_sums <- sapply(spl_indata, function(x) return(sum(x$ct)))
    if(all(sum(check_sums==0)!=0,tolower(tm)!="zir")) {
      warning("some data is all-zero: using zir to test")
      tm <- "zir"
    }

    if(tolower(tm)=="zir") {
      require(ZIR)
      require(dplyr)
      zir_p <- rep(NA,length(spl_indata))
      grp1s <- zir_p; grp2s <- zir_p; teststats <- zir_p
      for(i in 1:length(zir_p)) {
        grp1s[i] <- spl_indata[[compare_these[[i]][1]]]$add_col[1]
        grp2s[i] <- spl_indata[[compare_these[[i]][2]]]$add_col[2]
        zir_test <- ziw(x = spl_indata[[compare_these[[i]][1]]]$ct,
                        y = spl_indata[[compare_these[[i]][2]]]$ct, perm = FALSE)
        zir_p[i] <- zir_test$p.value
        teststats[i] <- zir_test$statistics
      }
      zir_p <- ifelse(is.na(zir_p),1,zir_p)
      teststats <- ifelse(is.na(teststats),0,teststats)
      zir_grp1 <- sapply(X = spl_indata, FUN = function(x) return(x[1]))
      test_res <- dplyr::as_tibble(data.frame('.y.' = 'ct',
                                              'group1' = grp1s,
                                              'group2' = grp2s,
                                              'n1' = 1,
                                              'n2' = 1,
                                              'statistic' = teststats,
                                              'p' = zir_p))
      if(all(length(zir_p)>1, apply_p.adjust)) {
        test_res <- rstatix::adjust_pvalue(data = test_res, p.col = "p", output.col = "p.adj")
        test_res <- rstatix::add_significance(data = test_res, p.col = "p.adj", output.col = "p.adj.signif")
      } else {
        test_res <- rstatix::add_significance(data = test_res, p.col = "p", output.col = "p.adj.signif")
      }
      attributes(test_res)$args$data <- indata[,c('ct','add_col')]
      attributes(test_res)$args$formula <- formula(ct ~ add_col)
    } else if(tm=="wilcox") {
      test_res <- rstatix::wilcox_test(data = indata, formula = ct ~ add_col, paired = FALSE, comparisons = compare_these)
    } else {
      stop("only wilcox or zir (zero-inflated wilcox variant) testing supported for now")
    }

    indata_range <- diff(range(indata$ct)) # smaller range = smaller step.increase; larger range = larger step.increase
    if(astep) {
      test_res <- rstatix::add_y_position(test = test_res, step.increase = indata_range*.1)
    } else {
      test_res <- rstatix::add_y_position(test = test_res, step.increase = 0)
    }

    plt <- ggplot()
    if(shq) {
      plt <- plt + geom_violin(data = indata, aes(x = add_col, y = ct, fill = test_cat), scale = "width", trim = TRUE, alpha = 0.7, draw_quantiles = TRUE) +
        ggrastr::geom_jitter_rast(data = indata_jitter, aes(x = add_col, y = ct), width = 0.2, height = 0, size = 1, alpha = 0.5) +
        # scale_fill_viridis_d() +
        stat_pvalue_manual(as.data.frame(test_res), label = "p.adj.signif", tip.length = 0.01, size = 7*texp) +
        theme_minimal() +
        ylab("normalized count") +
        theme(legend.position = "none",
              axis.title.y = element_text(size = 15*texp, face = "bold"),
              axis.text.x = element_text(size = 15*texp, face = "bold", angle = ifelse(rotx, 90, 0), hjust = 0.5),
              axis.title.x = element_blank())
    } else {
      plt <- plt + geom_violin(data = indata, aes(x = add_col, y = ct, fill = test_cat), scale = "width", trim = TRUE, alpha = 0.7) +
        ggrastr::geom_jitter_rast(data = indata_jitter, aes(x = add_col, y = ct), width = 0.2, height = 0, size = 1, alpha = 0.5) +
        # scale_fill_viridis_d() +
        stat_pvalue_manual(as.data.frame(test_res), label = "p.adj.signif", tip.length = 0.01, size = 7*texp) +
        theme_minimal() +
        ylab("normalized count") +
        theme(legend.position = "none",
              axis.title.y = element_text(size = 15*texp, face = "bold"),
              axis.text.x = element_text(size = 15*texp, face = "bold", angle = ifelse(rotx, 90, 0), hjust = 0.5),
              axis.title.x = element_blank())
    }
    # plt <- plt + geom_text(data = freqs, aes(x = add_col, y = max(indata$ct), label = freq), size = 5*texp,
    # nudge_y = ifelse(max(freqs$freq>=10), nudge_nz*1, nudge_nz))
    plt <- plt + geom_text(data = freqs, aes(x = add_col, y = min(indata$ct)-0.1, label = freq), size = 5*texp)
    if(!is.null(cond)) {
      plt <- plt + ggtitle(paste0("**",indata$gene[1],"** (",cond,")")) + theme(plot.title = ggtext::element_markdown(size = 18*texp, face = "bold", hjust = 0.5))
    } else {
      plt <- plt + ggtitle(indata$gene[1]) + theme(plot.title = element_text(size = 18*texp, face = "bold", hjust = 0.5))
    }
    return(plt)
  }
  vplots <- lapply(X = spl_mat, FUN = violin_internal, texp = text_expansion)
  return(vplots)
}


test_clusters_cat <- function(pid, clusters, condition, cat, stat_compares, y_axis_subset = "PBMC",
                              coord_stretch_factor = 0.1, text_size_factor = 0.8, color_map = NA,
                              dot_size = 0.5, bin_width = 0.5, point_size = 1, x_levels = NA,
                              subtract_background = FALSE, y_stretch_method = c("multiply","add")) {
  require(ggplot2)
  require(ggpubr)

  # testing
  # pid = py_meta$bid
  # clusters = py_meta$leiden2
  # condition = py_meta$Condition
  # cat = py_meta$progr_contr
  # y_axis_subset = "PBMC"
  # coord_stretch_factor = 0.1
  # text_size_factor = 0.7
  # color_map = NA
  # stat_compares = pc_compares
  # point_size = 2.5
  # subtract_background = TRUE

  cat_map <- data.frame(pid = pid, cat = cat)
  cat_map <- cat_map[!duplicated(cat_map$pid),]

  indata <- data.frame(pid = pid, clusters = clusters, condition = condition, cat = cat)
  indata$merge_id <- paste0(indata$pid,"_",indata$cat)
  uclus <- unique(indata$clusters)
  spl_uclus <- split(uclus,f=gsub("[0-9]+","",uclus))
  for(i in 1:length(spl_uclus)) {
    spl_uclus[[i]] <- spl_uclus[[i]][order(as.numeric(stringr::str_extract(string = spl_uclus[[i]], pattern = "[0-9]+")))]
  }
  uclus <- as.character(unlist(spl_uclus))
  ucondition <- unique(indata$condition)
  upid <- unique(indata$pid)
  freq_list <- vector("list", length = length(ucondition)); names(freq_list) <- ucondition
  freq_list_1 <- vector("list", length = length(uclus)); names(freq_list_1) <- uclus
  freq_matrix <- matrix(data = NA, nrow = length(upid), ncol = length(uclus))
  row.names(freq_matrix) <- upid; colnames(freq_matrix) <- uclus
  for(i in 1:length(freq_list)) {
    tmp_mat <- freq_matrix
    tmp_data <- indata[which(indata$condition==names(freq_list)[i]),]
    for(j in 1:nrow(freq_matrix)) {
      tmp_data_2 <- tmp_data[which(tmp_data$pid==row.names(freq_matrix)[j]),]
      for(k in 1:ncol(freq_matrix)) {
        tmp_mat[j,k] <- mean(tmp_data_2$clusters==colnames(freq_matrix)[k]) * 100
      }
    }
    freq_list[[i]] <- tmp_mat
  }
  if(subtract_background) {
    back_condition <- freq_list[["Media"]]
    freq_list <- lapply(X = freq_list, FUN = function(arg1, subtr_mat = back_condition){
      return(arg1 - subtr_mat)
    })
  }
  # for(i in 1:length(freq_list)) {
  #   row.names(freq_list[[i]]) <- paste0(names(freq_list)[i],"_",row.names(freq_list[[i]]))
  # }
  cl_list <- vector("list", length = length(uclus)); names(cl_list) <-  uclus
  for(i in 1:length(cl_list)) {
    cl_list[[i]] <- do.call(cbind, lapply(X = freq_list, FUN = function(arg1, cn = names(cl_list)[i]){
      return(arg1[,which(colnames(arg1)==cn)])
    }))
  }
  melted_list <- vector("list", length = length(cl_list)); names(melted_list) <- names(cl_list)
  ucat <- unique(cat_map$cat)
  for(i in 1:length(melted_list)) {
    melted_list[[i]] <- reshape2::melt(data = cl_list[[i]])
    colnames(melted_list[[i]]) <- c("pid","condition","frequency")
    melted_list[[i]]$cluster <- names(cl_list)[i]
    add_cat <- rep(NA,nrow(melted_list[[i]]))
    for(j in 1:length(ucat)) {
      add_cat[which(melted_list[[i]]$pid %in% cat_map$pid[which(cat_map$cat==ucat[j])])] <- ucat[j]
    }
    melted_list[[i]]$cat <- add_cat
    melted_list[[i]]$joined_cat <- paste0(melted_list[[i]]$condition,"_",melted_list[[i]]$cat)
    melted_list[[i]]$frequency <- ifelse(is.na(melted_list[[i]]$frequency),0,melted_list[[i]]$frequency)
  }

  iter_boxpl <- function(arg1, tsf = text_size_factor, my_compares = stat_compares, ptsize = point_size,
                         csf = coord_stretch_factor, yax_lab = y_axis_subset, over_col = color_map,
                         ysm = y_stretch_method, xlev = x_levels) {
    # testing
    # arg1 <- melted_list[[1]]
    # tsf = text_size_factor
    # my_compares = stat_compares
    # ptsize = point_size
    # csf = coord_stretch_factor
    # yax_lab = y_axis_subset
    # over_col = color_map

    arg1$joined_cat <- gsub("_","\n",arg1$joined_cat)
    if(!is.na(xlev[1])) {
      arg1$joined_cat <- factor(arg1$joined_cat, levels = xlev)
    }
    freq_range <- abs(diff(range(arg1$frequency)))
    # fraction_of_range <- max(arg1$frequency, na.rm = TRUE)/freq_range

    plt <- ggplot(data = arg1, mapping = aes(x = joined_cat, y = frequency, color = joined_cat, fill = joined_cat)) +
      geom_boxplot(fill = "#bfbfbf", lwd = 0.5, alpha = 0.4, width = 0.45, outlier.shape = NA) +
      # geom_dotplot(data = arg1, mapping = aes(fill = joined_cat), color = "black", stroke = 1.1,
      #              binaxis = "y", stackdir = "center", position = "dodge", binpositions = "all", binwidth = binw) +
      geom_jitter(data = arg1, width = 0.2, size = ptsize, alpha = 0.75, pch = 21, stroke = 0.25, color = "black") +
      stat_compare_means(paired = FALSE, method = "wilcox", comparisons = my_compares, size = 7.6*tsf, label.y = max(arg1$frequency))
    if(ysm=="multiply") {
      plt <- plt + coord_cartesian(ylim = c(min(arg1$frequency, na.rm = TRUE), max(arg1$frequency, na.rm = TRUE)*(1 + length(compare_these)*csf)))
    } else if(ysm=="add") {
      plt <- plt + coord_cartesian(ylim = c(min(arg1$frequency, na.rm = TRUE), max(arg1$frequency, na.rm = TRUE)+(freq_range*csf)))
    }
    if(!is.na(over_col)) {
      plt <- plt + scale_fill_manual(values = over_col) +
        scale_color_manual(values = over_col)
    }
    plt <- plt +
      ylab(paste0("% of ",yax_lab)) +
      ggtitle(paste0("cluster ",arg1[,"cluster"][1])) +
      theme_minimal() +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_text(size = 23*tsf, color = "black"),
            axis.text.x = element_text(size = 18*tsf, face = "bold", color = "black"),
            axis.text.y = element_text(size = 20*tsf),
            plot.title = element_text(size = 25*tsf, hjust = 0.5, face = "bold"),
            legend.position = "none")
    return(plt)
  }
  outplots <- lapply(X = melted_list, FUN = iter_boxpl, tsf = text_size_factor,
                     my_compares = stat_compares, csf = coord_stretch_factor,
                     yax_lab = y_axis_subset, over_col = color_map,
                     ptsize = point_size, ysm = y_stretch_method)
  return(outplots)
}

seurat_dge <- function(seurat_object,
                       dge_method = c("mast", "wilcox", "pseudobulk"),
                       assay = "RNA",
                       freq_expressed = 0.1,
                       fc_threshold = log2(1.5),
                       test_clusters = "all", # one or more clusters to test, or "all"
                       mast_lane = NULL,
                       cluster_column = "cell_type",
                       category_column = "age_group",
                       test_categories = c("younger","older"), # order matters: translates into (test_categories[1]/test_categories[2]) for DESeq2 pseudobulk
                       test_condition = "all",
                       condition_column = "condition",
                       pid_column = "pid",
                       pseudobulk_test_mode = c("cluster_identity","cluster_by_category","cluster_by_condition"),
                       return_all_pseudobulk = FALSE)
{
  suppressPackageStartupMessages({
    require(Seurat)
  })
  # testing
  # seurat_object = seu_small
  # dge_method = "mast"
  # assay = "RNA"
  # freq_expressed = 0.1
  # fc_threshold = log2(1.5)
  # test_clusters = "all"
  # cluster_column = "cell_type"
  # category_column = "age_group"
  # test_categories = c("younger","older")
  # test_condition = c("stim","media")
  # condition_column = "condition"
  # pid_column = "pid"
  # pseudobulk_test_mode = "cluster_identity"
  # mast_lane = NULL


  if(test_condition[1]!="all") {
    conditions <- test_condition[test_condition %in% seurat_object@meta.data[,condition_column]]
  } else {
    conditions <- unique(seurat_object@meta.data[,condition_column])
  }
  if(length(conditions)==0) {
    conditions <- NULL
    warning("None of the requested 'conditions' found in seurat_object@meta.data[,condition_column]. Continuing without subsetting on any condition(s).")
  } else {
    which_is_condition <- seurat_object@meta.data[,condition_column] %in% conditions
    seurat_object <- subset(x = seurat_object, cells = which(which_is_condition))
    print(paste0("seurat_object subset for: ",paste0(unique(seurat_object@meta.data[,condition_column]), collapse = ",")))
  }
  seu_conditions <- ifelse(test_condition[1]=="all","all",unique(seurat_object@meta.data[,condition_column]))
  if(length(seu_conditions)!=1) {
    seu_conditions <- paste0(seu_conditions, collapse = ".")
  }
  if(test_clusters[1]!="all") {
    test_clusters <- test_clusters[test_clusters %in% seurat_object@meta.data[,cluster_column]]
    if(length(test_clusters)==0) {
      annos <- NULL
      warning("None of 'test_clusters' found in seurat_object@meta.data[,cluster_column]. Continuing with test(s) on all cells.")
    } else {
      annos <- test_clusters
    }
  } else {
    annos <- unique(seurat_object@meta.data[,cluster_column])
  }
  if(is.null(test_categories)) {
    test_cats <- NULL
  } else if(test_categories[1]!="all") {
    test_categories <- test_categories[test_categories %in% seurat_object@meta.data[,category_column]]
    if(length(test_categories)==0) {
      test_cats <- NULL
      warning("no categories by 'test_categories' were found in seurat_object@meta.data[,category_column]")
    } else {
      test_cats <- test_categories
    }
  } else {
    test_categories <- unique(seurat_object@meta.data[,category_column])
    test_cats <- test_categories
  }
  if(all(is.null(annos), is.null(test_cats))) {
    stop("Nothing to compare")
  }

  if(tolower(dge_method)=="pseudobulk") {
    require(data.table)
    require(DESeq2)
    require(ashr)
    # if(length(test_cats)!=2) {
    #   stop("Test_categories must be length 2 when doing pseudobulk. If testing cluster vs rest, do c('in','out') mapped to what cluster is being tested. If testing ")
    # }
    capture_dir <- system.file(package = "seutools")
    Matrix::writeMM(obj = seurat_object@assays[[assay]]@layers$counts,
                    file = paste0(capture_dir,"/temp_files/__python_ct_matrix__.mtx"))
    write.csv(x = seurat_object@meta.data,
              file = paste0(capture_dir,"/temp_files/__python_obs_matrix__.csv"), row.names = FALSE)
    write.csv(x = data.frame(v1 = row.names(seurat_object)),
              file = paste0(capture_dir,"/temp_files/__python_feature_names__.csv"), row.names = FALSE)
    system(command = paste0("python ",
                            paste0(capture_dir,"/python/create_pseudobulk_profile.py")," ",
                            paste0(capture_dir,"/temp_files/__python_ct_matrix__.mtx")," ",
                            capture_dir,"/temp_files"," ",
                            paste0(capture_dir,"/temp_files/__python_obs_matrix__.csv")," ",
                            ifelse(pseudobulk_test_mode=="cluster_identity",cluster_column,category_column)," ",
                            paste0(capture_dir,"/temp_files/__python_feature_names__.csv")," ",
                            condition_column," ",
                            pid_column," ",
                            pseudobulk_test_mode))
    # ct_data <- read.csv(paste0(capture_dir,"/temp_files/__pseudobulk_sum_counts__.csv"), check.names = FALSE, header = FALSE)

    drop_novar <- function(arg1, rpat = gsub("_$","",id_pattern)) {
      if(ncol(arg1)>nrow(arg1)) {
        row.names(arg1) <- stringr::str_extract(string = row.names(arg1), pattern = rpat)
        check_var <- apply(X = arg1, MARGIN = 2, FUN = stats::var)
        which_novar <- which(check_var==0)
        if(length(which_novar)!=0) {
          arg1 <- arg1[,-which_novar]
        }
      } else {
        colnames(arg1) <- stringr::str_extract(string = colnames(arg1), pattern = rpat)
        check_var <- apply(X = arg1, MARGIN = 1, FUN = stats::var)
        which_novar <- which(check_var==0)
        if(length(which_novar)!=0) {
          arg1 <- arg1[-which_novar,]
        }
      }
      if(ncol(arg1)>nrow(arg1)) {
        arg1 <- t(arg1)
      }
      return(arg1)
    }
    drop_novar2 <- function(arg1) {
      if(ncol(arg1)>nrow(arg1)) {
        check_var <- apply(X = arg1, MARGIN = 2, FUN = stats::var)
        which_novar <- which(check_var==0)
        if(length(which_novar)!=0) {
          arg1 <- arg1[,-which_novar]
        }
      } else {
        check_var <- apply(X = arg1, MARGIN = 1, FUN = stats::var)
        which_novar <- which(check_var==0)
        if(length(which_novar)!=0) {
          arg1 <- arg1[-which_novar,]
        }
      }
      if(ncol(arg1)>nrow(arg1)) {
        arg1 <- t(arg1)
      }
      return(arg1)
    }

    if(pseudobulk_test_mode=="cluster_identity") {
      ct_spl <- vector("list", length = length(unique(seurat_object@meta.data[,cluster_column])))
      names(ct_spl) <- gsub("( |\\-|\\/)","_",unique(seurat_object@meta.data[,cluster_column]))
      meta_list <- vector("list", length = length(ct_spl)); names(meta_list) <- names(ct_spl)
      for(i in 1:length(ct_spl)) {
        ##
        ct_spl[[i]] <- as.data.frame(data.table::fread(paste0(capture_dir,"/temp_files/__pseudobulk_sum_counts_", names(ct_spl)[i], "__.csv"), check.names = FALSE, header = FALSE))
        if(file.exists(paste0(capture_dir,"/temp_files/__pseudobulk_sum_counts_", names(ct_spl)[i], "__.csv"))) {
          file.remove(paste0(capture_dir,"/temp_files/__pseudobulk_sum_counts_", names(ct_spl)[i], "__.csv"))
        }
        meta_list[[i]] <- read.csv(paste0(capture_dir,"/temp_files/__pseudobulk_obs_", names(ct_spl)[i], "__.csv"), check.names = FALSE, row.names = 1)
        meta_list[[i]]$cell_group <- row.names(meta_list[[i]])
        if(file.exists(paste0(capture_dir,"/temp_files/__pseudobulk_obs_", names(ct_spl)[i], "__.csv"))) {
          file.remove(paste0(capture_dir,"/temp_files/__pseudobulk_obs_", names(ct_spl)[i], "__.csv"))
        }
        obj_var <- read.csv(paste0(capture_dir,"/temp_files/__pseudobulk_var_", names(ct_spl)[i], "__.csv"))[,1]
        if(file.exists(paste0(capture_dir,"/temp_files/__pseudobulk_var_", names(ct_spl)[i], "__.csv"))) {
          file.remove(paste0(capture_dir,"/temp_files/__pseudobulk_var_", names(ct_spl)[i], "__.csv"))
        }
        row.names(ct_spl[[i]]) <- meta_list[[i]]$cell_group
        colnames(ct_spl[[i]]) <- obj_var
        ct_spl[[i]] <- drop_novar2(ct_spl[[i]])
        ##
        meta_list[[i]]$cluster_id <- meta_list[[i]][,cluster_column]
        meta_list[[i]]$sample_id <- meta_list[[i]][,pid_column]
        meta_list[[i]]$group_id <- factor(meta_list[[i]][,cluster_column], levels = c("in_cluster","out"))
        row.names(meta_list[[i]]) <- meta_list[[i]]$cell_group
        drop_indices <- which(!meta_list[[i]][,pid_column] %in% names(table(meta_list[[i]][,pid_column])[table(meta_list[[i]][,pid_column])==2]))
        if(length(drop_indices)!=0) {
          meta_list[[i]] <- meta_list[[i]][which(!1:nrow(meta_list[[i]]) %in% drop_indices),]
        }
        ct_spl[[i]] <- ct_spl[[i]][,which(colnames(ct_spl[[i]]) %in% meta_list[[i]]$cell_group)]
      }
    } else {
      ##
      ct_data <- as.data.frame(data.table::fread(paste0(capture_dir,"/temp_files/__pseudobulk_sum_counts__.csv"), check.names = FALSE, header = FALSE))
      if(file.exists(paste0(capture_dir,"/temp_files/__pseudobulk_sum_counts__.csv"))) {
        file.remove(paste0(capture_dir,"/temp_files/__pseudobulk_sum_counts__.csv"))
      }
      obj_obs <- read.csv(paste0(capture_dir,"/temp_files/__pseudobulk_obs__.csv"), check.names = FALSE, row.names = 1)
      obj_obs$cell_group <- row.names(obj_obs)
      if(file.exists(paste0(capture_dir,"/temp_files/__pseudobulk_obs__.csv"))) {
        file.remove(paste0(capture_dir,"/temp_files/__pseudobulk_obs__.csv"))
      }
      obj_var <- read.csv(paste0(capture_dir,"/temp_files/__pseudobulk_var__.csv"))[,1]
      if(file.exists(paste0(capture_dir,"/temp_files/__pseudobulk_var__.csv"))) {
        file.remove(paste0(capture_dir,"/temp_files/__pseudobulk_var__.csv"))
      }
      row.names(ct_data) <- obj_obs$cell_group
      colnames(ct_data) <- obj_var
      ##
      id_pattern <- paste0(gsub("-","\\\\-",paste0("(",paste0(unique(seurat_object@meta.data[,pid_column]), collapse = "|"),")")),"_")
      ct_spl <- split(x = ct_data, f = gsub(id_pattern,"",row.names(ct_data)))
      ct_spl <- lapply(X = ct_spl, FUN = drop_novar)
      meta_list <- vector("list", length = length(ct_spl)); names(meta_list) <- names(ct_spl)
      for(i in 1:length(meta_list)) {
        tmp_meta <- obj_obs
        # tmp_meta <- tmp_meta[which(tmp_meta[,pid_column] %in% colnames(ct_spl[[i]])),]
        tmp_meta <- tmp_meta[which(gsub(id_pattern,"",row.names(obj_obs)) %in% names(ct_spl)[i]),]
        if(mean(tmp_meta[,pid_column]==colnames(ct_spl[[i]]))!=1) {
          stop("pid mismatch (1)")
        }
        tmp_meta$cluster_id <- names(meta_list)[i]
        # tmp_meta$sample_id <- tmp_meta[,pid_column]
        tmp_meta$sample_id <- tmp_meta[,pid_column]
        tmp_meta$group_id <- factor(tmp_meta[,category_column], levels = test_categories)
        row.names(tmp_meta) <- tmp_meta$sample_id
        meta_list[[i]] <- tmp_meta
      }
    }

    deseq_input <- vector("list", length = length(ct_spl)); names(deseq_input) <- names(ct_spl)
    for(i in 1:length(deseq_input)) {
      deseq_input[[i]] <- list(ct_spl[[i]], meta_list[[i]])
    }

    # return(deseq_input) # check input is correct for "identity" testing

    use_adj_p <- TRUE # hard coding use adjusted p values; consider allowing use unadjusted for discovery
    padj_threshold <- 0.05 # hard coding 0.05; consider allowing to change
    do_deseq <- function(arg1, p_return_threshold = padj_threshold, stim = seu_conditions,
                       use_adj = use_adj_p) {
      # testing
      # arg1 <- deseq_input[[1]]
      # p_return_threshold = padj_threshold
      # stim = seu_conditions
      # use_adj = use_adj_p

      count_data <- arg1[[1]]
      meta_data <- arg1[[2]]

      test_cats_order <- levels(meta_data$group_id)

      count_data <- round(count_data)

      if(nrow(meta_data)==0) {
        fn_flag <- 0
        res <- "Not enough cells from either of the tested groups. Min cells is 10."
      } else {
        group_id <- levels(meta_data$group_id)
        num_g1 <- sum(meta_data$group_id==group_id[1])
        num_g2 <- sum(meta_data$group_id==group_id[2])
        if(any(num_g1<=1,num_g2<=1)) {
          fn_flag <- 0
          res <- paste0("[",num_g1,"] from group [",group_id[1],"] and [",num_g2,"] from group [",group_id[2],"] available for testing. Not enough to test. Min cells is 10.")
        } else {
          fn_flag <- 1
          dds <- DESeqDataSetFromMatrix(countData = count_data,
                                        colData = meta_data,
                                        design = ~ group_id)
          dds <- DESeq(dds)
          normalized_counts <- counts(dds, normalized = TRUE)
          contrast <- c("group_id", test_cats_order[1], test_cats_order[2])
          res <- results(dds,
                         contrast = contrast,
                         alpha = 0.05)
          res <- lfcShrink(dds,
                           contrast = contrast,
                           res=res,
                           type = "ashr")
          norm_ct <- DESeq2::counts(dds, normalized = TRUE)
          res <- data.frame(res); res <- res[!is.na(res$padj),]
          if(sum(is.na(res$padj))!=0) {
            res$padj <- p.adjust(p = res$pvalue, method = "fdr")
          }
          res$gene <- row.names(res)
          res$condition <- stim
          if(use_adj) {
            if(!is.null(p_return_threshold)) {
              res <- res[res$padj<=p_return_threshold,]
            }
          } else {
            if(!is.null(p_return_threshold)) {
              res <- res[res$pvalue<=p_return_threshold,]
            }
          }
        }
      }
      # if(nrow(res)!=0) {
      #   if(nrow(res)>40) {
      #     sig_genes <- res$gene[order(abs(res$log2FoldChange))][1:40]
      #   } else {
      #     sig_genes <- res$gene
      #   }
      #   sig_norm <- as.data.frame(norm_ct[which(row.names(norm_ct) %in% sig_genes),])
      #
      #   pheatm <- pheatmap(sig_norm,
      #                      cluster_rows = T,
      #                      show_rownames = T,
      #                      # annotation = data.frame(group_id = meta_data[, c("group_id")]),
      #                      annotation = meta_data[,c(cluster_column,"group_id")],
      #                      border_color = NA,
      #                      fontsize = 10,
      #                      scale = "row",
      #                      fontsize_row = 10,
      #                      height = 20)
      # } else {
      #   pheatm <- "no features plotted"
      # }
      if(fn_flag==1) {
        return(list(res = res, norm_ct = norm_ct, meta = meta_data))#, hm = pheatm))
      } else {
        return(list(res = res, meta = meta_data))#, hm = pheatm))
      }
    }
    deseq_out <- lapply(X = deseq_input, FUN = do_deseq, p_return_threshold = padj_threshold,
                        stim = seu_conditions, use_adj = use_adj_p)

    if(return_all_pseudobulk) {
      return(deseq_out)
    }

    drop_indices <- c()
    for(i in 1:length(deseq_out)) {
      if(nrow(deseq_out[[i]][[1]])!=0) {
        deseq_out[[i]][[1]]$cluster <- names(deseq_out)[i]
      } else {
        drop_indices <- append(drop_indices, i)
      }
    }
    if(length(drop_indices)==1) {
      deseq_out <- deseq_out[[-drop_indices]]
    } else if(length(drop_indices)>1) {
      deseq_out <- deseq_out[-drop_indices]
    }

    deseq_res <- lapply(X = deseq_out, FUN = function(arg1) return(arg1[[1]]))
    deseq_write <- do.call(rbind, deseq_res)
    return(deseq_write)
  }

  dge_outs <- vector("list", length = ifelse(is.null(conditions), 1, length(conditions)))
  if(is.null(conditions)) {
    names(dge_outs) <- "all"
  } else {
    names(dge_outs) <- conditions
  }
  dge_outs <- lapply(X = dge_outs, FUN = function(arg1, ann = annos){
    tmpv <- vector("list", length = ifelse(is.null(ann), 1, length(ann)))
    if(is.null(ann)) {
      names(tmpv) <- "all"
    } else {
      names(tmpv) <- ann
    }
    return(tmpv)
  })

  start_test <- Sys.time()
  for(i in 1:length(dge_outs)) {
    start_i <- Sys.time()
    if(names(dge_outs)[i]=="all") {
      subs1 <- seurat_object
    } else {
      seurat_object <- AddMetaData(object = seurat_object, metadata = seurat_object@meta.data[,condition_column], col.name = "c")
      subs1 <- subset(x = seurat_object, subset = c == names(dge_outs)[i])
    }
    for(j in 1:length(dge_outs[[i]])){
      start_j <- Sys.time()
      print(paste0("starting on [",names(dge_outs[[i]])[j],"] in [",names(dge_outs)[i],"] at ",Sys.time()))
      if(!is.null(annos)) {
        if(!is.null(test_cats)) {
          subs1 <- AddMetaData(object = subs1, metadata = subs1@meta.data[,cluster_column], col.name = "l")
          subs2 <- subset(x = subs1, subset = l == names(dge_outs[[i]])[j])
        } else {
          subs1 <- AddMetaData(object = subs1, metadata = ifelse(subs1@meta.data[,cluster_column]==names(dge_outs[[i]])[j], names(dge_outs[[i]])[j], "other"), col.name = "l")
          subs2 <- subs1
        }
      } else {
        subs2 <- subs1
      }
      if(!is.null(test_cats)) {
        ident1 <- test_categories[1]
        ident2 <- test_categories[2]
        Idents(subs2) <- subs2@meta.data[,category_column]
      } else {
        ident1 <- names(dge_outs[[i]])[j] # will be either "all" or a cluster
        if(is.null(ident1)) {
          stop("Nothing to test")
        }
        ident2 <- "other"
        Idents(subs2) <- ifelse(subs2@meta.data[,cluster_column]==ident1, ident1, "other")
      }
      print(paste0("for [",names(dge_outs[[i]])[j],"] in [",names(dge_outs)[i],"] ident.1 = ",ident1,", ident.2 = ",ident2))
      if(tolower(dge_method) == "mast") {
        suppressPackageStartupMessages({
          require(MAST)
          require(data.table)
          require(ggplot2)
        })
        test_idents <- c(ident1, ident2)
        print("dropping lowly expressed genes. Threshold >= 10% of cells.")
        keep_genes <- vector("list", length = length(test_idents)); names(keep_genes) <- test_idents
        num_unique_genes <- nrow(subs2)
        for(i in 1:2) {
          subs2_grp <- subset(x = subs2, cells = which(subs2@meta.data[,category_column]==test_idents[i]))
          bem <- subs2_grp@assays[[assay]]@layers$counts > 0
          percent_expr <- rowSums(bem) / ncol(bem); names(percent_expr) <- row.names(subs2)
          keep_genes[[i]] <- row.names(subs2)[percent_expr >= freq_expressed]
        }
        grp_genes <- unlist(keep_genes); grp_genes <- grp_genes[!is.na(grp_genes)]
        genes_to_keep <- names(table(grp_genes)[which(table(grp_genes)==2)])
        if(length(genes_to_keep)==0) {
          warning("no genes passed thresholding")
          next
        } else {
          print(paste0("number of lowly expressed genes dropped: ",num_unique_genes - length(genes_to_keep)," genes. ",length(genes_to_keep)," genes left for MAST testing."))
        }
        subs2 <- subset(x = subs2, features = genes_to_keep)
        subs2$category <- subs2@meta.data[,category_column]
        subs2$category <- ifelse(subs2$category==ident1,"Group1","Group2")
        Idents(subs2) <- subs2$category
        seu_as_sce <- as.SingleCellExperiment(subs2, assay = "RNA")
        seu_as_sce <- scuttle::logNormCounts(seu_as_sce, log = TRUE)
        cdr <- colSums(seu_as_sce@assays@data@listData[["logcounts"]]>0) # cellular detection rate; number of genes detected which is a well-known confounder (https://bioconductor.org/packages/release/bioc/vignettes/MAST/inst/doc/MAITAnalysis.html#22_Filtering)
        seu_as_sce$cngeneson <- scale(cdr)
        my_sca <- SceToSingleCellAssay(seu_as_sce, check_sanity = FALSE)
        if(!is.null(mast_lane)) {
          latv <- c("cngeneson",pid_column,mast_lane)
        } else {
          latv <- c("cngeneson",pid_column)
        }
        cond<-factor(colData(my_sca)$category)
        cond<-relevel(cond,"Group1")
        colData(my_sca)$category <- cond

        if (!is.null(pid_column)) {
          if (!pid_column %in% latv) {
            stop("Random effect variable (sample ID) should be included in latent variables! Specify sample ID using arg 'pid_column'")
          }
          latv <- latv[!latv %in% pid_column]
          fmla <- as.formula(object = paste0(" ~ category + ", paste(latv, collapse = "+"), glue::glue(" + (1|{pid_column})")))
          print(fmla)
          zlmCond <- MAST::zlm(formula = fmla,
                               sca = my_sca,
                               exprs_value = 'logcounts',
                               method="glmer",
                               ebayes=FALSE,
                               silent=T,
                               fitArgsD = list(nAGQ = 0),
                               strictConvergence = FALSE)
                               # fitArgsD = list(nAGQ = 1)
        } else {
          stop("Trying to run without mixed effect model but this is not supported. Sample ID must be specified with arg 'pid_column'")
        }
        summaryCond <- MAST::summary(object = zlmCond, doLRT = 'categoryGroup2')
        summaryDt <- summaryCond$datatable

        fcHurdle <- merge(summaryDt[contrast=='categoryGroup2' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                          summaryDt[contrast=='categoryGroup2' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients

        fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
        fcHurdleSig <- merge(fcHurdle[fdr<.05 & abs(coef)>fc_threshold], data.table::as.data.table(mcols(my_sca)), by='primerid')
        setorder(fcHurdleSig, fdr)
        mast_res <- as.data.frame(fcHurdleSig)

        if(nrow(fcHurdleSig)!=0) {
          lfcs <- MAST::logFC(zlmfit = zlmCond); lfc1 <- as.data.frame(lfcs$logFC[mast_res$primerid,]); lfc1$primerid <- row.names(lfc1)
          as.data.frame(lfcs$varLogFC[mast_res$primerid,])

          getlfcs <- as.data.frame(MAST::getLogFC(zlmfit = zlmCond)); getlfcs <- getlfcs[which(getlfcs$contrast!="cngeneson"),]
          row.names(getlfcs) <- getlfcs$primerid; getlfcs <- getlfcs[mast_res$primerid,]

          mast_res <- merge(x = mast_res, y = getlfcs, by = "primerid", sort = FALSE, all.x = TRUE)
          mast_res$contrast <- paste0(category_column,".",ident2)
          mast_res$cluster <- names(dge_outs[[i]])[j]

          entrez_to_plot <- fcHurdleSig[,primerid]
          flat_dat <- as(my_sca[entrez_to_plot,], 'data.table')
          flat_dat$category <- ifelse(flat_dat$category=="Group1",ident1,ident2)
          ggbase <- ggplot(flat_dat, aes(x=category, y=logcounts, color=category)) + geom_jitter() + facet_wrap(~primerid, scale='free_y')+ggtitle("DE Genes")
          mast_violins <- ggbase+geom_violin()

          dge_outs[[i]][[j]] <- list(res = mast_res, gene_plots = mast_violins)
        } else {
          dge_outs[[i]][[j]] <- list("no dge", list(fdr = 0.05, fc_threshold = fc_threshold))
        }

        ##### plotting stuff from mast; leaving for now; ggbase+geom_violin() may be particularly useful
        # entrez_to_plot <- fcHurdleSig[,primerid]
        # flat_dat <- as(my_sca[entrez_to_plot,], 'data.table')
        # ggbase <- ggplot(flat_dat, aes(x=category, y=logcounts, color=category)) + geom_jitter() + facet_wrap(~primerid, scale='free_y')+ggtitle("DE Genes")
        # ggbase+geom_violin()
        #
        # flat_dat[,lmPred:=lm(logcounts~cngeneson + category)$fitted, key=primerid]
        # ggbase + aes(x=cngeneson) + geom_line(aes(y=lmPred), lty=1) + xlab('Standardized Cellular Detection Rate')
        #
        # MM <- model.matrix(~category,unique(colData(my_sca)[,c("category"),drop=FALSE]))
        # rname_map <- colData(my_sca)$category; names(rname_map) <- colData(my_sca)$barcode
        # rownames(MM) <- rname_map[row.names(MM)]
        # predicted <- predict(zlmCond,modelmatrix=MM)
        #
        # predicted[, primerid:=as.character(primerid)]
        # predicted_sig <- merge(mcols(my_sca), predicted[primerid%in%entrez_to_plot], by='primerid')
        # predicted_sig <- as.data.table(predicted_sig)
        #
        # ggplot(predicted_sig)+aes(x=invlogit(etaD),y=muC,xse=seD,yse=seC,col=sample)+
        #   facet_wrap(~primerid,scales="free_y")+theme_linedraw()+
        #   geom_point(size=0.5)+scale_x_continuous("Proportion expression")+
        #   scale_y_continuous("Estimated Mean")+
        #   stat_ell(aes(x=etaD,y=muC),level=0.95, invert='x')
        #
        # # heatmap, exprs per cell, split by group
        # mat_to_plot <- assay(my_sca[entrez_to_plot,])
        # symbols_to_plot <- fcHurdleSig[,primerid]
        # rownames(mat_to_plot) <- symbols_to_plot
        # assay_as_matrix <- as.matrix(mat_to_plot)
        # NMF::aheatmap(assay_as_matrix,annCol=colData(my_sca)[,"category"],main="DE genes",
        #               col=rev(colorRampPalette(colors = RColorBrewer::brewer.pal(name="PiYG",n=10))(20)), labCol = NA)
        #
        # # gsea
        # boots <- bootVcov1(zlmCond, R = 50)
        # module <- "BTM"
        # min_gene_in_module <- 5
        # packageExt <- system.file("extdata", package='MAST')
        # module_file <- list.files(packageExt, pattern = module, full.names = TRUE)
        # gene_set <- getGmt(module_file)
        # gene_ids <- geneIds(gene_set)
        # gene_ids <- gene_ids[!names(gene_ids)%like%"TBA"&!names(gene_ids)%like%"B cell"]
        # sets_indices <- limma::ids2indices(gene_ids, mcols(sca)$symbolid)
        # # Only keep modules with at least min_gene_in_module
        # sets_indices <- sets_indices[sapply(sets_indices, length) >= min_gene_in_module]
        # gsea <- gseaAfterBoot(zlmCond, boots, sets_indices, CoefficientHypothesis("conditionStim"))
        # z_stat_comb <- summary(gsea, testType='normal')
        #####
      } else if(tolower(dge_method) == "wilcox") {
        suppressPackageStartupMessages({
          require(SeuratWrappers)
        })
        wilcox_res <- SeuratWrappers::RunPresto(object = subs2, assay = assay, ident.1 = ident1, ident.2 = ident2,
                                                test.use = "wilcox", only.pos = FALSE, min.pct = 0.01)
        # colnames(wilcox_res)[which(colnames(wilcox_res)=="pct.1")] <- gsub(" ","",paste0("pct.",ident1))
        # colnames(wilcox_res)[which(colnames(wilcox_res)=="pct.2")] <- gsub(" ","",paste0("pct.",ident2))
        wilcox_res$gene <- row.names(wilcox_res)
        wilcox_res$cluster <- ident1
        wilcox_res <- wilcox_res[wilcox_res$p_val_adj<0.05,]
        wilcox_res <- wilcox_res[which(abs(wilcox_res$avg_log2FC)>fc_threshold),]

        dge_outs[[i]][[j]] <- wilcox_res
      } else {
        stop("no supported 'dge_method' requested: supported algorithms are: 'wilcox', 'mast', or 'pseudobulk'")
      }
      print(paste0("[",names(dge_outs[[i]])[j],"] in [",names(dge_outs)[i],"] testing took ",round(as.numeric(difftime(Sys.time(), start_j, units = "mins")),2)," mins"))
    }
    print(paste0("[",names(dge_outs)[i],"] testing took ",round(as.numeric(difftime(Sys.time(), start_i, units = "mins")),2)," mins"))
  }
  print(paste0("total test time: ",round(as.numeric(difftime(Sys.time(), start_test, units = "hours")),3)," hours"))
  return(dge_outs)
}


seurat_correlate <- function(seurat_object, pid_column, condition = "Media", condition_column = "condition",
                             correlate_column = "time_to_progr", correlate_to_object, correlate_to_column,
                             correlate_to_pid_column, correlation_method = "spearman") {
  seurat_object <- ""
}
