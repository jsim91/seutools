heatmap_calculate <- function(seurat_obj, gene_set, set_name, clusters)
{
  require(ComplexHeatmap)
  require(FCSimple)
  require(ggplot2)

  # seurat_obj = sobj
  # gene_set = gene_set_list[[8]]
  # gene_set = gene_set_list[["Cytokine Signaling"]]
  # clusters = c(14,23,31)
  # set_name = "Cytokine Signaling"

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

tile_reduction <- function(seurat_obj, condition, cluster_numbers,
                           color_clusters, label_clusters = "all",
                           pt_alpha = 0.05, anno_text_size = 6, pt_size = 1, color_seed = 123,
                           outline_method = c("nudge","fontsize"), postfix_title_string = NA,
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
  # coords = py_rna_umap
  # condition = py_meta$Condition
  # cluster_numbers = type_cl
  # color_clusters = "all"
  # label_clusters = "all"
  # pt_alpha = 0.05
  # anno_text_size = 6
  # pt_size = 1
  # plot_order = c(1,2,3)
  # annotation_method = "repel"
  # color_seed = 123
  # override_color_aes = NA

  coords <- seurat_obj@reductions$

  plot_data <- data.frame(UMAP1 = coords$UMAP1, UMAP2 = coords$UMAP2,
                          cluster = cluster_numbers, condition = condition)
  unique_clus <- unique(plot_data$cluster)
  set.seed(color_seed)
  plot_data$cluster <- factor(x = plot_data$cluster, levels = sample(unique_clus,length(unique_clus),replace=F))
  xrange <- range(plot_data$UMAP1); yrange = range(plot_data$UMAP2)

  uclus <- unique(plot_data$cluster); uclus <- uclus[order(uclus,decreasing=F)]
  clusx <- rep(NA,length(uclus)); names(clusx) <- uclus; clusy <- clusx
  for(i in 1:length(clusx)) {
    clusx[i] <- median(plot_data$UMAP1[which(plot_data$cluster==names(clusx)[i])])
    clusy[i] <- median(plot_data$UMAP2[which(plot_data$cluster==names(clusx)[i])])
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
                       palpha = pt_alpha, anno_ts = anno_text_size, psize = pt_size,
                       plimx = xrange, plimy = yrange, amethod = annotation_method,
                       cseed = color_seed,text_omethod=outline_method, oca = override_color_aes,
                       flimx = force_xlim, flimy = force_ylim, pts = postfix_title_string,
                       fo = frameon)
  {
    # testing
    # input <- spl_data[[1]]
    # color_clus = color_clusters
    # xanno = clusx
    # yanno = clusy
    # palpha = pt_alpha
    # anno_ts = anno_text_size
    # psize = pt_size
    # plimx = xrange
    # plimy = yrange
    # amethod = annotation_method
    # cseed = color_seed
    # pts = postfix_title_string
    # oca = override_color_aes

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
        text_add <- text_add[-subs_rows,]

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
      plt <- ggplot(data = input, mapping = aes(x = UMAP1, y = UMAP2)) +
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
        plt <- ggplot(data = input, mapping = aes(x = UMAP1, y = UMAP2)) +
          geom_point(pch = 19, alpha = palpha, size = psize) + theme_void() +
          xlim(plimx) + ylim(plimy)
        plt <- plt + geom_point(data = foreground,
                                mapping = aes(x = UMAP1, y = UMAP2, color = cluster),
                                alpha = palpha, size = psize) +
          guides(color = guide_legend(override.aes = list(size = 6, alpha = 1))) +
          theme(legend.position = "none",
                legend.title = element_blank(),
                legend.text = element_text(size = 16))
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
                     mapping = aes(x = UMAP1, y = UMAP2, color = cluster),
                     alpha = palpha, size = psize) +
          guides(color = guide_legend(override.aes = list(size = 6, alpha = 1))) +
          theme(legend.position = "none",
                legend.title = element_blank(),
                legend.text = element_text(size = 16))
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
    # if(length(xanno)!=0) {
    #   plt <- plt + annotate("shadowtext", x = xanno, y = yanno, label = names(xanno),
    #                         aes(color = names(xanno)), size = anno_ts)
    # }
    if(amethod[1]!="none") {
      if(nrow(text_add)!=0) {
        if(amethod[1]=="shadowtext") {
          plt <- plt + annotate("shadowtext", x = color_text_add$UMAP1, y = color_text_add$UMAP2,
                                label = color_text_add$cluster, size = anno_ts)
        } else if(amethod[1]=="repel") {
          plt <- plt + ggrepel::geom_text_repel(data = color_text_add, mapping = aes(x = UMAP1, y = UMAP2, color = cluster, label = cluster),
                                                size = anno_ts, bg.color = "black", bg.r = 0.05, seed = 123)
        } else if(amethod[1]=="text") {
          plt <- plt + geom_text(data = color_text_add, mapping = aes(x = UMAP1, y = UMAP2, color = cluster, label = cluster),
                                 size = anno_ts)
        }
      } else if(nrow(color_text_add)!=0) {
        if(amethod[1]=="shadowtext") {
          plt <- plt + annotate("shadowtext", x = color_text_add$UMAP1, y = color_text_add$UMAP2,
                                label = color_text_add$cluster, size = anno_ts)
        } else if(amethod[1]=="repel") {
          plt <- plt + ggrepel::geom_text_repel(data = color_text_add, mapping = aes(x = UMAP1, y = UMAP2, color = cluster, label = cluster),
                                                size = anno_ts, bg.color = "black", bg.r = 0.05, seed = 123)
        } else if(amethod[1]=="text") {
          if(text_omethod[1]=="nudge") {
            plt <- plt + geom_text(data = color_text_add, mapping = aes(x = UMAP1, y = UMAP2, label = cluster),
                                   size = anno_ts, color = "black", nudge_x = 0.005, nudge_y = 0.005) +
              plt <- plt + geom_text(data = color_text_add, mapping = aes(x = UMAP1, y = UMAP2, color = cluster, label = cluster),
                                     size = anno_ts)
          } else if(text_omethod[1]=="fontsize") {
            plt <- plt + geom_text(data = color_text_add, mapping = aes(x = UMAP1, y = UMAP2, label = cluster),
                                   size = anno_ts*1.05, color = "black") +
              plt <- plt + geom_text(data = color_text_add, mapping = aes(x = UMAP1, y = UMAP2, color = cluster, label = cluster),
                                     size = anno_ts)
          }
        }
      }
    } else {
      plt <- plt + guides(color = guide_legend(override.aes = list(alpha = 1, stroke = 0.1, size = ifelse(psize<5,5,psize)))) +
        theme(legend.position = "bottom", legend.text = element_text(size = 16), legend.title = element_blank())
    }
    plt <- plt + ggtitle(ifelse(!is.na(pts), paste0(input$condition[1]," - ",pts), input$condition[1])) +
      theme(axis.text = element_blank(),
            axis.title = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 22))
    if(fo) {
      plt <- plt + theme_bw() + theme(axis.ticks = element_blank(),
                                      axis.title = element_blank(),
                                      axis.text = element_blank(),
                                      legend.title = element_blank(),
                                      panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank())
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

test_by_igra <- function(pid, condition, igra_status, test_clusters,
                         sum_clusters = FALSE,
                         test_method = c("sequential","cumulative"),
                         clusters, background_condition = "Media",
                         bar_colors = c("GRV" = "#2e7d32", "MTB300" = "#2962ff"),
                         y_axis_subset = "T/NK")
{
  require(ggplot2)
  require(ggpubr)
  require(scales)

  # testing
  # pid <- meta$soc_b_id
  # condition <- meta$soc_condition_id
  # igra_status <- data.frame(key = row.names(igra), value = igra$`IGRA Status`)
  # test_clusters <- c(14,23,31,34,21,28)
  # sum_clusters = TRUE
  # test_method <- "sequential"
  # clusters <- meta$seurat_clusters
  # background_condition <- "Media"
  # bar_colors = c("GRV" = "#2e7d32", "MTB300" = "#2962ff")
  # y_axis_subset = "T/NK"

  # testing
  # pid = meta$soc_b_id; condition = meta$soc_condition_id;
  # igra_status = data.frame(key = row.names(igra), value = igra$`IGRA Status`);
  # test_clusters = c(14,23,31); sum_clusters = TRUE;
  # test_method = c("sequential","cumulative");
  # clusters = meta$seurat_clusters; background_condition = "Media";
  # bar_colors = c("gRV" = "#2e7d32", "MTB300" = "#2962ff");
  # y_axis_subset = "T/NK"

  igra_short <- igra_status[which(igra_status$key %in% unique(pid)),]
  if(nrow(igra_short)==0) {
    stop("igra pids do not match or are not found in 'pid'")
  }

  upid <- unique(pid); upid <- upid[order(upid)]
  uclus <- unique(clusters); uclus <- uclus[order(uclus)]
  ucond <- unique(condition)

  freq_mat <- matrix(data = NA, nrow = length(upid), ncol = length(uclus))
  row.names(freq_mat) <- upid; colnames(freq_mat) <- uclus

  freq_by_condition <- vector("list", length = length(ucond))
  names(freq_by_condition) <- ucond

  query_data <- data.frame(pid = pid, condition = condition, cluster = clusters)

  for(i in 1:length(freq_by_condition)) {
    tmp_mat <- freq_mat
    tmp_query_data <- query_data[which(query_data$condition==names(freq_by_condition)[i]),]
    for(j in 1:nrow(tmp_mat)) {
      tmp_nums <- tmp_query_data$cluster[which(tmp_query_data$pid==row.names(tmp_mat)[j])]
      for(k in 1:ncol(tmp_mat)) {
        tmp_mat[j,k] <- mean(tmp_nums==as.numeric(colnames(tmp_mat)[k]))*100
      }
    }
    freq_by_condition[[i]] <- tmp_mat
  }
  which_background <- which(names(freq_by_condition)==background_condition)
  if(length(which_background)!=1) {
    stop("could not find one background condition from which to subtract stims")
  } else {
    background_data <- freq_by_condition[[which_background]]
    stim_data <- freq_by_condition[-which_background]
  }
  stim_data_subtr <- stim_data
  for(i in 1:length(stim_data_subtr)) {
    stim_data_subtr[[i]] <- stim_data_subtr[[i]] - background_data
  }
  for(i in 1:length(stim_data_subtr)) {
    stim_data_subtr[[i]] <- stim_data_subtr[[i]][,which(colnames(stim_data_subtr[[i]]) %in% as.character(test_clusters))]
  }
  if(sum_clusters) {
    for(i in 1:length(stim_data_subtr)) {
      new_cname <- paste0(colnames(stim_data_subtr[[i]]), collapse = "+")
      stim_data_subtr[[i]] <- data.frame(var1 = rowSums(stim_data_subtr[[i]]))
      colnames(stim_data_subtr[[i]]) <- new_cname
    }
  }

  # cluster_list <- vector("list", length = length(uclus)); names(cluster_list) <- uclus
  cluster_list <- vector("list", length = ncol(stim_data_subtr[[1]])); names(cluster_list) <- colnames(stim_data_subtr)
  tmp_df <- merge(x = data.frame(pid = row.names(freq_mat)),
                  y = igra_status, by.x = "pid", by.y = "key")

  for(i in 1:length(cluster_list)) {
    for(j in 1:length(stim_data_subtr)) {
      # if(j==1) {
      # tmp_df <- data.frame(var1 = stim_data_subtr[[j]][,i])
      if(j==1) {
        tmp_df2 <- merge(x = tmp_df, y = data.frame(mergevar = row.names(stim_data_subtr[[j]]),
                                                    freqval = stim_data_subtr[[j]][,i]), by.x = "pid", by.y = "mergevar")
        colnames(tmp_df2)[ncol(tmp_df2)] <- paste0(names(stim_data_subtr)[j],"_",colnames(stim_data_subtr[[j]])[i])
      } else {
        tmp_df2 <- merge(x = tmp_df2, y = data.frame(mergevar = row.names(stim_data_subtr[[j]]),
                                                     freqval = stim_data_subtr[[j]][,i]), by.x = "pid", by.y = "mergevar")
        colnames(tmp_df2)[ncol(tmp_df2)] <- paste0(names(stim_data_subtr)[j],"_",colnames(stim_data_subtr[[j]])[i])
      }

      # } else {
      # tmp_df <- cbind(tmp_df,data.frame(tmp_var = stim_data_subtr[[j]][,i]))
      # colnames(tmp_df)[ncol(tmp_df)] <- paste0('var',j)
      # }
    }
    # cluster_list[[i]] <- cbind(tmp_df, data.frame(cluster = rep(names(cluster_list)[i],nrow(stim_data_subtr[[1]])),
    #                                               pid = row.names(stim_data_subtr[[1]])))
    if(mean(tmp_df2$value=="Unknown")!=0) {
      tmp_df2 <- tmp_df2[-which(tmp_df2$value=="Unknown"),]
    }
    cluster_list[[i]] <- tmp_df2
    # colnames(cluster_list[[i]])[1:length(stim_data_subtr)] <- names(stim_data_subtr)
  }

  iterate_boxplot <- function(arg1, b_colors = bar_colors,
                              should_connect = FALSE,
                              # anchor_stim = backround_condition,
                              tp = tile_plots,
                              # dp = data_paired,
                              yax_lab = y_axis_subset)
  {
    require(reshape2)
    require(combinat)

    # testing
    # arg1 = cluster_list[[1]]; b_colors = bar_colors;
    # should_connect = FALSE;
    # # anchor_stim = backround_condition;
    # tp = TRUE
    # # dp = data_paired;
    # yax_lab = y_axis_subset

    arg1$group <- as.character(1:nrow(arg1))
    arg_melt <- reshape2::melt(data = arg1, value.name = "cluster")
    # colnames(arg_melt) <- c("cluster","pid","group","condition","frequency")
    colnames(arg_melt) <- c("pid","igra","group","condition","frequency")
    # compare_these <- combinat::permn(unique(arg_melt$condition))
    # for(i in 1:length(compare_these)) {
    #   compare_these[[i]] <- as.character(compare_these[[i]][1:2])
    # }
    # find_dupls <- rep(NA,length=length(compare_these))
    # for(i in 1:length(find_dupls)) {
    #   find_dupls[i] <- paste0(compare_these[[i]][order(compare_these[[i]])],collapse="")
    # }
    # rm_elements <- which(duplicated(find_dupls))
    # if(length(rm_elements)>0) {
    #   compare_these <- compare_these[-rm_elements]
    # }
    arg_melt$compare_group <- paste0(arg_melt$igra," ",gsub("_.+$","",arg_melt$condition))
    arg_melt$compare_group <- factor(x = arg_melt$compare_group,
                                     levels = c("Negative MTB300","Positive MTB300","Negative gRV","Positive gRV"))
    # compare_these <- c("Positive","Negative")
    compare_these <- list(c("Negative MTB300","Positive MTB300"), c("Negative gRV","Positive gRV"))

    plt <- ggplot(data = arg_melt, mapping = aes(x = compare_group, y = frequency, group = compare_group)) +
      geom_boxplot(fill = "#bfbfbf", lwd = 0.5, alpha = 0.4, width = 0.45) +
      geom_dotplot(data = arg_melt, fill = "grey", color = "black", stroke = 1,
                   binaxis = "y", stackdir = "center", position = "dodge", binpositions="all") + #, binwidth = 0.5) +
      scale_x_discrete(labels = c("Negative\nMTB300","Positive\nMTB300","Negative\ngRV","Positive\ngRV")) +
      stat_compare_means(paired = FALSE, method = "wilcox", comparisons = compare_these, size = 6) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
      # scale_fill_manual(values = b_colodata_pairedrs) +
      # scale_color_manual(values = b_colors) +
      ylab(paste0("% of ",yax_lab)) +
      # ggtitle(paste0("cluster ",internal_arg1[,"cluster"][1])) +
      ggtitle(paste0("cluster ",gsub("(MTB300|gRV)_", "", arg_melt$condition))) +
      theme_minimal() +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_text(size = 17, color = "black"),
            axis.text.x = element_text(size = 15, face = "bold", color = "black"),
            axis.text.y = element_text(size = 15),
            plot.title = element_text(size = 19, hjust = 0.5, face = "bold"),
            legend.position = "none")
    return(plt)
  }
  test_plots <- lapply(X = cluster_list, FUN = iterate_boxplot, b_colors = bin_colors)
  return(test_plots)
}

mean_count_hm <- function(exprs_matrix, cluster_numbers, include_clusters, gene_set, pid,
                          low_mid_high_cols = c("#DA29D9","black","#fff176"), show_other = TRUE,
                          scale_per_gene = TRUE, split_by_pid = TRUE, cluster_rows = FALSE,
                          cluster_annotation_ref = "none", cluster_annotation_color = NA,
                          prefix_cluster = TRUE, auto_order_genes = FALSE, get_legend = FALSE)
{
  require(ComplexHeatmap)
  require(grid)
  library(circlize)

  # testing
  # exprs_matrix = readRDS("J:/10x/TB7233/TNK_server_exports/TNK_RNA_exprs_data.rds")$rna_data;
  # cluster_numbers = meta$seurat_clusters
  # include_clusters = c(14,23,31,34,21,28,25)

  # exprs_matrix = rna_exprs_data
  # cluster_numbers = meta$seurat_clusters
  # include_clusters = c(14,23,31,34,21,28)
  # gene_set = unique(prio_genes)
  # pid = meta$ptid
  # low_mid_high_cols = c("#313695","#FFFFBF","#A50026")
  # scale_per_gene = TRUE
  # split_by_pid = TRUE
  # show_other = FALSE
  # cluster_rows = FALSE
  # prefix_cluster = TRUE
  # cluster_annotation_ref = list('CD4 T cell' = c(14,23,31), 'NK' = c(21,28), 'MAIT/gdT' = 34)
  # cluster_annotation_color = c("CD4 T cell" = "red", "MAIT/gdT" = "orange", "NK" = "blue")
  # auto_order_genes = TRUE
  # get_legend = FALSE

  # pathways_2 <- as.data.frame(readxl::read_xlsx(path = "J:/10x/TB7233/TNK7233_outs/pathways_and_genes.xlsx",
  #                                               sheet = 2))
  # gene_set <- strsplit(x = pathways_2$geneID[grep("R-HSA-168898",pathways_2$ID)], split = "/")[[1]]

  # gene_set = unique(prio_genes)
  # pid = meta$soc_b_id
  # low_mid_high_cols = c("#DA29D9","#000000","#fff176");
  # low_mid_high_cols = c(PurpleAndYellow()[1],"black",PurpleAndYellow()[length(PurpleAndYellow())]);
  # scale_per_gene = TRUE
  # split_by_pid = TRUE
  # cluster_annotation_ref <- list('CD4' = c(14,23,31), 'NK' = c(21,28), 'MAIT/gdT' = 34)
  # cluster_annotation_color <- c("CD4" = "red", "NK" = "blue", "MAIT/gdT" = "orange")
  # show_other = FALSE; prefix_cluster = FALSE

  # testing
  # exprs_matrix = readRDS("./TNK_server_exports/TNK_RNA_exprs_data.rds")$rna_data;
  # cluster_numbers =meta$seurat_clusters; include_clusters = "all";
  # gene_set = unique(prio_genes); scale_per_gene = TRUE

  # testing no cluster
  # exprs_matrix = readRDS("./TNK_server_exports/TNK_RNA_exprs_data.rds")$rna_data;
  # cluster_numbers = meta$seurat_clusters;
  # include_clusters = hm_clusters;
  # gene_set = top_genes_hm; pid = meta$soc_b_id;
  # low_mid_high_cols = c(Seurat::PurpleAndYellow()[1],"black",
  #                       Seurat::PurpleAndYellow()[length(Seurat::PurpleAndYellow())]);
  # scale_per_gene = TRUE; split_by_pid = TRUE; cluster_rows = FALSE
  # cluster_annotation_ref <- list('CD4' = c(14,23,31), 'NK' = c(21,28), 'MAIT/gdT' = 34)
  # cluster_annotation_color <- c("CD4" = "red", "NK" = "blue", "MAIT/gdT" = "orange")
  # show_other = FALSE; prefix_cluster = FALSE

  # testing
  # exprs_matrix = hm_data
  # cluster_numbers = hm_clusters
  # include_clusters <- gsub(pattern = "_", replacement = "-", x = names(rename_clus))
  # gene_set = unique(prio_genes); pid = hm_pid
  # low_mid_high_cols = c("#313695","#FFFFBF","#A50026")
  # scale_per_gene = TRUE; split_by_pid = TRUE; show_other = FALSE
  # cluster_rows = FALSE; prefix_cluster = TRUE
  # cluster_annotation_ref = list('CD4 T cell' = c("CD4-1", "CD4-2", "CD4-3"), 'NK' = c("NK-1","NK-2"), 'MAIT/gdT' = "MAIT/gdT")
  # cluster_annotation_color = c("CD4 T cell" = "red", "MAIT/gdT" = "orange", "NK" = "blue")
  # auto_order_genes = TRUE; get_legend = FALSE

  # include_clusters <- include_clusters[order(include_clusters)]
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
      if(prefix_cluster) {
        hm_groups <- paste0("cluster ",hm_groups)
        hm_groups <- c(hm_groups, rep(x = "other", times=length(upid)))
        hm_groups <- factor(hm_groups,levels=c(paste0("cluster ",include_clusters),"other"))
        if(!show_other) {
          other_ind <- which(hm_groups=="other")
          if(length(other_ind)!=0) {
            large_mat <- large_mat[,-other_ind]
            hm_groups <- hm_groups[-other_ind]
          }
        }
      } else {
        hm_groups <- c(hm_groups, rep(x = "other", times=length(upid)))
        hm_groups <- factor(hm_groups,levels=c(include_clusters,"other"))
        if(!show_other) {
          other_ind <- which(hm_groups=="other")
          if(length(other_ind)!=0) {
            large_mat <- large_mat[,-other_ind]
            hm_groups <- hm_groups[-other_ind]
          }
        }
      }
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
      if(prefix_cluster) {
        hm_groups <- paste0("cluster ",hm_groups)
        hm_groups <- factor(hm_groups,levels=paste0("cluster ",include_clusters))
      } else {
        hm_groups <- factor(hm_groups,levels=include_clusters)
      }
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
  if(split_by_pid) {
    col_fun = colorRamp2(breaks = c(min(range(large_mat)), mean(range(large_mat)), max(range(large_mat))),
                         colors = low_mid_high_cols)
    if(unlist(cluster_annotation_ref)[1]!="none") {
      type_anno <- rep(NA,length(hm_groups))
      for(i in 1:length(cluster_annotation_ref)) {
        if(prefix_cluster) {
          anno_pos_match <- which(hm_groups %in% paste0("cluster ",cluster_annotation_ref[[i]]))
        } else {
          anno_pos_match <- which(hm_groups %in% cluster_annotation_ref[[i]])
        }
        type_anno[anno_pos_match] <- names(cluster_annotation_ref)[i]
      }

      ha <- HeatmapAnnotation(bar = factor(type_anno, levels = names(cluster_annotation_color)),
                              col = list(bar = cluster_annotation_color), show_legend = FALSE,
                              show_annotation_name = FALSE, which = "column", simple_anno_size = unit(3, "mm"),
                              annotation_legend_param = list(legend_height=unit(3,"cm"),
                                                             grid_width=unit(0.6,"cm"),title_position="topleft",
                                                             labels_gp=gpar(fontsize=16),
                                                             title_gp=gpar(fontsize=0)))
      outhm <- ComplexHeatmap::Heatmap(matrix = large_mat, name = "GEX", column_split = hm_groups, # ComplexHeatmap doc 2.7.2 split single heatmap
                                       na_col = "black", col = col_fun,
                                       heatmap_legend_param=list(at=round(c(min(range(large_mat)),max(range(large_mat))),2),
                                                                 labels = c("low", "high"),
                                                                 legend_height=unit(3,"cm"),
                                                                 grid_width=unit(0.6,"cm"),title_position="topleft",
                                                                 labels_gp=gpar(fontsize=14),
                                                                 title_gp=gpar(fontsize=0,fontface="bold")),
                                       show_row_dend = FALSE, show_column_dend = FALSE, cluster_rows = cluster_rows,
                                       cluster_columns = FALSE, top_annotation = ha, show_heatmap_legend = FALSE,
                                       column_title_gp = gpar(fontsize=18), column_title_side = "top",
                                       show_column_names = FALSE, row_names_side = "left", cluster_column_slices = FALSE)
    } else {
      outhm <- ComplexHeatmap::Heatmap(matrix = large_mat, name = "GEX", column_split = hm_groups, # ComplexHeatmap doc 2.7.2 split single heatmap
                                       na_col = "black", col = col_fun,
                                       heatmap_legend_param=list(at=round(c(min(range(large_mat)),max(range(large_mat))),2),
                                                                 labels = c("low", "high"),
                                                                 legend_height=unit(3,"cm"),
                                                                 grid_width=unit(0.6,"cm"),title_position="topleft",
                                                                 labels_gp=gpar(fontsize=14),
                                                                 title_gp=gpar(fontsize=0,fontface="bold")),
                                       show_row_dend = FALSE, show_column_dend = FALSE, cluster_rows = cluster_rows,
                                       cluster_columns = FALSE, show_heatmap_legend = FALSE,
                                       column_title_gp = gpar(fontsize=18), column_title_side = "top",
                                       show_column_names = FALSE, row_names_side = "left", cluster_column_slices = FALSE)
    }
    if(get_legend) {
      scaled_leg <- Legend(at=round(c(min(range(large_mat)),max(range(large_mat))),2),
                           labels = c("low", "high"), col_fun = col_fun,
                           legend_height=unit(3,"cm"),
                           grid_width=unit(0.6,"cm"),title_position="topleft",
                           labels_gp=gpar(fontsize=14),
                           title_gp=gpar(fontsize=0,fontface="bold"))
      discrete_leg <- Legend(at = names(cluster_annotation_color), title = "none",
                             title_position="topleft", labels_gp=gpar(fontsize=14),
                             # title_gp=gpar(fontsize=14,fontface="bold"),
                             title_gp = gpar(fontsize=0, alpha = 0),
                             legend_gp = gpar(fill = cluster_annotation_color),
                             grid_height = unit(7, "mm"), grid_width = unit(0.6,"cm"))
      return(list(discrete_leg,scaled_leg))
    } else {
      return(list(hm = outhm,
                  tile_data = large_mat))
    }
  } else {
    outhm <- ComplexHeatmap::Heatmap(name = "feature\nmean\ncount", matrix = hm_matrix, cluster_rows = cluster_rows,
                                     heatmap_legend_param=list(at=round(c(min(range(hm_matrix)),max(range(hm_matrix))),2),
                                                               legend_height=unit(3,"cm"),
                                                               grid_width=unit(0.6,"cm"),title_position="topleft",
                                                               labels_gp=gpar(fontsize=14),
                                                               title_gp=gpar(fontsize=15,fontface="bold")),
                                     show_row_dend = FALSE, row_names_gp=gpar(fontsize=13,fontface="bold"),
                                     column_names_gp = gpar(fontsize=18,fontface="bold"), column_names_side = "top")
    if(get_legend) {
      scaled_leg <- Legend(at=round(c(min(range(large_mat)),max(range(large_mat))),2),
                           labels = c("low", "high"), col_fun = col_fun,
                           legend_height=unit(3,"cm"),
                           grid_width=unit(0.6,"cm"),title_position="topleft",
                           labels_gp=gpar(fontsize=14),
                           title_gp=gpar(fontsize=0,fontface="bold"))
      return(scaled_leg)
    } else {
      return(list(hm = outhm,
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
                                 shape_key = NULL)
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
                              shp = shape_key)
  {
    require(reshape2)
    require(combinat)
    require(ggplot2)
    require(ggpubr)

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
        pl <- pl + geom_point(mapping = aes(fill = variable, shape = shape_group), size = 4, color = "black") +
          scale_shape_manual(values = 21:25)
      } else {
        pl <- pl + geom_point(mapping = aes(fill = variable), pch = 21, size = 3, color = "black")
      }
    } else {
      pl <- pl + geom_dotplot(data = arg_melt, aes(fill = variable), color = "black", stroke = 1.1,
                              binaxis = "y", stackdir = "center", position = "dodge", binpositions="all")
    }
    if(dp) {
      pl <- pl +
        stat_compare_means(paired = TRUE, method = "wilcox", comparisons = compare_these, size = 7.6*tsf)
    } else {
      pl <- pl +
        stat_compare_means(paired = FALSE, method = "wilcox", comparisons = compare_these, size = 7.6*tsf)
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
      ggtitle(paste0("cluster ",arg_melt[,"cluster"][1])) +
      theme_minimal() +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_text(size = 23*tsf, color = "black"),
            axis.text.x = element_text(size = 22*tsf, face = "bold", color = "black"),
            axis.text.y = element_text(size = 20*tsf),
            plot.title = element_text(size = 25*tsf, hjust = 0.5, face = "bold"))
    if(is.null(shp[1])) {
      pl <- pl + theme(legend.position = "none")
    }
    return(pl)
  }
  test_plots <- lapply(X = cluster_list, FUN = iterate_boxplot, b_colors = bin_colors, dp = data_paired)
  return(test_plots)
}

seurat_size_bar <- function(seurat_object, pid_column = "pid", condition_column = "condition", cluster_column = "cell_type")
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

  plot_data <- data.frame(pid = seu@meta.data[,pid_column],
                          condition = seu@meta.data[,condition_column],
                          cluster = seu@meta.data[,cluster_column])
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
    # geom_bar(aes(fill = PID), position = "dodge", stat = "identity") +
    geom_bar(aes(fill = PID), stat = "identity") +
    theme_minimal() +
    ylab("number of cells") +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = 12),
          axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5),
          # axis.title = element_text(size = 15, face = "bold"))
          axis.title = element_blank())

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
                                     draw.ulim = F, draw.llim = F, ticks.linewidth = 0.5)) +
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


umap_by_value <- function(coords = meta[,c("UMAP1","UMAP2")], values, feature = "% TCR coverage",
                          tex = 1, point_alpha = 0.25, point_cex = 0.7, return_legend = FALSE,
                          plot_w_legend = TRUE, legx = 0.5, legy = 0.9, legend_orientation = "horizontal",
                          leg_height = unit(0.6,"cm"), leg_width = unit(3,"cm"))
{
  require(ggplot2)
  require(ggpubr)
  require(shadowtext)
  require(ggrastr)
  require(viridis)
  # testing
  # coords = meta_join[,c("UMAP1","UMAP2")]
  # values = coverage*100
  # feature = "% TCR coverage"
  # tex = 1
  # point_alpha = 0.05
  # point_cex = 0.3

  indata <- cbind(coords, data.frame(value = values))
  # subset for testing
  # indata <- indata[sample(x=1:nrow(indata),size=100000,replace=F),]

  # col_fun = colorRamp2(breaks = seq(from=floor(min(range(indata$value))),
  #                                   to=ceiling(max(range(indata$value))),
  #                                   length.out=5),
  #                      colors = viridis(n = 5))

  plt <- ggplot(data = indata, mapping = aes(x=UMAP1,y=UMAP2)) +
    geom_point_rast(aes(color=value),pch=19, alpha = point_alpha, cex = point_cex) +
    scale_color_viridis(option = "D", name = feature) +
    # guides(color = guide_colorbar(title.position = "top", title.hjust = 0.5, title = feature)) +
    guides(color = guide_colorbar(title.position="top", title.hjust = 0.5, title.vjust = 0.5,
                                  title = feature, barwidth = grid::unit(x = 0.44, units = "npc"),
                                  barheight = grid::unit(x = 6, units = "mm"), direction = "horizontal")) +
    theme_void() +
    theme(plot.title = element_blank(),
          legend.position = c(legx, legy),
          # legend.text = element_text(angle=45, size = 22*tex),
          legend.text = element_text(size = 22*tex),
          legend.title = element_text(size = 24*tex))
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
                    legend_height = leg_height, title = "% TCR coverage",
                    legend_width = leg_width, title_position = "topcenter",
                    labels_gp = gpar(fontsize=14),
                    title_gp = gpar(fontsize = 15, fontface = "bold"))
    } else if(legend_orientation=="vertical") {
      leg <- Legend(at = leg_lab,
                    labels = leg_lab, col_fun = col_fun,
                    legend_height = leg_height, title = "% TCR coverage",
                    grid_width = leg_width, title_position = "topleft",
                    labels_gp = gpar(fontsize=14),
                    title_gp = gpar(fontsize=15, fontface="bold"))
    }
    return(leg)
  } else {
    return(plt)
  }
}


feature_violin <- function(feature_counts, plot_features,
                           meta = py_meta,  meta_cat_column = "leiden2", # categorical column
                           plot_cat_features = "all", text_expansion = 1, nudge_nonzero = 0.35,
                           y_limit_expansion_factor = 0.5, condition = "Media", condition_cat = "Condition") {
  # testing
  # feature_counts = rna_ct
  # plot_features = clec_genes
  # meta = py_meta
  # meta_cat_column = "leiden2" # categorical column
  # plot_cat_features = "all"
  # text_expansion = 1
  # nudge_nonzero = 0.35
  # y_limit_expansion_factor = 0.5
  # condition = "Media"
  # condition_cat = "Condition"

  require(ggplot2)
  # require(dplyr)

  spl_mat2 <- vector("list", length = length(plot_features)); names(spl_mat2) <- plot_features
  for(i in 1:length(plot_features)) {
    spl_mat2[[i]] <- feature_counts[which(row.names(feature_counts)==plot_features[i]),]
  }
  spl_mat <- vector("list",length(spl_mat2)); names(spl_mat) <- names(spl_mat2)
  for(i in 1:length(spl_mat2)) {
    # spl_mat[[i]] <- data.frame(ct = spl_mat2[[i]], cat = factor(meta[,meta_cat_column]),
    spl_mat[[i]] <- data.frame(ct = spl_mat2[[i]], cat = meta[,meta_cat_column],
                               condition = meta[,condition_cat])
    if(!is.na(condition[1])) {
      spl_mat[[i]] <- spl_mat[[i]][which(spl_mat[[i]]$condition %in% condition),]
    }
    if(plot_cat_features[1]!="all") {
      spl_mat[[i]] <- spl_mat[[i]][which(spl_mat[[i]]$cat %in% plot_cat_features),]
      # spl_mat[[i]]$cat <- factor(as.character(spl_mat[[i]]$cat))
    }
    # cat_levels <- levels(spl_mat[[i]]$cat)
    # spl_levels <- split(x = cat_levels, f = gsub(pattern = "[0-9]+", replacement = "", x = cat_levels))
    # spl_levels <- lapply(X = spl_levels, FUN = function(arg1){
    #   get_nums <- as.numeric(stringr::str_extract(string = arg1, pattern = "[0-9]+"))
    #   return(arg1[order(get_nums)])
    # })
    # levels(spl_mat[[i]]$cat) <- unlist(spl_levels)
    spl_mat[[i]]$gene <- names(spl_mat)[i]
  }

  violin_internal <- function(indata, texp = text_expansion,
                              nudge_nz = nudge_nonzero, yle = y_limit_expansion_factor) {
    # testing
    # indata <- spl_mat[[1]]
    # texp = text_expansion
    # nudge_nz = nudge_nonzero
    # yle = y_limit_expansion_factor

    freqs <- data.frame(cat = unique(indata$cat), freq = rep(NA,length(unique(indata$cat))))
    for(i in 1:nrow(freqs)) {
      freqs$freq[i] <- round(mean(indata$ct[which(indata$cat==freqs$cat[i])]!=0)*100,3)
    }
    indata_jitter <- indata[which(indata$ct!=0),]

    expand_out <- max(indata$ct)*(yle/max(indata$ct) + 1)

    plt <- ggplot(data = indata, aes(x = cat, y = ct, fill = cat)) +
      geom_violin(scale = "width", trim = TRUE, alpha = 0.7) +
      geom_jitter(data = indata_jitter, width = 0.2, height = 0, size = 1, alpha = 0.5) +
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

feature_violin_test <- function(feature_counts, plot_features, meta, meta_cat_column = "leiden2", # categorical column
                                plot_cat_features = c("m0","m1","m2","m3","m4","m5"), # or "all"
                                text_expansion = 1, condition = "Media", condition_cat = "Condition",
                                test_cat = "progr_contr", stat_compares = c(), nudge_nonzero = 0.35,
                                y_limit_expansion_factor = 0.5, condition_in_title = FALSE,
                                stat_text_offset = 1) {
  # testing
  # feature_counts = rna_ct
  # meta = py_meta
  # plot_features = cc_trail_genes
  # meta_cat_column = "leiden2"
  # plot_cat_features = tnkclus
  # text_expansion = 1
  # condition = "Media"
  # condition_cat = "Condition"
  # test_cat = "progr_contr"
  # stat_compares = v5_compares
  # nudge_nonzero = 1
  # y_limit_expansion_factor = 1.1

  # feature_counts = rna_ct
  # plot_features = cc_trail_genes
  # meta = py_meta
  # meta_cat_column = "leiden2"
  # plot_cat_features = tnkclus
  # text_expansion = 1
  # condition = "Media"
  # condition_cat = "Condition"
  # test_cat = "progr_contr"
  # stat_compares = v4_compares
  # nudge_nonzero = 1
  # y_limit_expansion_factor = 1.1

  # testing
  # feature_counts = rna_ct
  # meta = py_meta
  # plot_features = gene_1
  # meta_cat_column = "annotation"
  # plot_cat_features = cat_1
  # text_expansion = 1
  # condition = "GRV"
  # condition_cat = "condition"
  # test_cat = "progr_contr"
  # stat_compares = compares_1
  # nudge_nonzero = 1
  # y_limit_expansion_factor = 2
  # condition_in_title = FALSE
  # stat_text_offset = 1


  require(ggplot2)

  # mat_ct <- as.matrix(feature_counts[which(row.names(feature_counts) %in% plot_features),])

  # spl_mat2 <- split(x = mat_ct, f = row.names(mat_ct))
  spl_mat2 <- vector("list", length = length(plot_features)); names(spl_mat2) <- plot_features
  for(i in 1:length(plot_features)) {
    spl_mat2[[i]] <- feature_counts[which(row.names(feature_counts)==plot_features[i]),]
  }
  spl_mat <- vector("list",length(spl_mat2)); names(spl_mat) <- names(spl_mat2)
  for(i in 1:length(spl_mat2)) {
    # spl_mat[[i]] <- data.frame(ct = spl_mat2[[i]], cat = factor(meta[,meta_cat_column]),
    spl_mat[[i]] <- data.frame(ct = spl_mat2[[i]], cat = meta[,meta_cat_column],
                               test_col = meta[,test_cat], condition = meta[,condition_cat])
    if(!is.na(condition[1])) {
      spl_mat[[i]] <- spl_mat[[i]][which(spl_mat[[i]]$condition %in% condition),]
    }
    if(plot_cat_features[1]!="all") {
      spl_mat[[i]] <- spl_mat[[i]][which(spl_mat[[i]]$cat %in% plot_cat_features),]
      # spl_mat[[i]]$cat <- factor(as.character(spl_mat[[i]]$cat))
    }
    # cat_levels <- levels(spl_mat[[i]]$cat)
    # spl_levels <- split(x = cat_levels, f = gsub(pattern = "[0-9]+", replacement = "", x = cat_levels))
    # spl_levels <- lapply(X = spl_levels, FUN = function(arg1){
    #   get_nums <- as.numeric(stringr::str_extract(string = arg1, pattern = "[0-9]+"))
    #   return(arg1[order(get_nums)])
    # })
    # levels(spl_mat[[i]]$cat) <- unlist(spl_levels)
    spl_mat[[i]]$gene <- names(spl_mat)[i]
    spl_mat[[i]]$compare_col <- paste0(spl_mat[[i]]$cat,"_",spl_mat[[i]]$test_col)
  }
  # spl_mat <- spl_mat[plot_features]
  # mean(spl_mat$TNFSF10$ct[which(spl_mat$TNFSF10$compare_col=="tnk17_controller")]!=0)

  violin_internal_test <- function(indata, texp = text_expansion, my_compares = stat_compares,
                                   nudge_nz = nudge_nonzero, yle = y_limit_expansion_factor,
                                   cit = condition_in_title, stat_off = stat_text_offset) {
    # testing
    indata <- spl_mat[[1]]
    texp = text_expansion
    my_compares = stat_compares
    nudge_nz = nudge_nonzero
    yle = y_limit_expansion_factor
    cit = condition_in_title
    stat_off = stat_text_offset

    rm_compares <- c()
    for(i in 1:length(my_compares)) {
      sum_nonzero <- rep(NA,length(my_compares[[i]])); names(sum_nonzero) <- my_compares[[i]]
      for(j in 1:length(my_compares[[i]])) {
        sum_nonzero[j] <- sum(indata$ct[which(indata$compare_col==my_compares[[i]][j])]!=0)
      }
      if(sum(sum_nonzero)==0) {
        rm_compares <- append(rm_compares, i)
      }
    }
    if(length(rm_compares)!=0) {
      my_compares <- my_compares[-rm_compares]
    }

    # freqs <- indata %>%
    #   group_by(compare_col) %>%
    #   summarise(freq = round(sum(ct > 0) / n() * 100, 3))
    freqs <- data.frame(compare_col = unique(indata$compare_col), freq = rep(NA,length(unique(indata$compare_col))))
    for(i in 1:nrow(freqs)) {
      freqs$freq[i] <- round(mean(indata$ct[which(indata$compare_col==freqs$compare_col[i])]!=0)*100,3)
    }

    # format_pvalue <- function(p) {
    #   if(is.na(p)) {
    #     return(NA)
    #   }
    #   if(p > 0.0001) {
    #     return(format(p, digits = 4, scientific = FALSE))
    #   } else {
    #     return(format(p, digits = 4, scientific = TRUE))
    #   }
    # }

    indata_jitter <- indata[which(indata$ct!=0),]

    expand_out <- max(indata$ct)*(yle/max(indata$ct) + 1)

    plt <- ggplot(data = indata, aes(x = compare_col, y = ct, fill = compare_col)) +
      geom_violin(scale = "width", trim = TRUE, alpha = 0.7) +
      geom_jitter(data = indata_jitter, width = 0.2, height = 0, size = 1, alpha = 0.5) +
      scale_fill_viridis_d() +
      coord_cartesian(ylim = c(0, max(indata$ct)))
    if(length(my_compares)!=0) {
      plt <- plt +
        stat_compare_means(paired = FALSE, method = "wilcox", comparisons = my_compares, size = 7.6*texp, label.y = max(indata$ct)*stat_off)
    }
    plt <- plt +
      coord_flip() +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 0),
            legend.position = "none") +
      geom_text(data = freqs, aes(x = compare_col, y = max(indata$ct), label = freq), size = 5*texp,
                # nudge_y = ifelse(max(freqs$freq)>=10, nudge_nz*1, nudge_nz)) +
                nudge_y = nudge_nz) +
      expand_limits(y = c(0, ifelse(max(freqs$freq)>=10, expand_out*1.025, expand_out))) +
      ylab("normalized count") +
      ggtitle(ifelse(cit, paste0(indata$gene[1]," - ", indata$condition[1]), indata$gene[1])) +
      theme(axis.text.y = element_text(size = 15*texp, face = "bold"),
            axis.title.y = element_blank(),
            axis.text.x = element_text(size = 14*texp),
            axis.title.x = element_text(size = 15*texp, face = "bold"),
            plot.title = element_text(size = 18*texp, face = "bold", hjust = 0.5))
    return(plt)
  }
  vplots <- lapply(X = spl_mat, FUN = violin_internal_test, texp = text_expansion,
                   my_compares = stat_compares, nudge_nz = nudge_nonzero,
                   yle = y_limit_expansion_factor, cit = condition_in_title)
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

scvi_volcano <- function(indata, plotx = c("lfc","fc"), gene_set = "none", text_size_factor = 1) #order_by = "fold", # gene_set one of c("none","all","top n") or 1d-array of genes to label // sort_by one of c("fold","bayes")
{
  library(ggplot2)
  library(ggrepel)
  library(stringr)
  # testing
  # indata <- dge_list$cl17
  # plotx = c("lfc","fc")
  # gene_set = "top5"
  # text_size_factor = 1
  # order_by = "fc"

  if(length(plotx)!=1) {
    plotx <- "lfc"
  }

  df_segment <- data.frame(
    xstart = c(0,min(min(indata$lfc_mean),0)),
    xend = c(0,max(indata$lfc_mean)),
    ystart = c(min(min(indata$bayes_factor),0),min(min(indata$bayes_factor),0)),
    yend = c(max(indata$bayes_factor),min(min(indata$bayes_factor),0))
  )

  if(length(gene_set)==1) {
    if(gene_set=="all") {
      labdata <- indata
    } else if(regexpr(pattern = "top", text = gene_set)!=-1) {
      topn <- as.numeric(stringr::str_extract(string = gene_set, pattern = "[0-9]+"))
      indata_n <- indata[which(indata$lfc_mean<0),]
      indata_p <- indata[which(indata$lfc_mean>0),]
      if(nrow(indata_n)<topn) {
        lows_lfc <- indata_n$gene[order(indata_n$lfc_mean, decreasing = FALSE)][1:nrow(indata_n)]
        lows_bf <- indata_n$gene[order(indata_n$bayes_factor, decreasing = TRUE)][1:nrow(indata_n)]
      } else {
        lows_lfc <- indata_n$gene[order(indata_n$lfc_mean, decreasing = FALSE)][1:topn]
        lows_bf <- indata_n$gene[order(indata_n$bayes_factor, decreasing = TRUE)][1:topn]
      }
      if(nrow(indata_p)<topn) {
        highs_lfc <- indata_p$gene[order(indata_p$lfc_mean, decreasing = TRUE)][1:nrow(indata_p)]
        highs_bf <- indata_p$gene[order(indata_p$bayes_factor, decreasing = TRUE)][1:nrow(indata_p)]
      } else {
        highs_lfc <- indata_p$gene[order(indata_p$lfc_mean, decreasing = TRUE)][1:topn]
        highs_bf <- indata_p$gene[order(indata_p$bayes_factor, decreasing = TRUE)][1:topn]
      }
      labdata <- indata[which(indata$gene %in% c(lows_lfc, lows_bf, highs_lfc, highs_bf)),]
    }
  } else {
    labdata <- indata[which(indata$gene %in% gene_set),]
  }

  plt <- ggplot() +
    geom_segment(data = df_segment, mapping = aes(x = xstart, y = ystart, xend = xend, yend = yend), lwd = 0.2) +
    geom_point(data = indata, mapping = aes(x = lfc_mean, y = bayes_factor), pch = 21, stroke = 0.3, color = "black", size = 1) +
    theme_minimal() +
    labs(x = "logFC", y = "bayes factor", title = indata$group1[1]) +
    theme(axis.title = element_text(size = 15*text_size_factor),
          axis.text = element_text(size = 12*text_size_factor),
          plot.title = element_text(size = 18*text_size_factor, hjust = 0.5))
  if(gene_set[1]!="none" && nrow(labdata)!=0) {
    plt <- plt +
      geom_point(data = labdata,
                 mapping = aes(x = lfc_mean, y = bayes_factor),
                 pch = 21, fill = "red", color = "black", stroke = 0.3, size = 2) +
      ggrepel::geom_label_repel(data = labdata,
                                mapping = aes(x = lfc_mean, y = bayes_factor,
                                              label = gene), seed = 123,
                                max.overlaps = 15, size = 4*text_size_factor,
                                verbose = FALSE, color = "red", fontface = "bold")
  }
  # plt
  return(plt)
}

seurat_dge <- function(seurat_object,
                       dge_method = c("MAST", "wilcox"),
                       freq_expressed = 0.1,
                       fc_threshold = log2(1.5),
                       test_per_cluster = TRUE,
                       test_clusters = "all",
                       cluster_column = "cell_type",
                       test_per_category = TRUE,
                       test_categories = c("younger","older"),
                       category_column = "age_group",
                       test_per_condition = TRUE, # or FALSE to test all cells without subsetting by condition
                       test_condition = "all", # only considered if test_per_condition, test_condition = "all" suggests
                       condition_column = "condition",
                       pid_column = "pid")
{
  require(Seurat)
  # testing
  # seurat_object = seu
  # seurat_object = seu_test_dge
  # # dge_method = "wilcox"
  # dge_method = "MAST"
  # freq_expressed = 0.1
  # fc_threshold = log2(1.5)
  # test_per_cluster = TRUE
  # # test_clusters = und_cl
  # test_clusters = c("Activ NK","Mono")
  # cluster_column = "cell_type"
  # test_per_category = FALSE
  # test_categories = c("younger","older")
  # category_column = "age_group"
  # test_per_condition = FALSE
  # test_condition = "all"
  # condition_column = "condition"
  # pid_column = "pid"

  if(test_per_condition) {
    conditions <- test_condition[test_condition %in% seurat_object@meta.data[,condition_column]]
    if(length(conditions)==0) {
      conditions <- NULL
      warning("None of the requested 'conditions' found in seurat_object@meta.data[,condition_column]. Continuing without subsetting on any condition(s).")
    }
  } else {
    conditions <- NULL
  }
  if(test_per_cluster) {
    test_clusters <- unique(test_clusters[test_clusters %in% seurat_object@meta.data[,cluster_column]])
    if(length(test_clusters)==0) {
      annos <- NULL
      warning("None of 'test_clusters' found in seurat_object@meta.data[,cluster_column]. Continuing with test(s) on all cells.")
    } else {
      annos <- test_clusters
    }
  } else {
    annos <- NULL
  }
  if(test_per_category) {
    test_cats <- test_categories[test_categories %in% seurat_object@meta.data[,category_column]]
    if(length(test_cats)!=2) {
      stop("After filtering for viable categories, 'test_categories' must be of length 2; which two categories in 'category_column' should be compared?")
    }
  } else {
    test_cats <- NULL
  }
  if(all(is.null(annos), is.null(test_cats))) {
    stop("Nothing to compare")
  }

  dge_outs <- vector("list", length = ifelse(is.null(conditions), 1, length(conditions))); names(dge_outs) <- ifelse(is.null(conditions), "all", conditions)
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
    start_loop <- Sys.time()
    if(names(dge_outs)[i]=="all") {
      subs1 <- seurat_object
    } else {
      seurat_object <- AddMetaData(object = seurat_object, metadata = seurat_object@meta.data[,condition_column], col.name = "c")
      subs1 <- subset(x = seurat_object, subset = c == names(dge_outs)[i])
    }
    for(j in 1:length(dge_outs[[i]])){
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
      if(tolower(dge_method) == "mast") {
        require(MAST)
        binary_expr_matrix <- subs2@assays$RNA@layers$counts > 0
        percent_expr <- rowSums(binary_expr_matrix) / ncol(binary_expr_matrix); names(percent_expr) <- row.names(subs2)
        genes_to_keep <- row.names(subs2)[percent_expr >= freq_expressed]
        if(sum(is.na(genes_to_keep))>0) {
          stop("one or more genes is NA")
        }
        if(length(genes_to_keep)==0) {
          next
        }
        subs2 <- subset(x = subs2, features = genes_to_keep)
        seu_as_sce <- as.SingleCellExperiment(subs2, assay = "RNA")
        logcounts(seu_as_sce) <- log2(counts(seu_as_sce) + 1)
        cdr <- colSums(seu_as_sce@assays@data@listData[["logcounts"]]>0) # cellular detection rate; number of genes detected which is a well-known confounder (https://bioconductor.org/packages/release/bioc/vignettes/MAST/inst/doc/MAITAnalysis.html#22_Filtering)
        seu_as_sce$cngeneson <- scale(cdr)
        subs2@assays$RNA@layers$data <- seu_as_sce@assays@data@listData[["logcounts"]]
        subs2@meta.data <- as.data.frame(seu_as_sce@colData)

        mast_res <- Seurat::FindMarkers(object = subs2, assay = "RNA", ident.1 = ident1, ident.2 = ident2,
                                        test.use = "MAST", only.pos = FALSE, latent.vars = c("cngeneson",pid_column))
        colnames(mast_res)[which(colnames(mast_res)=="pct.1")] <- gsub(" ","_",paste0("pct.",ident1))
        colnames(mast_res)[which(colnames(mast_res)=="pct.2")] <- gsub(" ","",paste0("pct.",ident2))
        mast_res$gene <- row.names(mast_res)
        mast_res <- mast_res[mast_res$p_val_adj<0.05,]
        mast_res <- mast_res[which(abs(mast_res$avg_log2FC)>fc_threshold),]

        dge_outs[[i]][[j]] <- mast_res
      } else if(tolower(dge_method) == "wilcox") {
        require(SeuratWrappers)
        wilcox_res <- SeuratWrappers::RunPresto(object = subs2, assay = "RNA", ident.1 = ident1, ident.2 = ident2,
                                                test.use = "wilcox", only.pos = FALSE, min.pct = 0.01)
        colnames(wilcox_res)[which(colnames(wilcox_res)=="pct.1")] <- gsub(" ","",paste0("pct.",ident1))
        colnames(wilcox_res)[which(colnames(wilcox_res)=="pct.2")] <- gsub(" ","",paste0("pct.",ident2))
        wilcox_res$gene <- row.names(wilcox_res)
        wilcox_res <- wilcox_res[wilcox_res$p_val_adj<0.05,]
        wilcox_res <- wilcox_res[which(abs(wilcox_res$avg_log2FC)>fc_threshold),]

        dge_outs[[i]][[j]] <- wilcox_res
      }
    }
    print(paste0("[",names(dge_outs[[i]])[j],"] in [",names(dge_outs)[i],"] testing took ",round(as.numeric(difftime(Sys.time(), start_loop, units = "mins")),2)," mins"))
  }
  print(paste0("total test time: ",round(as.numeric(difftime(Sys.time(), start_test, units = "hours")),3)," hours"))
  return(dge_outs)
}
