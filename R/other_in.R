seu_tile_plots <- function(plotlist, n_row = 2, n_col = 2, rm_legend = FALSE, common_leg = FALSE) {
  require(ggplot2)
  require(ggpubr)
  # testing
  # plotlist <- test_pc
  # n_row = 2
  # n_col = 2

  iseq <- seq(from = 1, to = length(plotlist), by = n_row * n_col)
  start_seq <- iseq
  end_seq <- iseq[2:length(iseq)]-1
  if(max(end_seq)<length(plotlist)) {
    end_seq <- append(end_seq, length(plotlist))
  }
  outplots <- vector("list", length = length(start_seq))
  for(i in 1:length(start_seq)) {
    outplots[[i]] <- ggpubr::ggarrange(plotlist = plotlist[start_seq[i]:end_seq[i]], nrow = n_row,
                                       ncol = n_col, legend = ifelse(rm_legend, "none","bottom"), common.legend = common_leg)
    print(start_seq[i]:end_seq[i])
  }
  return(outplots)
}


seu_plot_volcano <- function(dge_input, plot_clusters = "all",
                             gene_set = NA, prio_top_genes = 5,
                             pval_threshold = 1, table_height = 50,
                             fc_threshold = log2(1.5),
                             de_method = c("seurat_presto","pseudobulk_py","mast","edger"),
                             feature_gsub_pattern = "Hu\\.|Mu\\.",
                             include_gene_table = TRUE,
                             expand_y_axis_factor = 1) # input here is the seurat_dge output list
{
  suppressPackageStartupMessages({
    require(ggplot2)
    require(ggpubr)
    require(ggrepel)
    require(ggpmisc)
    require(scales)
  })
  options(scipen = 999)

  # testing
  # dge_input <- dge_pseudobulk
  # plot_clusters = "all"
  # gene_set = c("","","","","","","","","","","","","","")
  # prio_top_genes = 5
  # pval_threshold = 1
  # table_height = 50
  # fc_threshold = log2(1.5)
  # de_method = "pseudobulk_py"
  # feature_gsub_pattern = "Hu\\.|Mu\\."
  # include_gene_table = TRUE

  de_method <- tolower(de_method)
  if(length(de_method)!=1) {
    stop("use one of 'seurat_presto', 'pseudobulk_py', or 'mast' for 'de_method'; length of 'de_method' must be 1")
  }
  if(de_method=="seurat_presto") {
    logFC_colname <- "avg_log2FC"
    padj_colname <- "p_val_adj"
    cluster_colname <- "cluster"
    if(class(dge_input)=="list") {
      dge_input <- do.call(rbind, dge_input)
      if(!padj_colname %in% colnames(dge_input)) {
        stop("'padj_colname not found in dge_input")
      }
    }
    dge_input$p_val_adj_nlog10 <- -log10(dge_input[,padj_colname])
  } else if(de_method=="pseudobulk_py") {
    logFC_colname <- "log2FoldChange"
    padj_colname <- "padj"
    cluster_colname <- "cluster"
    if(class(dge_input)=="list") {
      dge_input <- do.call(rbind, dge_input)
      if(!padj_colname %in% colnames(dge_input)) {
        stop("'padj_colname not found in dge_input")
      }
    }
    dge_input$p_val_adj_nlog10 <- -log10(dge_input[,padj_colname])
  } else if(de_method=="edger") {
    logFC_colname <- "logFC"
    padj_colname <- "weighted_fdr"
    if(!padj_colname %in% colnames(dge_input)) {
        padj_colname <- "FDR"
    }
    cluster_colname <- "cluster"
    if (class(dge_input) == "list") {
      dge_input <- do.call(rbind, dge_input)
      if (!padj_colname %in% colnames(dge_input)) {
        stop("'padj_colname not found in dge_input")
      }
    }
    dge_input$p_val_adj_nlog10 <- -log10(dge_input[,padj_colname])
  }

  if(plot_clusters[1]!="all") {
    keep_cluster_rows <- which(dge_input[,cluster_colname] %in% plot_clusters)
    if(length(keep_cluster_rows)!=0) {
      dge_input <- dge_input[keep_cluster_rows,]
    }
  }
  dge_input$gene <- gsub(pattern = feature_gsub_pattern, replacement = "", x = dge_input$gene)
  dge_input$avg_FC <- 2^dge_input[,logFC_colname]
  dge_input[,padj_colname] <- ifelse(dge_input[,padj_colname]==0, min(dge_input[,padj_colname][which(dge_input[,padj_colname]!=0)]), dge_input[,padj_colname])
  which_rm1 <- which(dge_input[,padj_colname]>=pval_threshold)
  # dge_input$p_val_adj_nlog10 <- -log10(dge_input[,padj_colname])
  colnames(dge_input)[which(colnames(dge_input)==padj_colname)] <- "padj"
  which_rm2 <- which(abs(dge_input[,logFC_colname])<fc_threshold)
  which_rm <- union(which_rm1, which_rm2)
  if(length(which_rm)!=0) {
    dge_input <- dge_input[-which_rm,]
  }
  dge_input$p_val_adj_nlog10 <- ifelse(is.infinite(dge_input$p_val_adj_nlog10),
                                       max(dge_input$p_val_adj_nlog10[!is.infinite(dge_input$p_val_adj_nlog10)])+1,
                                       dge_input$p_val_adj_nlog10)
  colnames(dge_input)[which(colnames(dge_input)==logFC_colname)] <- "avg fold diff"
  dge_input$`avg fold diff` <- ifelse(dge_input$avg_FC<1,
                                      -1/dge_input$avg_FC,
                                      dge_input$avg_FC)
  if(nrow(dge_input)==0) { # if no dge, return empty plot
    # return("no dge")
    return(ggplot2::ggplot(data = cars, mapping = ggplot2::aes(x = dist, y = speed)) + ggplot2::theme_void())
  }
  avgdfc <- rep(NA,nrow(dge_input))
  for(i in 1:length(dge_input$`avg fold diff`)) {
    if(dge_input$`avg fold diff`[i]<0){
      avgdfc[i] <- -1*log(-1*dge_input$`avg fold diff`[i],2)
    } else {
      avgdfc[i] <- log(dge_input$`avg fold diff`[i],2)
    }
  }
  dge_input$avg_directional_log2FC <- avgdfc

  dge_split <- split(x = dge_input, f = dge_input[,cluster_colname])

  gset <- gene_set

  # if(all(prio_top_genes[1]!=0, !is.na(prio_top_genes[1]))) {
  #   dge_up <- dge_input[which(dge_input$`avg fold diff`>0),]
  #   if(nrow(dge_up)==0) {
  #     top_g <- c()
  #   } else {
  #     top_g <- dge_input$gene[order(dge_input$`avg fold diff`, decreasing = T)][1:prio_top_genes[1]]
  #   }
  #   dge_dn <- dge_input[which(dge_input$`avg fold diff`<0),]
  #   if(nrow(dge_dn)==0) {
  #     bottom_g <- c()
  #   } else {
  #     bottom_g <- dge_input$gene[order(dge_input$`avg fold diff`, decreasing = F)][1:prio_top_genes[1]]
  #   }
  #   gene_set <- c(top_g, bottom_g); gene_set <- gene_set[!is.na(gene_set)]
  # } else {
  #   gene_set <- c()
  # }
  # if(all(!is.na(gset[1]), length(gset)!=0)) {
  #   gene_set <- unique(c(gene_set, gset))
  # }

  tmp_volcano <- function(vol_in, prio_genes = gset, maxyval = max(dge_input$p_val_adj_nlog10),
                          y_thresh = pval_threshold, t_height = table_height, demethod = de_method,
                          igt = include_gene_table, topg = prio_top_genes, eyaf = expand_y_axis_factor)
  {
    # testing
    # vol_in = dge_split[[1]]
    # prio_genes = gset
    # maxyval = max(dge_input$p_val_adj_nlog10)
    # y_thresh = pval_threshold
    # t_height = table_height
    # demethod = de_method
    # igt = include_gene_table
    # topg = prio_top_genes

    # if(all(topg[1]!=0, !is.na(topg[1]))) {
    #   add_genes_1 <- vol_in$gene[order(vol_in$avg_directional_log2FC, decreasing = FALSE)][1:topg]
    #   add_genes_2 <- vol_in$gene[order(vol_in$avg_directional_log2FC, decreasing = TRUE)][1:topg]
    #   prio_genes <- c(prio_genes,add_genes_1,add_genes_2)
    # }

    minx_seg <- ifelse(min(vol_in$avg_directional_log2FC)<0, min(vol_in$avg_directional_log2FC), 0)
    maxx_seg <- ifelse(max(vol_in$avg_directional_log2FC)>0, max(vol_in$avg_directional_log2FC), 0)
    miny_seg <- ifelse(is.na(y_thresh[1]),0,y_thresh)
    maxy_seg <- max(vol_in$p_val_adj_nlog10)

    dfx_segment <- data.frame(x1 = minx_seg, y1 = ifelse(is.na(y_thresh[1]),0,y_thresh),
                              x2 = maxx_seg, y2 = ifelse(is.na(y_thresh[1]),0,y_thresh))
    dfy_segment <- data.frame(x1 = 0, y1 = miny_seg,
                              x2 = 0, y2 = maxy_seg)

    vol_in$pscaled <- scales::rescale(x = vol_in$p_val_adj_nlog10, to = c(0,max(abs(vol_in$avg_directional_log2FC))))
    vol_in$distance <- sqrt(vol_in$avg_directional_log2FC^2 + vol_in$pscaled^2)

    upgenes <- vol_in[vol_in$avg_directional_log2FC>0,]; gene_set_up <- upgenes$gene[order(upgenes$distance, decreasing = TRUE)]
    dngenes <- vol_in[vol_in$avg_directional_log2FC<0,]; gene_set_dn <- dngenes$gene[order(dngenes$distance, decreasing = TRUE)]

    lab_genes <- vol_in[vol_in$gene %in% c(gene_set_up[1:topg], gene_set_dn[1:topg]),]

    high_gene <- vol_in[vol_in$gene %in% c(gene_set_up[1:t_height]),]
    low_gene <- vol_in[vol_in$gene %in% c(gene_set_dn[1:t_height]),]

    pgene <- c(gene_set_up[1:topg], gene_set_dn[1:topg], prio_genes); pgene <- pgene[!is.na(pgene)]
    pgene <- pgene[pgene %in% vol_in$gene]
    vol_in <- vol_in[vol_in$p_val_adj_nlog10>y_thresh,] # new

    volc <- ggplot(data = vol_in, mapping = aes(x = avg_directional_log2FC, y = p_val_adj_nlog10)) +
      # ylim(maxy_seg,0) +
      geom_segment(data = dfx_segment, mapping = aes(x = x1, y = y1, xend = x2, yend = y2), color = "#757575") +
      geom_segment(data = dfy_segment, mapping = aes(x = x1, y = y1, xend = x2, yend = y2), color = "#757575")
    # if(all(!is.na(gene_set[1]), length(gene_set)!=0)) {
    if(length(pgene)!=0) {
      # if(prio_top_genes!=0) {
      # add_genes <- vol_in$gene[order(vol_in$avg_directional_log2FC, decreasing = TRUE)][1:prio_top_genes]
      # pgene <- c(pgene, add_genes)
      # }
      # prio_row <- which(vol_in$gene %in% pgene)
      # if(length(prio_row)!=0) {
      # if(length(pgene)!=0) {
      prio_data <- vol_in[which(vol_in$gene %in% pgene),]
      prio_data <- prio_data[prio_data$p_val_adj_nlog10>y_thresh,] # new
      volc <- volc + geom_point(pch = 21) +
        geom_point(data = prio_data,
                   mapping = aes(x = avg_directional_log2FC, y = p_val_adj_nlog10),
                   pch = 21, fill = "red", color = "black", stroke = 0.5, size = 3) +
        ggrepel::geom_label_repel(data = prio_data, min.segment.length = 0,
                                  mapping = aes(x = avg_directional_log2FC, y = p_val_adj_nlog10,
                                                label = gene), seed = 123,
                                  max.overlaps = 20, size = 6,
                                  verbose = FALSE, color = "red", fontface = "bold")
      # } else {
      #   volc <- volc + geom_point() + ggrepel::geom_label_repel(mapping = aes(label = gene), seed = 123,
      #                                                           max.overlaps = 15, min.segment.length = 0,
      #                                                           label.size = NA, fill = NA)
      # }
    } else {
      volc <- volc + geom_point() + ggrepel::geom_label_repel(mapping = aes(label = gene), seed = 123,
                                                              max.overlaps = 15, min.segment.length = 0,
                                                              label.size = NA, fill = NA)
    }
    # volc <- volc + ggtitle(paste0("cluster ",vol_in[,cluster_colname][1])) +
    volc <- volc + ggtitle(paste0(vol_in[,cluster_colname][1])) +
      labs(x = "log2 fold difference", y = "-log10(p-value)") +
      theme_minimal() +
      coord_cartesian(ylim = c(ifelse(is.na(y_thresh[1]),0.1,y_thresh),maxyval*eyaf)) +
      # scale_y_continuous(limits = c(ifelse(is.na(y_thresh[1]),0.1,y_thresh),maxyval*eyaf), trans='pseudo_log') +
      scale_y_continuous(trans='pseudo_log') +
      theme(plot.title = element_text(hjust = 0.5, size = 26),
            axis.title = element_text(size = 23),
            axis.text = element_text(size = 18))

    # if(demethod=="pseudobulk_py") {
    #   dn_direction <- vol_in[which(vol_in$stat<0),]
    #   low_gene <- dn_direction[order(dn_direction$stat,decreasing = FALSE),
    #                            c("stat","avg fold diff","padj","gene")]
    #   up_direction <- vol_in[which(vol_in$stat>0),]
    #   high_gene <- up_direction[order(up_direction$stat,decreasing = TRUE),
    #                             c("stat","avg fold diff","padj","gene")]
    # } else if(demethod=="seurat_presto") {
    #   up_df <- vol_in[which(vol_in$`avg fold diff`>0),]
    #   dn_df <- vol_in[which(vol_in$`avg fold diff`<0),]
    #   if(nrow(up_df)==0) {
    #     up_direction <- vol_in[0,]
    #   } else {
    #     up_direction <- vol_in[order(vol_in$`avg fold diff`, decreasing = TRUE),]
    #   }
    #   if(nrow(dn_df)==0) {
    #     dn_direction <- vol_in[0,]
    #   } else {
    #     dn_direction <- vol_in[order(vol_in$`avg fold diff`, decreasing = FALSE),]
    #   }
    #   low_gene <- dn_direction[,c("avg fold diff","padj","gene")]
    #   high_gene <- up_direction[,c("avg fold diff","padj","gene")]
    # }

    prep_table <- function(table_in) {
      # testing
      # table_in <- high_gene

      if(nrow(table_in)==0) {
        return(table_in)
      }

      row.names(table_in) <- table_in$gene; table_in <- table_in[,-which(colnames(table_in)=="gene")]
      # if(nrow(table_in)==0) {
      #   if(demet=="pseudobulk_py") {
      #     table_tmp <- data.frame('stat' = NA, 'avg fold diff' = NA, 'padj' = NA, check.names = FALSE)
      #     row.names(table_tmp) <- "none"
      #   } else if(demet=="seurat_pesto") {
      #     table_tmp <- data.frame('avg fold diff' = NA, 'padj' = NA, check.names = FALSE)
      #     row.names(table_tmp) <- "none"
      #   }
      #   return(table_tmp)
      # } else {
      #   if(demet=="pseudobulk_py") {
      #     table_in$`avg fold diff` <- round(table_in$`avg fold diff`,2)
      #     table_in$stat <- round(table_in$stat,2)
      #   } else if(demet=="seurat_presto") {
      #     table_in$`avg fold diff` <- round(table_in$`avg fold diff`,2)
      #   }
      # }

      convert_p <- function(pval) {
        if(pval==0) {
          return("0")
        }
        if(regexpr(pattern = "e", text = as.character(pval))!=-1) {
          p_exp <- gsub("e","",stringr::str_extract(string = as.character(pval), pattern = "e[0-9-]+$"))
          p_num <- stringr::str_extract(string = as.character(pval), pattern = "^[0-9.]+")
          p_num_rounded <- round(x = as.numeric(p_num), digits = 2) # this assumes exponential notation was used
          return(paste0(p_num_rounded,"e",p_exp))
        } else {
          # p_exp <- stringr::str_extract(string = as.character(pval), pattern = "[1-9]+$")
          # p_num <- gsub(pattern = "^0\\.", replacement = "", x = stringr::str_extract(string = as.character(pval), pattern = "^[0.]+"))
          # numzero <- nchar(p_num)
          # if(nchar(p_exp)>3) {
          #   num_pre <- substr(p_exp, start = 1, stop = 3)
          #   num_pre <- paste0(substr(num_pre,1,1),".",substr(num_pre,2,nchar(num_pre)))
          # } else {
          #   num_pre <- p_exp
          #   if(nchar(num_pre)!=1) {
          #     num_pre <- paste0(substr(num_pre,1,1),".",substr(num_pre,2,nchar(num_pre)))
          #   }
          # }
          return(formatC(pval, 3))
        }
      }
      table_in$padj <- sapply(X = table_in$padj, FUN = convert_p)
      table_out <- data.frame("FC" = ifelse(table_in$avg_directional_log2FC<0,(2^abs(table_in$avg_directional_log2FC))*-1,2^abs(table_in$avg_directional_log2FC)),
                              "padj" = table_in$padj, row.names = row.names(table_in))
      # if(nrow(table_in)>tht) {
      #   table_in <- table_in[1:tht,]
      # }
      #colnames(table_in)[2] <- "adj p-value"
      return(table_out)
    }

    low_gene_prep <- prep_table(table_in = low_gene)
    high_gene_prep <- prep_table(table_in = high_gene)

    # if(demethod=="pseudobulk_py") {
    #   low_gene_prep <- low_gene_prep[which(low_gene_prep$stat<0),]
    #   high_gene_prep <- high_gene_prep[which(high_gene_prep$stat>0),]
    # } else if(demethod=="seurat_presto") {
    #   low_gene_prep <- low_gene_prep[which(low_gene_prep$`avg fold diff`<0),]
    #   high_gene_prep <- high_gene_prep[which(high_gene_prep$`avg fold diff`>0),]
    # }

    blank_plt <- ggplot(data = cars, mapping = aes(x = speed, y = dist)) +
      geom_point(alpha = 0) + theme_void()
    if(nrow(high_gene_prep)==0) {
      high_tab <- blank_plt
    } else {
      high_tab <- blank_plt +
        annotate(geom='table',x=0.5,y=0.5,label=list(high_gene_prep),table.rownames=TRUE)
    }
    if(nrow(low_gene_prep)==0) {
      low_tab <- blank_plt
    } else {
      low_tab <- blank_plt +
        annotate(geom='table',x=0.5,y=0.5,label=list(low_gene_prep),table.rownames=TRUE)
    }

    if(igt) {
      arr_plot <- ggpubr::ggarrange(plotlist = list(volc, high_tab, low_tab),
                                    nrow = 1, widths = c(0.65,0.175,0.175))
      return(arr_plot)
    } else {
      return(volc)
    }
  }

  out_volc <- lapply(X = dge_split, FUN = tmp_volcano)

  return(out_volc)
}


seu_mast_venns <- function(mast_output, venn_colors = c("#F57336","#32671D","#3398DB","#7B1FA3"),
                       venn_groups = c("controller","progressor","progr less","progr more"),
                       repl_patterns = c("^c$","^p$","^l$","^m$")) { # '\n' will be paste0'd in front of repl_patterns, order matches venn_groups
  require(VennDiagram)
  #### mast_output should take the form of
  # "J:/10x/TB_sc/scbp2/r_out_report/mast_Media_dge_permuted_comparisons.rds"
  ####
  mast <- mast_output
  for(i in 1:length(mast)) {
    mast[[i]] <- append(mast[[i]], list(ct = names(mast)[i]))
  }
  venns <- lapply(X = mast, FUN = function(arg1, vc = venn_colors, vg = venn_groups, rp = repl_patterns) {
    gene_sets <- lapply(X = arg1[1:4], FUN = function(arg2) return(arg2$gene))
    cat_names <- names(gene_sets)
    repl_df <- data.frame(var1 = rp, var2 = vg)
    for(i in 1:length(cat_names)) {
      tmp_name <- strsplit(x = gsub("mast_","",cat_names[i]), split = "")[[1]]
      for(j in 1:nrow(repl_df)) {
        tmp_name[1] <- gsub(pattern = repl_df$var1[j], replacement = paste0(repl_df$var2[j],"\n"), x = tmp_name[1])
        tmp_name[2] <- gsub(pattern = repl_df$var1[j], replacement = repl_df$var2[j], x = tmp_name[2])
      }
      cat_names[i] <- paste0(tmp_name[1], tmp_name[2])
    }
    vd <- venn.diagram(x = gene_sets, filename = NULL, category.names = cat_names,
                       output = TRUE, fill = vc,
                       cat.cex = 1.25, main = arg1[[5]], main.cex = 2.5)
    return(vd)
  })
  return(venns)
  #### example of how to plot mast_venns output
  # pdf(paste0("J:/10x/TB_sc/scbp2/r_figure/venn/",stim_condition,"/",stim_condition,"_mast_dge_overlapping_genes.pdf"), width = 6, height = 6)
  # lapply(venns, function(gList) {
  #   grid.newpage()
  #   grid.draw(gList)
  # })
  # dev.off()
  ####
}

seu_mast_sets <- function() {
  packageExt <- system.file("extdata", package='seutools')
  # module_file <- list.files(path = packageExt, pattern = "gsea", full.names = TRUE)
  module_file <- list.files(path = "J:/seutools/seutools/inst/extdata", pattern = "gsea", full.names = TRUE)
  shortened_gmt <- sapply(X = module_file, FUN = function(arg1){
    string_out <- strsplit(x = arg1, split = "\\/")[[1]]
    neat_path <- gsub(pattern = "_modules_gsea\\.gmt", replacement = "", x = string_out[length(string_out)])
    return(neat_path)
  })
  data.frame(gene_set_db = as.character(shortened_gmt))
  print("gene set databases can be found here: https://www.gsea-msigdb.org/gsea/msigdb/human/genesets.jsp")
  print(paste0("use one or more of: ",paste0(as.character(shortened_gmt), collapse = ", ")))
}

seu_mast_gsea <- function(mast_dge_result, seu_mast_sets,
                          num_boots = 50, gs_min = 5, gs_max = Inf, gs_regex = NULL,
                          nthread = 1, prepare_plot_n = 20, verbose = TRUE) {
  # testing
  # mast_dge_result = mast_res[["ISG_Mono"]][["stim"]][["ISG_Mono"]]
  # seu_mast_sets = c("C2_CP_Reactome","C5_GO_BP",
  #                   "C5_GO_MF","H_Hallmark",
  #                   "C2_CP_KEGG_LEGACY","C2_CP_KEGG_MEDICUS",
  #                   "C7_ImmuneSigDB","C7_VAX")
  # num_boots = 8
  # gs_min = 3
  # prepare_plot_n = 20
  # gs_max = Inf
  # gs_regex = NULL
  # nthread = 2
  # verbose = TRUE

  if(sum(duplicated(mast_dge_result$raw_res$primerid))!=0) {
    stop(paste0('mast_dge_result$raw_res$primerid must not contain duplicates'))
  } else {
    row.names(mast_dge_result$raw_res) <- mast_dge_result$raw_res$primerid
  }
  positive_group <- gsub(pattern = "^.+\\.", replacement = "", x = mast_dge_result$raw_res$contrast[1])
  mast_zlmfit <- mast_dge_result[['zlmfit']]
  sca <- mast_dge_result[['sca']]

  packageExt <- system.file("extdata", package='seutools')

  if(verbose==TRUE) {
    if(num_boots<=50) {
      warning(paste0("'num_boots = ",num_boots,". MAST authors recommend >50 should be used if values will be reported. Less is fine for exploration."))
    }
    print(paste0("Users may specify their own gene sets by inserting gmt files here: ", packageExt," Inserted gmt file names must end with 'modules_gsea.gmt'. See existing gmt files for examples."))
  }

  db_pattern <- paste0("(",paste0(seu_mast_sets, collapse = "|"),")")

  if(nthread==1) {
    boots <- bootVcov1(zlmfit = mast_zlmfit, R = num_boots)
  } else {
    nthread <- round(nthread)
    makecl <- parallel::makePSOCKcluster(nthread)
    on.exit(parallel::stopCluster(makecl)) # Ensure the cluster is stopped
    boots <- pbootVcov1(zlmfit = mast_zlmfit, R = num_boots, cl = makecl)
  }
  # saveRDS(object = boots, file = "J:/seutools/seutools/test_script/mast_boots.rds")
  # boots <- readRDS(file = "J:/seutools/seutools/test_script/mast_boots.rds")
  # parallel::stopCluster(makecl)
  # lapply(makecl, function(node) {
  #   try(parallel::stopCluster(node), silent = TRUE)
  # })

  list_module_files <- list.files(packageExt, pattern = "modules_gsea\\.gmt", full.names = TRUE)
  module_files <- list_module_files[grep(pattern = db_pattern, x = list_module_files)]
  if(length(module_files)==0) {
    stop("No supported gene set databases found for testing. See args(seu_mast_gsea) for supported gene set databases or insert your own according to print messages given by 'verbose' = TRUE")
  }

  gsea_result <- vector("list", length = length(module_files))

  for(i in 1:length(module_files)) {
    module_file <- module_files[i]
    target_module <- stringr::str_extract(string = module_file, pattern = db_pattern)
    names(gsea_result)[i] <- target_module
    gene_set <- GSEABase::getGmt(module_file)
    gene_ids <- GSEABase::geneIds(gene_set)
    sets_indices <- limma::ids2indices(gene_ids, mcols(sca)$primerid)
    gs_min_thresh <- which(sapply(sets_indices, length) >= gs_min)
    gs_max_thresh <- which(sapply(sets_indices, length) <= gs_max)
    keep_gs_indices_length <- intersect(gs_min_thresh, gs_max_thresh)
    if(!is.null(gs_regex)) {
      keep_gs_indices_regex <- grep(pattern = gs_regex, x = names(sets_indices))
      keep_gs_indices <- intersect(keep_gs_indices_length, keep_gs_indices_regex)
    } else {
      keep_gs_indices <- keep_gs_indices_length
    }
    if(length(keep_gs_indices)==1) {
      sets_indices <- sets_indices[[keep_gs_indices]]
    } else if(length(keep_gs_indices)>1) {
      sets_indices <- sets_indices[keep_gs_indices]
    } else {
      stop("no gene sets to test with current 'gs_min', 'gs_max'")
    }
    geneset_meta <- vector('list', length = length(sets_indices)); names(geneset_meta) <- names(sets_indices)
    for(j in 1:length(sets_indices)) {
      genes_found <- mcols(sca)$primerid[sets_indices[[j]]]
      geneset_res <- mast_dge_result$raw_res[gene_ids[[j]],]
      geneset_res$primerid <- gene_ids[[j]]
      geneset_res$geneset <- names(sets_indices)[j]
      row.names(geneset_res) <- NULL
      geneset_meta[[j]] <- geneset_res
    }

    gsea <- gseaAfterBoot(mast_zlmfit, boots, sets_indices, CoefficientHypothesis("categoryGroup2"))
    gsea_table <- summary(gsea, testType='normal')
    gsea_table_top <- gsea_table[gsea_table$combined_Z>0,]; gsea_table_top <- gsea_table_top[order(gsea_table_top$combined_Z, decreasing = T),]
    if(nrow(gsea_table_top)>prepare_plot_n) {
      gsea_table_top <- gsea_table_top[1:prepare_plot_n,]
    }
    gsea_table_bottom <- gsea_table[gsea_table$combined_Z<0,]; gsea_table_bottom <- gsea_table_bottom[order(gsea_table_bottom$combined_Z, decreasing = F),]
    if(nrow(gsea_table_bottom)>prepare_plot_n) {
      gsea_table_bottom <- gsea_table_bottom[1:prepare_plot_n,]
    }
    gsea_top_bottom_table <- rbind(gsea_table_top, gsea_table_bottom[order(gsea_table_bottom$combined_Z, decreasing = T),])
    gseaTable_plot <- melt(gsea_top_bottom_table[,.(set, disc_Z, cont_Z, combined_Z)], id.vars='set')
    gseaTable_report <- as.data.frame(gsea_table)
    gseaTable_report$positive_score_group <- positive_group
    gsea_result[[i]] <- list(result = gseaTable_report,
                             melted_result = gseaTable_plot,
                             geneset_mast_res = geneset_meta)
  }

  return(gsea_result)
  # gseaTable <- melt(most_sigModules[,.(set, disc_Z, cont_Z, combined_Z)], id.vars='set')
  # ggplot(gseaTable, aes(y=set, x=variable, fill=value)) + geom_raster() +
  #   scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  #   theme_bw()
}


seu_pathfindr <- function(dge_dataframe, deseq_ct = NULL, deseq_meta = NULL,
                          pid_column = "pid", group_column = "group_id", group_str = c("group1","group2"),
                          gene_sets = "Reactome", nthread = 40, gs_min = 5, gs_max = 500) {
  require(pathfindR)

  if(class(dge_dataframe)=="list") {
    deseq_result <- deseq_output[[1]]
    deseq_ct <- deseq_output[[2]]
    pathR_cluster <- deseq_output[[1]]$cluster[1]
    pathR_meta <- deseq_output[[3]]
  } else {
    deseq_result <- dge_dataframe
  }

  if(any(is.null(deseq_ct),is.null(deseq_ct))) {
    stop("missing counts matrix; unable to score pathways")
  }

  gene_col <- grep(pattern = "gene", x = tolower(colnames(deseq_result)))
  logfc_col <- grep(pattern = "log2", x = tolower(colnames(deseq_result)))
  p_col <- grep(pattern = "padj", x = tolower(colnames(deseq_result)))

  # deseq_subset <- deseq_result[deseq_result$cluster==pathR_cluster,]
  path_input <- data.frame(Gene.symbol = deseq_result[,gene_col], logFC = deseq_result[,logfc_col], adj.P.Val = deseq_result[,p_col])
  output_df <- run_pathfindR(input = path_input,
                             gene_sets = gene_sets, min_gset_size = gs_min, max_gset_size = gs_max, enrichment_threshold = 0.05,
                             iterations = 25, n_processes = nthread, plot_enrichment_chart = FALSE)

  pathfindR_clustered <- try(cluster_enriched_terms(enrichment_res = output_df, plot_dend = FALSE,
                                                    plot_clusters_graph = FALSE, kappa_threshold = 0.35), silent = TRUE)
  if("try-error" %in% class(pathfindR_clustered)) {
    pathfindR_repr <- output_df
    upgene <- gsub(" ","",as.character(unlist(sapply(X = output_df$Up_regulated, FUN = strsplit, split = ","))))
    dngene <- gsub(" ","",as.character(unlist(sapply(X = output_df$Down_regulated, FUN = strsplit, split = ","))))
    pathgene <- unique(c(upgene, dngene))
    if(length(pathgene)==0) {
      return("no sig path genes")
    } else {
      normalized_counts <- matrix(deseq_ct[which(row.names(deseq_ct) %in% pathgene),], nrow = 1)
      colnames(normalized_counts) <- colnames(deseq_ct)
      row.names(normalized_counts) <- row.names(deseq_ct)[which(row.names(deseq_ct) %in% pathgene)]
    }
  } else {
    pathfindR_repr <- pathfindR_clustered[pathfindR_clustered$Status=="Representative",]
    upgene <- gsub(" ","",as.character(unlist(sapply(X = pathfindR_repr$Up_regulated, FUN = strsplit, split = ","))))
    dngene <- gsub(" ","",as.character(unlist(sapply(X = pathfindR_repr$Down_regulated, FUN = strsplit, split = ","))))
    pathgene <- unique(c(upgene, dngene))
    if(length(pathgene)==0) {
      return("no sig path genes")
    } else {
      normalized_counts <- deseq_ct[which(row.names(deseq_ct) %in% pathgene),]
    }
  }

  cases <- pathR_meta$bid[pathR_meta[,group_column]==group_str[1]]

  score_matrix <- score_terms(
    enrichment_table = pathfindR_repr,
    exp_mat = normalized_counts,
    cases = cases,
    plot_hmap = FALSE,
    use_description = TRUE, # default FALSE
    label_samples = FALSE, # default = TRUE
    case_title = group_str[1], # default = "Case"
    control_title = group_str[2], # default = "Control"
    low = "#35338e", # default = "green"
    mid = "#fffde4", # default = "black"
    high = "#b63d3d" # default = "red"
  )
#   return(list(path_scores = score_matrix, path_condition = deseq_result$condition[1],
#               path_cluster = pathR_meta$cluster_id[1], deseq_norm_counts = deseq_ct))
# }

  # scores <- lapply(X = deseq_out, FUN = pathR_scoring)

  # saveRDS(object = scores, file = paste0("/media/MPEdge16/TB_sc/sc/py/scbp2/py_out_obj/final_annotation_pseudobulk/",
  #                                        tolower(condition),"/",tolower(condition),"_deseq2_results_list.rds"))

  # scores_format <- lapply(X = scores, FUN = function(arg1){
  #   path_scores <- arg1[[1]]
  #   row.names(path_scores) <- gsub(" ","_",paste0(row.names(path_scores),".",arg1[[2]],".",arg1[[3]]))
  #   return(t(path_scores))
  # })

  # saveRDS(object = scores_format, file = paste0("/media/MPEdge16/TB_sc/sc/py/scbp2/py_out_obj/final_annotation_pseudobulk/",
  #                                               tolower(condition),"/",tolower(condition),"_pseudobulk_deseq2_pathfindR_scores.rds"))

  # deseq_norm_ct <- lapply(X = scores, FUN = function(arg1) return(arg1[[4]]))

  return(list(path_scores = score_matrix, path_condition = deseq_result$condition[1],
              path_cluster = pathR_meta$cluster_id[1], deseq_norm_counts = deseq_ct))
}
