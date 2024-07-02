tile_plots <- function(plotlist, n_row = 2, n_col = 2, rm_legend = FALSE) {
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
                                       ncol = n_col, legend = ifelse(rm_legend, "none","bottom"), common.legend = TRUE)
    print(start_seq[i]:end_seq[i])
  }
  return(outplots)
}
