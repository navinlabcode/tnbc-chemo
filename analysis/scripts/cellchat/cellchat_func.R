#' Compute and visualize the contribution of each ligand-receptor pair in the overall signaling pathways
#'
#' @param object CellChat object
#' @param signaling a signaling pathway name
#' @param sources.use a vector giving the index or the name of source cell groups
#' @param targets.use a vector giving the index or the name of target cell groups.
#' @param width the width of individual bar
#' @param vertex.receiver a numeric vector giving the index of the cell groups as targets in the first hierarchy plot
#' @param thresh threshold of the p-value for determining significant interaction
#' @param return.data whether return the data.frame consisting of the predicted L-R pairs and their contribution
#' @param x.rotation rotation of x-label
#' @param title the title of the plot
#' @param font.size font size of the text
#' @param font.size.title font size of the title
#' @param top show the top ligand-receptor pair
#' @importFrom dplyr select
#' @importFrom ggplot2 ggplot geom_bar aes coord_flip scale_x_discrete element_text theme ggtitle
#' @importFrom cowplot ggdraw draw_label plot_grid
#'
#' @return
#' @export
#' 
#' @note This function is different from `netAnalysis_contribution` in two ways: 1) no need to zoom into a specific pathway; 2) do not count the relative prob. Overall, this function makes the resulted contribution comparable across sources and targets. 
#' @examples
netAnalysis_contribution_allLR <- function(object, signaling = NULL, sources.use = NULL, targets.use = NULL,
                                           width = 0.1, vertex.receiver = NULL, thresh = 0.05, return.data = FALSE,
                                           x.rotation = 0, title = "Contribution of each L-R pair",
                                           font.size = 10, font.size.title = 10, top=10) {
  message('local func')
  if (is.null(signaling)) {signaling <- object@netP$pathways}
  pairLR <- searchPair(signaling = signaling, pairLR.use = object@LR$LRsig, key = "pathway_name", matching.exact = T, pair.only = T)
  pair.name.use = select(object@DB$interaction[rownames(pairLR),],"interaction_name_2")
  
  net <- object@net
  pairLR.use.name <- dimnames(net$prob)[[3]]
  pairLR.name <- intersect(rownames(pairLR), pairLR.use.name)
  pairLR <- pairLR[pairLR.name, ]
  prob <- net$prob
  pval <- net$pval
  
  prob[pval > thresh] <- 0
  
  if (!is.null(sources.use)) {
    if (is.character(sources.use)) {
      if (all(sources.use %in% dimnames(prob)[[1]])) {
        sources.use <- match(sources.use, dimnames(prob)[[1]])
      } else {
        stop("The input `sources.use` should be cell group names or a numerical vector!")
      }
    }
    idx.t <- setdiff(1:nrow(prob), sources.use)
    prob[idx.t, , ] <- 0
  }
  
  if (!is.null(targets.use)) {
    if (is.character(targets.use)) {
      if (all(targets.use %in% dimnames(prob)[[1]])) {
        targets.use <- match(targets.use, dimnames(prob)[[2]])
      } else {
        stop("The input `targets.use` should be cell group names or a numerical vector!")
      }
    }
    idx.t <- setdiff(1:nrow(prob), targets.use)
    prob[ ,idx.t, ] <- 0
  }
  
  if (length(pairLR.name) > 1) {
    pairLR.name.use <- pairLR.name[apply(prob[,,pairLR.name], 3, sum) != 0]
  } else {
    pairLR.name.use <- pairLR.name[sum(prob[,,pairLR.name]) != 0]
  }
  
  
  if (length(pairLR.name.use) == 0) {
    stop(paste0('There is no significant communication'))
  } else {
    pairLR <- pairLR[pairLR.name.use,]
  }
  
  prob <- prob[,,pairLR.name.use]
  
  if (length(dim(prob)) == 2) {
    prob <- replicate(1, prob, simplify="array")
    dimnames(prob)[3] <- pairLR.name.use
  }
  #prob <-(prob-min(prob))/(max(prob)-min(prob)) # this is the relative prob
  
  if (is.null(vertex.receiver)) {
    # pSum <- apply(prob, 3, function(pmat) {sum(pmat) / length(unique(c(sources.use, targets.use)))})
    pSum <- apply(prob, 3, function(pmat) { sum(pmat) / length(sum(pmat>0)) })
    # pSum.max <- sum(prob)
    # pSum <- pSum/pSum.max
    pSum[is.na(pSum)] <- 0
    y.lim <- max(pSum)
    
    pair.name <- unlist(dimnames(prob)[3])
    pair.name <- factor(pair.name, levels = unique(pair.name))
    if (!is.null(pairLR.name.use)) {
      pair.name <- pair.name.use[as.character(pair.name),1]
      pair.name <- factor(pair.name, levels = unique(pair.name))
    }
    mat <- pSum
    df1 <- data.frame(name = pair.name, contribution = mat)
    if(nrow(df1) < 10) {
      df2 <- data.frame(name = as.character(1:(10-nrow(df1))), contribution = rep(0, 10-nrow(df1)))
      df <- rbind(df1, df2)
    } else {
      df <- df1
    }
    df <- df[order(df$contribution, decreasing = TRUE), ]
    # df$name <- factor(df$name, levels = unique(df$name))
    df$name <- factor(df$name,levels=df$name[order(df$contribution, decreasing = TRUE)])
    df1$name <- factor(df1$name,levels=df1$name[order(df1$contribution, decreasing = TRUE)])
    gg <- ggplot(df, aes(x=name, y=contribution)) + geom_bar(stat="identity", width = 0.7) +
      theme_classic() + theme(axis.text.y = element_text(angle = x.rotation, hjust = 1,size=font.size, colour = 'black'), axis.text=element_text(size=font.size),
                              axis.title.y = element_text(size= font.size), axis.text.x = element_blank(), axis.ticks = element_blank()) +
      xlab("") + ylab("Relative contribution") + ylim(0,y.lim) + coord_flip() + theme(legend.position="none") +
      scale_x_discrete(limits = rev(levels(df$name)), labels = c(rep("", max(0, 10-nlevels(df1$name))),rev(levels(df1$name))))
    if (!is.null(title)) {
      gg <- gg + ggtitle(title)+ theme(plot.title = element_text(hjust = 0.5, size = font.size.title))
    }
    gg
    
  } else {
    pair.name <- factor(unlist(dimnames(prob)[3]), levels = unique(unlist(dimnames(prob)[3])))
    # show all the communications
    # pSum <- apply(prob, 3, function(pmat) {sum(pmat) / length(unique(c(sources.use, targets.use)))})
    pSum <- apply(prob, 3, sum)
    # pSum.max <- sum(prob)
    # pSum <- pSum/pSum.max
    pSum[is.na(pSum)] <- 0
    y.lim <- max(pSum)
    
    df<- data.frame(name = pair.name, contribution = pSum)
    gg <- ggplot(df, aes(x=name, y=contribution)) + geom_bar(stat="identity",width = 0.2) +
      theme_classic() + theme(axis.text=element_text(size=10),axis.text.x = element_text(angle = x.rotation, hjust = 1,size=8),
                              axis.title.y = element_text(size=10)) +
      xlab("") + ylab("Relative contribution") + ylim(0,y.lim)+ ggtitle("All")+ theme(plot.title = element_text(hjust = 0.5))#+
    
    # show the communications in Hierarchy1
    if (dim(prob)[3] > 1) {
      # pSum <- apply(prob[,vertex.receiver,], 3, function(pmat) {sum(pmat) / length(unique(c(sources.use, targets.use)))})
      pSum <- apply(prob[,vertex.receiver,], 3, sum)
    } else {
      # pSum <- sum(prob[,vertex.receiver,]) / length(unique(c(sources.use, targets.use)))
      pSum <- sum(prob[,vertex.receiver,])
    }
    
    # pSum <- pSum/pSum.max
    pSum[is.na(pSum)] <- 0
    
    df<- data.frame(name = pair.name, contribution = pSum)
    gg1 <- ggplot(df, aes(x=name, y=contribution)) + geom_bar(stat="identity",width = 0.2) +
      theme_classic() + theme(axis.text=element_text(size=10),axis.text.x = element_text(angle = x.rotation, hjust = 1,size=8), axis.title.y = element_text(size=10)) +
      xlab("") + ylab("Relative contribution") + ylim(0,y.lim)+ ggtitle("Hierarchy1") + theme(plot.title = element_text(hjust = 0.5))#+
    #scale_x_discrete(limits = c(0,1))
    
    # show the communications in Hierarchy2
    
    if (dim(prob)[3] > 1) {
      pSum <- apply(prob[,setdiff(1:dim(prob)[1],vertex.receiver),], 3, sum)
    } else {
      pSum <- sum(prob[,setdiff(1:dim(prob)[1],vertex.receiver),])
    }
    pSum <- pSum/pSum.max
    pSum[is.na(pSum)] <- 0
    
    df<- data.frame(name = pair.name, contribution = pSum)
    gg2 <- ggplot(df, aes(x=name, y=contribution)) + geom_bar(stat="identity", width=0.9) +
      theme_classic() + theme(axis.text=element_text(size=10),axis.text.x = element_text(angle = x.rotation, hjust = 1,size=8), axis.title.y = element_text(size=10)) +
      xlab("") + ylab("Relative contribution") + ylim(0,y.lim)+ ggtitle("Hierarchy2")+ theme(plot.title = element_text(hjust = 0.5))#+
    #scale_x_discrete(limits = c(0,1))
    title <- cowplot::ggdraw() + cowplot::draw_label(paste0("Contribution of each signaling"), fontface='bold', size = 10)
    gg.combined <- cowplot::plot_grid(gg, gg1, gg2, nrow = 1)
    gg.combined <- cowplot::plot_grid(title, gg.combined, ncol = 1, rel_heights=c(0.1, 1))
    gg <- gg.combined
    gg
  }
  if (return.data) {
    df <- subset(df, contribution > 0)
    return(list(LR.contribution = df, gg.obj = gg))
  } else {
    return(gg)
  }
}