
# custom function to plot DE genes between Trisomy and Normal by clusters
plot_de_volcano <- function(f, coef, title="Volcano Plot", maxq=0.1, minLFC=1, sigcol="red", sigLowFCcol="navyblue", showsub=TRUE, miRNA=FALSE) {
  if(miRNA) {
    x <- topTags(f, n=nrow(f$coefficients))$table
    colnames(x)[grep("FDR", colnames(x))] <- "adj.P.Val"
  } else {
    x <- topTable(f, coef=coef, n=nrow(f$coefficients))
  }
  
  x$logq <- -log10(x$adj.P.Val)
  minFDR <- -log10(maxq)
  sigs <- x$logq >= minFDR & abs(x$logFC) >= minLFC
  sigLowFC <- x$logq >= minFDR & abs(x$logFC) < minLFC
  
  cnt_plus <- sum(x$logq >= minFDR & x$logFC >= minLFC)
  cnt_neg <- sum(x$logq >= minFDR & x$logFC <= -minLFC)
  if(showsub) {
    subtitle <- paste0(sum(sigs), " DE genes: ", cnt_plus, " Pos and ", cnt_neg, " Neg",
                       " (FDR<", maxq, ", log2FC>", minLFC, ")")
  } else {
    subtitle <- NULL
  }
  
  p <- ggplot(x) +
    geom_point(aes(x=logFC, y=logq), size=0.5) + 
    geom_point(data=x[sigs,], aes(x=logFC, y=logq), size=0.5, col=sigcol) +
    geom_point(data=x[sigLowFC,], aes(x=logFC, y=logq), size=0.5, col=sigLowFCcol) +
    ggtitle(title, subtitle=subtitle) + xlab("log2 Fold Change") + ylab("log10 FDR") +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) + 
    geom_vline(xintercept=minLFC, col=sigcol, linetype="dashed", size=0.5) + 
    geom_vline(xintercept=-minLFC, col=sigcol, linetype="dashed", size=0.5) + 
    geom_hline(yintercept=minFDR, col=sigcol, linetype="dashed", size=0.5)
  
  return(p)
}

create_dt <- function(df,caption){
  datatable(df,
            caption = caption,
            extensions = 'Buttons',
            options = list(dom = 'Blfrtip',
                           buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                           lengthMenu = list(c(10,25,50,-1),
                                             c(10,25,50,"All"))))
}

# function to set a column as rownames and remove that column from the df; default is setting the first column as the rowname
setrowname =  function(df,col = 1){
  if (is.numeric(col)) {
    # print('it is a number')
    rownames(df) = df[,col]
    df = df[,-col]
  } else {
    # print('it is a string')
    rownames(df) = df[,which(names(df) %in% c(col))]
    df = df[,-which(names(df) %in% c(col))]
  }
  return(df)
} 

# function to get colors from a palatte
getcol = function(group,pal = "Darjeeling1",shift_col_palatte=0){
  groupcol = as.numeric(factor(group))
  groupcol = wes_palette(pal)[groupcol+shift_col_palatte]
  return(groupcol)
}

# Volcano plot
plot_de_volcano_scatter3 = function(df, maxq = 0.1, minLFC=1, showsub = T, sigcol = 'red', sigLowFCcol = 'navyblue', title = 'Volcano plot',
                                    axislim,leftarrowlim,rightarrowlim) {
  # calculate log10 q value. 
  df$qval_log10 = -log10(df$`adj.P.Val`)
  # generate summary statistics to show on caption
  minFDR <- -log10(maxq)
  sigs <- df$qval_log10 >= minFDR & abs(df$logFC) >= minLFC
  sigLowFC <- df$qval_log10 >= minFDR & abs(df$logFC) < minLFC
  
  cnt_plus <- sum(df$qval_log10 >= minFDR & df$logFC >= minLFC)
  cnt_neg <- sum(df$qval_log10 >= minFDR & df$logFC <= -minLFC)
  if(showsub) {
    subtitle <- paste0(sum(sigs), " DEGs: ", cnt_plus, " + and ", cnt_neg, " -",
                       " (q<", maxq, ", log2FC>", minLFC, ")")
  } else {
    subtitle <- NULL
  }
  
  p1 = ggplot(df) +
    geom_point(aes(x=logFC, y=qval_log10, color = logFC, size = qval_log10)) + 
    geom_point(data =  df %>%
                 subset(qval_log10 >= minFDR) %>%
                 dplyr::arrange(desc(logFC)) %>%
                 dplyr::slice(1:5),
               aes(x = logFC, y = qval_log10,
                   # fill = log2FoldChange,
                   size = qval_log10),
               shape = 21, show.legend = F, color = "#000000") +
    geom_text_repel(data =  df %>%
                      subset(qval_log10 >= minFDR) %>%
                      dplyr::arrange(desc(logFC)) %>%
                      dplyr::slice(1:5),
                    aes(x = logFC, y = qval_log10,label = geneID),
                    box.padding = 0.5,
                    nudge_x = 0.5,
                    nudge_y = 0.2,
                    segment.curvature = -0.1,
                    segment.ncp = 3,
                    # segment.angle = 10,
                    direction = "y", 
                    hjust = "left"
    ) + 
    geom_point(data =  df %>%
                 subset(qval_log10 >= minFDR) %>%
                 dplyr::arrange(logFC) %>%
                 dplyr::slice(1:5),
               aes(x = logFC, y = qval_log10,
                   # fill = log2FoldChange,
                   size = qval_log10),
               shape = 21, show.legend = F, color = "#000000") +
    geom_text_repel(data =  df %>%
                      subset(qval_log10 >= minFDR) %>%
                      dplyr::arrange(logFC) %>%
                      dplyr::slice(1:5),
                    aes(x = logFC, y = qval_log10,label = geneID),
                    box.padding = 0.5,
                    nudge_x = 0.5,
                    nudge_y = 0.2,
                    segment.curvature = -0.1,
                    segment.ncp = 3,
                    # segment.angle = 10,
                    direction = "y", 
                    hjust = "right"
    ) + 
    scale_color_gradientn(colours = c("#3288bd", "#66c2a5","#ffffbf", "#f46d43", "#9e0142"),
                          values = seq(0, 1, 0.2)) +
    scale_fill_gradientn(colours = c("#3288bd", "#66c2a5","#ffffbf", "#f46d43", "#9e0142"),
                         values = seq(0, 1, 0.2)) +
    ggtitle(title, subtitle=subtitle) + xlab("log2 FC") + ylab("-log10 q value") +
    labs(color = 'log2 FC', size = '-log10 q') +
    theme(panel.grid = element_blank(),
          legend.background = element_roundrect(color = "#808080", linetype = 1),
          axis.text = element_text(size = 13, color = "#000000"),
          axis.title = element_text(size = 15),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)) + 
    geom_vline(xintercept=minLFC, linetype=2) + 
    geom_vline(xintercept=-minLFC, linetype=2) + 
    geom_hline(yintercept=minFDR, linetype=4) +
    coord_cartesian(clip = "off") + 
    annotation_custom(
      grob = grid::segmentsGrob(
        y0 = unit(-10, "pt"),
        y1 = unit(-10, "pt"),
        arrow = arrow(angle = 45, length = unit(.2, "cm"), ends = "first"),
        gp = grid::gpar(lwd = 3, col = "#74add1")
      ), 
      xmin = leftarrowlim[1], 
      xmax = leftarrowlim[2],
      ymin = leftarrowlim[3],
      ymax = leftarrowlim[4]
    ) +
    annotation_custom(
      grob = grid::textGrob(
        label = "Down",
        gp = grid::gpar(col = "#74add1")
      ),
      xmin = leftarrowlim[1], 
      xmax = leftarrowlim[2],
      ymin = leftarrowlim[3],
      ymax = leftarrowlim[4]
    ) +
    annotation_custom(
      grob = grid::segmentsGrob(
        y0 = unit(-10, "pt"),
        y1 = unit(-10, "pt"),
        arrow = arrow(angle = 45, length = unit(.2, "cm"), ends = "last"),
        gp = grid::gpar(lwd = 3, col = "#d73027")
      ), 
      xmin = rightarrowlim[1], 
      xmax = rightarrowlim[2],
      ymin = rightarrowlim[3],
      ymax = rightarrowlim[4]
    ) +
    annotation_custom(
      grob = grid::textGrob(
        label = "Up",
        gp = grid::gpar(col = "#d73027")
      ),
      xmin = rightarrowlim[1], 
      xmax = rightarrowlim[2],
      ymin = rightarrowlim[3],
      ymax = rightarrowlim[4]
    ) 
  if (length(axislim) == 4) {
    p1 = p1 + xlim(axislim[1], axislim[2]) + ylim(axislim[3], axislim[4])
  }
  return(p1)
}


# select genes with satisifying conditions from the toptable
get_filtered_topTable <- function(theFit, theCoef=NULL, q=0.05, n=10, lFC=0.58, addGene=TRUE, miRNA=FALSE) {
  if(miRNA) {
    myDF <- topTags(theFit, n=nrow(theFit$coefficients))$table
    addGene=FALSE
    myDF <- myDF[myDF$FDR <= q,]
  } else {
    myDF <- topTable(theFit, coef=theCoef, number=nrow(theFit$coefficients))
    myDF <- myDF[myDF$adj.P.Val <= q,]
  }
  
  myDF <- myDF[abs(myDF$logFC) >= lFC,]
  if(addGene) { myDF <- insert_EnsemblGene(myDF) }
  myDF <- myDF[order(abs(myDF$logFC), decreasing=TRUE),]
  if (n != "all") {
    if (nrow(myDF) > n) {myDF <- myDF[1:n,]}
  }
  return(myDF)
}

# rename column names in place 
rename_colname_inplace = function(df,colname_vector){
  names(df)<-colname_vector
  return(df)
}

# set a column as rownames and remove that column from the df; default is setting the first column as the rowname
setrowname =  function(df,col = 1){
  if (is.numeric(col)) {
    # print('it is a number')
    rownames(df) = df[,col]
    df = df[,-col]
  } else {
    # print('it is a string')
    rownames(df) = df[,which(names(df) %in% c(col))]
    df = df[,-which(names(df) %in% c(col))]
  }
  return(df)
} 

# change a few column names in DESeq2 results to match to Limma results so that we can directly use the functions built for limma to plot volcano plots and such
rename_deseq2_results_to_limma = function(df){
  df <- df %>%
    dplyr::rename(logFC = log2FoldChange, P.Value = pvalue, adj.P.Val = padj)
  return(df)
}


# fGSEA wrapper function that parse the database, prepares input gene list, and run gsea function
run_gsea_v2 <- function(geneRank, qval = 0.05, refDB = msigDB,
                        enrichDB=c("KEGG", "REACTOME", "GO", "HALLMARK","IMMUNESIGDB","VAX")) {
  # # dummy variable for testing
  # geneRank = DEG.s.r.rank
  # qval = 0.05
  # refDB = msigDB
  # enrichDB = 'HALLMARK'
  
  # need to add the following line to the main script
  #msigbr to download all genes and their GO IDs #Ref:http://yulab-smu.top/biomedical-knowledge-mining-book/universal-api.html
  #msigDB <-  msigdbr::msigdbr(species = "Homo sapiens")
  
  # parse gene lists
  ## select the genes based on the desired gene list
  if (enrichDB == "KEGG") {
    db <- refDB %>% dplyr::filter(gs_subcat=="CP:KEGG") %>% 
      mutate(gs_name = str_remove(gs_name, "KEGG_")) %>% 
      dplyr::select(gs_name, gs_exact_source, gene_symbol) %>% #gene_symbol#ensembl_gene
      mutate(gs_name = str_replace_all(gs_name, "_", " "), gs_name = str_to_sentence(gs_name)) %>% 
      mutate(gs_name= paste(gs_exact_source, gs_name, sep = ":")) %>% 
      dplyr::select(!c(gs_exact_source))
    
  } else if (enrichDB == "REACTOME"){
    db <- refDB %>% dplyr::filter(gs_subcat=="CP:REACTOME") %>% 
      mutate(gs_name = str_remove(gs_name, "REACTOME_")) %>% 
      dplyr::select(gs_name, gs_exact_source, gene_symbol) %>% 
      mutate(gs_name = str_replace_all(gs_name, "_", " "), gs_name = str_to_sentence(gs_name)) %>% 
      mutate(gs_name= paste(gs_exact_source, gs_name, sep = ":")) %>% 
      dplyr::select(!c(gs_exact_source))
    
  }else if (enrichDB == "HALLMARK"){
    db <- refDB %>% dplyr::filter(gs_cat=="H") %>% 
      mutate(gs_name = str_remove(gs_name, "HALLMARK_")) %>% 
      dplyr::select(gs_name, gs_id, gene_symbol) %>% 
      mutate(gs_name = str_replace_all(gs_name, "_", " "), gs_name = str_to_sentence(gs_name)) %>% 
      mutate(gs_name= paste(gs_id, gs_name, sep = ":")) %>% 
      dplyr::select(!c(gs_id))
    
  }else if (enrichDB == "VAX"){
    db <- refDB %>% dplyr::filter(gs_subcat=="VAX") %>% 
      dplyr::select(gs_name, gs_id, gene_symbol) %>% 
      mutate(gs_name = str_replace_all(gs_name, "_", " "), gs_name = str_to_sentence(gs_name)) %>% 
      mutate(gs_name= paste(gs_id, gs_name, sep = ":")) %>% 
      dplyr::select(!c(gs_id))
    
  }else if (enrichDB == "IMMUNESIGDB"){
    db <- refDB %>% dplyr::filter(gs_subcat=="IMMUNESIGDB") %>% 
      dplyr::select(gs_name, gs_id, gene_symbol) %>% 
      mutate(gs_name = str_replace_all(gs_name, "_", " "), gs_name = str_to_sentence(gs_name)) %>% 
      mutate(gs_name= paste(gs_id, gs_name, sep = ":")) %>% 
      dplyr::select(!c(gs_id))
    
  } else if (enrichDB == "GO") {
    db <- refDB %>% dplyr::filter(gs_subcat=="GO:BP") %>%
      dplyr::select(gs_name, gs_exact_source, gene_symbol) %>% 
      mutate(gs_name = str_remove(gs_name, "GOBP_")) %>% 
      mutate(gs_name = str_replace_all(gs_name, "_", " "), gs_name = str_to_sentence(gs_name)) %>% 
      mutate(gs_name= paste(gs_exact_source, gs_name, sep = ":")) %>% dplyr::select(!c(gs_exact_source))
  }
  
  else{
    stop("Error: enrichDB must be one of 'KEGG', 'REACTOME', 'HALLMARK', 'VAX', 'IMMUNESIGDB', or 'GO'.")
  }
  
  # Group gene symbols by gs_name and remove duplicates (this will generate a list of pathway-named vectors of genes)
  geneSet.parsed <- lapply(split(db$gene_symbol, db$gs_name), unique)
  
  # Run GSEA using fgsea
  gsea.res = fgsea(pathways = geneSet.parsed,
                   stats = geneRank,
                   eps = 0)
  gsea.res = subset(gsea.res,padj<qval)
  gsea.res$direction = ifelse(gsea.res$NES>0,'Upregulation',
                              ifelse(gsea.res$NES<0,'Downregulation','0'))
  return(gsea.res)
}


# enrichment function that parse the database, prepares input gene list, and run clusterprofiler enricher function
geneEnrichments_v2 <- function(geneSet, bkgGenes = NULL, pval = 1, qval = 1, refDB = msigDB, refDBName = "org.Hs.eg.db",
                               enrichDB=c("KEGG", "REACTOME", "GO", "HALLMARK","IMMUNESIGDB","VAX"), targetedGO = NULL) {
  # need to add the following line to the main script
  #msigbr to download all genes and their GO IDs #Ref:http://yulab-smu.top/biomedical-knowledge-mining-book/universal-api.html
  #msigDB <-  msigdbr::msigdbr(species = "Homo sapiens")
  
  sigGenes <- unique(geneSet)
  backgroundGenes <- unique(bkgGenes)
  
  if (enrichDB == "KEGG") {
    kegg_db <- refDB %>% dplyr::filter(gs_subcat=="CP:KEGG") %>% 
      mutate(gs_name = str_remove(gs_name, "KEGG_")) %>% 
      dplyr::select(gs_name, gs_exact_source, gene_symbol) %>% #gene_symbol#ensembl_gene
      mutate(gs_name = str_replace_all(gs_name, "_", " "), gs_name = str_to_sentence(gs_name)) %>% 
      mutate(gs_name= paste(gs_exact_source, gs_name, sep = ":")) %>% 
      dplyr::select(!c(gs_exact_source))
    
    pathway <- clusterProfiler::enricher(sigGenes, TERM2GENE=kegg_db, universe = backgroundGenes, 
                                         pvalueCutoff = pval, qvalueCutoff = qval)
    
  } else if (enrichDB == "REACTOME"){
    reactome_db <- refDB %>% dplyr::filter(gs_subcat=="CP:REACTOME") %>% 
      mutate(gs_name = str_remove(gs_name, "REACTOME_")) %>% 
      dplyr::select(gs_name, gs_exact_source, gene_symbol) %>% 
      mutate(gs_name = str_replace_all(gs_name, "_", " "), gs_name = str_to_sentence(gs_name)) %>% 
      mutate(gs_name= paste(gs_exact_source, gs_name, sep = ":")) %>% 
      dplyr::select(!c(gs_exact_source))
    
    pathway <- clusterProfiler::enricher(sigGenes, TERM2GENE=reactome_db, universe = backgroundGenes, 
                                         pvalueCutoff = pval, qvalueCutoff = qval) 
    
  }else if (enrichDB == "HALLMARK"){
    hallmark_db <- refDB %>% dplyr::filter(gs_cat=="H") %>% 
      mutate(gs_name = str_remove(gs_name, "HALLMARK_")) %>% 
      dplyr::select(gs_name, gs_id, gene_symbol) %>% 
      mutate(gs_name = str_replace_all(gs_name, "_", " "), gs_name = str_to_sentence(gs_name)) %>% 
      mutate(gs_name= paste(gs_id, gs_name, sep = ":")) %>% 
      dplyr::select(!c(gs_id))
    
    pathway <- clusterProfiler::enricher(sigGenes, TERM2GENE=hallmark_db, universe = backgroundGenes, 
                                         pvalueCutoff = pval, qvalueCutoff = qval) 
    
  }else if (enrichDB == "VAX"){
    VAX_db <- refDB %>% dplyr::filter(gs_subcat=="VAX") %>% 
      dplyr::select(gs_name, gs_id, gene_symbol) %>% 
      mutate(gs_name = str_replace_all(gs_name, "_", " "), gs_name = str_to_sentence(gs_name)) %>% 
      mutate(gs_name= paste(gs_id, gs_name, sep = ":")) %>% 
      dplyr::select(!c(gs_id))
    
    pathway <- clusterProfiler::enricher(sigGenes, TERM2GENE=VAX_db, universe = backgroundGenes, 
                                         pvalueCutoff = pval, qvalueCutoff = qval) 
    
  }else if (enrichDB == "IMMUNESIGDB"){
    IMMUNESIGDB_db <- refDB %>% dplyr::filter(gs_subcat=="IMMUNESIGDB") %>% 
      dplyr::select(gs_name, gs_id, gene_symbol) %>% 
      mutate(gs_name = str_replace_all(gs_name, "_", " "), gs_name = str_to_sentence(gs_name)) %>% 
      mutate(gs_name= paste(gs_id, gs_name, sep = ":")) %>% 
      dplyr::select(!c(gs_id))
    
    pathway <- clusterProfiler::enricher(sigGenes, TERM2GENE=IMMUNESIGDB_db, universe = backgroundGenes, 
                                         pvalueCutoff = pval, qvalueCutoff = qval) 
    
  }else if (enrichDB == "IMMUNESIGDB_CD4Tcell_subset"){
    IMMUNESIGDB_db <- IMMUNESIGDB_CD4Tcell_subset %>% 
      dplyr::select(gs_name, gs_id, gene_symbol) %>% 
      mutate(gs_name = str_replace_all(gs_name, "_", " "), gs_name = str_to_sentence(gs_name)) %>% 
      mutate(gs_name= paste(gs_id, gs_name, sep = ":")) %>% 
      dplyr::select(!c(gs_id))
    
    pathway <- clusterProfiler::enricher(sigGenes, TERM2GENE=IMMUNESIGDB_db, universe = backgroundGenes, 
                                         pvalueCutoff = pval, qvalueCutoff = qval) 
    
  }
  else{
    if(!is.null(targetedGO)){ 
      pat <- c("GOBP_", "REACTOME_", "WP_")
      GO_db <- refDB  %>% dplyr::filter(str_detect(gs_name, paste(targetedGO, collapse = "|"))) %>% 
        mutate(gs_name = str_remove(gs_name, paste(pat, collapse = "|"))) %>% 
        mutate(gs_name = str_replace_all(gs_name, "_", " "), gs_name = str_to_sentence(gs_name)) %>% 
        mutate(gs_name= paste(gs_exact_source, gs_name, sep = ":")) %>%
        dplyr::select(gs_name, gene_symbol)
    }else{
      GO_db <- refDB %>% dplyr::filter(gs_subcat=="GO:BP") %>%
        dplyr::select(gs_name, gs_exact_source, gene_symbol) %>% 
        mutate(gs_name = str_remove(gs_name, "GOBP_")) %>% 
        mutate(gs_name = str_replace_all(gs_name, "_", " "), gs_name = str_to_sentence(gs_name)) %>% 
        mutate(gs_name= paste(gs_exact_source, gs_name, sep = ":")) %>% dplyr::select(!c(gs_exact_source))
    }
    pathway <- clusterProfiler::enricher(sigGenes, TERM2GENE=GO_db, universe = backgroundGenes, 
                                         pvalueCutoff = pval, qvalueCutoff = qval)
    
  }
  
  if(is.null(pathway)==TRUE){
    pathwayData <- as.data.frame(matrix(nrow=0, ncol = 3))
  }else{
    pathwayData <- pathway %>% as.data.frame(.) %>% arrange(desc(Count))
  }
  # 
  #   if(nrow(pathwayData) != 0) {
  #     pathwayData <- pathwayData %>%
  #       mutate(GeneRatio2 = unname(sapply(GeneRatio, function(x) eval(parse(text = x)))),
  #              BgRatio2 = unname(sapply(BgRatio, function(x) eval(parse(text = x)))), GeneEnrichment = GeneRatio2/BgRatio2) %>%
  #       dplyr::select(!c(GeneRatio2, BgRatio2)) %>% dplyr::filter(Count>=3)
  #   }
  return(pathwayData)
}
