# compares expression of genes for specified contrast, returns results of glmTreat function and TopTags 
# differentially expressed genes, prints summary and top genes
# can produce an MD and Volcano plot of data


#' Compare expression in two strains specified by contrats using edgeR package glmTreat
#'
#' @param fit fit object made using edgeR package
#' @param contrast_comp strains to compare, have to present in the design
#' @param design design specifying which libraries belong to which strain
#' @param genes data_frame containing genes data
#' @param plot_md should md plot be shown? defaut TRUE
#' @param plot_volcano should volcano plot be shown? default TRUE
#' @param complete_list should it return complete list of genes or only those differentially expressed? default FALSE
#' @param fold_change fold change used in glmTreat function
#' @param toptags_print how many genes should be printed?
#'
#' @return result of topTags function for all or differentially expressed genes
#' @export
#'
#' @examples
make_comparison <- function(fit, contrast_comp, design, genes, plot_md = TRUE, plot_volcano = TRUE,
                            complete_list = FALSE, fold_change = 1.5, toptags_print = 5){
  
  fit$genes <- dplyr::left_join(fit$genes, genes[,c(1,5,6,2)], by = c('genes' = 'gene'))
  
  contr <- limma::makeContrasts(contrasts = contrast_comp, levels = design)
  
  res <- edgeR::glmTreat(fit, contrast = contr, lfc = log2(fold_change))
  #print("Top 5 differentially expressed genes")
  #print(topTags(res, n = toptags_print))
  
  is.de <- edgeR::decideTestsDGE(res)
  
  x <- summary(is.de)
  print("Summarized differential expresion")
  print(summary(is.de))
  
  if(complete_list == TRUE){
    geny_de <- edgeR::topTags(res, n = sum(x))
  }else{
    geny_de <- edgeR::topTags(res, n = x[1] + x[3])
  }
  geny_de$table$contrast <- contrast
  
  if(plot_md == TRUE){
    limma::plotMD(res, status = is.de, values = c(1,-1), col = c('red', 'blue'), legend = 'topright')
  }
  
  if(plot_volcano == TRUE){
    
    geny_de_volcano <- edgeR::topTags(res, n = sum(x))
    res_plot <- geny_de_volcano$table
    
    top_genes <- res_plot$name[1:20]
    plot <- EnhancedVolcano::EnhancedVolcano(res_plot,
                                             lab = res_plot$name,
                                             selectLab = top_genes,
                                             x = 'logFC',
                                             y = 'FDR',
                                             #xlim = c(-4, 6),
                                             ylim = c(0, 10),
                                             pCutoff = 0.05,
                                             FCcutoff = 1.5, 
                                             title = contrast,
                                             subtitle = '',
                                             drawConnectors = TRUE
    )
    print(plot)
  }
  
  
  
  wyniki <- list(res, geny_de)
  
  return(wyniki)
  
}


#' will produce a data.frame with gene names, protein id and function and a note from genbank file
#'
#' @param genbank genbank file read in using genbankr package
#' @param add_names add names from uniprot file? default TRUE
#' @param uniprot_file uniprot file containing protein name etc.
#'
#' @return data frame containing gene id, locus_tag, product, protein_id and note from genbank file
#' @export
#'
#' @examples
get_transcripts_from_gb <- function(genbank, add_names = TRUE, uniprot_file){
  
  trans <- transcripts(genbank)
  
  res <- data.frame(gene = trans@elementMetadata@listData$gene,
                    locus_tag = trans@elementMetadata@listData$locus_tag,
                    product = trans@elementMetadata@listData$product,
                    protein_id = trans@elementMetadata@listData$protein_id,
                    note = trans@elementMetadata@listData$note,
                    stringsAsFactors = FALSE)
  
  if(add_names == TRUE){
    # load gene names from uniprot
    uniprot <- read_delim(uniprot_file, "\t", escape_double = FALSE, trim_ws = TRUE)
    colnames(uniprot) <- c('entry', 'entry_name', 'status', 'protein_name', 'gene_names', 'organism', 'length', 'gene')
    strsplit(uniprot$gene_names, ' ') -> y2
    names <- character(length = nrow(uniprot))
    for(i in 1:length(y2)){
      names[i] <- y2[[i]][1]
    }
    uniprot$name <- names
    
    uniprot <- uniprot[, c(8, 9, 4)]
    
    res %>% left_join(uniprot) -> res
  }
  
  return(res)
}

# 

#' Plots heatmap for n genes most differentially expressed in comparison specified by contrast, 
#' specific strains can be chosen in groups
#'
#' @param data_edger experiment data from edgeR analysis
#' @param comparison list containing table with pvalues for gene ordering
#' @param groups specific strains to show on heatmap
#' @param n number of genes shown on heatmap
#'
#' @return heatmap 
#' @export
#'
#' @examples
draw_heatmap <- function(data_edger, comparison, groups, n = 50){
  
  keep <- which(data_edger$samples$group %in% groups)
  
  keep_data <- data_edger[,keep]
  
  logCPM <- edgeR::cpm(keep_data, prior.count = 2, log = TRUE)
  rownames(logCPM) <- keep_data$genes$genes
  colnames(logCPM) <- paste(keep_data$samples$group, 1:3, sep = '-')
  
  #o <- order(comparison[[1]]$table$PValue)
  
  o <- order(comparison$table$PValue)
  
  logCPM <- logCPM[o[1:n],]
  
  logCPM <- t(scale(t(logCPM)))
  
  col.pan <- gplots::colorpanel(100, 'blue', 'white', 'red')
  
  p <- heatmap(logCPM, col = col.pan, margins = c(9,5), cexCol = 0.75, scale = 'none')
  
  return(p)
  
}

compare_strains <- function(list_res, grupy){
  list_de_genes <- list()
  for(i in 1:length(list_res)){
    de_genes <- list_res[[i]][[2]]@.Data[[1]]
    de_genes$grupa <- grupy[i]
    list_de_genes[[i]] <- de_genes
  }
  table_de_genes <- do.call(rbind, list_de_genes)
  
  table_de_genes %>% mutate(a = 0, b = 0, c = 0, d = 0, e = 0, f = 0, g = 0) -> table_de_genes
  
  colnames(table_de_genes)[(ncol(table_de_genes)-6):ncol(table_de_genes)] <- c(grupy[1], grupy[2], grupy[3],
                                                                               paste(grupy[1], grupy[2], sep = '_'),
                                                                               paste(grupy[1], grupy[3], sep = '_'),
                                                                               paste(grupy[2], grupy[3], sep = '_'),
                                                                               paste(grupy[1], grupy[2], grupy[3], sep = '_'))
  
  table_de_genes %>% dplyr::count(genes) -> table_counts
  genes_same <- filter(table_counts, n==3)$genes
  
  #table_de_genes[,17] <- ifelse(table_de_genes[,1] %in% genes_same, 1, 0)
  
  for(i in 1:nrow(table_de_genes)){
    count <- filter(table_counts, genes == table_de_genes$genes[i])$n
    
    if(count == 3){
      table_de_genes[i,18] <- 1
    }
    if(count == 1){
      if(table_de_genes$grupa[i] == grupy[1]){table_de_genes[i,12] <- 1}
      if(table_de_genes$grupa[i] == grupy[2]){table_de_genes[i,13] <- 1}
      if(table_de_genes$grupa[i] == grupy[3]){table_de_genes[i,14] <- 1}
    }
    if(count == 2){
      jakie_grupy <- filter(table_de_genes, genes == table_de_genes$genes[i])$grupa
      if(all(jakie_grupy == grupy[c(1,2)])){table_de_genes[i,15] <- 1}
      if(all(jakie_grupy == grupy[c(1,3)])){table_de_genes[i,16] <- 1}
      if(all(jakie_grupy == grupy[c(2,3)])){table_de_genes[i,17] <- 1}
    }
    
  }
  
  res <- table_de_genes[,12:18]  
  
  res %>% summarise_all(sum) -> res
  res[4:6] <- res[4:6]/2
  res[7] <- res[7]/3
  
  print(res)
  
  return(table_de_genes)
}


#' Reads cluster file prepared using clust program
#'
#' @param file clust result file
#' @param genes_list list of genes 
#'
#' @return data frame with all genes with assigned cluster or NA
#' @export
#'
#' @examples
read_cluster_file <- function(file, genes_list){
  
  data <- read_delim(file, 
                     "\t", escape_double = FALSE, trim_ws = TRUE, 
                     skip = 1)
  
  colnames(data) <- paste0('cluster_', 1:ncol(data))
  data %>% gather(key = 'cluster', value = 'gene') %>% filter(!is.na(gene)) -> data
  
  genes <- as.data.frame(genes_list, stringsAsFactors = FALSE)
  
  data %>% full_join(genes, by = c('gene' = 'genes_list')) %>% arrange(gene) -> data
  
  return(data)
}


plot_cluster_results <- function(clust, norm, genes, order=NULL, plot = TRUE){
  
  
  clust %>% filter(!is.na(cluster)) %>%
    left_join(genes) -> clust
  
  norm %>% filter(Genes %in% clust$gene) %>% 
    pivot_longer(-Genes, names_to = 'strain', 
                 values_to = 'nor_expression') ->
    data
  
  data %>% left_join(clust, by = c('Genes' = 'gene')) -> data
  
  if(!is.null(order)){
    data$strain <- factor(data$strain, levels = order)
  }
  
  p <- ggplot(data, aes(x = strain, y = nor_expression))
  p <- p + geom_jitter(alpha = 0.5)+ coord_flip()+
    facet_wrap(~cluster)+
    stat_summary(fun = 'mean', geom = 'point', color = 'red3')+
    theme_bw()
  if(plot){
    print(p)
  }
  return(list(p, data))
}

check_gene_cluster <- function(clust, genes, genes_list){
  
  clust_data <- read_cluster_file(clust, genes_list)
  
  #wynik <- subset(clust_data, gene %in% genes)
  
  wynik <- genes %in% na.omit(clust_data)$gene
  
  names(wynik) <- genes
  
  wynik2 <- clust_data$cluster[clust_data$gene %in% genes]
  names(wynik2) <- genes
  #wynik <- ifelse(wynik, clust_data$cluster[which(clust_data$gene == names(wynik))],  wynik)
  
  return(list(wynik, wynik2))
  
}


combine_rnaseq_kegg <- function(rnaseq_data, kegg_data, level = 'C', wide = TRUE, 
                                filter, direction = 'both', logFC_level = 1){
  
  if(wide == TRUE){
    rnaseq_data %>%
      gather(key = 'contrast', value = 'logFC', -genes, -name, -product) ->
      rnaseq_data
  }
  
  if(direction == 'down'){
    rnaseq_data %>% filter(logFC < 0) -> rnaseq_data
  }
  if(direction == 'up'){
    rnaseq_data %>% filter(logFC > 0) -> rnaseq_data
  }
  
  rnaseq_data %>% 
    left_join(kegg_data, by = c('genes' = 'gene')) %>%
    filter(contrast == filter, abs(logFC) >= logFC_level) -> rnaseq_data
  
  if(level == 'C'){ 
    kegg_data %>% 
      filter(grepl(x = gene, pattern = 'SCO[1-9]')) %>%
      group_by(B, C) %>%
      distinct(gene, .keep_all = TRUE) %>%
      summarize(n_kegg = n()) ->
      kegg_count_BC
    
    rnaseq_data %>%
      filter(!is.na(logFC)) %>%
      group_by(B, C) %>%
      distinct(genes, .keep_all = TRUE) %>%
      count() %>%
      left_join(kegg_count_BC) %>%
      mutate(percent = n/n_kegg) %>%
      arrange(-percent) -> result
    
  }
  if(level == 'B'){ 
    kegg_data %>% 
      filter(grepl(x = gene, pattern = 'SCO[1-9]')) %>%
      group_by(B) %>%
      distinct(gene, .keep_all = TRUE) %>%
      summarize(n_kegg = n()) ->
      kegg_count_B
    
    rnaseq_data %>%
      filter(!is.na(logFC)) %>%
      group_by(B) %>%
      distinct(genes, .keep_all = TRUE) %>%
      count() %>%
      left_join(kegg_count_B) %>%
      mutate(percent = n/n_kegg) %>%
      arrange(-percent) -> result
  }
  
  return(result)
}




plot_genes_coverage <- function(genes_positions, strains, sco_start, sco_end, ylim = NA, flank = 500,
                                chromosome = "NC_003888.3") {
  # start <-1305967
  # end <-1310596
  
  gene_start <- genes_positions %>% filter(name == sco_start)
  gene_end <- genes_positions %>% filter(name == sco_end)
  
  start_g <- gene_start$start
  end_g <- gene_end$end
  
  genes_all <- genes_positions %>% filter(start <= end_g + flank & start >= start_g - flank |
                                            end <= end_g + flank & end >= start_g - flank) %>%
    mutate(strand = ifelse(strand == '+', 1, -1))
  
  print(nrow(genes_all))
  
  region <- GRanges(Rle(c(chromosome), c(1)),  IRanges(start = start_g - flank, end = end_g + flank))
  
  print(region)
  
  x <- plot_bam_coverage2(bam.files = strains$file, region, lib.sizes = strains$libsize, 
                          names = strains$strains)
  
  if(is.na(ylim)){
    ylim <- max(x[[1]]$cov)
  }
  
  p <- x[[2]]
  
  p <- p + coord_cartesian(ylim = c(0, ylim))
  
  p_genes <- ggplot(genes_all, aes(xmin = start, xmax = end, fill = factor(strand), 
                                   forward = strand, label = name, y = '')) +
    geom_gene_arrow()+
    geom_gene_label()+theme_genes()+
    theme(legend.position = 'none')+
    ylab('')+
    coord_cartesian(xlim = c(start(region), end(region)))
  
  #print(p)
  
  print(p + p_genes + plot_layout(ncol = 1, heights = c(15, 1)))
  
}





make_table_rnaseq <- function(gene_i, data_all, genes_information) {
  gen <- data_all %>% filter(genes == gene_i)
  gen_info <- genes_information %>% filter(gene == gene_i)
  print(gen_info)
  
  gen <- gen[,c(1:9, 14:19)]
  
  gen %>% select(-1:-3) %>% gather() %>% separate(key, into = c('strain', 'time', 'NaCl'), sep = '_') %>%
    mutate(time = ifelse(is.na(NaCl), time, paste(time, NaCl, sep = '_'))) %>% spread(key = time, value = value) %>% select(1,3,5,4,6) -> gen_table
  
  gen_table[c(2,4,6),1:3] %>% left_join(gen_table[c(1,3,5), c(1,4,5)])  %>% mutate_if(is.numeric, round, digits = 2) -> gen_table
  
  return(gen_table)
  
  #kable(gen_table)
  
}





#' function reads from antismash files query genes with best subject matching genes
#'
#' @param file antismash result file in txt
#'
#' @return data frame with antismash results with gene information and blast data 
#' (% identity, blast score evalue) for each found similar gene 
#' @export
#'
#' @examples
read_antismash <- function(file){
  
  # read in raw file
  smash <- read_lines(file)
  
  # start of genes table for results
  line_query <- which(grepl('Table of genes, locations, strands and annotations of query cluster:', smash))
  
  # empty line separate the file into parts
  line_empty <- which(grepl('^$', smash))
  
  # end of result table
  line_query_last <- line_empty[line_empty>line_query][1]
  
  # prepare result table from the query genes
  query <- smash[(line_query+1):(line_query_last-1)]
  as.data.frame(query) %>% separate(query, sep = '\t', 
                                    into = c('gene', 'start', 'end', 'strand', 'annotation', 'x')) %>%
    select(-x)-> 
    result_table
  
  cluster <- sub('.*_', '', file)
  cluster <- sub('.txt', '', cluster)
  
  result_table$cluster <- cluster
  
  # find the best cluster (will have number 1)
  subject_name_line <- which(grepl('^1\\.', smash))[2]+1
  
  # if file is empty stop function
  if(is.na(subject_name_line)){
    cat('No significant cluster found \n')
    return(NULL)
  }
  
  # cluster product name
  subject_name <- sub('Source: ', '', smash[subject_name_line])
  result_table$subject <- subject_name
  
  # make table with annotated subject genes
  subject_table <- smash[(line_empty[line_empty > subject_name_line][1] + 2) : (line_empty[line_empty > subject_name_line][2]-1)]
  as.data.frame(subject_table) %>% separate(subject_table, sep = '\t', 
                                            into = c('gene_subject', 'gene_name', 'start_subject', 
                                                     'end_subject','strand_subject', 'annotation_subject')) %>%
    select(gene_subject, annotation_subject) -> subject_table
  
  # make table with query genes connected to subject genes
  subject_query_table <- smash[(line_empty[line_empty > subject_name_line][2] + 2) : (line_empty[line_empty > subject_name_line][3]-1)]
  as.data.frame(subject_query_table) %>% separate(subject_query_table, sep = '\t', 
                                                  into = c('gene', 'gene_subject', '%identity', 
                                                           'blast_score','%coverage', 'e-value')) -> subject_query_table
  # join all tables together, return result table
  result_table %>% left_join(subject_query_table) %>% 
    mutate(subject = ifelse(is.na(gene_subject), NA, subject)) %>%
    left_join(subject_table) -> result_table
  
  return(result_table)
  
}


plot_fitted_genes <- function(fitted, genes_positions, strains,
                              gene_start, gene_end, flank = 500, chromosome = "NC_003888.3"){
  
  g_start <- genes_positions %>% filter(name == gene_start)
  g_end <- genes_positions %>% filter(name == gene_end)
  
  start_g <- g_start$start
  end_g <- g_end$end
  
  genes_all <- genes_positions %>% filter(start <= end_g + flank & start >= start_g - flank |
                                            end <= end_g + flank & end >= start_g - flank) %>%
    mutate(strand = ifelse(strand == '+', 1, -1))
  
  fitted_genes <- fitted %>% filter(gene %in% genes_all$name, strain %in% strains) %>% 
    left_join(genes_all, by = c('gene' = 'name')) %>%
    mutate(name = factor(gene.y, levels = unique(gene.y)))
  
  fitted_genes %>% group_by(strain, name) %>%
    summarise(mean = mean(count),
              sd = sd(count),
              moe = sd) -> fitted_summary
  
  p <- ggplot(fitted_genes, aes(y = strain, x = count))
  p + 
    theme_bw()+
    stat_confidence_density(data = fitted_summary, aes(x = mean, y = strain, 
                                                       moe = moe,  fill = stat(ndensity)),
                            fill = 'grey50')+
    geom_point(position = position_dodge(width = 0.5))+
    facet_wrap(~name, ncol = 1)+
    xlab('normalized counts')
  
}

plot_rpkm_genes <- function(fitted, genes_positions, strains, fit,
                            gene_start, gene_end, flank = 0, chromosome = "NC_003888.3",
                            plot_type = 'gene_separate', log = FALSE, format = 'RPKM'){
  
  g_start <- genes_positions %>% dplyr::filter(name == gene_start)
  g_end <- genes_positions %>% dplyr::filter(name == gene_end)
  
  start_g <- g_start$start
  end_g <- g_end$end
  
  genes_all <- genes_positions %>% dplyr::filter(start <= end_g + flank & start >= start_g - flank |
                                                   end <= end_g + flank & end >= start_g - flank) %>%
    dplyr::mutate(strand = ifelse(strand == '+', 1, -1),
                  strand_plot = ifelse(strand == 1, TRUE, FALSE))
  
  if(format == 'RPKM'){
    
    fitted_genes <- fitted %>% dplyr::filter(gene %in% genes_all$name, strain %in% strains) %>% 
      dplyr::left_join(genes_all, by = c('gene' = 'name')) %>%
      dplyr::mutate(name = factor(gene.y, levels = unique(gene.y)))
    
  } else if(format == 'TPM'){
    
    counts <- as.data.frame(fit$counts)
    
    counts$gene <- rownames(counts)
    
    counts %>% dplyr::left_join(genes_positions, by = c('gene' = 'name')) %>% 
      dplyr::select(-gene.y) %>%
      tidyr::gather(key = 'strain', value = 'count', -gene, -start, -end, -width, -strand) %>%
      dplyr::group_by(strain)%>%
      dplyr::mutate(RPK = count/(width/1000),
                    scaling_factor = sum(RPK)/1000000,
                    count = RPK / scaling_factor,
                    strain = sub('_[123]', '', strain)) %>%
      dplyr::group_by(strain, gene) %>%
      dplyr::filter(gene %in% genes_all$name, strain %in% strains) %>%
      dplyr::summarise(count = mean(count),
                       start = mean(start),
                       end = mean(end)) -> fitted_genes
    
    
  }
  
  if(plot_type == 'gene_separate'){
    p <- ggplot2::ggplot(fitted_genes, ggplot2::aes(y = strain, x = count, yend = strain, xend = 0))
    p <- p + 
      ggplot2::theme_bw()+
      ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.5))+
      ggplot2::geom_segment(position = ggplot2::position_dodge(width = 0.5))+
      ggplot2::facet_wrap(~name, ncol = 1)+
      ggplot2::xlab('RPKM')
    if(log){
      p <- p + ggplot2::scale_x_log10()
    }
    
  }
  
  if(plot_type == 'gene_position'){
    p <- ggplot2::ggplot(fitted_genes, ggplot2::aes(x = (start+end)/2, 
                                                    y = count, 
                                                    xmax = (start+end)/2, 
                                                    ymin = 0,
                                                    ymax = count, 
                                                    xmin = (start+end)/2,
                                                    color = strain, 
                                                    fill = strain))
    p <- p+ 
      ggplot2::theme_bw()+
      ggplot2::geom_point(position = ggplot2::position_dodge(width = 350))+
      ggplot2::geom_linerange(position = ggplot2::position_dodge(width = 350))+
      ggplot2::xlab('Genome posiiton [bp]')+
      ggplot2::ylab('RPKM')+
      ggplot2::xlim(genes_all$start[1], genes_all$end[nrow(genes_all)])+
      ggplot2::scale_color_brewer(palette = 'Set1')+
      ggplot2::theme(text = ggplot2::element_text(size = 15),
                     legend.position = 'top')
    
    if(log){
      p <- p + ggplot2::scale_y_log10()
    }
    
    p_genes <- ggplot2::ggplot(genes_all, ggplot2::aes(xmin = start, xmax = end, fill = factor(strand), 
                                                       forward = strand_plot, label = gene, y = '')) +
      gggenes::geom_gene_arrow(arrowhead_height = grid::unit(9, 'mm'),
                               arrow_body_height = grid::unit(7, 'mm'))+
      gggenes::geom_gene_label(grow = FALSE, reflow = TRUE, height = grid::unit(2, 'cm'))+
      gggenes::theme_genes()+
      ggplot2::theme(legend.position = 'none', text = ggplot2::element_text(size = 12))+
      ggplot2::ylab('')
    
    #print(p)
    
    p <- p + p_genes + patchwork::plot_layout(ncol = 1, heights = c(10, 1))
    
  }
  
  return(p)
}


draw_heatmap_cluster <- function(fit, gene_start, gene_end, strains, genes_positions){
  
  fit_counts <- fit$fitted.values
  
  logCPM <- edgeR::cpm(fit_counts, prior.count = 2, log = TRUE)
  rownames(logCPM) <- fit$genes$genes
  
  logCPM <- t(scale(t(logCPM)))
  
  logCPM_table <- as.data.frame(logCPM)
  
  logCPM_table$gene <- rownames(logCPM_table)
  
  logCPM_table %>% tidyr::gather(strain, cpm, -gene) %>%  
    dplyr::mutate(replicate = sub('.*_', '', strain),
                  strain = sub('_[123]', '', strain)) -> logCPM_table
  
  g_start <- genes_positions %>% dplyr::filter(name == gene_start)
  g_end <- genes_positions %>% dplyr::filter(name == gene_end)
  
  start_g <- g_start$start
  end_g <- g_end$end
  
  genes_all <- genes_positions %>% dplyr::filter(start <= end_g & start >= start_g |
                                                   end <= end_g & end >= start_g) 
  
  cpm_genes <- logCPM_table %>% dplyr::filter(gene %in% genes_all$name, strain %in% strains) %>% 
    dplyr::group_by(strain, gene) %>% dplyr::summarise(cpm = mean(cpm)) 
  
  p <- ggplot2::ggplot(cpm_genes, ggplot2::aes(y = gene, strain, fill = cpm))
  
  p + ggplot2::geom_tile()+
    ggplot2::theme_minimal()+
    ggplot2::theme(axis.title = ggplot2::element_blank(),
                   axis.text.x.top = ggplot2::element_text(angle = 90))+
    ggplot2::scale_fill_distiller(type = 'div', palette = 5)+
    #scico::scale_fill_scico(palette = 'vikO')+ # tried scico palette, but I like the colorbrewer more
    ggplot2::scale_x_discrete(position = 'top')+
    ggplot2::theme(text = ggplot2::element_text(size = 15))
  
}


plot_volcano <- function(fit, contrast_comp, design, genes,
                         fold_change = 1.5, top_tag_plot = 20){
  
  fit$genes <- dplyr::left_join(fit$genes, genes[,c(1,5,6,2)], by = c('genes' = 'gene'))
  
  contr <- limma::makeContrasts(contrasts = contrast_comp, levels = design)
  
  res <- edgeR::glmTreat(fit, contrast = contr, lfc = log2(fold_change))
  
  is.de <- edgeR::decideTestsDGE(res)
  x <- summary(is.de)
  
  geny_de_volcano <- edgeR::topTags(res, n = sum(x))
  res_plot <- geny_de_volcano$table
  
  maks <- max(-log10(res_plot$FDR)) + 0.5
  
  top_genes <- res_plot$name[1:top_tag_plot]
  plot <- EnhancedVolcano::EnhancedVolcano(res_plot,
                                           lab = res_plot$name,
                                           selectLab = top_genes,
                                           x = 'logFC',
                                           y = 'FDR',
                                           #xlim = c(-4, 6),
                                           ylim = c(0, maks),
                                           pCutoff = 0.05,
                                           FCcutoff = fold_change, 
                                           title = contrast_comp,
                                           subtitle = '',
                                           drawConnectors = TRUE,
                                           labSize = 5
  )
  
  return(list(plot, res))
}


plot_logfc_genes <- function(logFC, genes_positions, contrasts,
                             gene_start, gene_end, flank = 0, chromosome = "NC_003888.3",
                             plot_type = 'gene_separate'){
  
  g_start <- genes_positions %>% dplyr::filter(name == gene_start)
  g_end <- genes_positions %>% dplyr::filter(name == gene_end)
  
  start_g <- g_start$start
  end_g <- g_end$end
  
  genes_all <- genes_positions %>% dplyr::filter(start <= end_g + flank & start >= start_g - flank |
                                                   end <= end_g + flank & end >= start_g - flank) %>%
    dplyr::mutate(strand = ifelse(strand == '+', 1, -1),
                  strand_plot = ifelse(strand == 1, TRUE, FALSE))
  
  #logFC <- data_adpA$data_all_notsig
  
  names(logFC)[[3]] <- 'locus_tag'
  
  logFC %>% dplyr::filter(locus_tag %in% genes_all$name, contrast %in% contrasts) -> logFC

    # logFC %>% tidyr::pivot_longer(cols = 4:ncol(logFC), names_to = 'contrast', values_to = 'logFC') -> logFC
  
  logFC %>% dplyr::left_join(genes_all, by = c('locus_tag' = 'name')) %>%
    dplyr::mutate(name = factor(gene.y, levels = unique(gene.y)),
                  significant = ifelse(FDR <= 0.05, TRUE, FALSE)) -> logFC
  
  p <- ggplot2::ggplot(logFC, ggplot2::aes(x = (start+end)/2, 
                                           y = logFC, 
                                           xmax = (start+end)/2, 
                                           ymin = 0,
                                           ymax = logFC, 
                                           xmin = (start+end)/2,
                                           color = contrast, 
                                           shape = significant))
  p <- p+ 
    ggplot2::theme_bw()+
    ggplot2::geom_point(position = ggplot2::position_dodge(width = 350), size = 3)+
    ggplot2::geom_linerange(position = ggplot2::position_dodge(width = 350))+
    ggplot2::xlab('Genome posiiton [bp]')+
    ggplot2::ylab('LogFC')+
    ggplot2::xlim(genes_all$start[1], genes_all$end[nrow(genes_all)])+
    #ggplot2::scale_color_brewer(palette = 'Set1')+
    ggplot2::theme(text = ggplot2::element_text(size = 15),
                   legend.position = 'top')+
    ggplot2::scale_shape_manual(values = c(1, 19))
  
  p_genes <- ggplot2::ggplot(genes_all, ggplot2::aes(xmin = start, xmax = end, fill = factor(strand), 
                                                     forward = strand_plot, label = gene, y = '')) +
    gggenes::geom_gene_arrow(arrowhead_height = grid::unit(9, 'mm'),
                             arrow_body_height = grid::unit(7, 'mm'))+
    gggenes::geom_gene_label(grow = FALSE, reflow = TRUE, height = grid::unit(2, 'cm'))+
    gggenes::theme_genes()+
    ggplot2::theme(legend.position = 'none', text = ggplot2::element_text(size = 12))+
    ggplot2::ylab('')
  
  #print(p)
  
  p <- p + p_genes + patchwork::plot_layout(ncol = 1, heights = c(10, 1))
  
  return(p)
}