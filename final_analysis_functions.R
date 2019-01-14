plot_ordered_heatmap <- function(hm_df = data.frame())  {
  temp_mat <- reshape2::acast(hm_df, cluster_id ~ GENENAME, value.var = 'prop')
  col_order <- colnames(temp_mat)[hclust(dist(t(temp_mat)))$order]
  row_order <- rownames(temp_mat)[hclust(dist(temp_mat))$order]
  ggplot2::ggplot({dplyr::mutate(hm_df, cluster_id = factor(cluster_id, levels = row_order, ordered = TRUE), GENENAME = factor(GENENAME, levels = col_order, ordered = TRUE))}, aes(x = GENENAME, y = cluster_id, fill = prop)) + geom_tile(color = 'white') + 
    scale_fill_gradient2(low = 'white', mid = 'deepskyblue', high = 'darkorange', midpoint = 0.5) + theme(axis.text.x = element_text(angle = -90)) + 
    theme(text = element_text(size = 10, vjust = 0.5), axis.text = element_text(size = 6, hjust = 0), legend.position = 'bottom') +
    labs(fill = 'Proportion', x = 'Gene Name', y = 'Cluster ID')
}

get_prop_mat_from_top_prop_dif <- function(sample_to_cluster_df = data.frame(), expression_mat = matrix(), top_n = 30)  {
  print(top_n)
  top_n_df <- get_enriched_genes_in_each_cluster(sample_to_cluster_df, expression_mat, top_n)
  ensembl_ids <- unique(top_n_df$ENSEMBL)
  prop_df <- get_prop_expressed_df(sample_to_cluster_df, expression_mat[,ensembl_ids])
  prop_mat <- reshape2::acast(prop_df, cluster_id ~ GENENAME, value.var = 'prop')
  return(prop_mat)
}

get_prop_expressed_df <- function(sample_to_cluster_df = data.frame(), expression_mat = matrix())  {
  cluster_counts  <- table(sample_to_cluster_df$cluster_id)
  cluster_ids <- unique(sample_to_cluster_df$cluster_id)
  expr_count_per_cluster_list <- lapply(cluster_ids, function(.cluster_id)  {
    sample_names_in_current_cluster <- dplyr::filter(sample_to_cluster_df, cluster_id == .cluster_id)$sample_name
    current_cluster_count <- apply(expression_mat[sample_names_in_current_cluster,], 2, function(.col)  {
      sum(.col > 0)
    })
    print("One cluster count completed")
    return(current_cluster_count)
  })
  #return(expr_count_per_cluster_list)
  cluster_count_mat <- do.call(rbind, expr_count_per_cluster_list)
  rownames(cluster_count_mat) <- cluster_ids
  gene_counter <- 1
  exp_genes_per_cluster_mat <- apply(cluster_count_mat, 2, function(.col)  {
    props <- sapply(1:length(.col), function(.ind)  {
      return(.col[.ind]/cluster_counts[cluster_ids[.ind]])
    })
    if(gene_counter %% 1000 == 0)  {
      print("1000 genes processed")
    }
    gene_counter <<- gene_counter + 1
    return(props)
  })
  props_df <- plyr::ldply(1:nrow(exp_genes_per_cluster_mat), function(.cluster_ind)  {
    .cluster <- exp_genes_per_cluster_mat[.cluster_ind,]
    #top_n_inds <- order((.cluster), decreasing = TRUE)[1:n_top_genes]
    return(data.frame(cluster_id = cluster_ids[.cluster_ind], ENSEMBL = colnames(exp_genes_per_cluster_mat), prop = .cluster, stringsAsFactors = FALSE))
  })
  ensembl_to_gene_name_df <- AnnotationDbi::select(EnsDb.Mmusculus.v79, keys=unique(props_df$ENSEMBL), keytype = 'GENEID', columns = c("GENENAME", "GENEID")) %>% dplyr::rename(ENSEMBL =GENEID)
  return(dplyr::left_join(props_df, ensembl_to_gene_name_df))
}

get_enriched_genes_in_each_cluster <- function(sample_to_cluster_df = data.frame(), expression_mat = matrix(), top_n = 30)  {
  
  if(top_n == 'all')  {
    n_top_genes <- ncol(expression_mat)
  }  else if(is.numeric(top_n) && top_n <= ncol(expression_mat))  {
    n_top_genes <- top_n
  }  else  {
    stop("Error!")
  }
  
  cluster_counts  <- table(sample_to_cluster_df$cluster_id)
  cluster_ids <- unique(sample_to_cluster_df$cluster_id)
  expr_count_per_cluster_list <- lapply(cluster_ids, function(.cluster_id)  {
    sample_names_in_current_cluster <- dplyr::filter(sample_to_cluster_df, cluster_id == .cluster_id)$sample_name
    current_cluster_count <- apply(expression_mat[sample_names_in_current_cluster,], 2, function(.col)  {
      sum(.col > 0)
    })
    print("One cluster count completed")
    return(current_cluster_count)
  })
  #return(expr_count_per_cluster_list)
  cluster_count_mat <- do.call(rbind, expr_count_per_cluster_list)
  rownames(cluster_count_mat) <- cluster_ids
  gene_counter <- 1
  exp_genes_per_cluster_mat <- apply(cluster_count_mat, 2, function(.col)  {
    prop_difs <- sapply(1:length(.col), function(.ind)  {
      return((.col[.ind]/cluster_counts[cluster_ids[.ind]]) - (sum(.col[-.ind])/sum(cluster_counts[names(cluster_counts) != cluster_ids[.ind]])))
    })
    if(gene_counter %% 1000 == 0)  {
      print("1000 genes processed")
    }
    gene_counter <<- gene_counter + 1
    return(prop_difs)
  })
  
  
  top_n_genes_per_cluster_df <- plyr::ldply(1:nrow(exp_genes_per_cluster_mat), function(.cluster_ind)  {
    .cluster <- exp_genes_per_cluster_mat[.cluster_ind,]
    top_n_inds <- order((.cluster), decreasing = TRUE)[1:n_top_genes]
    return(data.frame(cluster_id = cluster_ids[.cluster_ind], ENSEMBL = colnames(exp_genes_per_cluster_mat)[top_n_inds], prop_dif = .cluster[top_n_inds], stringsAsFactors = FALSE))
  })
  #colnames(top_n_genes_per_cluster_mat) <- paste("c", as.character(cluster_ids), sep = "")
  #top_n_genes_per_cluster_melted_df <- reshape2::melt(top_n_genes_per_cluster_mat)
  ensembl_to_gene_name_df <- select(EnsDb.Mmusculus.v79, keys=unique(top_n_genes_per_cluster_df$ENSEMBL), keytype = 'GENEID', columns = c("GENENAME", "GENEID")) %>% dplyr::rename(ENSEMBL =GENEID)
  return(dplyr::left_join(top_n_genes_per_cluster_df, ensembl_to_gene_name_df))
}

#creates a ggplot image of the tsne plot. The tsne_coordinates_mat matrix should be a sample_size x 2 matrix. The sample_name_to_cluster_id_df data.frame must have a 'sample_name' and a 'cluster_id' column.
get_clustered_scatterplot <- function(tsne_coordinates_mat, sample_name_to_cluster_id_df)  {
  temp_coord_df <- data.frame(tsne_coordinates_mat) %>% dplyr::rename(D1 = X1, D2 = X2) %>% dplyr::mutate(sample_name = sample_name_to_cluster_id_df$sample_name, sample_replicate = 'A', sample_date = sub("[A-Z]_.+", "", sample_name)) %>% dplyr::full_join(., sample_name_to_cluster_id_df)
  
  temp_coord_annos_for_fig_df <-  group_by(temp_coord_df, cluster_id) %>% dplyr::summarize(cluster_label_x = mean(D1, na.rm = TRUE), cluster_label_y = mean(D2, na.rm = TRUE))
  
  ggplot(temp_coord_df, aes(x = D1, y = D2, color = factor(cluster_id))) + geom_point() + geom_label(data = temp_coord_annos_for_fig_df, aes(x = cluster_label_x, y = cluster_label_y, label = as.character(cluster_id)))
}
  
get_treecuts <- function(expr_mat, deep_splits = 1:4)  {
  library(dynamicTreeCut)
  #On cluster:
  expr_mat_dist <- dist(expr_mat)
  expr_mat_hclust <- hclust(expr_mat_dist, method = "ward.D2")
  cut_list <- lapply(deep_splits, function(.split_int)  {
    return(cutreeDynamic(dendro = expr_mat_hclust, distM = as.matrix(expr_mat_dist), minClusterSize= 40, method = 'hybrid', deepSplit = .split_int, verbose = 4))
  })
  return(list(hclust = expr_mat_hclust, cut_list = cut_list, deep_split_ints = deep_splits))
}

get_dispersion <- function(log_normed_mat, nbins = 20, top= 1000)  {
  temp_mean_and_var <- apply(log_normed_mat, 2, function(.col)  {
    return(c(mean(.col), var(.col)))
  })
  temp_disp_mean_df <- data.frame(ensembl = colnames(temp_mean_and_var), mean_normed_gene_expr = temp_mean_and_var[1,], var_normed_gene_exp = temp_mean_and_var[2,], exp_bin = cut(temp_mean_and_var[1,], 20), stringsAsFactors = FALSE) %>% dplyr::mutate(gene_normed_exp_disp = var_normed_gene_exp/mean_normed_gene_expr)
  temp_disp_mean_df <- dplyr::left_join(temp_disp_mean_df, {dplyr::group_by(temp_disp_mean_df, exp_bin) %>% dplyr::summarize(mean_bin_expression = mean(mean_normed_gene_expr), mean_bin_dispersion = mean(gene_normed_exp_disp), sd_bin_dispersion = sd(gene_normed_exp_disp))}) %>% dplyr:: mutate(abs_normalized_bin_dispersion_deviation = abs((gene_normed_exp_disp - mean_bin_dispersion) / sd_bin_dispersion))
  temp_disp_mean_df <- dplyr::arrange(temp_disp_mean_df, desc(abs_normalized_bin_dispersion_deviation))
  topn_dispersed_mat <- log_normed_mat[,temp_disp_mean_df$ensembl[1:1000]]
  return(topn_dispersed_mat)
}

#the hm_df data.frame in get_ggplot_heatmap requires columns with names 'cluster_id', 'GENENAME', and 'mean_norm_expr'. The hm_df should be the mean of the normalized expression values of all the cells and genes included in the heatmap.
get_ggplot_heatmap <- function(col_order_char, row_order_char, hm_df = data.frame())  {
  temp_ggplot_df <- dplyr::group_by(hm_df, GENENAME) %>% dplyr::mutate(scaled_mean_norm_expr = (mean_norm_expr - mean(mean_norm_expr))/sd(mean_norm_expr)) %>% ungroup() %>% dplyr::mutate(GENENAME = factor(GENENAME, levels = col_order_char, ordered = TRUE), cluster_id = factor(cluster_id, levels =row_order_char, ordered = TRUE))
  print(head(temp_ggplot_df))
  ggplot(temp_ggplot_df, aes(x = GENENAME, y = cluster_id, fill = scaled_mean_norm_expr)) + geom_tile(color = 'white') + scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red') + theme(axis.text.x = element_text(angle = -90)) 
}

get_sample_name_to_cluster_id_df <- function(tree_splits_vec, expr_mat)  {
  data.frame(sample_name = rownames(expr_mat), cluster_id = tree_splits_vec, stringsAsFactors = FALSE)
}

#Returns the df of mean expression  of expressed genes that are enriched. The returned df 
get_mean_exp_by_cluster_with_genename_df <-function(enriched_mat, sample_name_to_cluster_id_df)  {
  enriched_mat_df <- data.frame(enriched_mat)
  enriched_mat_df$sample_name <- rownames(enriched_mat)
  #rm(giant_norpmt_gteq3500umi_lteq15000umi_noe14b_expressed_only_normed_mat)
  enriched_mat_melted_df <- reshape2::melt(enriched_mat_df, id.vars = 'sample_name', variable.name = 'GENEID', value.name = 'norm_expr') %>% dplyr::mutate(GENEID = as.character(GENEID)) %>% dplyr::left_join(., {sample_name_to_cluster_id_df %>% dplyr::select(sample_name, cluster_id)})# %>% dplyr::left_join(., enriched_mat_df)
  print(head(enriched_mat_melted_df))
  gene_id_df <- select(EnsDb.Mmusculus.v79, keytype = 'GENEID', keys = enriched_mat_melted_df$GENEID, columns = c('GENEID', 'GENENAME'))# %>% dplyr::rename(ENSEMBL = GENEID)
  enriched_mat_melted_df <- dplyr::left_join(enriched_mat_melted_df, gene_id_df) %>% dplyr::group_by(GENENAME, cluster_id) %>% dplyr::summarize(mean_norm_expr = mean(norm_expr)) %>% dplyr::ungroup()
  
  #Calculate the means per cluster
  #enriched_mat_scaled_expr_cluster_vs_genename_mat <-(reshape2::acast({dplyr::group_by(enriched_mat_melted_df, GENENAME, cluster_id) %>% dplyr::summarize(mean_norm_expr = mean(norm_expr))}, cluster_id ~ GENENAME))
  return(enriched_mat_melted_df)
}

get_annotation_specific_heatmap <- function(full_normed_mat = matrix(), ensembl_ids = character(), sample_name_to_cluster_id_df = data.frame(), row_order = NULL)  {
  temp_df <-
    get_mean_exp_by_cluster_with_genename_df(full_normed_mat[,unique(ensembl_ids)], sample_name_to_cluster_id_df = sample_name_to_cluster_id_df)
  temp_mat <- reshape2::acast(temp_df, cluster_id ~ GENENAME)
  temp_col_order <- colnames(temp_mat)[hclust(dist(scale(t(temp_mat))))$order]
  if(is.null(row_order))  {
    temp_cluster_order <- rownames(temp_mat)[hclust(dist(scale(temp_mat)))$order]    
  }  else  {
    temp_cluster_order <- row_order
  }

  get_ggplot_heatmap(temp_col_order, temp_cluster_order, temp_df)
}

remove_ribosomal_proteins_and_mitochondrial_genes_from_matrix <- function(cell_by_ensembl_mat = matrix(), species = 'hsapiens', rm_cell_cycle_genes = FALSE)  {
  library(biomaRt)
  if(!all(grepl("^ENS", colnames(cell_by_ensembl_mat))))  {
    stop("Column names should be ensembl ids and should start with 'ENS'")
  }
  if(species == 'hsapiens')  {
    ensembl <- biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl")
    useDataset("hsapiens_gene_ensembl", mart=ensembl)
    ensembl_gene_metadata_df <- biomaRt::select(ensembl, keys = colnames(cell_by_ensembl_mat), keytype = 'ensembl_gene_id', columns=c('ensembl_gene_id','chromosome_name', "hgnc_symbol", "description")) %>% dplyr::rename(ENSEMBL = ensembl_gene_id)
    ensembl_gene_metadata_MT_and_ribosomal_df <- dplyr::filter(ensembl_gene_metadata_df, grepl("^M?RP[SL]", hgnc_symbol) | (chromosome_name == 'MT'))
  }
  else if(species == 'mmusculus') {
    ensembl <- biomaRt::useMart("ensembl", dataset="mmusculus_gene_ensembl")
    useDataset("mmusculus_gene_ensembl", mart=ensembl)
    ensembl_gene_metadata_df <- biomaRt::select(ensembl, keys = colnames(cell_by_ensembl_mat), keytype = 'ensembl_gene_id', columns=c('ensembl_gene_id','chromosome_name', "mgi_symbol", "description")) %>% dplyr::rename(ENSEMBL = ensembl_gene_id)
    ensembl_gene_metadata_MT_and_ribosomal_df <- dplyr::filter(ensembl_gene_metadata_df, grepl("^M?RP[SL]", mgi_symbol, ignore.case = TRUE) | (chromosome_name == 'MT'))
  }
  else {
    warning("This is not a valid species option")
    stop()
  }
  #At this point we know we have either hsapiens or mmusculus species, so we go ahead and determine whether to remove the cell cycle genes
  cell_cycle_genes <- character()
  if(rm_cell_cycle_genes)  {
    cell_cycle_df <- read.table(dir(file.path(system.file(package = "rnaseqUtils"), 'extdata'), paste0(species, "_cell_cycle_genes.tsv"), full.names = TRUE), sep = '\t', header = TRUE, stringsAsFactors = FALSE)
    #current_gene_count <- ncol(ensembl_gene_metadata_MT_and_ribosomal_df)
    cell_cycle_genes <- cell_cycle_df$ENSEMBL
  }
  print(head(ensembl_gene_metadata_MT_and_ribosomal_df))
  ensembl_overlaps <- colnames(cell_by_ensembl_mat)[which(colnames(cell_by_ensembl_mat) %in% unique(c(ensembl_gene_metadata_MT_and_ribosomal_df$ENSEMBL, cell_cycle_genes)))]
  message(paste(length(ensembl_overlaps), 'of the', ncol(cell_by_ensembl_mat), 'ensembl ids are being removed'))
  return(cell_by_ensembl_mat[,!colnames(cell_by_ensembl_mat) %in% ensembl_overlaps])
}

get_df_from_named_char_vector <- function(char_vec = character(), col_names = character())  {
  temp_df <- as.data.frame(char_vec)
  names(temp_df) <- col_names[2]
  temp_df[col_names[1]] <- row.names(temp_df)
  temp_df <- temp_df[,2:1]
  row.names(temp_df) <- NULL
  return(temp_df)
}
