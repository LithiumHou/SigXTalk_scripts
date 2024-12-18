# Get the expression matrix from a seurat object
Get_Exp_Clu <- function(SeuratObj, clusterID = NULL, assay = "RNA", datatype = "data", cutoff = 0.05) {

  if (is.null(clusterID)) {
    subobject <- SeuratObj
  } else {
    subobject <- subset(SeuratObj, idents = clusterID)
  }
  if("SCT" %in% names(SeuratObj@assays)){
    Count_clu <- subobject@assays$SCT$counts
  }else if("RNA" %in% names(SeuratObj@assays)){
    Count_clu <- subobject@assays$RNA$counts
  }else{
    Count_clu <- subobject@assays$Spatial$counts
  }
  
  allgenes <- rownames(Count_clu)
  ncells <- ncol(Count_clu)
  expressedgenes <- rownames(Count_clu)[which(rowSums(Count_clu > 0) > cutoff * ncells)]
  if(assay %in% names(subobject@assays)){
    if(datatype == "data"){
      Exp_clu <- subobject@assays[[assay]]$data
    }else if(datatype == "counts"){
      Exp_clu <- subobject@assays[[assay]]$counts
    }else{
      Exp_clu <- subobject@assays[[assay]]$scale.data
    }
  }else if ("RNA" %in% names(subobject@assays)){
    if(datatype == "data"){
      Exp_clu <- subobject@assays$RNA$data
    }else if(datatype == "counts"){
      Exp_clu <- subobject@assays$RNA$counts
    }else{
      Exp_clu <- subobject@assays$RNA$scale.data
    }
  }else if ("SCT" %in% names(subobject@assays)){
    if(datatype == "data"){
      Exp_clu <- subobject@assays$SCT$data
    }else if(datatype == "counts"){
      Exp_clu <- subobject@assays$SCT$counts
    }else{
      Exp_clu <- subobject@assays$SCT$scale.data
    }
  }else if ("Spatial" %in% names(subobject@assays)){
    if(datatype == "data"){
      Exp_clu <- subobject@assays$Spatial$data
    }else if(datatype == "counts"){
      Exp_clu <- subobject@assays$Spatial$counts
    }else{
      Exp_clu <- subobject@assays$Spatial$scale.data
    }
  }

  # filter low-expressed genes
  expressedgenes <- intersect(expressedgenes,rownames(Exp_clu))
  Exp_clu <- Exp_clu[expressedgenes, ]
  Exp_clu <- as(Exp_clu, "sparseMatrix")
  gc()

  return(Exp_clu)
}

# Convert from a cellchat object or network prob 3-dim tensor to an LR prob matrix
Extract_LR_Prob <- function(result, source_type = NULL, target_type = NULL, cellchat_use = T, pv_threshold = 0.05) {
  if(is.null(source_type) & is.null(target_type)){
    stop("Must input the source or the target!\n")
  }
  if (cellchat_use) {
    probs <- result@net$prob
    pvals <- result@net$pval
    if (is.null(source_type)) { # consider all sender types
      LR_vector <- c()
      for (type in dimnames(probs)[[1]]) {
        LR_df <- data.frame(interaction = names(probs[type, target_type, ]),Weight = probs[type, target_type, ],pval = pvals[type, target_type, ])
        LR_df <- LR_df %>%
          dplyr::filter(Weight > 0) %>%
          filter(pval < pv_threshold)
        if (!rlang::is_empty(LR_df)) LR_vector <- rbind(LR_vector, LR_df)
      }
    } else {
      if (!is.null(target_type)) {
        LR_df <- data.frame(interaction = names(probs[source_type, target_type, ]),Weight = probs[source_type, target_type, ],pval = pvals[source_type, target_type, ])
        LR_df <- LR_df %>%
          dplyr::filter(Weight > 0) %>%
          filter(pval < pv_threshold)
        LR_vector <- LR_df
      } else {
        LR_vector <- c()
        for (type in dimnames(probs)[[2]]) {
          LR_df <- data.frame(interaction = names(probs[source_type, type, ]),Weight = probs[source_type, type, ],pval = pvals[source_type, type, ])
          LR_df <- LR_df %>%
            dplyr::filter(Weight > 0) %>%
            filter(pval < pv_threshold)
          if (!rlang::is_empty(LR_df)) LR_vector <- rbind(LR_vector, LR_df)
        }
      }
    }
  } else {
    if (is.null(source_type)) { # consider all sender types
      LR_vector <- c()
      for (type in dimnames(probs)[[1]]) {
        LR_df <- data.frame(interaction = names(probs[type, target_type, ]),Weight = probs[type, target_type, ])
        if (!rlang::is_empty(LR_df)) LR_vector <- rbind(LR_vector, LR_df)
      }
    } else { #consider all receiver types
      LR_vector <- c()
      if(!is.null(target_type)){
        LR_vector <- data.frame(interaction = names(probs[source_type, target_type, ]),Weight = probs[source_type, target_type, ])
        LR_vector <- LR_vector %>%
            dplyr::filter(Weight > 0)
      }else{
        for (type in dimnames(probs)[[2]]) {
          LR_df <- data.frame(interaction = names(probs[source_type, type, ]),Weight = probs[source_type, type, ],pval = pvals[source_type, type, ])
          LR_df <- LR_df %>%
            dplyr::filter(Weight > 0) %>%
            filter(pval < pv_threshold)
          if (!rlang::is_empty(LR_df)) LR_vector <- rbind(LR_vector, LR_df)
        }
      }
      
    }
  }
  LR_vector$Weight <- as.numeric(LR_vector$Weight)
  Pairprob_final <- c()

  for (i in 1:nrow(LR_vector)) {
    tempname <- strsplit(LR_vector$interaction[i], "_") %>% unlist()
    if (length(tempname) == 3) {
      temp_frame1 <- data.frame("From" = tempname[1], "To" = tempname[2], "Weight" = LR_vector$Weight[i] / 2)
      temp_frame2 <- data.frame("From" = tempname[1], "To" = tempname[3], "Weight" = LR_vector$Weight[i] / 2)
      Pairprob_final <- rbind(Pairprob_final, temp_frame1, temp_frame2)
    } else {
      temp_frame1 <- data.frame("From" = tempname[1], "To" = tempname[2], "Weight" = LR_vector$Weight[i])
      Pairprob_final <- rbind(Pairprob_final, temp_frame1)
    }
  }
  Pairprob_final <- as.data.frame(Pairprob_final)
  Pairprob_final$Weight <- as.numeric(Pairprob_final$Weight)
  Pairprob_sum <- aggregate(Pairprob_final$Weight, list(Pairprob_final$From, Pairprob_final$To), sum)
  colnames(Pairprob_sum) <- c("From", "To", "Weight")
  return(Pairprob_sum)
}

# Do Pseudo-bulk aggregation
Aggregate_Data <- function(SeuratObj, target_type, k_neigh = 50, atacbinary = TRUE, max_overlap = 0.8, reduction.name = NULL,
                           size_factor_normalize = TRUE, verbose = TRUE) {
  cell_sample <- matrix(0, nrow = 1, ncol = k_neigh)
  allclusters <- SeuratObj@active.ident %>%
    as.character() %>%
    unique()

  if ("SCT" %in% names(SeuratObj@assays)) {
    rna_new_all <- matrix(0, nrow = nrow(SeuratObj@assays$SCT@counts), ncol = 1)
  } else {
    rna_new_all <- matrix(0, nrow = nrow(SeuratObj@assays$RNA@counts), ncol = 1)
  }
  if ("peaks" %in% names(SeuratObj@assays)) {
    atac_new <- matrix(0, nrow = nrow(SeuratObj@assays$peaks@counts), ncol = 1)
  } else if ("ATAC" %in% names(SeuratObj@assays)) {
    atac_new <- matrix(0, nrow = nrow(SeuratObj@assays$ATAC@counts), ncol = 1)
  } else {
    atac_new <- NULL
  }

  for (i in 1:length(allclusters)) {
    message(paste0("Aggregating cluster ", allclusters[i]))
    subobject <- subset(SeuratObj, idents = allclusters[i])
    sub_index <- 1:ncol(subobject)

    if ("wnn.umap" %in% names(subobject@reductions)) {
      cell_coord_i <- subobject@reductions$wnn.umap@cell.embeddings
    } else {
      cell_coord_i <- subobject@reductions$umap.rna@cell.embeddings
    }

    sub_aggregated_data <- Aggregation_Single(subobject, cell_coord_i, k_neigh, atacbinary, max_overlap)

    rna_new_all <- cbind(rna_new_all, sub_aggregated_data$rna)

    if (!is.null(atac_new)) {
      atac_new <- cbind(atac_new, sub_aggregated_data$atac)
    }

    sub_cell_sample <- sub_aggregated_data$cell_sample
    if (ncol(sub_cell_sample) < k_neigh) {
      sub_cell_sample_new <- as.matrix(sub_cell_sample)
      sub_cell_sample_new <- cbind(sub_cell_sample_new, matrix(0, nrow = 1, ncol = k_neigh - ncol(sub_cell_sample_new)))
    } else {
      sub_cell_sample_new <- apply(sub_cell_sample, 2, function(x) {
        sub_index[x] # for each column return original index
      })
      sub_cell_sample_new <- as.data.frame(sub_cell_sample_new)
      sub_cell_sample_new <- as.matrix(sub_cell_sample_new)
    }
    cell_sample <- rbind(cell_sample, sub_cell_sample_new)
  }

  new_data <- list()
  new_data$rna <- rna_new_all[, -1]
  if (!is.null(atac_new)) {
    new_data$atac <- atac_new[, -1]
  }
  new_data$cell_sample <- cell_sample[-1, ]
  return(new_data)
}


Aggregation_Single <- function(object, cell_coord, k_neigh = 50, atacbinary = TRUE, max_overlap = 0.8) {
  new_data <- list()
  if (nrow(cell_coord) > k_neigh) {
    # Create a k-nearest neighbors map
    nn_map <- as.data.frame(FNN::knn.index(cell_coord,
      k = (k_neigh - 1)
    ))
    row.names(nn_map) <- row.names(cell_coord)
    nn_map$agg_cell <- 1:nrow(nn_map)
    good_choices <- 1:nrow(nn_map)

    message("Sample cells randomly.\n")

    # Sample cells randomly
    choice <- sample(1:length(good_choices), size = 1, replace = FALSE)
    chosen <- good_choices[choice]
    good_choices <- good_choices[good_choices != good_choices[choice]]

    it <- 0
    ## Slow (contain calculating of overlapping between cell groups)
    while (length(good_choices) > 0 & it < nrow(cell_coord) / ((1 - max_overlap) * k_neigh)) {
      it <- it + 1
      choice <- sample(1:length(good_choices), size = 1, replace = FALSE)
      new_chosen <- c(chosen, good_choices[choice])
      good_choices <- good_choices[good_choices != good_choices[choice]]
      cell_sample <- nn_map[new_chosen, ]

      # calculate overlapping between cell groups
      combs <- data.frame(1:(nrow(cell_sample) - 1), nrow(cell_sample))
      shared <- apply(combs, 1, function(x) { # Slow
        (k_neigh * 2) - length(unique(as.vector(as.matrix(cell_sample[x, ]))))
      })
      if (max(shared) < max_overlap * k_neigh) {
        chosen <- new_chosen
      }
    }

    # aggregating both scRNA-seq and scATAC-seq counts of cells within one group
    if ("SCT" %in% names(object@assays)) {
      rna_old <- as.matrix(object@assays$SCT@data)
    } else {
      rna_old <- as.matrix(object@assays$RNA@counts)
    }
    rna_mask <- sapply(seq_len(nrow(cell_sample)), function(x) seq_len(ncol(rna_old)) %in% cell_sample[x, , drop = FALSE])
    rna_mask <- Matrix::Matrix(rna_mask)
    rna_new <- rna_old %*% rna_mask %>% as.matrix()

    if ("peaks" %in% names(object@assays)) {
      atac_old <- object@assays$peaks@counts
    } else if ("ATAC" %in% names(SeuratObj@assays)) {
      atac_old <- object@assays$ATAC@counts
    } else {
      atac_old <- NULL
    }
    # binarize

    if (!is.null(atac_old)) {
      if (atacbinary) {
        atac_old <- atac_old > 0
      }

      atac_mask <- sapply(seq_len(nrow(cell_sample)), function(x) seq_len(ncol(atac_old)) %in% cell_sample[x, , drop = FALSE])
      atac_mask <- Matrix::Matrix(atac_mask)
      atac_new <- atac_old %*% atac_mask %>% as.matrix()
      new_data$atac <- atac_new
    }
  } else {
    if ("SCT" %in% names(object@assays)) {
      rna_old <- as.matrix(object@assays$SCT@counts)
    } else {
      rna_old <- as.matrix(object@assays$RNA@counts)
    }
    rna_new <- rowSums(rna_old) %>% as.matrix()

    if ("peaks" %in% names(object@assays)) {
      atac_old <- object@assays$peaks@counts
    } else if ("ATAC" %in% names(SeuratObj@assays)) {
      atac_old <- object@assays$ATAC@counts
    }
    # binarize
    if (!is.null(atac_old)) {
      if (atacbinary) {
        atac_old <- atac_old > 0
      }
      atac_new <- rowSums(atac_old) %>% as.matrix()
      new_data$atac <- atac_new
    }

    cell_sample <- as.data.frame(t(matrix(seq(from = 1, to = nrow(cell_coord)))))
  }

  new_data$rna <- rna_new
  new_data$cell_sample <- cell_sample
  return(new_data)
}

# Convert interaction matrices to From-To lists
Convert_Mat_To_Relationship <- function(input_mat, rm_zero = F) {
  tryCatch(
    {
      output_df <- tidyfst::mat_df(input_mat)
      colnames(output_df) <- c("From", "To", "Weight")
      output_df$Weight <- as.numeric(output_df$Weight)
      return(output_df)
    },
    error = function(e) {
      message("Warning: package 'tidyfst' is not installed, using other methods...\n")
      output_df <- c()
      from_size <- nrow(input_mat)
      to_size <- ncol(input_mat)
      for (i in 1:from_size) {
        tempfrom <- rownames(input_mat)[i]
        tempdata <- cbind(rep(tempfrom, to_size), colnames(input_mat), input_mat[i, ])
        tempdata <- as.data.frame(tempdata)
        colnames(tempdata) <- c("From", "To", "Weight")
        if (rm_zero) {
          tempdata <- tempdata[tempdata$Weight != 0, ]
        }
        output_df <- rbind(output_df, tempdata)
      }
      output_df$Weight <- as.numeric(output_df$Weight)
      rownames(output_df) <- NULL
      return(output_df)
    }
  )
}


Filter_DB <- function(DB, allgenes) {
  require(rlang)
  DB_new <- DB[, 1:2]
  colnames(DB_new) <- c("From", "To")
  selfregu <- which(DB[, 1] == DB[, 2])
  if (!is_empty(selfregu)) {
    DB_new <- DB_new[-selfregu, ]
  }
  DB_new <- DB_new[DB_new$From %in% allgenes, ]
  DB_new <- DB_new[DB_new$To %in% allgenes, ]

  return(DB_new)
}

# Run the fisher test
Fisher_Test <- function(subset1,subset2,background){
  subset1 <- intersect(subset1, background)
  subset2 <- intersect(subset2, background)
  a=intersect(subset1,subset2)
  la <- length(a)
  b=setdiff(subset2,a)
  lb <- length(b)
  c=setdiff(subset1,a)
  lc <- length(c)
  abc <- union(a,b)
  abc <- union(abc,c)
  d <- setdiff(background,abc)
  ld <- length(d)
  matrix=matrix(c(la,lc,lb,ld),nrow=2)
  fisher.test(matrix,alternative="greater")$p.value
}

# Run the Chi-square test
Chisq_Test <- function(subset1,subset2,background){
  subset1 <- intersect(subset1, background)
  subset2 <- intersect(subset2, background)
  a=intersect(subset1,subset2)
  la <- length(a)
  b=setdiff(subset2,a)
  lb <- length(b)
  c=setdiff(subset1,a)
  lc <- length(c)
  abc <- union(a,b)
  abc <- union(abc,c)
  d <- setdiff(background,abc)
  ld <- length(d)
  matrix=matrix(c(la,lc,lb,ld),nrow=2)
  chisq.test(matrix)$p.value
}
# Calculate the overlap of genes
FindOverlap <- function(exp,TF,TG, type = "overlap"){
  
  # exp: genes as rows, cells as columns
  if(type == "overlap"){
    cell1 <- which(exp[TF,]>0)
    cell2 <- which(exp[TG,]>0)
    bothcell <- intersect(cell1,cell2)
    return(length(bothcell)/length(cell2))
  }else {
    return(cor(exp[TF,],exp[TG,],method = "spearman"))
  }
  
}


Calculate_Coexp <- function(DB, exp, method = 'spearman', Ncores){
  
  # method: overlap, spearman, pearson
  colnames(DB)[1:2] <- c("from","to") 
  DB_list <- split(DB[,1:2], f = DB$to)
  targets <- names(DB_list)

  if(Ncores == 1){
    coexps <- lapply(targets, function(x){
      regulators <- DB_list[[x]]$from
      if(method == "overlap"){
        cell1 <- which(exp[x,]>0)
        coexp <- sapply(regulators,function(y){
          cell2 <- which(exp[y,]>0)
          length(intersect(cell1,cell2))/length(cell2)
        })
      }else{
        coexp <- sapply(regulators,function(y){
          cor(exp[x,]+1e-5,exp[y,]+1e-5,method = method)
        })
      }
      data.frame(from = regulators, to = x, coexp = coexp %>% as.numeric() %>% abs())
      
    })
  }else{
    library(parallel)
    cl <- makeCluster(Ncores, type = "FORK")
    # doParallel::registerDoParallel(cl) 
    clusterEvalQ(cl, library(dplyr))
    clusterExport(cl, varlist = c("exp", "DB_list", "method"), envir=environment())
    coexps <- parLapply(cl, targets, function(x){
      regulators <- DB_list[[x]]$from
      if(method == "overlap"){
        cell1 <- which(exp[x,]>0)
        coexp <- sapply(regulators,function(y){
          cell2 <- which(exp[y,]>0)
          length(intersect(cell1,cell2))/length(cell2)
        })
      }else{
        coexp <- sapply(regulators,function(y){
          cor(exp[x,]+1e-5,exp[y,]+1e-5,method = method)
        })
      }
      data.frame(from = regulators, to = x, coexp = coexp %>% as.numeric())
      
    })      
    stopCluster(cl)
  }

  DB_coexp <- bind_rows(coexps)
  
  return(DB_coexp)
}

myPageRank <- function(gg, regulators, targets) {
  # Modified from scSeqComm
  # Convert dataframe to igraph object

  # Modified from scSeqComm
  association <- data.frame(regulator = character(0), target = character(0), PPR = numeric(0))
  nod <- sapply(igraph::V(gg)$name, function(z) strsplit(strsplit(z,":")[[1]][2],",")) #nodes of graph and their gene components
  if(is.na(nod[[1]])[1]){
      nod <- sapply(igraph::V(gg)$name, function(z) strsplit(strsplit(z,":")[[1]][1],",")) #nodes of graph and their gene components
  }
  r <- nodes_find(regulators, gg, nod)
  tf <- nodes_find(targets, gg, nod)

  if (length(r) == 0 | length(tf) == 0) return(association) # pathway does not contains at least one R and one TF

  # for each receptor
  for (j in 1:length(r)){
    # node reachable by all nodes associated to the current receptor
    tf_connected <- lapply(seq_along(r[[j]]), function(i) {
      return(is_connected_to(r[[j]][i], gg))
    })

    E <- rep(0, times = length(igraph::V(gg))) # preference vector for pagerank algorithm

    # computation of PERSONALIZED PAGERANK setting receptor node as seed node
    ppr <- single_receptor_ppr_wrapper(r[[j]], E, gg, delta = 0.85)
    colnames(ppr) <- igraph::V(gg)$name

    # for each transcription factor
    for (k in 1:length(tf)) {

      # a TF can be contained in more than one node
      for (t in 1:length(tf[[k]])) {
        ans <- c()

        # if multiple receptor node (a receptorr R1 contained in more than on one): the Persornalized PageRank algoritm is computed taking all receptor nodes as seed nodes
        if (length(r[[j]]) > 1) {
          # which receptor node is connected to the current tf?
          recp <- unlist(lapply(tf_connected, function(i) tf[[k]][t] %in% i))

          if (any(recp)) {
            # recomputation of pageranks setting the receptor nodes connected to tf as seed nodes
            ppr <- single_receptor_ppr_wrapper(r[[j]][recp], E, gg, delta = 0.85)
            # if(any(ppr<0)) print(gg$name)

            colnames(ppr) <- igraph::V(gg)$name

            # result
            ans <- c(ans, ppr[1, tf[[k]]][t])
          }
        } else {
          # the receptor is contained in only one node

          if (tf[[k]][t] %in% unlist(tf_connected)) # tf node is reachable by receptor node?
            {
              # result
              ans <- c(ans, ppr[1, tf[[k]]][t])
            }
        }
      } # for multiple tf

      # save results

      if (length(ans) != 0) {
        association <- rbind(association, data.frame(regulator = names(r)[j], target = names(tf)[k], PPR = mean(ans)))
      }
    } # for tf
  } # for r

  association <- association[order(-association$PPR), ]
  
  return(association)
}


# Do single time PPR algorithm
single_receptor_ppr_wrapper <- function(r,E,graph,delta){

  # set receptor as seed in preference vector E
  E[which(igraph::V(graph) %in% r,arr.ind = T)] <- 1

  #Personalized PageRank algorithm
  ppr <- igraph::page_rank(graph, algo = c("prpack"), vids = igraph::V(graph),directed = TRUE, damping = delta, personalized = E, weights = NULL) %>% .$vector
  ppr_matrix = matrix(unlist(ppr), ncol = length(E), byrow = TRUE)
  return(ppr_matrix)
}

nodes_find <- function(x,gg,nod){

  # Modified from scseqcomm

  ans <- list()
  for(j in 1:length(x)) {
    rec <- strsplit(x[[j]],",")[[1]]
    X <- c()

    #multi-subunit regulators
    for (r in rec)
    {
      #find nodes that contain current gene
      redok <- vapply(seq_along(nod),function(i) any(nod[[i]]==r), logical(1))

      if (any(redok))
      {
        X <- unique(c(X,igraph::V(gg)[igraph::V(gg)$name %in% names(nod[redok])]))  #ID node
      }
    }

    #save results
    if (length(X)>0)
    {
      t <- list(X)
      names(t) <- paste(x[j])
      ans <- c(ans,t)
    }


  }

  return(ans)
}

#function that finds all nodes reachable by the given node
is_connected_to <- function(r,gg)
{
  R_connection <- igraph::subcomponent(graph = gg, v = r, mode = "out") # NOTE: node "r" is reachable by itself, so it is included in the output
  R_connection <- sort(R_connection[R_connection!=r]) # remove node "r" and sort the vertex
  return(R_connection)
}


Multi_PPR <- function(igraphs, regulators, targets){
  
  R_TF_association <- data.frame(receptor = character(0), tf = character(0), tf_PPR = numeric(0))

  for (i in 1:length(igraphs))
  {
    gg <- igraphs[[i]] #graph
    #print(gg$name)
    nod <- sapply(igraph::V(gg)$name, function(x) strsplit(strsplit(x,":")[[1]][2],",")) #nodes of graph and their gene components
    #find nodes containing regulators and TRANSCRIPTION FACTORS
    r <- nodes_find(regulators,gg,nod)
    tf <- nodes_find(targets,gg,nod)

    if(length(r)==0 | length(tf)==0) next #pathway does not contains at least one R and one TF

    for (j in 1:length(r)) #for each receptor
    {

      #node reachable by all nodes associated to the current receptor
      tf_connected <- lapply(seq_along(r[[j]]), function(i) return(is_connected_to(r[[j]][i],gg)))


      E = rep(0,times=length(igraph::V(gg))) #preference vector for pagerank algorithm

      #computation of PERSONALIZED PAGERANK setting receptor node as seed node
      ppr <- single_receptor_ppr_wrapper(r[[j]],E,gg,delta=0.85)
      colnames(ppr) = igraph::V(gg)$name

      for (k in 1:length(tf)) #for each transcription factor
      {
        for (t in 1:length(tf[[k]])) #a TF can be contained in more than one node
        {
          ans <- c()

          #if multiple receptor node (a receptorr R1 contained in more than on one): the Persornalized PageRank algoritm is computed taking all receptor nodes as seed nodes
          if (length(r[[j]])>1)
          {
            #which receptor node is connected to the current tf?
            recp <- unlist(lapply(tf_connected, function(i) tf[[k]][t] %in% i))

            if(any(recp)){

              #recomputation of pageranks setting the receptor nodes connected to tf as seed nodes
              ppr <- single_receptor_ppr_wrapper(r[[j]][recp],E,gg,delta=0.85)
              #if(any(ppr<0)) print(gg$name)

              colnames(ppr) = igraph::V(gg)$name

              #result
              ans <- c(ans,ppr[1,tf[[k]]][t])
            }


          } else {
            #the receptor is contained in only one node

            if(tf[[k]][t] %in% unlist(tf_connected)) #tf node is reachable by receptor node?
            {
              #result
              ans <- c(ans,ppr[1,tf[[k]]][t])
            }
          }


        }#for multiple tf


        #save results

        if(length(ans) != 0)
        {
          R_TF_association <- rbind(R_TF_association,data.frame(receptor=names(r)[j],pathway=gg$name,tf=names(tf)[k],tf_PPR=mean(ans)))
        }

      }#for tf
    }#for r
  }#for igraphs

  return(R_TF_association)
}


Complete_groups <- function(df) {

  cnames <- colnames(df)
  # Ensure column names are treated as symbols
  x_col <- as.name(cnames[1])
  y_col <- as.name(cnames[2])
  group_col <- as.name(cnames[3])
  
  groups <- unique(df[,3])
  # Identify unique values of x in each group
  x_in_A <- unique(df[[x_col]][df[[group_col]] == groups[1]])
  x_in_B <- unique(df[[x_col]][df[[group_col]] == groups[2]])
  
  # Determine missing x values in each group
  missing_in_B <- setdiff(x_in_A, x_in_B)
  missing_in_A <- setdiff(x_in_B, x_in_A)
  
  # Create new rows with y = 0 for missing x values in each group
  new_rows_B <- data.frame(x = missing_in_B, y = 0, group = groups[2])
  new_rows_A <- data.frame(x = missing_in_A, y = 0, group = groups[1])
  
  # Rename columns to match the input dataframe's column names
  names(new_rows_B) <- names(df)
  names(new_rows_A) <- names(df)
  
  # Combine the original dataframe with the new rows
  df_complete <- rbind(df, new_rows_A, new_rows_B)
  
  # Sort the resulting dataframe for readability (optional)
  df_complete <- df_complete[order(df_complete[[group_col]], df_complete[[x_col]]), ]
  
  return(df_complete)
}
