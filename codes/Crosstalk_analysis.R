Count_Crosstalk <- function(CC_results, KeyGenes = NULL, verbose = T, datatype = "SSC") {
  if (datatype == "SSC") {
    if (is.null(KeyGenes)) {
      KeyGenes <- CC_results$SSC %>% unique()
    }
    KeySSC <- KeyGenes
    SSC_used <- intersect(KeySSC, CC_results$SSC) %>% unique()
    if (!length(SSC_used)) {
      stop("The input SSCs have no crosstalk!")
    }
    Count_Xmodule <- array(0, length(SSC_used))
    for (i in 1:length(SSC_used)) {
      SSC <- SSC_used[i]
      temp_df <- CC_results[CC_results$SSC == SSC, ]
      nXTalk <- nrow(temp_df)
      if (nXTalk > 1) {
        if (verbose) {
          message(paste0("There are ", nXTalk, "crosstalk pathways for SSC:", SSC, "!\n"))
        }
        Count_Xmodule[i] <- nXTalk
      } else {
        if (verbose) {
          message(paste0("There is no crosstalk for SSC:", SSC, "!\n"))
        }
        Count_Xmodule[i] <- 0
      }
    }
    names(Count_Xmodule) <- SSC_used
  } else if (datatype == "Target") {
    if (is.null(KeyGenes)) KeyGenes <- CC_results$Target %>% unique()
    KeyTG <- KeyGenes
    TG_used <- intersect(KeyTG, CC_results$Target) %>% unique()
    if (!length(TG_used)) {
      stop("The input TGs have no crosstalk!")
    }
    Count_Xmodule <- array(0, length(TG_used))
    for (i in 1:length(TG_used)) {
      tg <- TG_used[i]
      temp_df <- CC_results[CC_results$Target == tg, ]
      nXTalk <- nrow(temp_df)
      if (nXTalk > 1) {
        if (verbose) {
          message(paste0("There are", nXTalk, " crosstalk pathways for Target:", tg, "!\n"))
        }
        Count_Xmodule[i] <- nXTalk
      } else {
        if (verbose) {
          message(paste0("There is no crosstalk for Target:", tg, "!\n"))
        }
        Count_Xmodule[i] <- 0
      }
    }
    names(Count_Xmodule) <- TG_used
  }else{
    stop("Must specify the type of gene to be counted (SSC or Target)!")
  }

  return(Count_Xmodule)
}


# Calculate the Rec->TG Causality using the copula sum

Aggregate_Causality <- function(CC_results, sum_type, data_type = "Target") {
  if (data_type == "Target") {
    CC_df <- CC_results[, c("Receptor", "Target", "Weight")]
    CC_df$Weight <- 0.999 * CC_df$Weight / max(CC_df$Weight) # Normalize
    CC_aggregated <- aggregate(CC_df$Weight, list(CC_df$Rec, CC_df$Target), function(x) Nonlinear_Sum(x, type = sum_type))
    colnames(CC_aggregated) <- c("Receptor", "Target", "Weight")
    CC_sorted <- CC_aggregated[order(CC_aggregated$Rec, -CC_aggregated$Weight), ]
  } else if (data_type == "SSC") {
    CC_df <- CC_results[, c("Receptor", "SSC", "Weight")]
    CC_df$Weight <- 0.999 * CC_df$Weight / max(CC_df$Weight) # Normalize
    CC_aggregated <- aggregate(CC_df$Weight, list(CC_df$Receptor, CC_df$TF), function(x) Nonlinear_Sum(x, type = sum_type))
    colnames(CC_aggregated) <- c("Receptor", "TF", "Weight")
    CC_sorted <- CC_aggregated[order(CC_aggregated$Receptor, -CC_aggregated$Weight), ]
  }

  return(CC_sorted)
}


# Obtain a FIDELITY matrix for a fixed target gene.

Calculate_Fidelity_Matrix <- function(CC_results, TG_used) {

  colnames(CC_results) <- c("Receptor","SSC","Target","Weight")

  if (!TG_used %in% CC_results$Target) {
    stop("The input gene must be a target gene!")
    return(NULL)
  }
  if (length(TG_used) != 1) {
    stop("The input gene must be a target gene!")
    return(NULL)
  }
  CC_used <- CC_results %>% filter(Target == TG_used)
  CC_mat <- df_mat(CC_used, row = Receptor, col = SSC, value = Weight)
  CC_mat[is.na(CC_mat)] <- 0
  Fid_SSC <- sweep(CC_mat, 1, rowSums(CC_mat), FUN = "/") # for a given Receptor-Target, the fidelity for each SSC
  Fid_Receptor <- sweep(CC_mat, 2, colSums(CC_mat), FUN = "/") # for a given SSC-Target, the fidelity for each rec
  Fid_all <- CC_mat / sum(CC_mat) # fidelity for each pathway Receptor-SSC-Target

  Fid_results <- list(
    SSC = Fid_SSC,
    Receptor = Fid_Receptor,
    all = Fid_all
  )
  return(Fid_results)
}

Calculate_Specificity_Matrix <- function(CC_results, KeyRec) {

  if (!KeyRec %in% CC_results$Receptor) {
    message("The input gene must be a receptor!")
    return(NULL)
  }
  if (length(KeyRec) != 1) {
    stop("The input gene must be a receptor!")
    return(NULL)
  }
  CC_used <- CC_results %>% filter(Receptor == KeyRec)
  CC_used <- CC_used[!duplicated(CC_used[,1:3]),]
  CC_mat <- tidyfst::df_mat(CC_used, row = Target, col = SSC, value = Weight)
  CC_mat[is.na(CC_mat)] <- 0

  Spe_SSC <- sweep(CC_mat, 1, rowSums(CC_mat), FUN = "/") # for a given Receptor-Target, the specificity for each SSC
  Spe_Target <- sweep(CC_mat, 2, colSums(CC_mat), FUN = "/") # for a given Receptor-SSC, the specificity for each Target
  Spe_all <- CC_mat / sum(CC_mat) # specificity for each pathway Receptor-SSC-Target

  Spe_results <- list(
    SSC = Spe_SSC,
    Target = Spe_Target,
    all = Spe_all
  )
  return(Spe_results)
}

# Calculate the L-R-T pathway fidelity/specificity
Calculate_Pathway_Fidelity <- function(CC_results, KeyTG, KeySSC) {

  if (length(KeyRec) + length(KeySSC) != 2) stop("Must input one Receptor and one SSC! \n")
  CC_used <- CC_results %>%
    filter(Target == KeyTG) %>%
    filter(SSC == KeySSC)
  Fid_pathway <- CC_used$Weight / sum(CC_used$Weight)
  names(Fid_pathway) <- CC_used$Target
  return(Fid_pathway)
}

Calculate_Pathway_Specificity <- function(CC_results, KeyRec, KeySSC) {

  if (length(KeyRec) + length(KeySSC) != 2) stop("Must input one Receptor and one SSC! \n")
  CC_used <- CC_results %>%
    filter(Receptor == KeyRec) %>%
    filter(SSC == KeySSC)
  Spe_pathway <- CC_used$Weight / sum(CC_used$Weight)
  names(Spe_pathway) <- CC_used$Target
  return(Spe_pathway)
}


# Calculate the pairwise pathway fidelity/specificity
Calculate_Pairwise <- function(CC_pair_results, KeyGene = NULL, type = "Fid"){

  if(type == "Fid" || type == "Fidelity"){
    if(is.null(KeyGene)){
      KeyGene <- CC_pair_results$Target
    }else{
      KeyGene <- intersect(KeyGene, CC_pair_results$Target)
    }
    if(length(KeyGene) == 0){
      stop("The input gene is not a target gene!\n")
    }else{
      CC_used <- filter(CC_pair_results, Target %in% KeyGene)
      temp <- CC_used %>% group_by(Target) %>% mutate(Fidelity = Weight / sum(Weight))
      temp <- temp[,c("Receptor","Target","Fidelity")]
    }
  }
  if(type == "Spe" || type == "Specificity"){
    if(is.null(KeyGene)){
      KeyGene <- CC_pair_results$Receptor
    }else{
      KeyGene <- intersect(KeyGene, CC_pair_results$Receptor)
    }
    if(length(KeyGene) == 0){
      stop("The input gene is not a Receptor!\n")
    }else{
      CC_used <- filter(CC_pair_results, Receptor %in% KeyGene)
      temp <- CC_used %>% group_by(Receptor) %>% mutate(Specificity = Weight / sum(Weight))
      temp <- temp[,c("Receptor","Target","Specificity")]
    }
  }
  return(temp)
}


Compare_Multitype_CC <- function(CC_pair_list, type = "Fid") {

  alltypes <- names(CC_pair_list)
  Rec_list <- lapply(CC_pair_list, function(x) {
    x$Receptor %>% unique()
  })
  names(Rec_list) <- alltypes

  TG_list <- lapply(CC_pair_list, function(x) {
    x$Target %>% unique()
  })
  names(TG_list) <- alltypes

  common_Recs <- Reduce(intersect, Rec_list)
  common_TGs <- Reduce(intersect, TG_list)

  if (type == "Fid") {
    Fid_all <- lapply(alltypes, function(x) {
      CC_pair_res <- CC_pair_list[[x]]
      CC_pair_used <- CC_pair_res %>% filter(Target %in% common_TGs)
      CC_pair_mat <- df_mat(CC_pair_used, row = Receptor, col = Target, value = Weight)
      CC_pair_mat[is.na(CC_pair_mat)] <- 0
      Fid_pair_mat <- CC_pair_mat / colSums(CC_pair_mat)
      Fid_pair_df <- mat_df(Fid_pair_mat)
      colnames(Fid_pair_df) <- c("Receptor", "Target", "Fidelity")
      Fid_pair_df$celltype <- x
      Fid_pair_df
    })
    Fid_df <- do.call(rbind, Fid_all)

    return(Fid_df)
  }
  if (type == "Spe") {
    Spe_all <- lapply(alltypes, function(x) {
      CC_pair_res <- CC_pair_list[[x]]
      CC_pair_used <- CC_pair_res %>% filter(Target %in% common_TGs)
      CC_pair_mat <- df_mat(CC_pair_used, row = Receptor, col = Target, value = Weight)
      CC_pair_mat[is.na(CC_pair_mat)] <- 0
      Spe_pair_mat <- CC_pair_mat / rowSums(CC_pair_mat)
      Spe_pair_df <- mat_df(Spe_pair_mat)
      colnames(Spe_pair_df) <- c("Receptor", "Target", "Specificity")
      Spe_pair_df$celltype <- x
      Spe_pair_df
    })
    Spe_df <- do.call(rbind, Spe_all)

    return(Spe_df)
  }
}


XT_Subclustering <- function(SeuratObj, pair_results, KeyRecs, target_type, topk = 10) {
  require(flexclust)

  Tar_weight <- pair_results %>% filter(Receptor %in% KeyRecs)
  Tar_weight <- aggregate(Tar_weight$Weight, list(Tar_weight[,2]), sum)
  colnames(Tar_weight) <- c("Tar", "Weight")
  Tar_weight <- Tar_weight[order(Tar_weight$Weight, decreasing = T), ]
  topk <- min(topk, nrow(Tar_weight))
  Tar_weight <- Tar_weight[1:topk,]
  SeuratObj$celltype <- Idents(SeuratObj)
  subobject <- subset(SeuratObj, celltype == target_type)
  Exp_Tars <- subobject@assays$SCT$counts[Tar_weight$Tar, ]
  ws <- Tar_weight$Weight
  ws <- ws/sum(ws)

  cluster_results <- cclust(Exp_Tars %>% t(), k = 2, dist = "manhattan", method = "hardcl", weights = ws)
  cluster_center <- cluster_results@centers
  cluster_df <- data.frame(cell = colnames(Exp_Tars), clusterID = cluster_results@second, interacted = NA)
  norm1 <- sum(cluster_center[1, ])
  norm2 <- sum(cluster_center[2, ])

  if (norm1 < norm2) {
    cluster_df[which(cluster_df$clusterID == 1), 3] <- FALSE
    cluster_df[which(cluster_df$clusterID == 2), 3] <- TRUE
  } else {
    cluster_df[which(cluster_df$clusterID == 1), 3] <- TRUE
    cluster_df[which(cluster_df$clusterID == 2), 3] <- FALSE
  }

  return(cluster_df)
}
