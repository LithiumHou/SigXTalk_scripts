##### The methods used for benchmark

############## Run in parallel with 10 cores
## CytoTalk ## Run in R 4.1.0 version
############## Filter the spots with the number of cells less than 3
CytoTalk_function <- function(SeuratObj, receiver = NULL, species = "human") {
  suppressMessages(library(CytoTalk))
  suppressMessages(library(Seurat))
  suppressMessages(library(tidyverse))
  suppressMessages(library(igraph))
  set.seed(123)

  # create conda environment
  if (F) {
    reticulate::conda_create(envname = "r_reticulate_CytoTalk", python_version = "3.8") # Create a new Conda environment to facilitate the Python module installation.
    reticulate::conda_install(envname = "r_reticulate_CytoTalk", "pybind11") # Install two necessary Python modules for correctly compiling and using the "pcst_fast" Python module.
    reticulate::conda_install(envname = "r_reticulate_CytoTalk", "numpy")
    reticulate::conda_install(envname = "r_reticulate_CytoTalk", "git+https://github.com/fraenkel-lab/pcst_fast.git", pip = TRUE) # To install the "pcst_fast" module.
  }
  if (species == "human") {
    pcg <- CytoTalk::pcg_human
    lrp <- CytoTalk::lrp_human
  } else {
    pcg <- CytoTalk::pcg_mouse
    lrp <- CytoTalk::lrp_mouse
  }

  counts <- SeuratObj@assays$SCT$data
  meta <- Idents(SeuratObj) %>% as.character()
  senders <- setdiff(unique(meta), receiver)

  lst.sc <- CytoTalk:::new_named_list(counts, meta)
  reticulate::use_condaenv("r_reticulate_CytoTalk")
  result <- list()

  Ncores <- parallel::detectCores() - 2
  for (sender in senders) {
    res <- CytoTalk::run_cytotalk(lst.sc, sender, receiver,
      cores = Ncores,
      cutoff_a = 0.05, cutoff_b = 0.05,
      pcg = pcg, lrp = lrp
    )
    cp <- paste(sender, receiver, sep = "_")
    result[[cp]] <- res
  }
  result_tmp <- lapply(names(result), function(cp){

    sender <-  gsub("_.*", "",cp)
    tmp <- list()
    res <- result[[cp]]

    if (!is.null(res$pathways)) {
      network <- res$pcst$final_network
      network <- network[!(network$node1_type == sender & network$node2_type == sender), ]
      network <- network[!(network$node1_type == receiver & network$node2_type == sender), ]
      network$node1 <- toupper(network$node1)
      network$node2 <- toupper(network$node2)
      network$node1 <- gsub("ORF", "orf", network$node1)
      network$node2 <- gsub("ORF", "orf", network$node2)
      LR <- network[which(network$node1_type == sender & network$node2_type == receiver), ]
      LR <- paste(LR$node1, LR$node2, sep = "_")
      Ligand <- network[which(network$node1_type == sender & network$node2_type == receiver), "node1"] %>%
        unique()
      Receptor <- network[which(network$node1_type == sender & network$node2_type == receiver), "node2"] %>%
        unique()
      Target <- c(
        network[which(network$node1_type == receiver & network$node2_type == receiver), "node1"],
        network[which(network$node1_type == receiver & network$node2_type == receiver), "node2"]
      ) %>%
        unique()
      Target <- setdiff(Target, c(Ligand, Receptor))
      tmp <- list(
        Edge = network,
        LR = LR,
        Ligand = Ligand,
        Receptor = Receptor,
        Target = Target
      )
    } else {
      tmp <- NA
    }
    tmp
  })

  names(result_tmp) <- senders
  result_tmp[which(is.na(result_tmp))] <- NULL

  edge_list <- lapply(result_tmp, function(res) {
    res <- res$Edge[, c(1, 2)]
    colnames(res) <- c("from", "to")
    res
  })
  edge_df <- do.call(rbind, edge_list)
  edge_df <- distinct(edge_df, from, to)
  intranetwork <- graph_from_edgelist(as.matrix(edge_df), directed = FALSE)

  Receptor <- lapply(result_tmp, function(res) {
    res$Receptor
  }) %>%
    unlist() %>%
    unique()
  Target <- lapply(result_tmp, function(res) {
    res$Target
  }) %>%
    unlist() %>%
    unique()

  distance <- distances(intranetwork,
    v = Receptor,
    to = Target
  )
  distance <- reshape2::melt(distance)
  distance <- distance[!is.infinite(distance$value), ]
  colnames(distance) <- c("regulon", "target", "value")

  return(distance)
}


########################
## NicheNet-databases ## Nichenet using pure databases
########################
NicheNet_database_function <- function(SeuratObj, geneset, receiver = NULL) {
  suppressMessages(library(nichenetr))
  suppressMessages(library(Seurat))
  suppressMessages(library(tidyverse))
  set.seed(123)

  # load database

  ligand_target_matrix <- readRDS("/home/jiawen/myMLnet/benchmark/ESICCC/nichenet/ligand_target_matrix.rds")
  lr_network <- readRDS("/home/jiawen/myMLnet/benchmark/ESICCC/nichenet/lr_network.rds")
  lr_network <- lr_network[, 1:2]
  all_types <- Idents(SeuratObj) %>%
    as.array() %>%
    unique()
  senders <- setdiff(all_types, receiver)

  # define the sender genes
  expressed_genes_sender <- lapply(senders, function(ct) {
    expressed_genes_sender_ct <- get_exp_genes(ct, SeuratObj, pct = 0.05)
    expressed_genes_sender_ct
  })
  expressed_genes_sender <- expressed_genes_sender %>%
    unlist() %>%
    unique()

  # define the receiver genes
  expressed_genes_receiver <- get_exp_genes(receiver, SeuratObj, pct = 0.05)
  geneset <- intersect(geneset, expressed_genes_receiver)
  # define the potential ligands
  ligands <- lr_network %>%
    pull(from) %>%
    unique()
  receptors <- lr_network %>%
    pull(to) %>%
    unique()

  expressed_ligands <- intersect(ligands, expressed_genes_sender)
  expressed_receptors <- intersect(receptors, expressed_genes_receiver)
  potential_ligands <- lr_network %>%
    filter(from %in% expressed_ligands & to %in% expressed_receptors) %>%
    pull(from) %>%
    unique()

  potential_targets <- intersect(geneset, rownames(ligand_target_matrix))
  lt_filtered <- ligand_target_matrix[potential_targets, potential_ligands]
  ltf_df <- tidyfst::mat_df(lt_filtered %>% t())
  colnames(ltf_df) <- c("regulon", "target", "value")
  result <- Convert_LT_to_RT(ltf_df, lr_network, potential_targets)

  return(result)
}

##########################
## NicheNet-computation ## Nichenet using databases and computation
##########################
NicheNet_computational_function <- function(SeuratObj, geneset, receiver = NULL) {

  suppressMessages(library(nichenetr))
  suppressMessages(library(Seurat))
  suppressMessages(library(tidyverse))
  set.seed(123)

  # load database

  ligand_target_matrix <- readRDS("/home/jiawen/myMLnet/benchmark/ESICCC/nichenet/ligand_target_matrix.rds")
  lr_network <- readRDS("/home/jiawen/myMLnet/benchmark/ESICCC/nichenet/lr_network.rds")
  weighted_networks <- readRDS("/home/jiawen/myMLnet/benchmark/ESICCC/nichenet/weighted_networks.rds")
  weighted_networks_lr <- weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from, to), by = c("from", "to"))

  all_types <- Idents(SeuratObj) %>%
    as.array() %>%
    unique()
  senders <- setdiff(all_types, receiver)

  # define the sender genes
  expressed_genes_sender <- lapply(senders[which(senders != receiver)], function(ct) {
    expressed_genes_sender_ct <- get_exp_genes(ct, SeuratObj, pct = 0.05)
    expressed_genes_sender_ct
  })
  expressed_genes_sender <- expressed_genes_sender %>%
    unlist() %>%
    unique()
  # define the potential ligands
  ligands <- lr_network %>%
    pull(from) %>%
    unique()
  receptors <- lr_network %>%
    pull(to) %>%
    unique()

  expressed_genes_receiver <- get_exp_genes(receiver, SeuratObj, pct = 0.05)
  expressed_ligands <- intersect(ligands, expressed_genes_sender)
  expressed_receptors <- intersect(receptors, expressed_genes_receiver)

  potential_ligands <- lr_network %>%
    filter(from %in% expressed_ligands & to %in% expressed_receptors) %>%
    pull(from) %>%
    unique()

  # define the receiver genes
  background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

  geneset <- geneset %>% .[. %in% rownames(ligand_target_matrix)]

  # define the potential ligands
  ligands <- lr_network %>%
    pull(from) %>%
    unique()
  receptors <- lr_network %>%
    pull(to) %>%
    unique()

  expressed_ligands <- intersect(ligands, expressed_genes_sender)
  expressed_receptors <- intersect(receptors, expressed_genes_receiver)

  potential_ligands <- lr_network %>%
    filter(from %in% expressed_ligands & to %in% expressed_receptors) %>%
    pull(from) %>%
    unique()

  # Perform NicheNet ligand activity analysis
  ligand_activities <- predict_ligand_activities(
    geneset = geneset,
    background_expressed_genes = background_expressed_genes,
    ligand_target_matrix = ligand_target_matrix,
    potential_ligands = potential_ligands
  )
  ligand_activities <- ligand_activities %>%
    arrange(-pearson) %>%
    mutate(rank = rank(desc(pearson)))

  # Get all the ligands for downstreams analysis
  best_upstream_ligands <- ligand_activities %>%
    arrange(-pearson) %>%
    pull(test_ligand) %>%
    unique()
  head(best_upstream_ligands)

  # Get the active ligand-target links
  active_ligand_target_links_df <- best_upstream_ligands %>%
    lapply(get_weighted_ligand_target_links,
      geneset = geneset,
      ligand_target_matrix = ligand_target_matrix,
      n = nrow(ligand_target_matrix)
    ) %>%
    bind_rows() %>%
    drop_na()
  colnames(active_ligand_target_links_df) <- c("regulon", "target", "value")
  # active_ligand_target_links_df$type <- "ligand"
  result <- Convert_LT_to_RT(active_ligand_target_links_df, lr_network[, 1:2], expressed_genes_receiver)
  return(result)
}

##########################
## SigXTalk-computation ## SigXTalk using databases and computation
##########################
SigXTalk_computational_function <- function(SeuratObj, geneset, receiver = NULL) {

  gc()
  library(Seurat)
  library(dplyr)
  library(CellChat)
  library(MASS)
  library(tibble)
  library(ggplot2)
  library(stringr)
  options(stringsAsFactors = FALSE)

  source("/home/jiawen/myMLnet/codes/MLNconstruct.R")
  source("/home/jiawen/myMLnet/codes/Numericals.R")
  source("/home/jiawen/myMLnet/codes/Utilities.R")
  source("/home/jiawen/myMLnet/codes/Crosstalk_analysis.R")

  DefaultAssay(SeuratObj) <- "SCT"
  SeuratObj$celltype <- Idents(SeuratObj)
  allgenes <- rownames(SeuratObj)
  cell_anno <- data.frame(cell = colnames(SeuratObj), cluster = SeuratObj$celltype)
  LRDB <- readRDS("/home/jiawen/myMLnet/pathways/Nichenet/LR_human.rds") %>% distinct(from, to, .keep_all = T) %>% Filter_DB(allgenes)
  RecTFDB <- readRDS("/home/jiawen/myMLnet/pathways/Nichenet/RTF_human.rds") %>% distinct(from, to, .keep_all = T) %>% Filter_DB(allgenes)
  # RecTFDB <- readRDS("/home/jiawen/myMLnet/pathways/KEGG/RTF_human.rds")
  TFTGDB <- readRDS("/home/jiawen/myMLnet/pathways/Nichenet/TFT_human.rds") %>% distinct(from, to, .keep_all = T) %>% Filter_DB(allgenes)
  LR_original <- Infer_CCI(SeuratObj, cell_anno, LRDB = LRDB, cellchat_output = T, use_spatial = F, db_use = "human")

  input_dir <- "/home/jiawen/myMLnet/pythoncodes/inputs"

  Prepare_Input_New(SeuratObj, receiver, geneset, LR_original, RecTFDB, TFTGDB, data_dir = input_dir,
      assay = "RNA", datatype = "scale.data", imputation = F,exp_threshold = 0.1,
      CCC_threshold = 0.05, Fisher_threshold = 1)
  gc()
  args.project <- "benchmark"
  args.target <- receiver
  conda_python <- "/home/jiawen/anaconda3/envs/SigXTalk/bin/python"

  system2(conda_python, args = c("/home/jiawen/myMLnet/pythoncodes/main_new.py", paste("--project",shQuote(args.project)), paste("--target_type",args.target)))

  output_dir <- '/home/jiawen/myMLnet/results/'
  filen <- paste0(output_dir, 'benchmark/pathways_', receiver,'.csv')
  RTFTG_results <- read.csv(filen, header = T)
  RTFTG_results <- RTFTG_results[,c(1,2,4)]
  RTFTG_results <- filter(RTFTG_results, pred_label > 0.75)
  RTFTG_results <- RTFTG_results[!duplicated(RTFTG_results[,1:2]),]

  return(RTFTG_results)
}

##### Utilities
# Convert ligand-target prediction to receptor-target prediction
Convert_LT_to_RT <- function(result, LRDB, allgenes) {
  colnames(result) <- c("Ligand", "TG", "value")
  colnames(LRDB) <- c("Ligand", "Receptor")

  LRDB <- filter(LRDB, Ligand %in% result$Ligand)
  LRDB <- filter(LRDB, Receptor %in% allgenes)

  merged_LR <- merge(result, LRDB, by = "Ligand")
  merged_LR <- merged_LR %>%
    group_by(Receptor) %>%
    mutate(value = value / n()) %>%
    ungroup()

  result_new <- merged_LR %>% select(Receptor, TG, value)

  return(result_new)
}

get_exp_genes <- function(ident, seurat_obj, pct = 0.1, assay_oi = NULL) {
  requireNamespace("Seurat")
  requireNamespace("dplyr")
  if (!"RNA" %in% names(seurat_obj@assays)) {
    if ("Spatial" %in% names(seurat_obj@assays)) {
      if (class(seurat_obj@assays$Spatial$data) != "matrix" &
        class(seurat_obj@assays$Spatial$data) != "dgCMatrix") {
        warning("Spatial Seurat object should contain a matrix of normalized expression data. Check 'seurat_obj@assays$Spatial$data' for default or 'seurat_obj@assays$SCT$data' for when the single-cell transform pipeline was applied")
      }
      if (sum(dim(seurat_obj@assays$Spatial$data)) == 0) {
        stop("Seurat object should contain normalized expression data (numeric matrix). Check 'seurat_obj@assays$Spatial$data'")
      }
    }
  } else {
    if (class(seurat_obj@assays$RNA$data) != "matrix" & class(seurat_obj@assays$RNA$data) !=
      "dgCMatrix") {
      warning("Seurat object should contain a matrix of normalized expression data. Check 'seurat_obj@assays$RNA$data' for default or 'seurat_obj@assays$integrated$data' for integrated data or seurat_obj@assays$SCT$data for when the single-cell transform pipeline was applied")
    }
    if ("integrated" %in% names(seurat_obj@assays)) {
      if (sum(dim(seurat_obj@assays$RNA$data)) == 0 & sum(dim(seurat_obj@assays$integrated$data)) ==
        0) {
        stop("Seurat object should contain normalized expression data (numeric matrix). Check 'seurat_obj@assays$RNA$data' for default or 'seurat_obj@assays$integrated$data' for integrated data")
      }
    } else if ("SCT" %in% names(seurat_obj@assays)) {
      if (sum(dim(seurat_obj@assays$RNA$data)) == 0 & sum(dim(seurat_obj@assays$SCT$data)) ==
        0) {
        stop("Seurat object should contain normalized expression data (numeric matrix). Check 'seurat_obj@assays$RNA$data' for default or 'seurat_obj@assays$SCT$data' for data corrected via SCT")
      }
    } else {
      if (sum(dim(seurat_obj@assays$RNA$data)) == 0) {
        stop("Seurat object should contain normalized expression data (numeric matrix). Check 'seurat_obj@assays$RNA$data'")
      }
    }
  }
  if (sum(ident %in% unique(Idents(seurat_obj))) != length(ident)) {
    stop("One or more provided cell clusters is not part of the 'Idents' of your Seurat object")
  }
  if (!is.null(assay_oi)) {
    if (!assay_oi %in% Seurat::Assays(seurat_obj)) {
      stop("assay_oi should be an assay of your Seurat object")
    }
  }
  cells_oi <- Idents(seurat_obj) %>%
    .[Idents(seurat_obj) %in%
      ident] %>%
    names()
  if (!is.null(assay_oi)) {
    cells_oi_in_matrix <- intersect(
      colnames(seurat_obj[[assay_oi]]$data),
      cells_oi
    )
    exprs_mat <- seurat_obj[[assay_oi]]$data %>% .[, cells_oi_in_matrix]
  } else {
    if ("integrated" %in% names(seurat_obj@assays)) {
      warning("Seurat object is result from the Seurat integration workflow. The expressed genes are now defined based on the integrated slot. You can change this via the assay_oi parameter of the get_exp_genes() functions. Recommended assays: RNA or SCT")
      cells_oi_in_matrix <- intersect(
        colnames(seurat_obj@assays$integrated$data),
        cells_oi
      )
      if (length(cells_oi_in_matrix) != length(cells_oi)) {
        stop("Not all cells of interest are in your expression matrix (seurat_obj@assays$integrated$data). Please check that the expression matrix contains cells in columns and genes in rows.")
      }
      exprs_mat <- seurat_obj@assays$integrated$data %>%
        .[, cells_oi_in_matrix]
    } else if ("SCT" %in% names(seurat_obj@assays) & !"Spatial" %in%
      names(seurat_obj@assays)) {
      # warning("Seurat object is result from the Seurat single-cell transform workflow. The expressed genes are defined based on the SCT slot. You can change this via the assay_oi parameter of the get_exp_genes() functions. Recommended assays: RNA or SCT")
      cells_oi_in_matrix <- intersect(
        colnames(seurat_obj@assays$SCT$data),
        cells_oi
      )
      if (length(cells_oi_in_matrix) != length(cells_oi)) {
        stop("Not all cells of interest are in your expression matrix (seurat_obj@assays$SCT$data). Please check that the expression matrix contains cells in columns and genes in rows.")
      }
      exprs_mat <- seurat_obj@assays$SCT$data %>% .[, cells_oi_in_matrix]
    } else if ("Spatial" %in% names(seurat_obj@assays) & !"SCT" %in%
      names(seurat_obj@assays)) {
      warning("Seurat object is result from the Seurat spatial object. The expressed genes are defined based on the Spatial slot. If the spatial data is spot-based (mixture of cells) and not single-cell resolution, we recommend against directly using nichenetr on spot-based data (because you want to look at cell-cell interactions, and not at spot-spot interactions! ;-) )")
      cells_oi_in_matrix <- intersect(
        colnames(seurat_obj@assays$Spatial$data),
        cells_oi
      )
      if (length(cells_oi_in_matrix) != length(cells_oi)) {
        stop("Not all cells of interest are in your expression matrix (seurat_obj@assays$Spatial$data). Please check that the expression matrix contains cells in columns and genes in rows.")
      }
      exprs_mat <- seurat_obj@assays$Spatial$data %>% .[
        ,
        cells_oi_in_matrix
      ]
    } else if ("Spatial" %in% names(seurat_obj@assays) & "SCT" %in%
      names(seurat_obj@assays)) {
      warning("Seurat object is result from the Seurat spatial object, followed by the SCT workflow. If the spatial data is spot-based (mixture of cells) and not single-cell resolution, we recommend against directly using nichenetr on spot-based data (because you want to look at cell-cell interactions, and not at spot-spot interactions! The expressed genes are defined based on the SCT slot, but this can be changed via the assay_oi parameter.")
      cells_oi_in_matrix <- intersect(
        colnames(seurat_obj@assays$SCT$data),
        cells_oi
      )
      if (length(cells_oi_in_matrix) != length(cells_oi)) {
        stop("Not all cells of interest are in your expression matrix (seurat_obj@assays$Spatial$data). Please check that the expression matrix contains cells in columns and genes in rows.")
      }
      exprs_mat <- seurat_obj@assays$SCT$data %>% .[, cells_oi_in_matrix]
    } else {
      if (sum(cells_oi %in% colnames(seurat_obj@assays$RNA$data)) ==
        0) {
        stop("None of the cells are in colnames of 'seurat_obj@assays$RNA$data'. The expression matrix should contain cells in columns and genes in rows.")
      }
      cells_oi_in_matrix <- intersect(
        colnames(seurat_obj@assays$RNA$data),
        cells_oi
      )
      if (length(cells_oi_in_matrix) != length(cells_oi)) {
        stop("Not all cells of interest are in your expression matrix (seurat_obj@assays$RNA$data). Please check that the expression matrix contains cells in columns and genes in rows.")
      }
      exprs_mat <- seurat_obj@assays$RNA$data %>% .[, cells_oi_in_matrix]
    }
  }
  n_cells_oi_in_matrix <- length(cells_oi_in_matrix)
  if (n_cells_oi_in_matrix < 5000) {
    genes <- exprs_mat %>%
      apply(1, function(x) {
        sum(x > 0) / n_cells_oi_in_matrix
      }) %>%
      .[. >= pct] %>%
      names()
  } else {
    splits <- split(1:nrow(exprs_mat), ceiling(seq_along(1:nrow(exprs_mat)) / 100))
    genes <- splits %>%
      lapply(function(genes_indices, exprs,
                      pct, n_cells_oi_in_matrix) {
        begin_i <- genes_indices[1]
        end_i <- genes_indices[length(genes_indices)]
        exprs <- exprs[begin_i:end_i, , drop = FALSE]
        genes <- exprs %>%
          apply(1, function(x) {
            sum(x > 0) / n_cells_oi_in_matrix
          }) %>%
          .[. >= pct] %>%
          names()
      }, exprs_mat, pct, n_cells_oi_in_matrix) %>%
      unlist() %>%
      unname()
  }
  return(genes)
}


