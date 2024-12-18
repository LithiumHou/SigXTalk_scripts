### Compare the results with crosstalk TF
# Given two receptors, find the shared TF as ground truth
FindSharedTF_single <- function(df, receptor1, receptor2, allgenes) {
  colnames(df)[1:2] <- c("From", "To")
  TFset1 <- df %>%
    filter(From == receptor1) %>%
    select(To)
  TFset1 <- TFset1 %>%
    unlist() %>%
    unique() %>%
    intersect(allgenes)

  TFset2 <- df %>%
    filter(From == receptor2) %>%
    select(To)
  TFset2 <- TFset2 %>%
    unlist() %>%
    unique() %>%
    intersect(allgenes)
  
  if(length(TFset1)*length(TFset2) == 0){
    return(data.frame(TF = character(), value = numeric()))
  }
  tfc <- intersect(TFset1, TFset2)
  if(length(tfc) > 0){
    TF_pos <- data.frame(TF = tfc, value = 1)
  }else{
    TF_pos <- data.frame(TF = character(), value = numeric())
  }
  tf1 <- setdiff(TFset1, TFset2)
  tf2 <- c(tf1, setdiff(TFset2, TFset1))
  if(length(tf2) > 0){
    TF_neg <- data.frame(TF = tf2, value = 0)
  }else{
    TF_neg <- data.frame(TF = character(), value = numeric())
  }
  
  return(rbind(TF_pos, TF_neg) %>% as.data.frame())
}

# Given multiple receptors, find the shared TF as ground truth based on all expressed genes
FindSharedTF_truth <- function(truth_df, receptors, allgenes) {

    receptors <- intersect(receptors, allgenes)

    if(length(receptors) < 2){ # No crosstalk detected
      return(list())
    }

    # Find the 2-combinations of receptors
    rec_comb <- combn(receptors, 2)
    rec_comb <- t(rec_comb) %>% as.data.frame()
    rec_comb$pair <- paste0(rec_comb[, 1], "_", rec_comb[, 2])

    shared_tfs <- lapply(1:nrow(rec_comb), function(i) {
        x <- rec_comb[i, ] %>% unlist()
        tryCatch({FindSharedTF_single(truth_df, x[1], x[2], allgenes)},
        error = function(e) {message("Error in element ", i, ": ", e$message)
        return(NULL)
        })
    })
    names(shared_tfs) <- rec_comb$pair

    return(shared_tfs)
}

# Identify the SSC from predictions
FindSharedTF_pred <- function(result, receptors, allgenes) {

    colnames(result) <- c("from", "to", "value")

    receptors_result <- result$from %>%
        as.character() %>%
        unique()

    receptors <- intersect(receptors, receptors_result)

    if(length(receptors) < 2){ # No crosstalk detected
      return(list())
    }

    rec_comb <- combn(receptors, 2)

    rec_comb <- t(rec_comb) %>% as.data.frame()
    rec_comb$pair <- paste0(rec_comb[, 1], "_", rec_comb[, 2])

    shared_tfs <- lapply(1:nrow(rec_comb), function(i) {
        x <- rec_comb[i, ] %>% unlist()
        tryCatch({FindSharedTF_single(result, x[1], x[2], allgenes)},
        error = function(e) {message("Error in element ", i, ": ", e$message)
        return(NULL)
        })
    })

    names(shared_tfs) <- rec_comb$pair

    return(shared_tfs)
}

##### Evaluate the methods

AUC_single <- function(tf_prediction, tf_truth) {

  tf_merged <- merge(tf_prediction, tf_truth, by = "ID")
  true_labels <- tf_merged$value
  predicted_probs <- tf_merged$value_prod

  pr <- PRROC::pr.curve(true_labels, predicted_probs)
  prc <- pr$auc.integral
  rc <- PRROC::roc.curve(true_labels, predicted_probs)
  roc <- rc$auc

  return(c(roc, prc))
}

Evaluate_single <- function(tf_prediction, tf_truth, alltf, method = "Coverage") {

  real_tf <- tf_truth %>%
    filter(value == 1) %>%
    select(TF) %>%
    unlist() %>%
    unique()
  pred_tf <- tf_prediction %>%
    filter(value > 0) %>%
    select(TF) %>%
    unlist() %>% 
    unique()

  if(method == "Coverage"){
    commons <- intersect(real_tf, pred_tf)
    temp <- length(commons) / length(real_tf)
  }else if(method == "Fisher"){
    temp <- Fisher_Test(subset1 = real_tf, subset2 = pred_tf, background = alltf)
  }else if(method == "Chisq"){
    temp <- Chisq_Test(subset1 = real_tf, subset2 = pred_tf, background = alltf)
  }
  return(temp)

}


Evaluate_method <- function(result, truth_df, receptors, allgenes, threshold = 0.05) {

    alltf <- union(truth_df[, 2] %>% unique(), result[, 2] %>% unique())

    result <- filter(result, value > max(result$value)*threshold)
    predicted_SSC_list <- FindSharedTF_pred(result, receptors, allgenes)
    rec_pairs <- names(predicted_SSC_list)

    true_SSC_list <- FindSharedTF_truth(truth_df, receptors, allgenes)
    rec_pairs2 <- names(true_SSC_list)

    if(length(predicted_SSC_list)*length(true_SSC_list) == 0){
      metric_df <- data.frame(Coverage = 0, Fisher = 0, Chisqaure = 0)
      return(metric_df)
    }

    rec_pairs <- intersect(rec_pairs, rec_pairs2)
    metrics <- lapply(rec_pairs, function(x) {
        tf_prediction <- predicted_SSC_list[[x]]
        tf_truth <- true_SSC_list[[x]]
        if (is.null(tf_prediction) | is.null(tf_truth)) {
          temp <- c(Coverage = numeric(), Fisher = numeric(), Chisqaure = numeric())
        } else {
        
        cover <- Evaluate_single(tf_prediction, tf_truth, alltf, method = "Coverage")
        fis <- Evaluate_single(tf_prediction, tf_truth, alltf, method = "Fisher")
        chis <- Evaluate_single(tf_prediction, tf_truth, alltf, method = "Chisq")
        temp <- c(Coverage = cover, Fisher = -log10(fis), Chisqaure = -log10(chis))
        }
        temp
    })

    metric_df <- do.call(rbind, metrics)
    metric_df <- na.omit(metric_df)
    if(nrow(metric_df) == 0){
      return(data.frame(Coverage = 0, Fisher = 0, Chisqaure = 0))
    }
    return(metric_df %>% as.data.frame())
}
