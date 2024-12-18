GRN_Output_For_SERGIO <- function(GRN,no_bins = 1,workdir = NULL){

  GRN1 <- as.data.frame(GRN)
  colnames(GRN1)[1:2] <- c("From","To")
  if(is.null(GRN$Weight)){
    GRN1$Weight <- -4+runif(nrow(GRN1),min=0,max=10)
  }
  if(is.null(GRN$ks)){
    GRN1$ks = 2.0
  }

  regulators <- setdiff(GRN1$From,GRN1$To) %>% unique()
  reg_info <- data.frame(reg = regulators)
  for(i in 1:no_bins){
    temp <- runif(length(regulators),min = 0,max = 5)
    reg_info <- cbind(reg_info,temp)
  }

  filename1 <- paste0(workdir,"Interaction_test.txt")
  filename2 <- paste0(workdir,"reg_test.txt")
  targets <- unique(GRN1$To)

  rm1 <- file.remove(filename1)
  rm2 <- file.remove(filename2)

  for(target in targets){

    us_regu <- GRN1$From[which(GRN1$To == target)]
    no_us <- length(us_regu)
    if(no_us > 0){
      weight_regu <- GRN1$Weight[which(GRN1$To == target)]
      ks_regu <- GRN1$ks[which(GRN1$To == target)]
      temp <- c(target,no_us,us_regu,weight_regu,ks_regu)
      temp <- t(as.matrix(temp))
      write.table(temp,file = filename1,append = T,quote = F,sep = ",",row.names = F,col.names = F)
    }
  }
  write.table(reg_info, file = filename2,quote = F,sep = ",",row.names = F,col.names = F)

  no_allgenes <- c(GRN1$From,GRN1$To) %>% unique() %>% length()
  returns <- paste0("Successfully generate a GRN containing ",no_allgenes," genes!\n")
  cat(returns)

  return(GRN1)
}

Generate_Random_MLnet <- function(no_rec,no_tf,no_tg,p1,p2){
  results <- c()
  maxnodes <- no_rec+no_tf+no_tg
  for(i in 1:no_rec){
    count <- 0
    for(j in 1:no_tf){
      rnd <- runif(1)
      if(rnd<p1){
        results <- rbind(results,c(i-1,j+no_rec-1))
        count <- count+1
      }
    }
    if(!count) results <- rbind(results,c(i-1,no_rec))
  }
  for(i in 1:no_tf){
    count <- 0
    for(j in 1:no_tg){
      rnd <- runif(1)
      if(rnd<p2){
        results <- rbind(results,c(i+no_rec-1,j+no_rec+no_tf-1))
        count <- count+1
      }
    }
    if(!count) results <- rbind(results,c(i+no_rec-1,no_rec+no_tf))
  }
  results <- as.data.frame(results)
  allnums <- as.vector(results) %>% unlist()
  newnums <- allnums
  names(newnums) <- NULL
  for(i in (maxnodes-1):0){
    if(!i %in% allnums){
      newnums[allnums>i] <- newnums[allnums>i]-1
    }
  }
  results_new <- data.frame(From = newnums[1:(length(newnums)/2)], To = newnums[(length(newnums)/2+1):length(newnums)])
  #alltfs <- setdiff(results_new$From,1:no_rec) %>% unique()

  #GRN <- list(network = results_new,
  #            Rec = 1:no_rec,
  #            TF = alltfs,
  #            TG = setdiff(results_new$To,alltfs) %>% unique())

  return(results_new)
}