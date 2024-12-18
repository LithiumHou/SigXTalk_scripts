# Calculate the normalized multi information between two data series
Calculate_NMI <- function(data1, data2, data3 = NULL, method = "emp"){
  require(infotheo)
  data1 <- discretize(data1)
  data2 <- discretize(data2)
  if(is.null(data3)){
    MI_value <- mutinformation(data1, data2, method)
  }else{
    data3 <- discretize(data3)
    MI_value <- condinformation(data1, data2, S = data3, method)-mutinformation(data1, data2, method)
  }
  MI_normalized <- MI_value/min(entropy(data1),entropy(data2))

  return(MI_normalized)
}

# Calculate the average value of a vector
Trimean <- function(x,nonzero = T) {
  if(nonzero){
    x[x == 0] <- NA
  }
  mean(stats::quantile(x, probs = c(0.25, 0.50, 0.50, 0.75), na.rm = T))
}

# Calculate the average expression of a gene
Calculate_Average_Expression <- function(targene, Expmat, nonzero = T)
{

  Ave_exp <- lapply(targene, function(gene) {
    if(!gene %in% rownames(Expmat)){
      0
    }else{
      Exp_gene <- Expmat[gene,]
      Trimean(Exp_gene, nonzero)
    }

  }) %>% as.data.frame()

  rownames(Ave_exp) <- "Exp"
  colnames(Ave_exp) <- targene
  return(Ave_exp)
}


# Calculate the distance matrix between cells
Calculate_Distance <- function(cellnames, cell_coords, p = 1){
  cellid <- cellnames[which(cellnames %in% rownames(cell_coords))]
  cell_coords_use <- cell_coords[cellid,]
  # Use the Manhattan distance
  distance_cells <- dist(cell_coords_use, method = "euclidean", diag = 1, upper = 1, p) %>% as.matrix()
  rownames(distance_cells) <- cellnames
  colnames(distance_cells) <- cellnames
  return(distance_cells)
}

# Calculate the copula sum
Nonlinear_Sum <- function(x, type = "Copula"){
  # input should be a numerical vector
  # sum type: copula, prob, sum
  if(!is.numeric(x)){
    stop("The input should be a numerical vector!\n")
  }
  x_len <- length(x)
  if(x_len < 2){
    if(x_len == 1){
      return(x)
    }else{
      stop("The input should at least contain 1 number!\n ")
    }
  }

  if(x_len == 2){
    if(type == "Copula"){
      return((x[1]+x[2]-2*x[1]*x[2])/(1-x[1]*x[2]))
    }
    else{
      return(1-(1-x[1])*(1-x[2]))
    }
  }else{
    if(type == "Copula"){
      Csum <- (x[1]+x[2]-2*x[1]*x[2])/(1-x[1]*x[2])
      for(i in 2:(x_len-1)){
        temp <- (Csum+x[i+1]-2*Csum*x[i+1])/(1-Csum*x[i+1])
        Csum <- temp
      }
    }else if(type == "prob"){
      Csum <- (1-(1-x[1])*(1-x[2]))
      for(i in 2:(x_len-1)){
        temp <- (1-(1-Csum)*(1-x[i+1]))
        Csum <- temp
      }
    }else if(type == "sum"){
      Csum <- sum(x)
    }

    return(Csum)
  }
}