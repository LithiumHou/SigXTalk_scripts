myPageRank <- function(receptors, tfactors, igraphs, allgenes) {

  recs <- intersect(receptors, allgenes)
  tfs <- intersect(tfactors, allgenes)  

  R_TF_association <- data.frame(receptor = character(0), pathway = character(0), tf = character(0), tf_PPR = numeric(0))

  for (i in 1:length(igraphs))
  {
    gg <- igraphs[[i]] #graph
    #print(gg$name)

    #find nodes containing RECEPTORS and TRANSCRIPTION FACTORS
    r <- nodes_find(recs,gg)
    tf <- nodes_find(tfs,gg)

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
    message(i)
  }#for igraphs

  return(R_TF_association)
}

nodes_find <- function (x, gg) {
    nod <- sapply(igraph::V(gg)$name, function(x) strsplit(strsplit(x, 
        ":")[[1]][2], ","))
    ans <- list()
    for (j in 1:length(x)) {
        rec <- strsplit(x[[j]], ",")[[1]]
        X <- c()
        for (r in rec) {
            redok <- vapply(seq_along(nod), function(i) any(nod[[i]] == 
                r), logical(1))
            if (any(redok)) {
                X <- unique(c(X, igraph::V(gg)[igraph::V(gg)$name %in% 
                  names(nod[redok])]))
            }
        }
        if (length(X) > 0) {
            t <- list(X)
            names(t) <- paste(x[j])
            ans <- c(ans, t)
        }
    }
    return(ans)
}

single_receptor_ppr_wrapper <- function(r,E,gg,delta){

    E[which(igraph::V(gg) %in% r, arr.ind = T)] <- 1
    ppr <- igraph::page_rank(gg, algo = c("prpack"), vids = igraph::V(gg), 
        directed = TRUE, damping = delta, personalized = E, weights = NULL) %>% 
        .$vector
    ppr_matrix = matrix(unlist(ppr), ncol = length(E), byrow = TRUE)
    return(ppr_matrix)
}

is_connected_to <- function(r,gg) {

  R_connection <- igraph::subcomponent(graph = gg, v = r, mode = "out")
  R_connection <- sort(R_connection[R_connection != r])
  return(R_connection)
  
}
