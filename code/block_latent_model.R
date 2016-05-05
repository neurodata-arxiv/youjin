block_latent_network <- function(popn, ratio, rho){
  # popn : number of nodes
  n <- rep(0, length(ratio))
  for(j in 1:length(ratio)){
    n[j] <- as.integer(popn*ratio[j] / sum(ratio))
  }
  n[length(ratio)] <- popn - sum(n[1:(length(ratio) - 1)])
  sample <- c()
  for(j in 1:length(ratio)){
    Sigma <- matrix(c(1,rho[j], rho[j],1),2,2) # Covariance matrix
    sample <- rbind(sample, mvrnorm(n = n[j], mu = c(0,0), Sigma ) )# Generate popn's samples (X,Y)
  }
  
  
  # Create network graph with latent variable dependence structure
  prob <- matrix(0, popn, popn) # initiate a probability (of sharing a tie) matrix 
  adj <- matrix(0, popn, popn) # initiate an adjacency matrix
  
  for (i in 1:popn) {
    for (j in i:popn) {
      if (i == j) { 
        prob[i,j] = 0.0
      }
      else if (i != j) {
        dif <- abs(sample[i,1] - sample[j,1])
        
        if (dif < 0.01) dif = 0.01
        
        if (dif < 0.50) {
          val = dif
        }else {
          val = -20*dif
        }        
        prob[i,j] <- exp(val) / (1 + exp(val))
      }
      adj[i,j] <- rbinom(1, 1 , prob[i,j]) 
      adj[j,i] <- adj[i,j]
    } 
  }
  
  
  
  G <- ngraph(adj)
  G <- connect_ngraph(G)
  V(G)$outcome <- sample[,2]
  
  return(G)
}