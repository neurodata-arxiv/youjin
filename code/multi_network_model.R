multi_network_model <- function(popn, nnet, rho){
  # popn : number of nodes
  # nnet : number of networks
  sample <- list()
  G <- list()
  
  for(k in 1:nnet){
    Sigma <- matrix(c(1,rho, rho,1),2,2) # Covariance matrix
    sample[[k]] <- mvrnorm(n = popn, mu = c(0,0), Sigma ) # Generate popn's samples (X,Y)
  }
  
  for(k in 1:nnet){
    # Create network graph with latent variable dependence structure
    prob <- matrix(0, popn, popn) # initiate a probability (of sharing a tie) matrix 
    adj <- matrix(0, popn, popn) # initiate an adjacency matrix
    
    for (i in 1:popn) {
      for (j in i:popn) {
        if (i == j) { 
          prob[i,j] = 0.0
        }
        else if (i != j) {
          dif <- abs(sample[[k]][i,1] - sample[[k]][j,1])
          
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
    
    G[[k]] <- ngraph(adj)
    G[[k]] <- connect_ngraph(G[[k]])
    V(G[[k]])$outcome <- sample[[k]][,2]
  }
  return(G)
}

