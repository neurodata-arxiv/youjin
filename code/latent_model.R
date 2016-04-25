library(igraph)
library(MASS)
library(mvrtn)
library(Matrix)

simdata <- function(popn, rho){
 # popn : number of nodes
 # rho : correlation coefficient between the outcomes and latent variables
  #popn <- 500
  #rho = 0.1
  Sigma <- matrix(c(1,rho, rho,1),2,2) # Covariance matrix
  sample <- mvrnorm(n = popn, mu = c(0,0), Sigma ) # Generate popn's samples (X,Y)
  
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
        
        if (dif < 0.10) {
          val = 2 / (dif)
        }else if(dif < 0.50){
          val = dif
        }else if(dif < 1.00){
          val = - dif
        }else if(dif < 1.50){
          val = - 2*dif
        }else{
          val = -10*dif 
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