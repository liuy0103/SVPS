# Functions ---------------------------------------------------------------

library(RSpectra)
library(ScorePlus)
library(rmatio)
library(igraph)
library(randnet)
library(irlba)
# library(NetSurv)
library(igraphdata)

Stepwise <- function(A, Kmax = 15, clus='SCORE', thres){
  d_i <- rowSums(A)
  n <- dim(A)[1]
  
  for(m in 1:Kmax){
    if(m == 1){
      label <- rep(1, n)
      Pi_hat <- as.matrix(rep(1, n))
    }else{
      if(clus == 'SCORE'){
        label <- SCORE(A, m)$labels
      }else{
        label <- reg.SSP(A, K = m, tau = 0.25, lap = TRUE)$cluster
      }
      
      Pi_hat <- matrix(0, n, m)
      for(i in 1:m){
        Pi_hat[label == i, i] <- 1
      }
    }
    
    theta_hat = rep(0, n)
    for(k in 1:m){
      num1 <- pmax(d_i[label==k], 1)
      num2 <- max(1, sqrt(sum(A[label==k, label==k])))
      denom <- max(1, sum(d_i[label==k]))
      theta_hat[label==k] <- num1 * num2 / denom 
    }
    
    if(m == 1){
      B_hat <- matrix(1, m, m)
    }else{
      B_hat <- matrix(0, m, m)
      for(k in 1:(m-1)){
        for(l in (k+1):m){
          denom1 <- max(1, sqrt(sum(A[label==k, label==k])))
          denom2 <- max(1, sqrt(sum(A[label==l, label==l])))
          num <- max(1, sum(A[label==k, label==l]))
          B_hat[k, l] <- num / (denom1 * denom2)
          B_hat[l, k] <- B_hat[k, l]
        }
      }
      diag(B_hat) = 1
    }
    
    M_hat <- t(theta_hat * (Pi_hat %*% B_hat %*% t(Pi_hat)) ) * theta_hat
    V_hat <- M_hat 

    mat <- V_hat
    
    x <- 1 / sqrt(rowSums(mat))
    sinkhorn <- x
    
    while(max(abs(rowSums(mat) - 1)) > 1e-8){
      for(i in 1:20){
        mat <- t(x * mat) * x
        x <- 1 / sqrt(rowSums(mat))
        sinkhorn <- sinkhorn * x
      }
    }
    
    sinkhorn <- sqrt(n) * sinkhorn
    
    A_inv_hat <- t(sqrt(sinkhorn) * A) * sqrt(sinkhorn)
    A_inv_hat_eig <- eigs(A_inv_hat, k=Kmax+5, which = "LM")$values
    A_inv_hat_eig <- A_inv_hat_eig[order(abs(A_inv_hat_eig), decreasing = T)]
    if(abs(A_inv_hat_eig[m+1]) <= thres*sqrt(n)){
      return(m)
      break
    }
  }
}



CBIC_ICL <- function(A, Kmax = 15, clus='SCORE', lambda = 1) 
{
  d_i <- rowSums(A)
  n <- dim(A)[1]
  
  cbic.result <- rep(0, Kmax)
  icl.result <- rep(0, Kmax)
  for (m in 1:Kmax) {
    if(m == 1){
      label <- rep(1, n)
      Pi_hat <- as.matrix(rep(1, n))
    }else{
      if(clus == 'SCORE'){
        label <- SCORE(A, m)$labels
      }else{
        label <- reg.SSP(A, K = m, tau = 0.25, lap = TRUE)$cluster
      }
      
      Pi_hat <- matrix(0, n, m)
      for(i in 1:m){
        Pi_hat[label == i, i] <- 1
      }
    }
    
    n_k = rep(0, m)
    omega_hat = rep(0, n)
    for(a in 1:m){
      n_k[a] = max(1, sum(label == a))
      num = pmax(d_i[label==a], 1)
      denom = max(1, sum(d_i[label==a]))
      omega_hat[label==a] = n_k[a] * num / denom 
    }
    
    M = matrix(0, m, m)
    n_kl = matrix(0, m, m)
    B_hat = matrix(0, m, m)
    for(a in 1:m){
      for(b in 1:m){
        if(a == b){
          M[a, a] = max(1, sum(A[label==a, label==a]) / 2 )
          n_kl[a, a] = max(1, n_k[a] * (n_k[a] - 1) / 2)
          B_hat[a, a] = M[a, a] / n_kl[a, a]
        }else{
          M[a, b] = max(1, sum(A[label==a, label==b]))
          n_kl[a, b] = max(1, n_k[a] * n_k[b])
          B_hat[a, b] = M[a, b] / n_kl[a, b]
        }
      }
    }
    
    loglike = sum(d_i * log(omega_hat))
    for(a in 1:m){
      for(b in a:m){
        loglike = loglike + M[a,b] * log(B_hat[a,b]) - n_kl[a, b] * B_hat[a, b]
      }
    }
    cbic.result[m] <- loglike - (lambda*n*log(m) + m*(m+1)*log(n)/2)
    
    icl.result[m] <- loglike - sum(n_k*log(n/n_k)) - m*(m+2)*log(n)/2
  }
  return(list(cbic = which.max(cbic.result), icl = which.max(icl.result)))
}



# SVPS estimating the number of communities in the weighted adjacency with regularization ----------------------------------------------------------------

thres <- 2.05

real_data <- read.graph("dataset/Realdataset/lesmis.gml", format = "gml")
A <- as.matrix(as_adjacency_matrix(real_data, attr = 'value'))
n <- nrow(A)

# Regularized weighted adjacency; adjust tau to get different regularization levels
tau = 0.5
A = A + tau

# Estimated number of communities with SVPS using SCORE and RSC as clustering method, respectively
Stepwise(A, Kmax=10, clus='SCORE', thres=thres)
Stepwise(A, Kmax=10, clus='RSC', thres=thres)



# CBIC and ICL estimating the number of communities in the weighted adjacency ----------------------------------------------------------------

thres <- 2.05

real_data <- read.graph("dataset/Realdataset/lesmis.gml", format = "gml")
A <- as.matrix(as_adjacency_matrix(real_data, attr = 'value'))
n <- nrow(A)

# Estimated number of communities with CBIC and ICL using SCORE and RSC as clustering method, respectively
CBIC_ICL(A, Kmax=10, clus='SCORE')
CBIC_ICL(A, Kmax=10, clus='RSC')


















