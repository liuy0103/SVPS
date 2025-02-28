# Other methods and spectral clustering with SCORE and RSC ---------------------------------------------------------------------

library(RSpectra)
library(irlba)
library(ggplot2)
library(reshape)

delta <- 10^(-6)

# Calculates number of communities by CBIC and ICL
CBIC_ICL <- function(A, Kmax = 10, clus = 'SCORE', lambda = 1, tau = 0.25) 
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
        Pi_hat[label==i, i] <- 1
      }
    }
    
    n_k = rep(0, m)
    theta_hat = rep(0, n)
    for(k in 1:m){
      n_k[k] = max(1, sum(label==k))
      
      num1 = pmax(d_i[label==k], 1)
      num2 = max(1, sqrt(sum(A[label==k, label==k])))
      denom = max(1, sum(d_i[label==k]))
      theta_hat[label==k] = num1 * num2 / denom 
    }
    
    if(m == 1){
      B_hat = matrix(1, m, m)
    }else{
      B_hat = matrix(0, m, m)
      for(k in 1:(m-1)){
        for(l in (k+1):m){
          denom1 = max(1, sqrt(sum(A[label == k, label == k])))
          denom2 = max(1, sqrt(sum(A[label == l, label == l])))
          num = max(1, sum(A[label == k, label == l]))
          B_hat[k, l] = num / (denom1 * denom2)
          B_hat[l, k] = B_hat[k, l]
        }
      }
      diag(B_hat) = 1
    }
    
    # Calculating test statistics for CBIC, ICL
    M_hat = t(theta_hat * (Pi_hat %*% B_hat %*% t(Pi_hat)) ) * theta_hat
    
    P_hat = pmin(1-delta, (M_hat/n_nb+1)^(-1))
    loglike = sum(upperTriangle(n_nb*log(P_hat) + A*log(1-P_hat) + log(factorial(A+n_nb-1)), diag = T))
    cbic.result[m] <- loglike - (lambda*n*log(m) + m*(m+1)*log(n)/2)
    icl.result[m] <- loglike - sum(n_k*log(n/n_k)) - m*(m+2)*log(n)/2
  }
  return(list(cbic = which.max(cbic.result), 
              icl = which.max(icl.result)))
}



# Regularized Spectral Clustering (RSC)
reg.SSP <- function (A, K, tau = 1, lap = FALSE, nstart = 100, iter.max = 1000) 
{
  avg.d <- mean(colSums(A))
  A.tau <- A + tau * avg.d/nrow(A)
  if (!lap) {
    SVD <- irlba(A.tau, nu = K, nv = K)
    V <- SVD$v[, 1:K]
    V.norm <- apply(V, 1, function(x) sqrt(sum(x^2)))
    V.normalized <- diag(1/V.norm) %*% V
  }
  else {
    d.tau <- colSums(A.tau)
    L.tau <- diag(1/sqrt(d.tau)) %*% A.tau %*% diag(1/sqrt(d.tau))
    SVD <- irlba(L.tau, nu = K, nv = K)
    V <- SVD$v[, 1:K]
    V.norm <- apply(V, 1, function(x) sqrt(sum(x^2)))
    V.normalized <- diag(1/V.norm) %*% V
  }
  km <- kmeans(V.normalized, centers = K, nstart = nstart, iter.max = iter.max)
  return(list(cluster = km$cluster, 
              loss = km$tot.withinss))
}



# Spectral clustering with SCORE
SCORE <- function (A, K, threshold = NULL, nstart = 100, iter.max = 1000) 
{
  if (!is.matrix(A)) {
    A = as.matrix(A)
    if (any(A < 0)) {
      stop("Entries of adjacency matrix A should be nonegative!")
    }
  }
  if (any(A != t(A))) {
    stop("Aadjacency matrix A is not symmetric!")
  }
  eig.out = RSpectra::eigs(A, k = K)
  ev = eig.out$vectors[, 1:K]
  R = ev[, 2:K]/matrix(rep(ev[, 1], K - 1), ncol = K - 1, byrow = F)
  R = as.matrix(R, nrow = nrow(ev))
  if (is.null(threshold)) {
    threshold = log(nrow(A))
  }
  R[R > threshold] = threshold
  R[R < -threshold] = -threshold
  labels.hat = kmeans(R, K, nstart = nstart, iter.max = iter.max)$cluster
  return(list(R = R, 
              labels = labels.hat, 
              eig.values = eig.out$values[1:K], 
              eig.vectors = ev))
}





# Functions for SVPS and DCSBM generation -------------------------------------------------------

library(randnet)

upperTriangle <- function(x, diag = FALSE, byrow = FALSE){
  if (byrow) {
    t(x)[rev(upper.tri(x, diag = diag))]
  }
  else x[upper.tri(x, diag = diag)]
}
"upperTriangle<-" <- function(x, diag = FALSE, byrow = FALSE, value){
  x[upper.tri(x, diag = diag)] = value
  x
}

DCSBM <- function(K, theta_i, rho, r_val, n_vec){
  B <- matrix(1, K, K)
  diag(B) <- 1+r_val
  Pi <- matrix(0, n, K)
  membership <- rep(1:K, times = n_vec)
  membership <- sample(membership)
  Pi[cbind(1:n, membership)] <- 1
  M <- rho * diag(theta_i) %*% Pi %*% B %*% t(Pi) %*% diag(theta_i)
  return(M)
}

A_sim <- function(M, n_nb = 5){
  A_temp = matrix(0,n,n)
  upperTriangle(A_temp) = rnbinom(n*(n-1)/2, size = n_nb, prob = 1 - upperTriangle(M)/n_nb)
  A = A_temp + t(A_temp)
  diag(A) = rnbinom(n, size = n_nb, prob = 1 - diag(M)/n_nb)
  return(A)
}

theta_gen <- function(theta_min, theta_max){
  theta_i <- runif(n, min = theta_min, max = theta_max)
  theta_final <- rep(0, n)
  for(i in 1:n){
    theta_final[i] <- sample(c(theta_i[i], 0.5, 1.5), size = 1, replace = T, 
                             prob = c(0.8, 0.1, 0.1))
  }
  return(theta_final)
}

Stepwise <- function(A, Kmax = 10, clus = 'SCORE', thres, tau = 0.25){
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
        label <- reg.SSP(A, K = m, tau = tau, lap = TRUE)$cluster
      }
      Pi_hat <- matrix(0, n, m)
      for(i in 1:m){
        Pi_hat[label == i, i] <- 1
      }
    }
    
    theta_hat = rep(0, n)
    for(k in 1:m){
      num1 = pmax(d_i[label==k], 1)
      num2 = max(1, sqrt(sum(A[label==k, label==k])))
      denom = max(1, sum(d_i[label==k]))
      theta_hat[label==k] = num1 * num2 / denom 
    }
    
    if(m == 1){
      B_hat = matrix(1, m, m)
    }else{
      B_hat = matrix(0, m, m)
      for(k in 1:(m-1)){
        for(l in (k+1):m){
          denom1 = max(1, sqrt(sum(A[label==k, label==k])))
          denom2 = max(1, sqrt(sum(A[label==l, label==l])))
          num = max(1, sum(A[label==k, label==l]))
          B_hat[k, l] = num / (denom1 * denom2)
          B_hat[l, k] = B_hat[k, l]
        }
      }
      diag(B_hat) = 1
    }
    
    M_hat <- t(theta_hat * (Pi_hat %*% B_hat %*% t(Pi_hat)) ) * theta_hat
    # V_hat <- M_hat / (1 - M_hat / n_nb)
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





# Simulation 3 with plot ----------------------------------------------------------

library(parallel)

parallel_train <- function(idx){
  set.seed(1000 + idx)
  theta_i <- theta_gen(0.6, 1.4)
  M <- DCSBM(K, theta_i, rho, r_val, n_vec)
  A <- A_sim(M, n_nb = n_nb)
  
  K_hat_est <- Stepwise(A, Kmax = 10, clus = 'SCORE', thres = thres)
  cbic_icl_est <- CBIC_ICL(A, Kmax = 10, clus='SCORE')
  
  K_hat_est_2 <- Stepwise(A, Kmax = 10, clus = 'RSC', thres = thres, tau = tau)
  cbic_icl_est_2 <- CBIC_ICL(A, Kmax = 10, clus = 'RSC', tau = tau)
  
  return(list(
    K_hat = K_hat_est,
    bic = cbic_icl_est$cbic,
    icl = cbic_icl_est$icl,
    mod = cbic_icl_est$modularity,
    K_hat_2 = K_hat_est_2,
    bic_2 = cbic_icl_est_2$cbic,
    icl_2 = cbic_icl_est_2$icl,
    mod_2 = cbic_icl_est_2$modularity
  ))
}



set.seed(16)

# Threshold for sequential testing in SVPS
thres <- 2.1

# Number of iterations
ite_max <- 100

# Weighted DCSBM generation parameters
n_nb <- 5 # Number of trials in the binomial distribution

r_val <- 4 # Unbalancedness, which is r in Section 5.1
rho <- 0.04 # Sparsity, which is rho in Section 5.1

n_vec_whole <- c(50, 100, 150, 50, 100, 150) # Community size
tau <- 0.25 # Parameter for RSC


rho_K_hat <- c()
K_vec <- seq(2, 6)

for(K in K_vec){
  n_vec <- n_vec_whole[1:K]
  n <- sum(n_vec)
  
  K_hat_vec <- c() 
  bic_vec <- c()
  icl_vec <- c()
  mod_vec <- c()
  
  K_hat_vec_2 <- c() 
  bic_vec_2 <- c()
  icl_vec_2 <- c()
  mod_vec_2 <- c()
  
  ite_cnt <- seq(1, ite_max)
  results <- mclapply(ite_cnt, parallel_train, mc.cores = detectCores() - 1)
  for (i in 1:ite_max) {
    K_hat_vec[i] <- results[[i]]$K_hat
    bic_vec[i] <- results[[i]]$bic
    icl_vec[i] <- results[[i]]$icl
    K_hat_vec_2[i] <- results[[i]]$K_hat_2
    bic_vec_2[i] <- results[[i]]$bic_2
    icl_vec_2[i] <- results[[i]]$icl_2
  }
  
  print(K_hat_vec)
  rho_K_hat <- cbind(rho_K_hat, c(sum(K_hat_vec == K), sum(bic_vec == K), sum(icl_vec == K), 
                                  sum(K_hat_vec_2 == K), sum(bic_vec_2 == K), 
                                  sum(icl_vec_2 == K)) / ite_max)
}

all_methods <- c('SVPS_score','CBIC_score','ICL_score','SVPS_rsc','CBIC_rsc','ICL_rsc')

plot_data <- data.frame(t(rho_K_hat))
colnames(plot_data) <- all_methods

df <- data.frame(K = K_vec, plot_data)
df <- melt(df, id.vars = 'K')

ggplot(df, aes(x = K, y = value, linetype = variable, shape = variable, color = variable)) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 4) +
  scale_colour_manual(name = element_blank(),
                      labels = all_methods,
                      values = c("blue","red","green","black","orange","purple")) +
  scale_shape_manual(name = element_blank(),
                     labels = all_methods,
                     values = c(0, 1, 2, 3, 4, 5)) +
  scale_linetype_manual(name = element_blank(),
                        labels = all_methods,
                        values = c('solid','twodash','12345678','longdash','dashed','dotted')) +
  xlab('K') + ylab('Accuracy Rate') + 
  scale_y_continuous(n.breaks=6, limits=c(0, 1)) +
  theme_light() + 
  theme(plot.title = element_text(size = 25, hjust = 0.5),
        text = element_text(size = 25),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text = element_text(size = 20),
        legend.title = element_blank(),
        legend.text = element_text(size = 18),
        legend.justification = c(1, 0), legend.position = c(0.3, 0.05))

ggsave("plots_new/negBinom_r=4.png", width = 200, height = 140, units = 'mm')

save.image(file = "workspace_new/negBinom_r=4.RData")

