# functions used for GGM demo
library("R.utils")
library("glasso")
library("flare")

# param: 
#      S: sample cov mat based on data
#      invSigma: estimated precision mat by glasso; path mat by clime or tiger
#      N: the number of observations
#      use_tol: take tolerance into consideration to decide the adjaceny mat based on invSigma
# return: the optimal glasso object satisfying the same zero pattern with invSigma
net_est <- function(S, invSigma, N, use_tol=TRUE, model="glasso")
{
    if(model == "glasso"){
        # find adj based on invSigma
        if(use_tol){
            adj_mat <- ifelse(abs(invSigma) > tol, 1, 0)
        }else{
            adj_mat <- ifelse(invSigma != 0, 1, 0)
        }
    }else{
        adj_mat = invSigma # path in clime and tiger
    }
  diag(adj_mat) <- 1
  # Compute zeroes:
  zeroes <- which(adj_mat==0, arr.ind=TRUE)
  # fit
  net <- NULL
  # fit
  if(nrow(zeroes)>0){
    net <- glasso(S, 0, zero = zeroes, nobs = N, thr = tol/100, approx = FALSE, trace = 0, penalize.diagonal=FALSE)
    # net2 <- glasso(S, 1e12*(1-adj_mat),nobs=N, penalize.diagonal=FALSE)
  }else{
    net <- glasso(S,0, nobs = N, thr = tol/100, approx = FALSE, trace = 0, penalize.diagonal=FALSE)
  }
  return(net)
}





# criteria
# alll criteria use loglik calculated by glasso

# param: 
#      S: sample cov mat based on data
#      K: estimated precision mat by glasso
#      L: loglik
#      E: the number of edges based on K
#      use_tol: take tolerance into consideration to decide the adjaceny mat based on invSigma
# return: AIC score
AIC_cal <- function(S,K,L,N,E,use_tol=TRUE,countDiagonal=FALSE)
{
  if (missing(E)){
    if(use_tol){
      E <- sum(abs(K[lower.tri(K,diag=countDiagonal)]) > tol)
    }else{
      E <- sum(K[lower.tri(K,diag=countDiagonal)] != 0)
    }
  }
  -2 * L + 2 * E
}

# param: 
#      S: sample cov mat based on data
#      K: estimated precision mat by glasso
#      L: loglik
#      N: sample size
#      gamma: 0 for BIC; >0 for EBIC
#      E: the number of edges based on K
#      use_tol: take tolerance into consideration to decide the adjaceny mat based on invSigma
# return: BIC/EBIC score
EBIC_cal <- function(S,K,L,N,gamma = 1,E,use_tol=TRUE,countDiagonal=FALSE)
{
  if (missing(E)){
    if(use_tol){
      E <- sum(abs(K[lower.tri(K,diag=countDiagonal)]) > tol)
    }else{
      E <- sum(K[lower.tri(K,diag=countDiagonal)] != 0)
    }
  }
  p <- nrow(K)
  
  # return EBIC:
  -2 * L + E * log(N) + 4 * E * gamma * log(p)  
}

CV_partition <- function(N, nfold){
  ntest = floor(N/nfold)
  ntrain = N - ntest
  # ith col is data used for train/test for ith fold
  train_mat = matrix(NA, nrow = ntrain, ncol = nfold)
  test_mat = matrix(NA, nrow = ntest, ncol = nfold)
  # fill in data
  N_seq = c(1:N)
  data_idx = sample(N)
  for(i in 1:nfold){
    test_idx = ((i-1)*ntest+1):(i*ntest)
    train_idx = N_seq[!N_seq %in% test_idx]
    test_mat[,i] = data_idx[test_idx]
    train_mat[,i] = data_idx[train_idx]
  }
  cv_data_list = list(train_mat=train_mat, test_mat=test_mat)
  return(cv_data_list)
}

CV <- function(N, nfold, data, fixed_lam, tol=1e-3){
    Nlam = length(fixed_lam)
    cv_data_list = CV_partition(N, nfold)
    cv_loss <- sapply(1:nfold, function(i){
        data.train = data[cv_data_list$train_mat[,i],]
        g_path = glasso::glassopath(s = cov(data.train), rholist = fixed_lam, thr = tol/100, approx = FALSE, trace = 0, penalize.diagonal = FALSE)
        data.test = data[cv_data_list$test_mat[,i],]
        # max loglik
        loss  <- sapply(1:Nlam, function(j){
            S = cov(data.test)*(1-1/nrow(data.test))
            K = g_path$wi[,,j]
            d = det(K)
            (N/2)*(log(d) - sum(diag(S%*%K)))
        })
    # cv_loss2[i,] = loss
    return(loss)
    })
    cv_loss.mean = apply(cv_loss, 1, mean)
    cv_loss.sd = apply(cv_loss, 1, sd)
    return(list(cv.mean=cv_loss.mean, cv.sd=cv_loss.sd))
}

# param: 
#      idx: between 1 to Nlam
#      adj_mat_True: the true adjacency matrix based on theta_0
#      use_tol: take tolerance into consideration to decide the adjaceny mat based on invSigma
# return: the number of mismatch between true adj_mat and adj_mat based on the optimal lasso with lam = fixed_lam[idx]
hd_loss <- function(idx, est, adj_mat_True, use_tol=TRUE, tol=1e-3){
    invSigma <- est[,idx]$wi
    if(use_tol){
        adj_mat <- ifelse(abs(invSigma) > tol, 1, 0)
    }else{
        adj_mat <- ifelse(invSigma !=0, 1, 0)
    }
    diag(adj_mat) = 1
    sum(abs(adj_mat_True - adj_mat))/2
}

# param: 
#      idx: between 1 to Nlam
#      use_tol: take tolerance into consideration to decide the adjaceny mat based on invSigma
# return: False Positive Rate(FDR)
FDR <- function(idx, est, adj_mat_True, use_tol=TRUE, tol=1e-3,countDiagonal=FALSE, file="output.txt"){
    K <- est[,idx]$wi
    if(use_tol){
        adj_mat <- ifelse(abs(K) > tol, 1, 0)
    }else{
        adj_mat <- ifelse(K !=0, 1, 0)
    }
    diag(adj_mat) <- 1
    true_neg <- which(adj_mat_True == 0, arr.ind = TRUE)
    true_pos <- which(adj_mat_True == 1, arr.ind = TRUE)
    # should be 0, but marked 1 (FP)
    FP = sum(adj_mat[true_neg])
    # should be 1, and marked 1(TP)
    TP = sum(adj_mat[true_pos])-nrow(adj_mat_True)
    if(FP+TP <= 0){
        cat(paste0("\nAt lambda idx = ", idx, ",  estimated adj mat is all zeros; invalid FDR.\n"), file=file,append=TRUE)
        return(0)
    }
    return(FP/(FP+TP))
}

# param: 
#      idx: between 1 to Nlam
#      use_tol: take tolerance into consideration to decide the adjaceny mat based on invSigma
# return: Positive Selection Rate (the number of edges in estimated adj_mat/total number of edges of true adj_mat)
TPR <- function(idx, est, adj_mat_True, use_tol=TRUE, tol=1e-3,countDiagonal=FALSE){
    K <- est[,idx]$wi
    if(use_tol){
        adj_mat <- ifelse(abs(K) > tol, 1, 0)
    }else{
        adj_mat <- ifelse(K !=0, 1, 0)
    }
    diag(adj_mat) <- 1
    true_pos <- which(adj_mat_True == 1, arr.ind = TRUE)
    # should be 1, and marked 1(TP)
    TP = sum(adj_mat[true_pos])-nrow(adj_mat_True)
    # all pos
    pos_sum = sum(adj_mat_True[true_pos])-nrow(adj_mat_True)
    #return(sum(adj_mat[true_pos])/total_edge_num)
    return(TP/pos_sum)
}

# param: 
#      idx: between 1 to Nlam
#      theta_0: true precision mat
# return: squared difference between estimated precision mat and true precision mat
mse_loss <- function(idx, est, theta_0){
    invSigma <- est[,idx]$wi
    sum(theta_0 - invSigma)^2
}

cal_num_edge <- function(K, total_edge_num, tol=1e-3){
  adj_mat_K <- ifelse(abs(K) > tol, 1, 0)
  diag(adj_mat_K) = 1
  return(sum(adj_mat_K) - nrow(K))
}
