library("glmnet")

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Initialize variables
n <- NULL
p <- NULL
nt <- NULL
gtype <- NULL

# Parse command line arguments
for (arg in args) {
    split_arg <- strsplit(arg, "=")[[1]]
    switch(split_arg[1],
            "--n" = {n <- as.integer(split_arg[2])},
            "--p" = {p <- as.integer(split_arg[2])},
            "--nt" = {nt <- as.integer(split_arg[2])},
            "--m" = {model <- split_arg[2]},
            "--g" = {gtype <- split_arg[2]})
}

# Print the arguments
cat(n,"\n")
cat(p, "\n")
cat(nt, "\n")
cat(gtype, "\n")
cat(model, "\n")

# LOAD DATA and File
cur_dir = getwd()
data_path  <- paste0(cur_dir, "/Data/")
g_file <- paste0(data_path,gtype,"_", p, "_", nt, ".rds")
lam_file <- paste0(data_path,gtype,"_",p, "_", nt,"_lamALL.rds")
cat(g_file,"\n")
cat(lam_file,"\n")

# GET TEMP SAVE
temp_save_dir <- paste0(cur_dir,"/temp_est_save/")
temp_name <- paste0("est_", gtype,"_",model,"_",n,"_",p,"_",nt)
filename <- paste0(temp_save_dir, temp_name,".txt")
data_filename <- paste0(temp_save_dir, temp_name,".rds")
cat("SAVE TO :", data_filename,"\n")


cat("Loading ", g_file, lam_file, "............ ")
G <- readRDS(g_file)
lambda.all_algo <- readRDS(lam_file)
# use glasso for lambda choices
lambda.all <- lambda.all_algo[["glasso"]]
cat("DONE!\n")


NS <- function(pall, data, fixed_lam){
    adj_mat_pred_n <- list()
    Nlam <- length(fixed_lam)
    for(lam_idx in 1:Nlam){
        adj_mat_pred_n[[lam_idx]] <- adj_mat_template
    }

    for(p_idx in 1:length(pall)){
        X_temp <- data[,-p_idx]
        Y_temp <- data[,p_idx]
        X_idx <- pall[-p_idx]
        # cat("p ", p_idx, X_idx,"\n")
        fit <- glmnet(x = X_temp, y = Y_temp, alpha = 1, lambda = rev(fixed_lam), intercept = 0)
        coefs <- coef(fit, s = rev(fixed_lam))
        ne.pred <- matrix(ifelse(coefs[-1,] != 0, 1, 0), ncol=length(fixed_lam))
        idx_2_fill <- which(ne.pred != 0, arr.ind = TRUE)
        if(nrow(idx_2_fill)>0){
            for(j in 1:nrow(idx_2_fill)){
                rlam_idx <- idx_2_fill[j,2]
                lam_idx <- which(fixed_lam == rev(fixed_lam)[rlam_idx])
                x_idx_temp <- X_idx[idx_2_fill[j,1]]
                adj_mat_temp <- adj_mat_pred_n[[lam_idx]]
                adj_mat_temp[p_idx, x_idx_temp] <- 1
                adj_mat_pred_n[[lam_idx]] <- adj_mat_temp 
            }
        }
    }
    for(lam_idx in 1:Nlam){
        adj_mat_temp <- adj_mat_pred_n[[lam_idx]]
        adj_mat_temp <- ifelse(adj_mat_temp+t(adj_mat_temp) == 2, 1,0)
        diag(adj_mat_temp) <- 1
        adj_mat_pred_n[[lam_idx]] <- adj_mat_temp
    }
    return(adj_mat_pred_n)
}



pall <- c(1:p)

n_list = as.integer(rownames(lambda.all))
n_idx <- which(n_list == n)
adj_mat_True = G$adj_mat_True
adj_mat_template <- matrix(rep(0, length(adj_mat_True)), nrow = nrow(adj_mat_True), ncol = ncol(adj_mat_True))
theta_0 = G$theta_0
data <- G$data[1:n,]
S = cov(data)
fixed_lam = lambda.all[n_idx,]
cat(paste0("\nN=", n, ",  P=",p), file=filename,append=TRUE)
cat("\nLmabdas: ", fixed_lam, file=filename,append=TRUE)


ns_path <- NS(pall = pall, data = data, fixed_lam = fixed_lam)
output <- list(est_path=ns_path, n=n, p=p, gtype=gtype, model=model, fixed_lam=fixed_lam, S=S, data=data, adj_mat_True=adj_mat_True, theta_0=theta_0)
cat("Saving est for ", gtype, " via NS with p, n = ", p, n, "in ", filename, "............")
saveRDS(output, file = data_filename)
cat("DONE\n")