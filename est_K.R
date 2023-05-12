# IMPORT
library("glasso")
library("flare")
source("/home/zlyu/R_exp/GGM2/GGM_functions.R")


# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Initialize variables
n <- NULL
p <- NULL
r <- 1 # Default to rerun
model <- NULL
gtype <- NULL

# Parse command line arguments
for (arg in args) {
    split_arg <- strsplit(arg, "=")[[1]]
    switch(split_arg[1],
            "--n" = {n <- as.integer(split_arg[2])},
            "--p" = {p <- as.integer(split_arg[2])},
            "--m" = {model <- split_arg[2]},
            "--g" = {gtype <- split_arg[2]})
}


arg_len <- sapply(c(n,p,gtype,model), length)
na_check <- sapply(c(n,p,gtype,model), is.na)
if(sum(arg_len) < 4 || any(na_check)){
    warning("Incomplete args!")
    stop()
}


# Print the arguments
cat(n,"\n")
cat(p, "\n")
cat(model, "\n")
cat(gtype, "\n")

# LOAD DATA and File
cur_dir = getwd()
# data_path  <- paste0(cur_dir, "/Data/")
# g_file <- paste0(data_path,gtype,"_", p, ".rds")
# lam_file <- paste0(data_path,gtype,"_",p,"_lamALL.rds")
# cat(g_file,"\n")
# cat(lam_file,"\n")

# GET TEMP SAVE
temp_save_dir <- paste0(cur_dir,"/temp_est_save/")
temp_name <- paste0("est_", gtype,"_",model,"_",n,"_",p)
filename <- paste0(temp_save_dir, temp_name,".txt")
data_filename <- paste0(temp_save_dir, temp_name,".rds")
if(file.exists(filename) && file.exists(data_filename)){
    est_path.all <- readRDS(data_filename)
    est_path <- est_path.all$est_path
    if(is.null(est_path)){
        cat("\nZero pattern exploration reached elapsed Time Limit!\n",file=filename,append=TRUE)
        stop()
    }
    fixed_lam <- est_path.all$fixed_lam
    S <- est_path.all$S
    data <- est_path.all$data
    adj_mat_True <- est_path.all$adj_mat_True
    theta_0 <- est_path.all$theta_0
    cat("\nSuccessfully loaded!\n")
}else{
    warning("\nIncomplete files! Run with same args and --r=1 for est_path.R\n")
    stop()
}

# SET PARAM
use_tol = TRUE
tol = 1e-3
new_data_filename <- paste0(temp_save_dir, "K_est_", gtype,"_",model,"_",n,"_",p,".rds")

est <- NULL
if(model == "glasso"){
    # a list of optimal lasso based on zero patterns of wi in g_path
    est <- sapply(seq_along(fixed_lam), function(i){
        invSigma <- est_path$wi[,,i]
        net_est(S = S, invSigma = invSigma, N = n, use_tol = use_tol)
    })
}else if(model == "clime" || model == "tiger"){
    # a list of optimal lasso based on zero patterns of wi in g_path
    est <- sapply(seq_along(fixed_lam), function(i){
        i_r = which(est_path$lambda == fixed_lam[i])
        invSigma <- est_path$path[[i_r]]
        net_est(S = S, invSigma = invSigma, N = n, use_tol = use_tol, model = "nonglasso")
    })
}

output <- list(est=est, n=n, p=p, gtype=gtype, model=model, fixed_lam=fixed_lam, S=S, data=data, adj_mat_True=adj_mat_True, theta_0=theta_0)

cat("Saving K_est for ", gtype, " via ", model, " with p, n = ", p, n, "in ", new_data_filename, "............")
saveRDS(output, file = new_data_filename)
# cat("DONE\n")

