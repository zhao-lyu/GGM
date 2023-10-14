# IMPORT
library("glasso")
library("flare")
source(paste0(getwd(), "/GGM_functions.R"))


# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Initialize variables
n <- NULL
p <- NULL
nt <- NULL
model <- NULL
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


arg_len <- sapply(c(n,p,gtype,model), length)
na_check <- sapply(c(n,p,gtype,model), is.na)
if(sum(arg_len) < 4 || any(na_check)){
    warning("Incomplete args!")
    stop()
}


# Print the arguments
cat(n,"\n")
cat(p, "\n")
cat(nt, "\n")
cat(model, "\n")
cat(gtype, "\n")

# LOAD DATA and File
cur_dir = getwd()

# GET TEMP SAVE
temp_save_dir <- paste0(cur_dir,"/temp_est_save/")
temp_name <- paste0("est_", gtype,"_",model,"_",n,"_",p,"_",nt)
filename <- paste0(temp_save_dir, temp_name,".txt")
data_filename <- paste0(temp_save_dir, temp_name,".rds")
cat(filename,"\n")
cat(data_filename,"\n")
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
new_data_filename <- paste0(temp_save_dir, "K_est_", gtype,"_",model,"_",n,"_",p,"_",nt,".rds")
cat(new_data_filename,"\n")

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
}else if(model == "NS"){
    # a list of optimal lasso based on zero patterns of wi in g_path
    est <- sapply(seq_along(fixed_lam), function(i){
        invSigma <- est_path[[i]]
        net_est(S = S, invSigma = invSigma, N = n, use_tol = use_tol, model = "NS")
    })
}

output <- list(est=est, n=n, p=p, gtype=gtype, model=model, fixed_lam=fixed_lam, S=S, data=data, adj_mat_True=adj_mat_True, theta_0=theta_0)

cat("Saving K_est for ", gtype, " via ", model, " with p, n = ", p, n, "in ", new_data_filename, "............")
saveRDS(output, file = new_data_filename)
cat("DONE\n")

