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
            "--r" = {r <- as.integer(split_arg[2])},
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
data_path  <- paste0(cur_dir, "/Data/")
g_file <- paste0(data_path,gtype,"_", p, ".rds")
lam_file <- paste0(data_path,gtype,"_",p,"_lamALL.rds")
cat(g_file,"\n")
cat(lam_file,"\n")

# GET TEMP SAVE
temp_save_dir <- paste0(cur_dir,"/temp_est_save/")
temp_name <- paste0("est_", gtype,"_",model,"_",n,"_",p)
filename <- paste0(temp_save_dir, temp_name,".txt")
data_filename <- paste0(temp_save_dir, temp_name,".rds")
cat("SAVE TO :", filename,"\n")

# NEED TO RERUN
if(r == TRUE){
    if(file.exists(filename)){
        file.remove(filename)
    }
    if(file.exists(data_filename)){
        file.remove(data_filename)
    }
}else if(!file.exists(filename) || !file.exists(data_filename)){
    warning("\nAt least one file is missing! Please set r=1 to rerun")
    
    if(file.exists(filename)){
        file.remove(filename)
    }
    if(file.exists(data_filename)){
        file.remove(data_filename)
    }
    stop()
}else{
    message("\nFILE ALREADY EXISTS!")
    stop()
}

# SET PARAM
use_tol = TRUE
tol = 1e-3
check_lammax = TRUE

cat("Loading ", g_file, lam_file, "............ ")
G <- readRDS(g_file)
lambda.all_algo <- readRDS(lam_file)
lambda.all <- lambda.all_algo[[model]]
cat("DONE!\n")

# get G data
n_list = as.integer(rownames(lambda.all))
n_idx <- which(n_list == n)
adj_mat_True = G$adj_mat_True
theta_0 = G$theta_0
data <- G$data[1:n,]
S = cov(data)
fixed_lam = lambda.all[n_idx,]
Nlam = length(fixed_lam)
cat(paste0("\nN=", n, ",  P=",p), file=filename,append=TRUE)
cat("\nLmabdas: ", fixed_lam, file=filename,append=TRUE)

# check max(fixed_lam)
if(check_lammax){
    if(model == "glasso"){
        glasso_lam_max <- glasso(s = S, rho = max(fixed_lam), thr = tol/100, approx = FALSE, trace = 0, penalize.diagonal = FALSE)
        if(use_tol){
            cat("\nGlasso: #zeros with lambda max = ", cal_num_edge(glasso_lam_max$wi, total_edge_num, tol = tol), file=filename,append=TRUE)
        }else{
            cat("\nGlasso: #zeros with lambda max = ", cal_num_edge(glasso_lam_max$wi, total_edge_num, tol = 0), file=filename,append=TRUE)
        }
    }else if(model == "clime"){
        clime_lam_max <- sugm(data = data, lambda = max(fixed_lam), method = "clime", sym ="and", prec = tol/100)
        cat("\nClime: #zeros with lambda max = ", sum(clime_lam_max$path[[1]]), file=filename,append=TRUE)
    }else if(model == "tiger"){
        tiger_lam_max <- sugm(data = data, lambda = max(fixed_lam), method = "tiger", sym ="and", prec = tol/100)
        cat("\nTiger: #zeros with lambda max = ", sum(tiger_lam_max$path[[1]]), file=filename,append=TRUE)
    }
}

# get model
est_path <- NULL
if(model == "glasso"){
    g_path <- glassopath(s = S, rholist = fixed_lam, thr = tol/100, approx = FALSE, trace = 0, penalize.diagonal = FALSE)
    # CHECK SYMMETRY for glasso
    # CHECK ALL GLASSO ESTIMATION wi PD (CAN IGNORE PRINT OUT, AND CHECK FOR WARNING MSG IN NEXT CELL)
    Sym_mat_theta <- matrix(0, Nlam)
    Sym_mat_wtol <- matrix(0, Nlam)
    Sym_mat_wotol <- matrix(0, Nlam)
    for(i in 1:Nlam){
        theta_hat <- g_path$wi[,,i]
        adj_mat1 <- ifelse(abs(theta_hat) > tol, 1, 0)
        adj_mat2 <- ifelse(theta_hat != 0, 1, 0)
        diag(adj_mat1) <- 1
        diag(adj_mat2) <- 1
        if(max(abs(theta_hat - t(theta_hat)))>tol){
            Sym_mat_theta[i] = 1
        }
        if(!isSymmetric(adj_mat1)){
            Sym_mat_wtol[i] = 1
        }
        if(!isSymmetric(adj_mat2)){
            Sym_mat_wotol[i] = 1
        }
    }

    cat("\nSYMMETRY CHECK:\n", file=filename,append=TRUE)
    if(length(which(Sym_mat_theta == 1))!=0){
        warning("Please check [Sym_mat_theta] to see where wi is NOT SYMMETRIC even with tolenrance\nCan decrease thr to calculate g_path...\n")
    }else{
        cat("With tolerance, all wi are symmetric\n", file=filename,append=TRUE)
    }

    if(length(which(Sym_mat_wtol == 1))!=0){
        warning("Please check [Sym_mat_wtol] to see where adjacency matrix with with tolerance is NOT SYMMETRIC\n", file=filename,append=TRUE)
    }
    if(length(which(Sym_mat_wotol == 1))!=0){
        warning("Please check [Sym_mat_wtol] to see where adjacency matrix with without tolerance is NOT SYMMETRIC\n", file=filename,append=TRUE)
    }
    est_path = g_path
}else if(model == "clime"){
    # est_path <- wrapper_with_timeout(time.limit, "clime", data, S, fixed_lam, tol)
    est_path <- sugm(data = data, lambda = rev(fixed_lam), method = "clime", sym ="and", prec = tol/100)
}else if(model == "tiger"){
    # est_path <- wrapper_with_timeout(time.limit, "tiger", data, S, fixed_lam, tol)
    est_path <- sugm(data = data, lambda = rev(fixed_lam), method = "tiger", sym ="and", prec = tol/100)
}else{
    stop("\nINVALID MODEL\n")
}

output <- list(est_path=est_path, n=n, p=p, gtype=gtype, model=model, fixed_lam=fixed_lam, S=S, data=data, adj_mat_True=adj_mat_True, theta_0=theta_0)

cat("Saving est for ", gtype, " via ", model, " with p, n = ", p, n, "in ", data_filename, "............")
saveRDS(output, file = data_filename)
cat("DONE\n")

