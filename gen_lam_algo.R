library(utils)

cal_num_edge <- function(K, total_edge_num, tol){
    adj_mat_K <- ifelse(abs(K) > tol, 1, 0)
    diag(adj_mat_K) = 1
    return(sum(adj_mat_K) - nrow(K))
}

cal_lam <- function(data,total_edge_num,model="glasso",Nlam = 100,factor = 0.95,edge_ratio = 2,tol=1e-3){
    library("glasso")
    S <- cov(data)
    n <- nrow(data)
    p <- ncol(data)
    lam.max_candidate <- sapply(seq(p), function(i){
        norm(t(data[,-i])%*%data[,i], "I")/n
    })
    lambda.max = max(lam.max_candidate)
    # Nlam = 100
    nlam = 1
    lam.min_candidate <- rep(0,Nlam)
    lam_temp = lambda.max
    if(model == "glasso"){
        g_temp = glasso(s = S, rho = lam_temp, thr = tol/100, approx = FALSE, trace = 0, penalize.diagonal = FALSE)
        while (cal_num_edge(g_temp$wi, total_edge_num, tol) < edge_ratio * total_edge_num && nlam <= Nlam) {
            lam.min_candidate[nlam] = lam_temp
            lam_temp = lam_temp*factor
            g_temp = glasso(s = S, rho = lam_temp, thr = tol/100, approx = FALSE, trace = 0, penalize.diagonal = FALSE)
            nlam = nlam + 1
        }
    }else if(model == "clime"){
        library(flare)
        factor <- 0.9
        c_temp = sugm(data = data, lambda = lam_temp, method = "clime", sym ="and", prec = tol/100)
        while (sum(c_temp$path[[1]]) < edge_ratio * total_edge_num && nlam <= Nlam) {
            lam.min_candidate[nlam] = lam_temp
            lam_temp = lam_temp*factor
            c_temp = sugm(data = data, lambda = lam_temp, method = "clime", sym ="and", prec = tol/100)
            nlam = nlam + 1
        }
    }else if(model == "tiger"){
        library(flare)
        factor <- 0.9
        t_temp = sugm(data = data, lambda = lam_temp, method = "tiger", sym ="and", prec = tol/100)
        while (sum(t_temp$path[[1]]) < edge_ratio * total_edge_num && nlam <= Nlam) {
            lam.min_candidate[nlam] = lam_temp
            lam_temp = lam_temp*factor
            t_temp = sugm(data = data, lambda = lam_temp, method = "clime", sym ="and", prec = tol/100)
            nlam = nlam + 1
        }
    }
    lambda.min = min(lam.min_candidate[which(lam.min_candidate!= 0)])
    fixed_lam <- exp(seq(log(lambda.min), log(lambda.max), length.out = Nlam))
    return(fixed_lam)
}

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

nt <- NULL
gtype <- NULL

# Parse command line arguments
for (arg in args) {
    split_arg <- strsplit(arg, "=")[[1]]
    switch(split_arg[1],
            "--nt" = {nt <- as.integer(split_arg[2])},
            "--g" = {gtype <- split_arg[2]})
}

cat("nt", nt,"\n")
cat("gtype", gtype,"\n")

# get cur dir
cur_dir <- getwd()
file_path <- paste0(cur_dir, "/Data/")

# get all data files
file_list <- list.files(file_path, pattern = paste0("_", nt, "\\.rds$"))
file_list <- file_list[grepl(gtype, file_list)]
G_file <- file_list[!grepl("lam", file_list)]

# # # DELETE BELOW TWO LINES AFTER TESTING.......................
# args = commandArgs(TRUE)
# gtype_str <- unlist(args)
# G_file <- G_file[grepl(gtype_str, G_file)]

file_list <- paste0(file_path, G_file)
# cat(file_list)



n_val <- Sys.getenv('N_VAL')
n_list <- as.integer(strsplit(n_val, ",")[[1]])
# n_list <- c(10,20,50,100,200,500,800,1000)
Nlam <- 100



for(file_name in file_list){
    G <- readRDS(file_name)
    data = G$data
    p = ncol(data)
    graph = G$gtype
    total_edge_num = G$total_edge_num
    # lam_file <- paste0(file_path, graph,"_lam.rds")
    lam_file <- paste0(strsplit(file_name,".rds")[1],"_lamALL.rds")
    # cat(lam_file,"\n")
    lams_g <- t(sapply(n_list, function(n){
        subdata <- data[1:n,]
        cal_lam(data=subdata, total_edge_num=total_edge_num,Nlam=Nlam,model="glasso")
    }))
    lams_c <- t(sapply(n_list, function(n){
        subdata <- data[1:n,]
        cal_lam(data=subdata, total_edge_num=total_edge_num,Nlam=Nlam,model="clime")
    }))
    lams_t <- t(sapply(n_list, function(n){
        subdata <- data[1:n,]
        cal_lam(data=subdata, total_edge_num=total_edge_num,Nlam=Nlam,model="tiger")
    }))
    lams_g = matrix(lams_g, nrow=length(n_list), ncol=Nlam)
    rownames(lams_g) <- n_list
    lams_c = matrix(lams_c, nrow=length(n_list), ncol=Nlam)
    rownames(lams_c) <- n_list
    lams_t = matrix(lams_t, nrow=length(n_list), ncol=Nlam)
    rownames(lams_t) <- n_list
    lams <- list(glasso = lams_g, clime = lams_c, tiger = lams_t)
    cat("Saving Lambda for ", graph, " with dim = ", p, "in ", lam_file, "............ ")
    saveRDS(lams, file = lam_file)
    cat("DONE\n")
}
