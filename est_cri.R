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
dir_name = paste0(cur_dir, "/", gtype, "/", model,"_", gtype, "/")
# GET TEMP SAVE
temp_save_dir <- paste0(cur_dir,"/temp_est_save/")
temp_name <- paste0("est_", gtype,"_",model,"_",n,"_",p,"_",nt)
K_temp_name <- paste0("K_est_", gtype,"_",model,"_",n,"_",p,"_",nt)
filename <- paste0(temp_save_dir, temp_name,".txt")
data_filename <- paste0(temp_save_dir, K_temp_name,".rds")
cat(filename,"\n")
cat(data_filename,"\n")
if(file.exists(filename) && file.exists(data_filename)){
    est_K.all <- readRDS(data_filename)
    est <- est_K.all$est
    if(is.null(est)){
        cat("\nOpt K exploration reached elapsed Time Limit!\n",file=filename,append=TRUE)
        stop()
    }
    fixed_lam <- est_K.all$fixed_lam
    Nlam <- length(fixed_lam)
    S <- est_K.all$S
    data <- est_K.all$data
    adj_mat_True <- est_K.all$adj_mat_True
    theta_0 <- est_K.all$theta_0
    cat("\nSuccessfully loaded!\n")
}else{
    warning("\nIncomplete files! Run with same args and --r=1 for est_path.R\n")
    stop()
}

# SET PARAM
use_tol = TRUE
tol = 1e-3
nfold = 5
gen_pdf <- FALSE

K_est <- sapply(seq_along(fixed_lam), function(i){
    g <- est[,i]
    g$wi
})
dim(K_est) <- c(p,p,Nlam)
        
L_est <- sapply(seq_along(fixed_lam), function(i){
    g <- est[,i]
    g$loglik
})

EBIC_val <- sapply(seq_along(fixed_lam), function(i){
    K <- K_est[,,i]
    L <- L_est[i]
    EBIC_cal(S = S, K = K, L = L, N = n, gamma = 1, countDiagonal = FALSE)
})
BIC_val <- sapply(seq_along(fixed_lam), function(i){
    K <- K_est[,,i]
    L <- L_est[i]
    EBIC_cal(S = S, K = K, L = L, N = n, gamma = 0, countDiagonal = FALSE)
})

Half_BIC_val <- sapply(seq_along(fixed_lam), function(i){
    K <- K_est[,,i]
    L <- L_est[i]
    EBIC_cal(S = S, K = K, L = L, N = n, gamma = 0.5, countDiagonal = FALSE)
})

AIC_val <- sapply(seq_along(fixed_lam), function(i){
    K <- K_est[,,i]
    L <- L_est[i]
    AIC_cal(S = S, K = K, L = L, N = n)
})

CV_result = CV(n, nfold, data, fixed_lam)

BIC.opt_idx <- which.min(BIC_val)
Half_BIC_val.opt_idx <- which.min(Half_BIC_val)
EBIC.opt_idx <- which.min(EBIC_val)
LOGLIK.opt_idx <- which.max(L_est)
AIC.opt_idx <- which.min(AIC_val)
CV.opt_idx = which.max(CV_result$cv.mean)

hd_err <- sapply(seq_along(fixed_lam), function(i){
    hd_loss(i, est, adj_mat_True)
})
SHD.opt_idx <- which.min(hd_err)

cat("\n-------------IDX-------------\nLOG: ", LOGLIK.opt_idx, "\nAIC: ", AIC.opt_idx, "\nCV: ", CV.opt_idx, "\nBIC: ", BIC.opt_idx, "\n0.5 BIC: ",Half_BIC_val.opt_idx, "\nEBIC: ", EBIC.opt_idx, "\nSHD: ", SHD.opt_idx, file=filename,append=TRUE)
# corresponding lambda
cri_idx <- c(LOGLIK.opt_idx, AIC.opt_idx, CV.opt_idx, BIC.opt_idx, Half_BIC_val.opt_idx, EBIC.opt_idx, SHD.opt_idx)
cat("\n-------------VALUE-------------\nLOG: ", fixed_lam[LOGLIK.opt_idx], "\nAIC: ", fixed_lam[AIC.opt_idx], "\nCV: ", fixed_lam[CV.opt_idx], "\nBIC: ", fixed_lam[BIC.opt_idx], "\n0.5 BIC: ",fixed_lam[Half_BIC_val.opt_idx], "\nEBIC: ", fixed_lam[EBIC.opt_idx], "\nSHD: ", fixed_lam[SHD.opt_idx], file=filename,append=TRUE)

mse_err <- sapply(seq_along(fixed_lam), function(i){
    mse_loss(i, est, theta_0)
})
hd_cri <- sapply(cri_idx, function(i){hd_loss(i, est, adj_mat_True)})
mse_cri <-  sapply(cri_idx, function(i){mse_loss(i, est, theta_0)})
tpr_percentage <- sapply(cri_idx, function(i){TPR(i, est, adj_mat_True)})
fdr_percentage <- sapply(cri_idx, function(i){FDR(i, est, adj_mat_True, file=filename)})

cat("\nTPR: ", tpr_percentage, file=filename,append=TRUE)
cat("\nFDR: ", fdr_percentage, file=filename,append=TRUE)
cat("\nSHD: ", hd_cri, file=filename,append=TRUE)
cat("\nSHD loss: ", hd_err, file=filename,append=TRUE)
cat("\n--------------------------------------------\n", file=filename,append=TRUE)

# if generate pdf
if(gen_pdf){
    cat("\nSaving pdf of SHD and MSE............")
    pdf_file_name=paste0(dir_name, gtype, "_", model, "_P=", P, "_N=", N, ".pdf")
    pdf(file=paste0(pdf_file_name, collapse =".pdf"), width = 10, height = 4)
    par(mar = c(6,1,5,1), mfrow = c(1,2))
    plot(fixed_lam, hd_err, pch=4)
    points(fixed_lam[which.min(hd_err)], min(hd_err), pch = 24, cex=2, col="blue", bg="red", lwd=2)
    abline(v=fixed_lam[cri_idx], col=c("blue", "red", "green", "brown2" ,"orange", "cadetblue1", "blueviolet"),lty=c(1,5,6,3,2,4,3))
    title(paste0("hd: ", paste0(hd_cri, collapse = " ")))
    legend("topright", legend = c("lam_LOGLIK", "lam_AIC", "lam_CV", "lam_BIC", "lam_BIC.5","lam_EBIC", "lam_SHD"), col=c("blue", "red", "green", "brown2","orange", "cadetblue1", "blueviolet"), lty=c(1,5,6,3,2,4,3), text.font=4, cex=0.8, box.lty=0)
    mtext(paste0("PSR: ", paste0(sprintf("%.2f", psr_percentage), collapse = ",")), side = 3, line = 1, adj = 0)
    mtext(paste0("FDR: ", paste0(sprintf("%.2f", fdr_percentage), collapse = ",")), side = 3, adj = 0)

    plot(fixed_lam, mse_err, pch=4)
    points(fixed_lam[which.min(mse_err)], min(mse_err), pch = 24, cex=2, col="blue", bg="red", lwd=2)
    abline(v=fixed_lam[cri_idx], col=c("blue", "red", "green", "brown2" ,"orange", "cadetblue1", "blueviolet"),lty=c(1,5,6,3,2,4,3))
    title(paste0("mse: ", paste0(sprintf("%.0f", mse_cri), collapse = ","), " min: ", min(mse_err)))
    legend("topright", legend = c("lam_LOGLIK", "lam_AIC", "lam_CV", "lam_BIC", "lam_BIC.5","lam_EBIC"), col=c("blue", "red", "green", "brown2","orange", "cadetblue1", "blueviolet"), lty=c(1,5,6,3,2,4,3), text.font=4, cex=0.8, box.lty=0)
    dev.off()
}

