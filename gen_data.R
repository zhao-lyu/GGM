# Generate data based on graph type and P

# generator func
generator <- function(N,P,g,ss){
    set.seed(ss)
    cat("ss: ", ss,"\n")
    # band and sf by flare
    if(g == "band"){
        library('flare')
        G <- sugm.generator(n=N, d=P, graph=g, v=0.3,u=0,g=1,seed=ss)
        theta_0 = G$omega
        sigma_0 = G$sigma
        data = G$data
        S = G$sigmahat
        adj_mat_True = G$theta
    }else if(g == "sf"){
        library('flare')
        G <- sugm.generator(n=N, d=P, graph="scale-free",seed=ss)
        theta_0 = G$omega
        sigma_0 = G$sigma
        data = G$data
        S = G$sigmahat
        adj_mat_True = G$theta
    }
    else if(g == "er"){
        library("igraph")
        gnp_m  = P-1
        G_renyi = erdos.renyi.game(n = P, p.or.m = gnp_m, type = "gnm", directed = FALSE, loops = FALSE)
        adj_mat_True = as_adjacency_matrix(graph = G_renyi, type = "both", sparse = TRUE)
        diag(adj_mat_True) = 0
        v = 0.3
        u = 0
        theta_0 = adj_mat_True*v
        diag(theta_0) = abs(min(eigen(theta_0)$values)) + 0.1 + u
        sigma_0 = cov2cor(solve(theta_0))
        theta_0 = solve(sigma_0)

    }else if(g == "knn"){
        library("mstknnclust")
        library("stats")
        library("igraph")
        # random data
        x <- matrix(runif(P*P, min = -10, max = 10), nrow=P, ncol=P)
        d <- base::as.matrix(stats::dist(x, method="euclidean"))
        cg <- generate.complete.graph(1:nrow(x),d)
        G_knn <- generate.knn(cg, suggested.k=2)
        adj_mat_True = as_adjacency_matrix(graph = G_knn$knn.graph, type = "both", sparse = TRUE)
        adj_mat_True[which(as.matrix(adj_mat_True) > 0, arr.ind = TRUE)] <- 1
        diag(adj_mat_True) = 0
        v = 0.3
        u = 0
        theta_0 = adj_mat_True*v
        diag(theta_0) = abs(min(eigen(theta_0)$values)) + 0.1 + u
        sigma_0 = cov2cor(solve(theta_0))
        theta_0 = solve(sigma_0)
    }
    # CONSTRUCT DIST FOR DATA
    library("MASS")
    # mean
    mu <- rep(0, P)
    # data: N \times P; each row is an observation
    data <- mvrnorm(N, mu, sigma_0)
    # sample covariance matrix
    S <- cov(data)
    diag(adj_mat_True) = 1 
    total_edge_num = sum(adj_mat_True) - P

    return(list(theta_0=theta_0, sigma_0=sigma_0, data=data, S=S, adj_mat_True=adj_mat_True, total_edge_num=total_edge_num, gtype=g))
}

args <- commandArgs(trailingOnly = TRUE)

# Parse the command line argument
nt <- args[1]
nt <- as.integer(nt)
cat("random seed set to", nt, "\n")

# get args for graph type and P
p_val <- Sys.getenv('P_VAL')
n_val <- Sys.getenv('N_VAL')
g_val <- Sys.getenv('G_VAL')
dist.Gtype <- strsplit(g_val,",")[[1]]
dist.p <- as.integer(strsplit(p_val, ",")[[1]])
n <- max(as.integer(strsplit(n_val, ",")[[1]]))

if(n<=5000){
    n=5000
}

cat(dist.p,"\n",dist.Gtype, "\n")
cat("Generating", n, "sample data............\n")

# save dir
cur_dir = getwd()
# TODO: ADD CHECK EXISTING FOLDER
save_main_dir = paste0(cur_dir, "/", "Data/")
cat(save_main_dir)

for(i in 1:length(dist.Gtype)){
    graph <- dist.Gtype[i]
    save_dir = paste0(save_main_dir,graph,"_")
    for(p_idx in 1:length(dist.p)){
        P <- as.integer(dist.p[p_idx])
        G <- generator(N = n, P = P, g = graph, ss=nt)
        save_dir_p = paste0(save_dir, P, "_", nt, ".rds")
        cat("Saving dist for ", graph, " with dim = ", P, "in ", save_dir_p, "............")
        saveRDS(G, file = save_dir_p)
        cat("DONE\n")
    }
}

