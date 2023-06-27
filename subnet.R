rm(list = ls())
library(glmnet)
library(parallel)
library(deSolve)
library(orthopolynom)
library(pbapply)

#metadata <- read.csv("F:/IBD/1 preprocess/Metadata.csv")
load("F:/IBD/2 fit/df1_result.Rdata")
load("F:/IBD/2 fit/df2_result.Rdata")
load("F:/IBD/2 fit/df3_result.Rdata")

trans1 <- data.frame(df1_result[["trans_data"]])
trans2 <- data.frame(df2_result[["trans_data"]])
trans3 <- data.frame(df3_result[["trans_data"]])
fit1 <- data.frame(df1_result[["power_fit"]])
fit2 <- data.frame(df2_result[["power_fit"]])
fit3 <- data.frame(df3_result[["power_fit"]])
par1 <- df1_result[["power_par"]]
par2 <- df2_result[["power_par"]]
par3 <- df3_result[["power_par"]]
meta_index1 <- as.numeric(df1_result$Time)
meta_index2 <- as.numeric(df2_result$Time)
meta_index3 <- as.numeric(df3_result$Time)

load("F:/IBD/3 cluster/Clu_results_AR1_9.Rdata")
cluster <- Clu_results12[["clustered_df1"]]$cluster
table(cluster)
m1 = which(cluster==2)
m2 = which(cluster==3)
m3 = which(cluster==4)
m4 = which(cluster==5)
m5 = which(cluster==7)
m6 = which(cluster==9)
m7 = which(cluster==10)
m8 = which(cluster==11)
m9 = which(cluster==12)



#lasso
get_interaction <- function(data,col){
  set.seed(2023)
  n <- nrow(data)
  clean_data <- data
  gene_list <- list()
  m <- clean_data[,col]
  M <- clean_data[,-col]
  x_matrix <- M
  x_matrix <- as.matrix(x_matrix)
  #vec <- sapply(1:length(M[1,]),function(c)cor(m,M[,c]))
  #x_matrix <- M[,which( vec %in% -sort(-vec)[1:(n/log(n))] )]
  #x_matrix <- as.matrix(x_matrix)
  name <- colnames(clean_data)
  
  ridge1_cv <- cv.glmnet(x = x_matrix, y = m,type.measure = "mse",family="gaussian",nfold = 10,alpha = 0)
  best_ridge_coef <- abs(as.numeric(coef(ridge1_cv, s = ridge1_cv$lambda.min))[-1])
  
  fit_res <- cv.glmnet(x = x_matrix, y = m,type.measure = "mse", family="gaussian",
                       nfold = 10,alpha = 1,
                       #penalty.factor = 1/best_ridge_coef,
                       keep = TRUE)
  best_alasso_coef1 <- coef(fit_res, s = fit_res$lambda.1se)
  
  return_obj = list(ind.name = name[col],
                    dep.name = best_alasso_coef1@Dimnames[[1]][best_alasso_coef1@i[-1]+1],
                    coefficient = best_alasso_coef1@x[-1])
  
  return(return_obj)
}
relationship1 <- pblapply(1:length(m1),function(c) get_interaction(t(fit1[m1,]),c))


get_legendre_matrix <- function(x,legendre_order){
  legendre_coef <- legendre.polynomials(n = legendre_order, normalized=F)
  legendre_matrix <- as.matrix(as.data.frame(polynomial.values(
    polynomials = legendre_coef, x = scaleX(x, u = -1, v = 1))))
  colnames(legendre_matrix) <- paste0("legendre_",0:legendre_order)
  return(legendre_matrix[,2:(legendre_order+1)])
}
get_legendre_par <- function(y,legendre_order,x){
  #lm_method
  legendre_par <- as.numeric(coef(lm(y~get_legendre_matrix(x,legendre_order))))
  return(legendre_par)
}
legendre_fit <- function(par,x){
  legendre_order = length(par)
  fit <- sapply(1:length(par), function(c)
    par[c]*legendre.polynomials(n=legendre_order, normalized=F)[[c]])
  legendre_fit <- as.matrix(as.data.frame(polynomial.values(
    polynomials = fit, x = scaleX(x, u = -1, v = 1))))
  x_interpolation <- rowSums(legendre_fit)
  return(x_interpolation)
}

qdODEmod <- function(Time, State, Pars, power_par) {
  nn = length(Pars)
  ind_effect = paste0("alpha","*",names(State)[1])
  dep_effect = sapply(2:nn, function(c) paste0(paste0("beta",c-1),"*",names(State)[c]))
  dep_effect = paste0(dep_effect, collapse = "+")
  all_effect = paste0(ind_effect, "+", dep_effect)
  expr = parse(text = all_effect)
  
  with(as.list(c(State, Pars)), {
    dx = eval(expr)
    dy <- power_par[,1]*power_par[,2]*Time^(power_par[,2]-1)
    dind = alpha*x
    for(i in c(1:(nn-1))){
      tmp = paste0(paste0("beta",i),"*",paste0("y",i))
      expr2 = parse(text = tmp)
      assign(paste0("ddep",i),eval(expr2))
    }
    return(list(c(dx, dy, dind, mget(paste0("ddep",1:(nn-1))))))
  })
}
qdODE_ls <- function(pars, data, Time, power_par,real_data){
  n = length(pars)
  power_par = as.matrix(power_par)
  if (n==2) {
    Pars = c(alpha = pars[1], beta1 = pars[2:n])
    power_par = t(power_par)
    State = c(x=data[1,1],y1 = matrix(data[-1,1], nrow = n-1, ncol=1),ind = as.numeric(real_data[1]),dep = rep(0,n-1))
  } else{
    Pars = c(alpha = pars[1], beta = pars[2:n])
    State = c(x=data[1,1],y = matrix(data[-1,1], nrow = n-1, ncol=1),ind = as.numeric(real_data[1]),dep = rep(0,n-1))
  }
  out = as.data.frame(ode(func = qdODEmod, y = State, parms = Pars,
                          times = Time, power_par = power_par))
  X = as.numeric(real_data)
  fit = as.numeric(out[,2])
  penalty = sum((out$ind[out$ind<0])^2)
  sse = crossprod(X-fit) + penalty 
  return(sse)
}
qdODE_fit <- function(pars, data, Time, power_par, real_data, LOP_order = 2, new_time = NULL, n_expand = 100){
  n = length(pars)
  if (n==2) {
    Pars = c(alpha = pars[1], beta1 = pars[2:n])
    power_par = t(power_par)
    State = c(x=data[1,1],y1 = matrix(data[-1,1], nrow = n-1, ncol=1),ind = as.numeric(real_data[1]),dep = rep(0,n-1))
  } else{
    Pars = c(alpha = pars[1], beta = pars[2:n])
    State = c(x=data[1,1],y = matrix(data[-1,1], nrow = n-1, ncol=1),ind = as.numeric(real_data[1]),dep = rep(0,n-1))
  }
  out = as.data.frame(ode(func = qdODEmod, y = State, parms = Pars,
                          times = Time, power_par = power_par))
  out2 = data.frame(x = out[,1], y = real_data, y.fit = out[,2],
                    ind = out[,(n+2)], dep = out[,(n+3):(ncol(out))])
  colnames(out2)[4:ncol(out2)] = c(rownames(data)[1], rownames(data)[2:n])
  rownames(out2) = NULL
  
  all_LOP_par = sapply(2:ncol(out2),function(c)get_legendre_par(out2[,c], LOP_order, out2$x))
  
  if (is.null(new_time)) {
    time2 = seq(min(Time), max(Time), length = n_expand)
    out3 = apply(all_LOP_par, 2, legendre_fit, x = time2)
    out3 = cbind(time2, out3)
  } else{
    out3 = apply(all_LOP_par, 2, legendre_fit, x = new_time)
    out3 = cbind(new_time, out3)
  }
  colnames(out3) = colnames(out2)
  result = list(fit = out2,
                predict = data.frame(out3),
                LOP_par = all_LOP_par)
  return(result)
}

qdODE_all <- function(result, relationship, Time = Time, power_par = power_par,real_data = real_data, i, init_pars = 1, LOP_order = 2, method = "ls",
                      new_time = NULL, n_expand = 100, maxit = 1e5){
  
  variable = c(relationship[[i]]$ind.name, relationship[[i]]$dep.name)
  data = result[variable,]
  real_data = as.numeric(real_data[variable[1],])
  
  if (length(variable)<=1) {
    qdODE.est = NA
    result = NA
    return.obj <- append(result, list(ODE.value = NA,
                                      parameters = NA))
  } else{
    power_par = power_par[variable,][-1,]
    n = nrow(data)
    pars_int = c(init_pars,relationship[[i]]$coefficient)
    if (method == "ls") {
      qdODE.est <- optim(pars_int, qdODE_ls, data = data, Time = Time, power_par = power_par,real_data = real_data,
                         method = "Nelder-Mead",
                         #method = "BFGS",
                         control = list(trace = TRUE, maxit = maxit))
      
      result <- qdODE_fit(pars = qdODE.est$par,
                          data = data,
                          power_par = power_par,
                          Time = Time,
                          real_data = real_data)
      return.obj <- append(result, list(ODE.value = qdODE.est$value,
                                        parameters = qdODE.est$par))
    } else{
      qdODE.est <- optim(pars_int, qdODE_ls, data = data, Time = Time, power_par = power_par,real_data = real_data,
                         method = "Nelder-Mead",
                         control = list(trace = TRUE, maxit = maxit))
      
      result <- qdODE_fit(pars = qdODE.est$par,
                          data = data,
                          power_par = power_par,
                          Time = Time,
                          real_data = real_data)
      return.obj <- append(result, list(ODE.value = qdODE.est$value,
                                        parameters = qdODE.est$par))
    }
  }
  return(return.obj)
}
qdODE_parallel <- function(result,relationship,Time,power_par,real_data, reduction = FALSE, thread = 2, maxit = 1e4){
  
  core.number <- thread
  cl <- makeCluster(getOption("cl.cores", core.number))
  clusterEvalQ(cl, {require(orthopolynom)})
  clusterEvalQ(cl, {require(deSolve)})
  clusterExport(cl, c("qdODEmod", "qdODE_ls", "qdODE_fit", "qdODE_all","get_legendre_matrix",
                      "get_legendre_par","legendre_fit","result","relationship","Time","power_par","real_data","maxit"), envir=environment())
  result = pblapply(1:nrow(result),function(c) qdODE_all(result = result,
                                                         relationship = relationship,
                                                         Time = Time,
                                                         power_par = power_par,
                                                         real_data = real_data,
                                                         i = c,
                                                         maxit = maxit), cl = cl)
  stopCluster(cl)
  names(result) = rownames(result)
  names(relationship) = rownames(result)
  return_obj <- list(ode_result = result,
                     relationship = relationship)
  return(return_obj)
}


net1 = qdODE_parallel(result = fit1[m1,], relationship = relationship1, Time = meta_index1, 
                      power_par = par1[m1,], real_data = trans1, maxit = 1e4)
