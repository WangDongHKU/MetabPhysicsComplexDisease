rm(list = ls())
library(glmnet)
library(parallel)
library(deSolve)
library(orthopolynom)
library(pbapply)
library(ggplot2)
library(patchwork)
library(cowplot)

metadata <- read.csv("F:/IBD/1 preprocess/Metadata.csv")
load("F:/IBD/1 preprocess/data_nonIBD.Rdata")
load("F:/IBD/1 preprocess/data_CD.Rdata")
load("F:/IBD/1 preprocess/data_UC.Rdata")

#power fit

power_equation <- function(x, power_par){ t(sapply(1:nrow(power_par),function(c) power_par[c,1]*x^power_par[c,2] ) )}
power_equation_base <- function(x, y){
  x <- as.numeric(x)
  y <- as.numeric(y)
  min_value = min(y[y!=0])
  
  set.seed(2023)
  lmFit <- lm( log( y + runif(1, min = 0, max = min_value))  ~ log(x))
  coefs <- coef(lmFit)
  a <- exp(coefs[1])
  b <- coefs[2]
  
  model <- try(nls(y~a*x^b,start = list(a = a, b = b),
                   control = nls.control(maxiter = 1e3, minFactor = 1e-200)))
  if( 'try-error' %in% class(model)) {
    result = NULL
  }
  else{
    result = model
  }
  return(result)
}
power_equation_all <- function(x,y, maxit=1e2){
  result <- power_equation_base(x,y)
  iter <- 1
  while( is.null(result) && iter <= maxit) {
    iter <- iter + 1
    try(result <- power_equation_base(x,y))
  }
  return(result)
}
power_equation_fit <- function(data, trans = log10, thread = 2) {
  data = data[,order(colSums(data))]
  if ( is.null(trans)) {  
    X = colSums(data)
    trans_data = data
  } else{
    X = trans(colSums(data)+1)
    trans_data = trans(data+1)
  }
  colnames(trans_data) = X
  
  core.number <- thread
  cl <- makeCluster(getOption("cl.cores", core.number))
  clusterExport(cl, c("power_equation_all", "power_equation_base", "trans_data", "X"), envir = environment())
  all_model = parLapply(cl = cl, 1:nrow(data), function(c) power_equation_all(X, trans_data[c,]))
  stopCluster(cl)
  
  
  names(all_model) = rownames(data)
  no = which(sapply(all_model, length)>=1)
  all_model2 = all_model[no]
  data2 = data[no,]
  trans_data2 = trans_data[no,]
  
  new_x = X
  power_par = t(vapply(all_model2, coef, FUN.VALUE = numeric(2), USE.NAMES = TRUE))
  power_fit = t(vapply(all_model2, predict, newdata = data.frame(x=new_x),
                       FUN.VALUE = numeric(length(new_x)), USE.NAMES = TRUE))
  
  colnames(power_fit) = new_x
  result = list(original_data = data2, trans_data = trans_data2,
                power_par = power_par, power_fit = power_fit,
                Time = X)
  return(result)
}

df1_result = power_equation_fit(data = df_CD, trans = log10, thread = 2)
df2_result = power_equation_fit(data = df_UC, trans = log10, thread = 2)
df3_result = power_equation_fit(data = df_nonIBD, trans = log10, thread = 2)

#save(df3_result,file = "F:/IBD/2 fit/df3_result.Rdata")


df1 <- data.frame(cbind(as.numeric(df1_result$Time),t(df1_result[["trans_data"]])))
linedata1 <- data.frame(cbind(as.numeric(df1_result$Time),t(df1_result[["power_fit"]])))
df2 <- data.frame(cbind(as.numeric(df2_result$Time),t(df2_result[["trans_data"]])))
linedata2 <- data.frame(cbind(as.numeric(df2_result$Time),t(df2_result[["power_fit"]])))
df3 <- data.frame(cbind(as.numeric(df3_result$Time),t(df3_result[["trans_data"]])))
linedata3 <- data.frame(cbind(as.numeric(df3_result$Time),t(df3_result[["power_fit"]])))


get_p <- function(i,data_df,data_line){
  p <- ggplot()+
    geom_point(data = data_df,mapping = aes(as.numeric(data_df[,1]),as.numeric(data_df[,i+1])),
               color="#DB7344",size=5,shape=1,alpha=1)+
    geom_line(data = data_line,mapping = aes(as.numeric(data_line[,1]),as.numeric(data_line[,i+1])),
              color="#276A6C",linewidth=2)+
    theme_bw()+theme(panel.grid =element_blank())+
    labs(title = as.vector(rownames(df1_result$original_data)[i]))+
    theme(plot.title = element_text(vjust = -6)) + 
    theme(plot.title = element_text(hjust = 0.5))+
    xlab('')+ylab('')+
    theme(plot.margin = unit(c(0,0,0,0),"lines"))
  
  return(p)
}


