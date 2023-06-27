
#--data for cluster
load("F:/IBD/2 fit/df1_result.Rdata")
load("F:/IBD/2 fit/df2_result.Rdata")
load("F:/IBD/2 fit/df3_result.Rdata")

df <- cbind(df1_result[["power_fit"]],df2_result[["power_fit"]],df3_result[["power_fit"]])
meta_index1 <- as.numeric(df1_result$Time)
meta_index2 <- as.numeric(df2_result$Time)
meta_index3 <- as.numeric(df3_result$Time)

power_equation <- function(par,x){
  y <- par[1]*x^par[2]
  y
}
f <- function(y,times){
  tmp <- c(0.5,0.5)
  par_est <- function(par,x){
    sum( (y - power_equation(par,x))^2 )
  }
  r <- optim(tmp,par_est,x=as.numeric(times), method = "Nelder-Mead")
  return(r$par)
}
get_init_par <- function(data,k){
  set.seed(2023)
  #get initial pars based on k-means
  init_cluster <- kmeans(data,centers = k,iter.max = 1000)
  cuM <- init_cluster$centers 
  cusd <- diag(cov(df))
  init_curve_par <- cbind(t(sapply(1:k, function(c)f(cuM[c,1:265],times = meta_index1))),
                          t(sapply(1:k, function(c)f(cuM[c,266:411],times = meta_index2))),
                          t(sapply(1:k, function(c)f(cuM[c,412:546],times = meta_index3))))
                          
  init_SAD_par <- c(mean(cusd[1:265]),0.4,mean(cusd[266:411]),0.4,mean(cusd[412:546]),0.4)
  
  init_pro <- table(init_cluster$cluster)/nrow(data)
  return_object <- list(init_SAD_par,init_curve_par,init_pro)
  names(return_object)<-c("init_SAD_par","init_curve_par","init_pro")
  return(return_object)
}
get_cluster <- function(data,k,input){
  Delta <- 10000; iter <- 1; itermax <- 100;
  
  get_biSAD1 <- function(par){
    n1=265;n2=146;n3=135
    
    AR1_get_matrix <- function(par,n){
      var_cov<-matrix(0,n,n)
      
      for (i in 1:n){
        for (j in 1:n){
          var_cov[i,j]<- par[2]^abs(j-i)*par[1]
        }
      }
      return(var_cov)
}
    sig1 <- get_SAD1_covmatrix(par[1:2],n1)
    sig2 <- get_SAD1_covmatrix(par[3:4],n2)
    sig3 <- get_SAD1_covmatrix(par[5:6],n3)
    
    sig12 <- array(0, dim=c(n1,n2))
    sig21 <- array(0, dim=c(n2,n1))
    sigma1 <- cbind(sig1,sig12)
    sigma2 <- cbind(sig21,sig2)
    sigma12 <- rbind(sigma1,sigma2)
    
    sig123 <- array(0, dim=c(n1+n2,n3))
    sig321 <- array(0, dim=c(n3,n1+n2))
    sigma1 <- cbind(sigma12,sig123)
    sigma2 <- cbind(sig321,sig3)
    sigma123 <- rbind(sigma1,sigma2)
    
    return(sigma123)

    return(sigma12)
  } 
  power_equation <- function(par,x){
    y <- par[1]*x^par[2]
    y
  }
  mle <- function(par,data,prob){
    par1 <- par[1:6]
    par2 <- matrix(par[-c(1:6)],nrow = k)
    miu <- t( sapply(1:k, function(c)c(power_equation(par2[c,1:2],meta_index1),
                                       power_equation(par2[c,3:4],meta_index2),
                                       power_equation(par2[c,5:6],meta_index3))))
    temp_S <- sapply(1:k,function(c)dmvnorm(data,
                                            miu[c,],
                                            get_biSAD1(par1))*prob[c])

    LL <- sum(-log(rowSums(temp_S)))
    
    return(LL)
  }
  
  cat(paste0("Start Clu Calculation ","\n","Cluster_number=",k))
  while ( Delta > 10 && iter <= itermax ) {
    # initiation
    if(iter == 1){
      init_SAD_par <- input[[1]]
      init_curve_par <- input[[2]]
      pro <- input[[3]]
    }
    #E step, calculate the posterior probability
    old_par <- c(init_SAD_par,init_curve_par)
    LL_mem <- mle(old_par,data,pro)
    
    miu <- t( sapply(1:k, function(c)c(power_equation(init_curve_par[c,1:2],meta_index1),
                                       power_equation(init_curve_par[c,3:4],meta_index2),
                                       power_equation(init_curve_par[c,5:6],meta_index3))))
    mvn.c <- sapply(1:k, function(c) dmvnorm(data,
                                             miu[c,],
                                             get_biSAD1(init_SAD_par))*pro[c] )
    omega <- mvn.c/rowSums(mvn.c)
    
    #M step, calculate parameters
    pro <- colSums(omega)/sum(omega)
    
    new_par <- try(optim(old_par, mle, data=data, prob=pro, method = "Nelder-Mead"))
    if ('try-error' %in% class(new_par))
      break
    L_Value <- new_par$value
    init_SAD_par <- new_par$par[1:6]
    init_curve_par <- matrix(new_par$par[-c(1:6)],nrow = k)
    Delta <- abs(L_Value-LL_mem)
    
    cat("iter=",iter,"LL=",L_Value,'\n')
    iter <- iter+1; LL_mem <- L_Value
  } 
  
  BIC <- 2*(L_Value)+log(nrow(data))*length(old_par)
  
  cluster <- apply(omega,1,which.max)
  clustered_df1 <- data.frame(row.names(data),data[,1:265],cluster)
  clustered_df2 <- data.frame(row.names(data),data[,266:411],cluster)
  clustered_df3 <- data.frame(row.names(data),data[,412:546],cluster)
  
  get_plot1 <- function(clustered_data){
    colnames(clustered_data) <- c("marker",1:265,"cluster")
    long_df <- melt(clustered_data,c("marker","cluster"))
    colnames(long_df) <- c("marker","cluster","time","effect")
    p <-  ggplot()+geom_line(long_df,mapping=aes(as.numeric(as.character(time)),effect,group=marker,
                                                 colour= as.character(cluster)),alpha=1)+
      facet_wrap(long_df$cluster,scales = "fixed")+ 
      theme(legend.position="none") + xlab("Time")+ylab("generic_effect")
    return(p)
  }
  get_plot2 <- function(clustered_data){
    colnames(clustered_data) <- c("marker",266:411,"cluster")
    long_df <- melt(clustered_data,c("marker","cluster"))
    colnames(long_df) <- c("marker","cluster","time","effect")
    p <-  ggplot()+geom_line(long_df,mapping=aes(as.numeric(as.character(time)),effect,group=marker,
                                                 colour= as.character(cluster)),alpha=1)+
      facet_wrap(long_df$cluster,scales = "fixed")+ 
      theme(legend.position="none") + xlab("Time")+ylab("generic_effect")
    return(p)
  }
  get_plot3 <- function(clustered_data){
    colnames(clustered_data) <- c("marker",412:546,"cluster")
    long_df <- melt(clustered_data,c("marker","cluster"))
    colnames(long_df) <- c("marker","cluster","time","effect")
    p <-  ggplot()+geom_line(long_df,mapping=aes(as.numeric(as.character(time)),effect,group=marker,
                                                 colour= as.character(cluster)),alpha=1)+
      facet_wrap(long_df$cluster,scales = "fixed")+ 
      theme(legend.position="none") + xlab("Time")+ylab("generic_effect")
    return(p)
  }
  p1 <- get_plot1(clustered_df1);  p2 <- get_plot2(clustered_df2);  p3 <- get_plot3(clustered_df3)
  
  return_object <- list(init_SAD_par,init_curve_par,pro,LL_mem,BIC,
                        clustered_df1,clustered_df2,clustered_df3,p1,p2,p3) 
  cat("Finish biFunClu Calculation")
  names(return_object)<-c("SAD_par", "curve_par", "pro", "LL","BIC", 
                          "clustered_df1","clustered_df2","clustered_df3","plot1","plot2","plot3")
  return(return_object)
}

Clu_intial_pars <- get_init_par(data=df,k = 2)
Clu_results12 <- get_cluster(data=df,k = 2,input = Clu_intial_pars)