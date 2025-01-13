library(deSolve)
library(rstan)
library(ggplot2)
library(psych)
library(bayesplot)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


# Parameters is a list with "params" and "initial conditions"
build_data <- function(myseed=0, noise_level=0.01, species=3, tmax=0, n_points=25, 
                       folder=".", plot_this=TRUE, parameters=list()) {
  
  # Lotka-Volterra equations for "K" species
  lotka_volterra <- function(t, initial_condition, parameters) {
    K <- length(initial_condition)
    x <- initial_condition
    r <- parameters[1:K]
    b <- t(matrix(parameters[-(1:K)],K,K))
    dxdt <- rep(0,K)
    for(i in 1:K) {
      dxdt[i] <- r[i]*x[i]
      for(j in 1:K) {
        dxdt[i] <- dxdt[i] + b[i,j]*x[i]*x[j]
      }
    }
    return(list(dxdt))
  }
  if(myseed ==0) myseed <- sample.int(100000,1)
  set.seed(myseed)
  
  if(length(parameters)==0) {
    var_names <- c(); 
    for(i in 1:species) var_names <- c(var_names,sprintf("r%d",i))
    for(i in 1:species) {
      for(j in 1:species) {
        var_names <- c(var_names,sprintf("a%d%d",i,j))
      }
    }
    
    flag <- FALSE
    while(!flag) {
      r0 <- rexp(species,.5)
      a0 <- matrix(rnorm(species^2),species,species)
      for(i in 1:species) a0[i,i]  <- -rexp(1,.5)
      params <- c(r0,as.vector(t(a0)))
      names(params) <- var_names
      
      xss <- solve(a0,-r0)
      d0 <- diag(xss)
      eigenvalues <- Re(eigen(d0 %*% a0)$values)
      # print(eigenvalues <- Re(eigen(d0 %*% a0)$values))
      flag <- all(eigenvalues <=0)
    }
    
    
    # Initial initial_condition
    to_eval <- "initial_condition <- c("
    for(i in 1:(species-1)){
      to_eval <- paste0(to_eval,sprintf("X%d = rexp(1,10), ",i))
      
    }
    to_eval <- paste0(to_eval,sprintf("X%d = rexp(1,10))",species))
    eval(parse(text=to_eval))
    # print(initial_condition)
  }
  else {
    cat("Species:",species,"\n")
    species <- parameters$K
    cat("Species:",species,"\n")
    myseed <- parameters$myseed
    params <- parameters$params
    initial_condition <- parameters$initial_condition
    r0 <- params[1:species]
    a0 <- t(matrix(params[-(1:species)],species,species))
    xss <- solve(a0,-r0)
    d0 <- diag(xss)
    print(eigenvalues <- Re(eigen(d0 %*% a0)$values))
  }
  
  
  # Time points (estimate from lowest eigenvalue)
  if(tmax==0) tmax <- -20/min(eigenvalues)
  times <- seq(0,tmax, length=n_points)
  
  # Solve differential equations
  out <- ode(y = initial_condition, times = times, func = lotka_volterra, parms = params)
  
  # add noise
  which(is.na(out))
  for(i in 2:(species+1)) out[,i] <- rlnorm(length(out[,i]),log(1e-6+abs(out[,i])),noise_level)
  
  
  # Plot
  if(plot_this) {
    x11("",12,8)
    par(mfrow=c(1,1))
    par(mar=c(5,5,3,3))
    # colors <- c("#B2001D", "#ED8B45", "#998EC3", "#007C7C", "#800080", "#505050")
    colors <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF")
    
    symbols <- rep(c(15,17,19,18),3)
    matplot(out[,1],out[,-1],type='b',pch=symbols,cex=1.5,cex.lab=1.7,cex.axis=1.4,
            xlab="Time",ylab="Populations",lwd=2,lty=1:3,col=colors)
    
    lgnd <- "c(expression(x[1])"
    for(idx in 2:species) lgnd <- paste0(lgnd,sprintf(",expression(x[%d])",idx))
    lgnd <- parse(text=paste0(lgnd,")"))
    # legend('topright',legend=c(expression(x[1]),expression(x[2]),expression(x[3])),
    legend('topright',legend=eval(lgnd),
           pch=symbols[1:species],col=colors[1:species],cex=2,bg="white",lty=1:species,lwd=3)
    
    
    dev.copy2pdf(file=sprintf("%s/simulation_%d_%f.pdf",folder,myseed,noise_level))
    Sys.sleep(1)
    dev.off()
  }
  
  output_list <- list(K=species,out=out,eigenvalues=eigenvalues,noise_level=noise_level,
                      myseed=myseed,params=params,initial_condition=initial_condition,
                      times=times, lotka_volterra=lotka_volterra)
}
########################################################################################
build_data_obs <- function(myseed=0, noise_level=0.01, species=7, obs=127, tmax=0, n_points=25, 
                       folder=".", plot_this=TRUE, parameters=list()) {
  
  # Lotka-Volterra equations for "K" species
  lotka_volterra <- function(t, initial_condition, parameters) {
    K <- length(initial_condition)
    x <- initial_condition
    r <- parameters[1:K]
    b <- t(matrix(parameters[-(1:K)],K,K))
    dxdt <- rep(0,K)
    for(i in 1:K) {
      dxdt[i] <- r[i]*x[i]
      for(j in 1:K) {
        dxdt[i] <- dxdt[i] + b[i,j]*x[i]*x[j]
      }
    }
    return(list(dxdt))
  }
  if(myseed ==0) myseed <- sample.int(100000,1)
  set.seed(myseed)
  
  if(length(parameters)==0) {
    var_names <- c(); 
    for(i in 1:species) var_names <- c(var_names,sprintf("r%d",i))
    for(i in 1:species) {
      for(j in 1:species) {
        var_names <- c(var_names,sprintf("a%d%d",i,j))
      }
    }
    
    flag <- FALSE
    while(!flag) {
      r0 <- rexp(species,.5)
      a0 <- matrix(rnorm(species^2),species,species)
      for(i in 1:species) a0[i,i]  <- -rexp(1,.5)
      params <- c(r0,as.vector(t(a0)))
      names(params) <- var_names
      
      xss <- solve(a0,-r0)
      d0 <- diag(xss)
      eigenvalues <- Re(eigen(d0 %*% a0)$values)
      # print(eigenvalues <- Re(eigen(d0 %*% a0)$values))
      flag <- all(eigenvalues <=0)
    }
    
    
    # Initial initial_condition
    to_eval <- "initial_condition <- c("
    for(i in 1:(species-1)){
      to_eval <- paste0(to_eval,sprintf("X%d = rexp(1,10), ",i))
      
    }
    to_eval <- paste0(to_eval,sprintf("X%d = rexp(1,10))",species))
    eval(parse(text=to_eval))
    # print(initial_condition)
  }
  else {
    cat("Species:",species,"\n")
    species <- parameters$K
    cat("Species:",species,"\n")
    myseed <- parameters$myseed
    params <- parameters$params
    obs_bin <- to_binary(obs,species)
    initial_condition <- parameters$initial_condition*obs_bin
    cat("###############\n",initial_condition,"\n######\n")
    
    
    r0 <- params[1:species]
    a0 <- t(matrix(params[-(1:species)],species,species))
    xss <- solve(a0,-r0)
    d0 <- diag(xss)
    print(eigenvalues <- Re(eigen(d0 %*% a0)$values))
  }
  
  
  # Time points (estimate from lowest eigenvalue)
  if(tmax==0) tmax <- -20/min(eigenvalues)
  times <- seq(0,tmax, length=n_points)
  initial_condition <- initial_condition*obs
  # Solve differential equations
  out <- ode(y = initial_condition, times = times, func = lotka_volterra, parms = params)
  
  # add noise
  which(is.na(out))
  for(i in 2:(species+1)) out[,i] <- rlnorm(length(out[,i]),log(1e-6+abs(out[,i])),noise_level)
  
  
  # Plot
  if(plot_this) {
    x11("",12,8)
    par(mfrow=c(1,1))
    par(mar=c(5,5,3,3))
    # colors <- c("#B2001D", "#ED8B45", "#998EC3", "#007C7C", "#800080", "#505050")
    colors <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF")
    
    symbols <- rep(c(15,17,19,18),3)
    matplot(out[,1],out[,-1],type='b',pch=symbols,cex=1.5,cex.lab=1.7,cex.axis=1.4,
            xlab="Time",ylab="Populations",lwd=2,lty=1:3,col=colors)
    
    lgnd <- "c(expression(x[1])"
    for(idx in 2:species) lgnd <- paste0(lgnd,sprintf(",expression(x[%d])",idx))
    lgnd <- parse(text=paste0(lgnd,")"))
    # legend('topright',legend=c(expression(x[1]),expression(x[2]),expression(x[3])),
    legend('topright',legend=eval(lgnd),
           pch=symbols[1:species],col=colors[1:species],cex=2,bg="white",lty=1:species,lwd=3)
    
    
    dev.copy2pdf(file=sprintf("%s/simulation_%d_%f.pdf",folder,myseed,noise_level))
    Sys.sleep(1)
    dev.off()
  }
  
  output_list <- list(K=species,out=out,eigenvalues=eigenvalues,noise_level=noise_level,
                      myseed=myseed,params=params,initial_condition=initial_condition,
                      times=times, lotka_volterra=lotka_volterra)
}
########################################################################################

############################ Main call to stan #########################################
fit_stan_rk4 <- function(data, model="lotka_volterra_rk_K_dimensions.stan", 
                         folder=".", n_chains=3, n_iter=4000, warmup=3000) {
  # Prepare data for Stan model
  N = nrow(data$out)
  Y = data$out[,-1]
  T = data$out[,1]
  K = data$K
  stan_data <- list(N=N,Y=Y,T=T,K=K)
  
  recover <- list(stan_data=stan_data,data=data)
  saveRDS(recover,sprintf("%s/recover_files/recover_%d_%f.rda",folder,data$myseed,data$noise_level))
  
  # Fit model
  fit <- stan(file = model, data = stan_data, iter = n_iter, warmup = warmup,chains = n_chains)
  
  return(fit)
}
########################################################################################


######################### Built-it stan plotting ########################################
make_stan_plots <- function(fit,data,folder=".") {
  x11("",12,8)
  pars <- c(); 
  color_scheme_set("teal")
  for(i in 1:data$K) for(j in 1:data$K) pars <- c(pars,sprintf("b[%d,%d]",i,j))
  print(rstan::stan_plot(fit,pars=pars))
  dev.copy2pdf(file=sprintf("%s/plotfit_%d_%f.pdf",folder,data$myseed,data$noise_level))
  print(rstan::traceplot(fit,pars=pars))
  dev.copy2pdf(file=sprintf("%s/traceplot_%d_%f.pdf",folder,data$myseed,data$noise_level))
  # rstan::stan_hist(fit,pars=pars)
  print(bayesplot::mcmc_hist_by_chain(fit,pars=pars))
  Sys.sleep(0.5)
  dev.copy2pdf(file=sprintf("%s/hist_%d_%f.pdf",folder,data$myseed,data$noise_level))
  dev.off()
}
#########################################################################################


####### Home-made plotting ###############################################
plot_analysis <- function(fit, data, folder=".", n_chains=3, chain=0, nsamples=500) {
  species <- data$K
  post <- extract(fit)  # Posterior
  n_samples <- dim(post$r)[1]/n_chains
  if(chain!=0) {
    post$r <- post$r[((chain-1)*n_samples+1):(chain*n_samples),]
    post$a <- post$b[((chain-1)*n_samples+1):(chain*n_samples),,]
  }
  
  x11("",12,8)
  K <- data$K  
  par(mfrow=c(K+1,K))
  par(mar=c(5,5,3,3))
  
  for(i in 1:K) { 
    error <- abs(signif(100*(mean(post$r[,i])-data$params[i])/data$params[i],2))
    error <- ecdf(post$r[,i])(data$params[i])
    n.wrong <- length(which(sign(post$r[,i])!=sign(data$params[i])))
    n <- length(post$r[,i])
    hist(post$r[,i],border='white',col='#47A173',xlab=parse(text=paste0("r[",i,"]")),
         # main=sprintf("True: %.3f. Error: %.0f%%",data$params[i],error),cex.lab=1.3,cex.axis=1.1)
         # main=sprintf("True: %.3f. Quantile: %.3f%%",data$params[i],error),cex.lab=1.3,cex.axis=1.1)
         # main=sprintf("P(wrong sign)=%.0f%% Quantile: %.3f%%",n.wrong/n*100,error),cex.lab=1.3,cex.axis=1.1)
         main=sprintf("Quantile: %.3f%%",error),cex.lab=1.3,cex.axis=1.1)
    abline(v=data$params[i])
  }
  for(i in 1:K) {
    for(j in 1:K) {
      error <- abs(signif(100*(mean(post$b[,i,j])-data$params[K*i+j])/data$params[K*i+j],2))
      error <- ecdf(post$b[,i,j])(data$params[K*i+j])
      n.wrong <- length(which(sign(post$b[,i,j])!=sign(data$params[K*i+j])))
      n <- length(post$r[,i])
      h <- hist(post$b[,i,j],plot=F)
      hist(post$b[,i,j],border='white',
           xlab=parse(text = paste0('b[', i,j,"]")),
           # main=sprintf("True: %.2f. Error: %.0f%%",data$params[K*i+j],error),
           # main=sprintf("True: %.2f. Quantile: %.3f%%",data$params[K*i+j],error),
           main=sprintf("P(wrong sign)=%.0f%% Quantile: %.3f%%",n.wrong/n*100,error),
           col=ifelse(h$mids*data$params[K*i+j]<0,'#8A4A88','#47A173'),
           xlim=c(min(data$params[K*i+j],post$b[,i,j]),max(post$b[,i,j],data$params[K*i+j])),
           probability = TRUE,ylab="Posterior",cex.lab=1.3,cex.axis=1.1)
      abline(v=data$params[K*i+j],lwd=2,col='darkred')
    }
  }
  
  
  dev.copy2pdf(file=sprintf("%s/posteriors_%d_%f.pdf",folder,data$myseed,data$noise_level))
  
  
  ## Plot posterior predictive checks
  # colors <- c("#B2001D", "#ED8B45", "#998EC3", "#007C7C", "#800080", "#505050")
  colors <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF")

  symbols <- rep(c(15,17,19,18),3)
  
  par(mfrow=c(1,1))
  
  simulation <- with(data,
                     ode(y = initial_condition, times = times, func = lotka_volterra, parms = params)
  )
  
  par(mar=c(5,5,3,3))
  matplot(simulation[,1],simulation[,-1],type='p',pch=symbols,cex=1.5,cex.lab=1.7,cex.axis=1.4,
          xlab="Time",ylab="Populations",col=colors,
          main=sprintf("Seed:%d Noise level: %.4f", data$myseed, data$noise_level))
  
  
  nsamples <- nsamples
  # for(val in sample(seq(1:length(post$r[,1])),nsamples)) {
  for(val in seq(length(post$r[,1])-nsamples,length(post$r[,1]))) {
    params_sample <- c(post$r[val,],as.vector(t(post$b[val,,])))
    names(params_sample) <- names(data$params)
    simulation <- with(data,
                       ode(y = initial_condition, times = times, func = lotka_volterra, parms = params_sample)
    )
    matlines(simulation[,1],simulation[,-1],type='l',col=colors,lwd=0.2,lty=3)
    
    
  }
  # Plot median solution
  am <- c(); for(i in 1:K) for(j in 1:K) am <- c(am,median(post$b[,i,j]))
  rm <- apply(post$r,2,median)
  params_sample <- c(rm,am)
  print("##################################################")
  print(params_sample)
  print("##################################################")
  simulation <- with(data,
                     ode(y = initial_condition, times = times, 
                         func = lotka_volterra, parms = params_sample)
  )
  # matlines(simulation[,1],simulation[,-1],type='l',col='black',lwd=2,lty=rep(1,K))
  
  # matpoints(simulation[,1],simulation[,-1],col='black',lwd=2,pch=c(15,17,19))
  lgnd <- "c(expression(x[1])"
  for(idx in 2:species) lgnd <- paste0(lgnd,sprintf(",expression(x[%d])",idx))
  # lgnd <- parse(text=paste0(lgnd,",'Median')"))
  lgnd <- parse(text=paste0(lgnd,")"))
  # legend('topright',legend=c(expression(x[1]),expression(x[2]),expression(x[3])),
  legend('topright',legend=eval(lgnd),
         pch=c(symbols[1:species]),col=c(colors[1:species]),cex=2,bg="white",lty=c(1:species),lwd=3)
         # pch=c(symbols[1:species],-1),col=c(colors[1:species],'black'),cex=2,bg="white",lty=c(1:species,1),lwd=3)
  
  # legend('topright',legend=c(expression(x[1]),expression(x[2]),expression(x[3]),"Median"),
         # pch=c(15,17,19,-1),cex=2,bg="white",lty=1:3,lwd=3,col=c(colors[1:K],'black'))
  # dev.copy2pdf(file=sprintf("stan_lotka/ppc_%d_%f.pdf",data$myseed,data$noise_level))
  dev.copy2pdf(file=sprintf("%s/ppc_%d_%f.pdf",folder,data$myseed,data$noise_level), useDingbats = FALSE)
  
  ### pairs
  pairs.panels(post$b)
  dev.copy2pdf(file=sprintf("%s/panels_%d_%f.pdf",folder,data$myseed,data$noise_level))
  
  dev.off()
  
  return(post)
}
########################################################################################
