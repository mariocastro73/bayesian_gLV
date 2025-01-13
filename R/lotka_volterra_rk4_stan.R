library(deSolve)
library(rstan)
library(ggplot2)
library(psych)


rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

myseed <- sample.int(1000,1)
set.seed(myseed)

var_names <- c("r1",  "r2",  "r3",  "a11", "a12", "a13", "a21", "a22", "a23", "a31", "a32", "a33")
flag <- FALSE

while(!flag) {
  r0 <- rexp(3,.5)
  a0 <- matrix(rnorm(9),3,3)
  for(i in 1:3) a0[i,i]  <- -rexp(1,.5)
  params <- c(r0,as.vector(t(a0)))
  names(params) <- var_names
  
  xss <- solve(a0,-r0)
  d0 <- diag(xss)
  print(eigenvalues <- Re(eigen(d0 %*% a0)$values))
  flag <- all(eigenvalues <=0)
}

# Lotka-Volterra equations
lotka_volterra <- function(t, initial_condition, parameters) {
  with(as.list(c(initial_condition, parameters)), {
    dX1 <- r1*X1 + a11*X1*X1 + a12*X1*X2 + a13*X1*X3
    dX2 <- r2*X2 + a21*X1*X2 + a22*X2*X2 + a23*X2*X3
    dX3 <- r3*X3 + a31*X1*X3 + a32*X2*X3 + a33*X3*X3
    return(list(c(dX1, dX2, dX3)))
  })
}

# Initial initial_condition
initial_condition <- c(X1 = rexp(1,10), X2 = rexp(1,10), X3 = rexp(1,10))

# Time points
tmax <- 10
tmax <- -20/min(eigenvalues)
times <- seq(0,tmax, length=25)

# Solve differential equations

out <- ode(y = initial_condition, times = times, func = lotka_volterra, parms = params)
noise_level <- 0.01
for(i in 2:4) {
  out[,i] <- rlnorm(length(out[,i]),log(out[,i]),noise_level)
}
# x11("",12,8)
par(mfrow=c(1,1))
par(mar=c(5,5,3,3))
matplot(out[,1],out[,-1],type='b',pch=c(15,17,19),cex=1.5,cex.lab=1.7,cex.axis=1.4,
        xlab="Time",ylab="Populations",lwd=2,lty=1:3)


legend('topright',legend=c(expression(x[1]),expression(x[2]),expression(x[3])),
       pch=c(15,17,19),col=1:3,cex=2,bg="white",lty=1:3,lwd=3)


dev.copy2pdf(file=sprintf("stan_lotka/ppd_%d_%f.pdf",myseed,noise_level))

print(myseed)
print(eigenvalues)



# Prepare data for Stan model
N = nrow(out)
Y = out[,-1]
T = out[,1]
stan_data <- list(N=N,Y=Y,T=T)
stan_data <- list(N=N,Y=Y,T=T,K=3)

recover <- list(stan_data=stan_data,initial_condition=initial_condition,
                eigenvalues=eigenvalues,params=params,a0=a0,r0=r0,out=out)

saveRDS(recover,sprintf("recover_%d.rda",myseed))
# list2env(readRDS("recover_473.rda"), envir = .GlobalEnv)


####################################################################################
# Fit model
n_iter <- 4000
warmup <- 3000
# fit <- stan(file = "lotka_volterra_rk.stan", data = stan_data, iter = n_iter, warmup = warmup,chains = 3)
fit <- stan(file = "lotka_volterra_rk_K_dimensions.stan", data = stan_data, iter = n_iter, warmup = warmup,chains = 3)
dev.off()

x11("",12,8)
# Print the results
print(fit)
pars <- c(); for(i in 1:3) for(j in 1:3) pars <- c(pars,sprintf("a[%d,%d]",i,j))
plot(fit,pars=pars)

post <- extract(fit)

####################################################################################
# x11("",12,8)
par(mfrow=c(4,3))
par(mar=c(5,5,3,3))

for(i in 1:3) { 
  error <- abs(signif(100*(mean(post$r[,i])-params[i])/params[i],2))
  hist(post$r[,i],border='white',col='#47A173',xlab=parse(text=paste0("r[",i,"]")),
       main=sprintf("True: %.3f. Error: %.0f%%",params[i],error),cex.lab=1.3,cex.axis=1.1)
  abline(v=params[i])
}
for(i in 1:3) {
  for(j in 1:3) {
    error <- abs(signif(100*(mean(post$a[,i,j])-params[3*i+j])/params[3*i+j],2))
    h <- hist(post$a[,i,j],plot=F)
    hist(post$a[,i,j],border='white',
         xlab=parse(text = paste0('a[', i,j,"]")),
         main=sprintf("True: %.2f. Error: %.0f%%",params[3*i+j],error),
         col=ifelse(h$mids*params[3*i+j]<0,'#8A4A88','#47A173'),
         xlim=c(min(params[3*i+j],post$a[,i,j]),max(post$a[,i,j],params[3*i+j])),
         probability = TRUE,ylab="Posterior",cex.lab=1.3,cex.axis=1.1)
    abline(v=params[3*i+j],lwd=2,col='darkred')
  }
}


dev.copy2pdf(file=sprintf("stan_lotka/posteriors_%d_%f.pdf",myseed,noise_level))

pairs.panels(post$a)
dev.copy2pdf(file=sprintf("stan_lotka/panels_%d_%f.pdf",myseed,noise_level))

## Plot posterior predictive checks
par(mfrow=c(1,1))

out <- ode(y = initial_condition, times = times, func = lotka_volterra, parms = params)

par(mar=c(5,5,3,3))
matplot(out[,1],out[,-1],type='p',pch=c(15,17,19),cex=1.5,cex.lab=1.7,cex.axis=1.4,
        xlab="Time",ylab="Populations",main=sprintf("Seed:%d Noise level: %.4f",myseed,noise_level))


nsamples <- 500
for(val in sample(seq(1:length(post$r[,1])),nsamples)) {
  params_sample <- c(post$r[val,],as.vector(t(post$a[val,,])))
  names(params_sample) <- names(params)
  out <- ode(y = initial_condition, times = times, func = lotka_volterra, parms = params_sample)
  
  matlines(out[,1],out[,-1],type='l')
  
}

legend('topright',legend=c(expression(x[1]),expression(x[2]),expression(x[3])),
       pch=c(15,17,19),col=1:3,cex=2,bg="white",lty=1:3,lwd=3)
dev.copy2pdf(file=sprintf("stan_lotka/ppd_%d_%f.pdf",myseed,noise_level))

dev.off()
