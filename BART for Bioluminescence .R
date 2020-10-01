library(rpart)
library(MASS)
library(mvtnorm)
library(BayesTree)
library(gam)

tmpdir <- tempdir()
url <- "http://www.highstat.com/Book2/ZuurDataMixedModelling.zip"
file <- basename(url)
download.file(url, file)
unzip(file, exdir = tmpdir)
file <- paste(tmpdir, "/ISIT.txt", sep = "")
ISIT <- read.table(file, header = TRUE)
Sources16 <- ISIT$Sources[ISIT$Station == 16]
Depth16 <- ISIT$SampleDepth[ISIT$Station == 16]

Sources16 <- Sources16[order(Depth16)]
Depth16 <- sort(Depth16)
plot(Depth16, Sources16, las = 1, ylim = c(0, 65), col = rgb(0, 0, 0, 0.25),
     pch = 19)

y <- Sources16
x <- Depth16

#Using Bart Package

m1 <- gam(y~s(x))

set.seed(99)
bartFit <- bart(x,y,ndpost=10000,x.test=seq(0,5000,by=10),ntree=1,power=6
                ,base=100)
plot(bartFit)
plot(x,bartFit$yhat.train.mean)

par(mfcol=c(1,1)) # plot layout
par(mar=c(5.1,4.7,2.1,2.1),oma=c(5.1,3.1,2.1,2.1)/1000) 
par(pch=19) # plotting character
cols <- c("#F8766D","#00BFC4","#7CAE00","#C77CFF") # color palette

plot(x,y,las=1,col=rgb(0,0,0,0.25),xlab=expression(italic(x[i])),
     ylab=expression(italic(y[i])),cex.lab=1.7,bty="n" )
points(seq(0,5000,by=10),bartFit$yhat.test.mean,typ="l",lwd=3)

#Confdence Interval
newx <- seq(min(x), max(x), length.out=100)
preds <- predict(bartFit, newdata = data.frame(x=newx), 
                 interval = 'confidence')


lwr.CI <- apply(bartFit$yhat.test,2,FUN=quantile,prob=c(0.025))
upper.CI <- apply(bartFit$yhat.test,2,FUN=quantile,prob=c(0.975))
polygon(c(seq(0,5000,by=10),rev(seq(0,5000,by=10))),c(lwr.CI,rev(upper.CI)),col=rgb(0.5,0.5,0.5,0.3),border=NA)
legend(x = 0, y = 10, cex = 1.3, legend = c(expression("E("*y[new]*"|"*y[obs]*")"),"95% CI"), bty = "n",
       lty = 1, lwd = 2, col = c("black",rgb(0.5,0.5,0.5,0.5)))


#Code from Scratch

lhs <- function(x,phi)ifelse(x<phi,phi-x,0)
rhs <- function(x,phi)ifelse(x<phi,0,x-phi)

m.draws <- 2000  #Number of MCMC samples to draw
samples <- as.data.frame(matrix(,m.draws+1,5))
colnames(samples) <- c("beta0","beta1","beta2","sigma2","phi")
samples[1,] <- c(rep(NA,3),1,2500)  #Starting values for sigma2 & phi
ifsplitsamples<-as.data.frame(matrix(,m.draws+1,5))

accept <- as.data.frame(matrix(0,m.draws+1,1))


sigma2.beta <- 100 #Prior variance for beta
q <- 2.01; r <- 1 #Inverse gamma prior with E() = r/(q-1) and Var()= r^2/((q-1)^2(q-2))
a <- 0; b <- 5000 #Uniform prior for phi

phi.tune <- 100
alpha<-0.8 #Probability of split
psplit<-alpha

set.seed(4403)

for(k in 1:m.draws){
  
  sigma2.e <- samples[k,4]
  phi <- samples[k,5]
  
  if(rbinom(1,1,psplit)==1){ #Determining whether to split or not
    
    X <- model.matrix(~lhs(x,phi)+rhs(x,phi))
    
    #Sample beta
    C <- solve(t(X)%*%X + diag(1/sigma2.beta,dim(X)[2]))
    c <- t(X)%*%y
    beta <-  mvrnorm(1,C%*%c,sigma2.e*C)
    
    #Sample sigma2.e
    sigma2.e <- 1/rgamma(1,q+length(y)/2,r+t(y-X%*%beta)%*%(y-X%*%beta)/2)
    
    #Sample phi
    phi.star <- rnorm(1,phi,phi.tune)
    #phi.star <- sample(x,1)
    if(phi.star > a & phi.star < b){
    #if(rbinom(1,1,psplit)==1){
      X.star <- model.matrix(~lhs(x,phi.star)+rhs(x,phi.star))
      mh1 <- dmvnorm(y,X.star%*%beta,sigma2.e*diag(1,length(y)),log=TRUE) +dnorm(beta, 0, sigma2.beta^0.5,log=T)+ dunif(phi.star,a,b,log=TRUE)
      mh2 <- dmvnorm(y,X%*%beta,sigma2.e*diag(1,length(y)),log=TRUE) + dnorm(beta, 0, sigma2.beta^0.5,log=T)+ dunif(phi,a,b,log=TRUE)
      R <- min(1,exp(mh1 - mh2))
      if (R > runif(1)){phi <- phi.star;accept[k+1,1] <- 1}}
    
    ifsplitsamples[k+1,] <- c(beta,sigma2.e,phi)
    
  }
  
  #Save samples
  samples[k+1,] <- c(beta,sigma2.e,phi)
  if (k%%1000 == 0)print(k)
}

mean(accept[-1,])

par(mfrow=c(5,1))
plot(samples[,1],typ="l",xlab="k",ylab=expression(beta[0]))
plot(samples[,2],typ="l",xlab="k",ylab=expression(beta[1]))
plot(samples[,3],typ="l",xlab="k",ylab=expression(beta[2]))
plot(samples[,4],typ="l",xlab="k",ylab=expression(sigma2[epsilon]))
plot(samples[,5],typ="l",xlab="k",ylab=expression(phi))

par(mfrow=c(5,1))
plot(ifsplitsamples[,1],typ="l",xlab="k",ylab=expression(beta[0]))
plot(ifsplitsamples[,2],typ="l",xlab="k",ylab=expression(beta[1]))
plot(ifsplitsamples[,3],typ="l",xlab="k",ylab=expression(beta[2]))
plot(ifsplitsamples[,4],typ="l",xlab="k",ylab=expression(sigma2[epsilon]))
plot(ifsplitsamples[,5],typ="l",xlab="k",ylab=expression(phi))


par(mfrow=c(1,1))
hist(ifsplitsamples[,5],col="grey",main="Posterior distribution of breakpoint")

