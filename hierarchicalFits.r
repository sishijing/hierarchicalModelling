####compare and try different approaches to fitting Gaussian
##-Gaussian hierarchical models
###fully bayesian, post-hoc modelling and so on
###approach 1: fully bayesian
###fully Bayesian fitting with Gibbs sampler, dont fit the
###group variance, only input the data and true values of group
##means, it returns the MSE
##prior on tau^2 is proportional to (tau^2)^(-1/2)
fb.no.var <-function(TSS=10000,len=1000,Data,Mu,sigma){
  I <-  dim(Data)[2] ##number of columns, i.e., groups
  SSpg <- dim(Data)[1] ##sample size per group
  fb.mu <- matrix(0,TSS,I)
  group.mean <- colMeans(Data)
  fb.mu[1,]<- group.mean ##initial values
  fb.gamma <- NULL
  fb.tausq <- NULL
  fb.gamma[1] <- mean(group.mean)
  fb.tausq[1] <- sum((group.mean-fb.gamma[1])^2)/(I-1)
  for(ind in 2:TSS){
    mu.var <- 1/(SSpg/sigma^2+1/fb.tausq[ind-1])
    mu.mean <- (SSpg*group.mean/sigma^2+fb.gamma[ind-1]/fb.tausq[ind-1])*mu.var
    fb.mu[ind,] <- rnorm(I,mu.mean,sqrt(mu.var))
    tausq.scale <- sum((fb.mu[ind,]-mean(fb.mu[ind,]))^2)
    fb.tausq[ind] <- rigamma(1,(I-2)/2,tausq.scale/2)
    fb.gamma[ind] <- rnorm(1,mean(fb.mu[ind,]),sqrt(fb.tausq[ind-1]/I))
  }
  begn <- (len+1)
  fb.mu.est <-colMeans(fb.mu[begn:TSS,])
  fb.mse <- sum((fb.mu.est-Mu)^2)
  fb.gamma.est <-mean(fb.gamma[begn:TSS])
  fb.tausq.est <- mean(fb.tausq[begn:TSS])
  lst <- list("fb.est"=fb.mu.est,"fb.mse"=fb.mse,"fb.mu"=fb.mu[begn:TSS,],
              "fb.gamma"=fb.gamma.est,"fb.tausq"=fb.tausq.est,
              "fb.gamma.sam"=fb.gamma[begn:TSS],"fb.tausq.sam"=fb.tausq[begn:TSS])
  return(lst)
}
###individual fitting with known variance for data
###sigma is the measurement error, standard deviation of
##the observed mass
indv.no.var <-function(SS=10000,Data,sigma){
  I <-  dim(Data)[2] ##number of columns, i.e., groups
  SSpg <- dim(Data)[1] ##sample size per group
  group.mean <- colMeans(Data)
  indv.mu <- matrix(0,SS,I) ##individual fitting
  for(i in 1:I){indv.mu[,i] <- rnorm(SS,group.mean[i],sigma/sqrt(SSpg))}
  lst <- list("indv.mu"=indv.mu)
  return(lst)
}
###post-hoc hierarchical analysis
## use the individual analysis result as data

###Empirical Bayes fitting
eb.no.var <-function(Data,Mu,sigma){
  I <-  dim(Data)[2] ##number of columns, i.e., groups
  SSpg <- dim(Data)[1] ##sample size per group
  group.mean <- colMeans(Data)
  eb.gamma <- mean(Data) #optimal gamma
  eb.tausq <- max(var(c(Data))-sigma^2,10^(-15)) 
  eb.mu.var <- 1/(1/sigma^2+1/eb.tausq)
  eb.mu.est <- (group.mean/sigma^2+eb.gamma/eb.tausq)*eb.mu.var
  eb.mse <- sum((eb.mu.est-Mu)^2)
  lst <- list("eb.est"=eb.mu.est,"eb.mse"=eb.mse,"eb.var"=eb.mu.var,
              "eb.gamma"=eb.gamma,"eb.tausq"=eb.tausq)
  return(lst)
}
####individual fitting approx. with approximation
indvFit.no.var <-function(SS= 10000,Data,Mu,sigma, len =1000){
  I <-  dim(Data)[2] ##number of columns, i.e., groups
  SSpg <- dim(Data)[1] ##sample size per group
  group.mean <- colMeans(Data)
  indvFit.mu <- matrix(0,SS,I) ##individual fitting mu
  indvFit.gamma <- NULL
  indvFit.tausq <- NULL
  for(i in 1:SS){
    indvFit.mu[i,] <-rnorm(I,group.mean,sigma/sqrt(SSpg))
    indvFit.tausq[i] <- rigamma(1,I/2-1,(I-1)*var(indvFit.mu[i,])/2)
    indvFit.gamma[i] <- rnorm(1,mean(group.mean),sigma/sqrt(I))
  }
  begn <- SS-len+1
  indvFit.mu.est <- colMeans(indvFit.mu[begn:SS,])
  indvFit.mse <- sum((indvFit.mu.est-Mu)^2)
  indvFit.gamma.est <- mean(indvFit.gamma[begn:SS])
  indvFit.tausq.est <- mean(indvFit.tausq[begn:SS])
  lst <- list("indvFit.est"=indvFit.mu.est,"indvFit.mu"=indvFit.mu,
              "indvFit.mse"=indvFit.mse,
              "indvFit.gamma"=indvFit.gamma.est,"indvFit.tausq"=indvFit.tausq.est,
              "indvFit.gamma.sam"=indvFit.gamma,"indvFit.sam"=indvFit.tausq)
  return(lst)
}

####The 6th method--approximation to fully bayesian
##try to approximate fully Bayes from individual fitting
##use metropolis-hasting algorithm to correct the Gibbs sampler
app.no.var <-function(SS=10000,len=5000,indv.sam,Data,Mu,sigma){
  I <-  dim(Data)[2] ##number of columns, i.e., groups
  SSpg <- dim(Data)[1] ##sample size per group
  TSS <- dim(indv.sam)[1] ##total sample size
  group.mean <- colMeans(Data)
  app.mu <- matrix(0,SS,I) ##sample for mu in this method
  app.mu[1,] <- indv.sam[1,]
  app.gamma <- NULL
  app.tausq <- NULL
  app.tausq.scale <-sum((app.mu[1,]-mean(app.mu[1,]))^2)
  app.tausq[1] <- rigamma(1,(I-2)/2,app.tausq.scale/2)
  app.gamma[1] <- rnorm(1,mean(app.mu[1,]),sqrt(app.tausq[1]/I)) #optimal gamma
  app.acp <- matrix(0,SS,I) ##acceptance matrix
    app.acp[1,] <- rep(1,I)
  for(aind in 2:SS){
    ran.num <- sample.int(TSS,I,replace=T) ##random numbers
    app.prop <- NULL
    for(i in 1:I){
      app.prop[i] <- indv.sam[ran.num[i],i] ##proposal
    }
    app.prob <- dnorm(app.prop,app.gamma[aind-1],sqrt(app.tausq[aind-1]))/
      dnorm(app.mu[aind-1,],app.gamma[aind-1],sqrt(app.tausq[aind-1]))
    app.unif <- runif(I)
    temp.ind<-I(app.unif<=app.prob) ##temporary indication
    ind.dist <- ran.num-app.acp[aind-1,]
    app.acp[aind,] <- app.acp[aind-1,] + ind.dist*temp.ind
    app.delta <-app.prop - app.mu[aind-1,] ##difference between proposal and previous
    app.mu[aind,] <- app.mu[aind-1,]+app.delta*temp.ind
    app.tausq.scale <- sum((app.mu[aind,]-mean(app.mu[aind,]))^2)
    app.tausq[aind] <- rigamma(1,(I-2)/2,app.tausq.scale/2)
    app.gamma[aind] <- rnorm(1,mean(app.mu[aind,]),sqrt(app.tausq[aind]/I))
  }
  begn <- len+1
  app.mu.est <- colMeans(app.mu[begn:SS,]) ##estimator of mu
  app.mse <- sum((app.mu.est-Mu)^2)
  app.gamma.est <- mean(app.gamma[begn:SS])
  app.tausq.est <- mean(app.tausq[begn:SS])
  lst <- list("app.est"=app.mu.est,"app.mse"=app.mse,"app.mu.sam"=app.mu[begn:SS,],
              "acp.ind"=app.acp,"app.gamma"=app.gamma.est,"app.tausq"=app.tausq.est,
              "app.gamma.sam"=app.gamma[begn:SS],"app.tausq.sam"=app.tausq[begn:SS])
  return(lst)
  }

####partial data proposal algorithm
##use metropolis-hasting algorithm to correct the Gibbs sampler
###perc parameter set the percentage of data used in the proposal
pd.no.var <-function(SS=10000,TSS= 10000,len=1000,perc=.5,Data,Mu,sigma){
  I <-  dim(Data)[2] ##number of columns, i.e., groups
  SSpg <- dim(Data)[1] ##sample size per group
  pSize<-floor(perc*SSpg) ##size of sample used in the proposal rule
  ##total sample size
  prop.sd <- sigma/sqrt(pSize)
  prop.mn <- colMeans(Data[1:pSize,]) ##mean vector of proposal dist.
  like.sd <- sigma/sqrt(SSpg) ##likelihood sd
  like.mn <- colMeans(Data) ##likelihood mean vector
  indv.sam <- matrix(0,TSS,I) ##individual fittings
  for(i in 1:I){indv.sam[,i] <- rnorm(TSS,prop.mn[i],prop.sd)}
  group.mean <- like.mn
  pd.mu <- matrix(0,SS,I) ##sample for mu in this method
  pd.mu[1,] <- indv.sam[1,]
  pd.gamma <- NULL
  pd.tausq <- NULL
  pd.tausq.scale <-sum((pd.mu[1,]-mean(pd.mu[1,]))^2)
  pd.tausq[1] <- rigamma(1,(I-2)/2,pd.tausq.scale/2)
  pd.gamma[1] <- rnorm(1,mean(pd.mu[1,]),sqrt(pd.tausq[1]/I)) #optimal gamma
  pd.acp <- matrix(0,SS,I) ##acceptance matrix
  pd.acp[1,] <- rep(1,I)
  for(aind in 2:SS){
    ran.num <- sample.int(TSS,I,replace=T) ##random numbers
    pd.prop <- NULL
    for(i in 1:I){
      pd.prop[i] <- indv.sam[ran.num[i],i] ##proposal
    }
    pd.prob <- dnorm(pd.prop,pd.gamma[aind-1],sqrt(pd.tausq[aind-1]))*
      dnorm(pd.prop,like.mn,like.sd)*dnorm(pd.mu[aind-1,],prop.mn,prop.sd)/
      dnorm(pd.mu[aind-1,],pd.gamma[aind-1],sqrt(pd.tausq[aind-1]))/
      dnorm(pd.mu[aind-1,],like.mn,like.sd)/dnorm(pd.prop,prop.mn,prop.sd)
    pd.unif <- runif(I)
    temp.ind<-I(pd.unif<=pd.prob) ##temporary indication
    ind.dist <- ran.num-pd.acp[aind-1,]
    pd.acp[aind,] <- pd.acp[aind-1,] + ind.dist*temp.ind
    pd.delta <-pd.prop - pd.mu[aind-1,] ##difference between proposal and previous
    pd.mu[aind,] <- pd.mu[aind-1,]+pd.delta*temp.ind
    pd.tausq.scale <- sum((pd.mu[aind,]-mean(pd.mu[aind,]))^2)
    pd.tausq[aind] <- rigamma(1,(I-2)/2,pd.tausq.scale/2)
    pd.gamma[aind] <- rnorm(1,mean(pd.mu[aind,]),sqrt(pd.tausq[aind]/I))
  }
  begn <- len+1
  pd.mu.est <- colMeans(pd.mu[begn:SS,]) ##estimator of mu
  pd.mse <- sum((pd.mu.est-Mu)^2)
  pd.gamma.est <- mean(pd.gamma[begn:SS])
  pd.tausq.est <- mean(pd.tausq[begn:SS])
  lst <- list("pd.est"=pd.mu.est,"pd.mse"=pd.mse,"pd.mu.sam"=pd.mu[begn:SS,],
              "acp.ind"=pd.acp[begn:SS,],"pd.gamma"=pd.gamma.est,"pd.tausq"=pd.tausq.est,
              "pd.gamma.sam"=pd.gamma[begn:SS],"pd.tausq.sam"=pd.tausq[begn:SS])
  return(lst)
}

##Metropolis-Hastings algorithm within Gibbs sampler
##use shifted individual as the proposal to draw the 
###simulations
strt <- Sys.time()
library(pscl) ##library packages with rigamma in it
##true values
gamma<-0.5
tau<-.1
lvar.ratio <- 1 #sigma^2/tau^2
llvar.ratio <- length(lvar.ratio)
I<-10 ##number of groups
for(ratio in 1:length(lvar.ratio)){
  SSpg <- 60 ##sample size per group
  sigma <- 1000
  #set.seed(12345) ##set random seed
  reps <- 1 ##number of replicates
  for(rep in 1:reps){
    Mu<-rnorm(I,gamma,tau) ##group means
    Data<-matrix(0,SSpg,I) ##each column for one group
    for(j in 1:I){
      Data[,j]<-rnorm(SSpg,Mu[j],sigma) ##data
    }
    ###different ways to fit the hierarchical model
    fb.fit <- fb.no.var(TSS=100000,len=1000,Data,Mu,sigma)
    pd.fit<-pd.no.var(SS=100000,TSS=50000,len=1000,perc=.5,Data,Mu,sigma)
    indvFit.fit <-indvFit.no.var(SS=50000,Data,Mu,sigma,len=1000)
    indv.sam <- indvFit.fit$indvFit.mu
    #post.fit <-post.no.var(TSS=10000,len=1000,Data,Mu,sigma)
    #post.temp[rep,1] <- post.fit$post.mse
    #post.temp[rep,3] <- (post.fit$post.tausq-tau^2)^2
    #post.temp[rep,2] <- (post.fit$post.gamma-gamma)^2
  #  eb.fit <-eb.no.var(Data,Mu,sigma)
   # eb.temp[rep,1] <- eb.fit$eb.mse
    #eb.temp[rep,2] <- (eb.fit$eb.gamma-gamma)^2
    #eb.temp[rep,3] <- (eb.fit$eb.tausq-tau^2)^2
    app.fit <-app.no.var(SS=100000,len=1000,indv.sam,Data,Mu,sigma)
    #app.temp[rep,1] <- app.fit$app.mse
    #app.temp[rep,2] <- (app.fit$app.gamma-gamma)^2
    #app.temp[rep,3] <- (app.fit$app.tausq-tau^2)^2
    #scl.fit<-scale.no.var(len=1000,indv.sam,Data,Mu,sigma)
    #app.acp.temp[rep,] <- app.fit$acprat
    #filter.fit <-filter.no.var(TSS=10000,len=1000,indv.sam,Data,Mu,sigma)
    #filter.temp[rep,1] <- filter.fit$filter.mse
    #filter.temp[rep,2] <- (filter.fit$filter.gamma-gamma)^2
    #filter.temp[rep,3] <- (filter.fit$filter.tausq-tau^2)^2
  }
  fb.res[ratio,] <- colMeans(fb.temp)
  indv.res[ratio,] <- colMeans(indv.temp)
  indvFit.res[ratio,] <- colMeans(indvFit.temp)
  eb.res[ratio,] <- colMeans(eb.temp)
  app.res[ratio,] <- colMeans(app.temp)
  app.acprat[ratio,] <- colMeans(app.acp.temp)
  filter.res[ratio,] <- colMeans(filter.temp)
  ###take square roots
  fb.sqrt <-sqrt(fb.res)
  indv.sqrt <- sqrt(indv.res)
  indvFit.sqrt <- sqrt(indvFit.res)
  eb.sqrt <- sqrt(eb.res)
  app.sqrt <- sqrt(app.res)
  filter.sqrt <- sqrt(filter.res)
  time.period <- Sys.time()-strt
  save(time.period,fb.sqrt,indv.sqrt,indvFit.sqrt,eb.sqrt,
       app.sqrt,app.acprat,filter.sqrt,file="newfitHier4Feb16.Rdata")
  }
####plots to compare the partial data proposal and fully bayesian
pdf(file = "partialDataT01S1000.pdf",width = 6,height = 7)
par(mfrow=c(2,2),mar=c(3,3,1,1),oma=c(0,0,3,0),tck=0.02,
    mgp=c(1.4,0.2,0),cex=0.6)
qqplot(pd.fit$pd.gamma.sam,fb.fit$fb.gamma.sam,xlab="gamma from partial data Algo.",ylab="gamma from FB")
 abline(0,1,col="red")
 qqplot(sqrt(pd.fit$pd.tausq.sam),sqrt(fb.fit$fb.tausq.sam),
        xlab="tau from partial data Algo.",ylab="tau from FB")
 abline(0,1,col="red")
 qqplot(app.fit$app.gamma.sam,fb.fit$fb.gamma.sam,xlab="gamma from MH correction",ylab="gamma from FB")
 abline(0,1,col="red")
 qqplot(sqrt(app.fit$app.tausq.sam),sqrt(fb.fit$fb.tausq.sam),
        xlab="tau from MH correction",ylab="tau from FB")
 abline(0,1,col="red")
 mtext(paste(paste("tau= ",tau),paste(",sigma=",sigma),sep=" "), outer = TRUE, cex = 1.5)
 dev.off()
###compare different method
fbmu <-fb.fit$fb.mu
mhmu <-app.fit$app.mu
shmu <-shift.fit$shift.mu
nshmu <-nshift.fit$shift.mu
 qqplot(fbmu[,9],mhmu[,9])
 abline(0,1,col="red")
 qqplot(fbmu[,9],shmu[,9])
 abline(0,1,col="red")
 qqplot(fbmu[,9],nshmu[,9])
 abline(0,1,col="red")
 ###Try shifted MH correction
 NoG <- 3 ##number of groups 
 gamma0<-0.2
 tau0<-.2
 Mu <- rnorm(NoG,gamma0,tau0)
 sigma=2
 Data <- matrix(rnorm(NoG,Mu,sigma),1,NoG)
 indvFit.fit <-indvFit.no.var(Data,Mu,sigma,len=1000)
 indv.sam <- indvFit.fit$indvFit.mu
 fb.fit <- fb.no.var(TSS=10000,len=1000,Data,Mu,sigma)
 shift.fit <-shift.no.var(TSS=10000,len=1000,Data,Mu,sigma)
 app.fit <-app.no.var(TSS=10000,len=1000,Data,Mu,sigma)
 nshift.fit <-nshift.no.var(TSS=10000,len=1000,indv.sam,Data,Mu,sigma)
 filter.fit<-filter.no.var(TSS=10000,len=1000,indv.sam,Data,Mu,sigma)
 nshift.acprat <- nshift.fit$acprat
 shift.acprat <- shift.fit$acprat
 app.acprat <- app.fit$acprat
 nsmu <-nshift.fit$shift.mu
 smu <- shift.fit$shift.mu
 amu <- app.fit$app.mu
 fmu <- filter.fit$filter.mu
 fbmu <- fb.fit$fb.mu
 nsgam <- nshift.fit$shift.gamma.sam
 sgam <- shift.fit$shift.gamma.sam
 agam <- app.fit$app.gamma.sam
 fgam <- filter.fit$filter.gamma.sam
 fbgam <- fb.fit$fb.gamma.sam
 nstau <- nshift.fit$shift.tausq.sam
 stau <- shift.fit$shift.tausq.sam
 atau <- app.fit$app.tausq.sam
 ftau <- filter.fit$filter.tausq.sam
 fbtau <- fb.fit$fb.tausq.sam
 
 i =2 ###i-th mu
 qqplot(nsmu[,i],fbmu[,i])
 abline(0,1,col='red')
 qqplot(smu[,i],fbmu[,i])
 abline(0,1,col='red')
 qqplot(fmu[,i],fbmu[,i])
 abline(0,1,col='red')
 qqplot(amu[,i],fbmu[,i])
 abline(0,1,col='red')
 qqplot(nsgam,fbgam)
 abline(0,1,col='red')
 qqplot(sgam,fbgam)
 abline(0,1,col='red')
 qqplot(fgam,fbgam)
 abline(0,1,col='red')
 qqplot(agam,fbgam)
 abline(0,1,col='red')
 qqplot(nstau,fbtau)
 abline(0,1,col='red')
 qqplot(stau,fbtau)
 abline(0,1,col='red')
 qqplot(ftau,fbtau)
 abline(0,1,col='red')
 qqplot(atau,fbtau)
 abline(0,1,col='red')
 ####comparison with fully bayesian results
 app.gam <- app.fit$app.gamma.sam
 app.tausq <- app.fit$app.tausq.sam
 app.mu <- app.fit$app.mu
 scl.mu <-scl.fit$scl.mu
 scl.gam <- scl.fit$scl.gamma.sam
 scl.tausq <- scl.fit$scl.tausq.sam
 fb.mu <- fb.fit$fb.mu
 fb.gam <- fb.fit$fb.gamma.sam
 fb.tausq <- fb.fit$fb.tausq.sam
 sir.gam <- filter.fit$filter.gamma.sam
 sir.tausq <- filter.fit$filter.tausq.sam
 sir.mu <- filter.fit$filter.mu
 qqplot(scl.gam,fb.gam)
 abline(0,1,col="red")
 qqplot(app.gam,fb.gam)
 abline(0,1,col="red")
 qqplot(app.tausq,fb.tausq)
 abline(0,1,col="red")
 qqplot(scl.tausq,fb.tausq)
 abline(0,1,col="red")
 qqplot(sir.tausq,fb.tausq)
 abline(0,1,col="red")
 
 qqplot(scl.mu[,1],fb.mu[,1])
 abline(0,1,col="red")
 qqplot(scl.mu[,1],fb.mu[,1])
 abline(0,1,col="red")
 qqplot(scl.mu[,2],fb.mu[,2])
 abline(0,1,col="red")
 qqplot(scl.mu[,2],fb.mu[,2])
 abline(0,1,col="red")
 qqplot(scl.mu[,3],fb.mu[,3])
 abline(0,1,col="red")
 qqplot(app.mu[,5],fb.mu[,5])
 abline(0,1,col="red")
 qqplot(app.mu[,2],fb.mu[,2])
 abline(0,1,col="red")
 qqplot(app.mu[,3],fb.mu[,3])
 abline(0,1,col="red")
 qqplot(sir.mu[,2],fb.mu[,2])
 abline(0,1,col="red")
 qqplot(sir.mu[,3],fb.mu[,3])
 abline(0,1,col="red")
 qqplot(sir.mu[,4],fb.mu[,4])
 abline(0,1,col="red")
 qqplot(sir.mu[,5],fb.mu[,5])
 abline(0,1,col="red")
 ###fully bayesian marginal posterior for mu and its normality
 fbm1 <- fb.mu[,3]
 fbm1.mn <- mean(fbm1)
 fbm1.var <- var(fbm1)
 ran1 <- range(fbm1)
 xran1 <- seq(ran1[1],ran1[2],by=0.05)
 hist(fbm1,freq=F,breaks=50)
 lines(xran1,dnorm(xran1,fbm1.mn,sqrt(fbm1.var)),col="red")
 