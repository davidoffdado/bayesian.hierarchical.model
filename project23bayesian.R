library(lme4)
library(rstudioapi)
library(ggplot2)
library(dplyr)
library(broom)
library(gridExtra)
library(SimDesign)
library(coda)
library(sbgcop)

data(sleepstudy)
dati=sleepstudy
set.seed(1)



#' preparation of data
N=NULL
Y=X=list()
m=length(table(dati$Subject))
id=levels(dati$Subject)
for (j in 1:m) {
  Y[[j]]=dati[dati[,3]==id[j], 1]
  N[j]=sum(dati[,3]==id[j])
  xj=dati[dati[,3]==id[j], 2]
  xj2=xj^2
  X[[j]]=cbind(rep(1,N[j]), xj, xj2)
}


#' fitting a pooled ols model 
ols.pooled=lm(Reaction~Days, dati) 
pooled.pred=predict(ols.pooled, Days=dati$Days)
subject.col <- c("#FF0000", "#00FF00", "#0000FF", "#FFFF00", "#FF00FF", "#00FFFF", "#000000", "#FFFFFF", "#FFA500", "#800080", "#008000", "#000080", "#FFC0CB", "#800000", "#808080", "#FFFFF0", "#F0FFF0", "#F0FFFF")
plot(dati$Days, dati$Reaction, type = "n", xlab="Days", ylab="Reaction")#, col=subject.col[dati$Subject])
lines(dati$Days, pooled.pred)
summary(ols.pooled)

#' plotting variances
ggplot(ols.pooled, aes(x=.fitted, y=.resid))+
  geom_point()+
  geom_abline(intercept=0, slope=0, color="red", linewidth=1)

#' fitting a linear OLS model for each subject
#' fitting a quadratic OLS model for each subject
#' obtaining BETA.LS[1]=intercept, BETA.LS[2]=slope and S2.LS=varinance for each subject
S2.LS=BETA.LS=NULL
fit.ols=R2.LS=list()
days.cont=seq(min(xj), max(xj), by=0.1)
react.pred.ols=list()
for(j in 1:m) { #'togliendo l'intercetta le linee partono da 2 circa
  fit.ols[[j]]<-lm(Y[[j]]~xj, as.data.frame(X[[j]])) #regressione lineare per ogni scuola usando xj
  BETA.LS<-rbind(BETA.LS,c(fit.ols[[j]]$coef)) #la regressione lineare produce due coefficienti, uno relativo all'intera scuola quindi il punteggio medio della scuola (intercetta), uno relativo al coefficiente dello status socio-economico di quella scuola (slope)
  S2.LS<-c(S2.LS, summary(fit.ols[[j]])$sigma^2) #varianza per ogni singola scuola
  R2.LS[j]=summary(fit.ols[[j]])$r.squared
  react.pred.ols[[j]]=predict(fit.ols[[j]], list(xj=days.cont))
} 

subject.col <- c("#FF0000", "#00FF00", "#0000FF", "#FFFF00", "#FF00FF", "#00FFFF", "#000000", "#FFFFFF", "#FFA500", "#800080", "#008000", "#000080", "#FFC0CB", "#800000", "#808080", "#FFFFF0", "#F0FFF0", "#F0FFFF")
plot(dati$Days, dati$Reaction, type="n", xlab = "Days", ylab = "Reaction")  #,col=subject.col[dati$Subject])
for(j in 1:m) {lines(days.cont, react.pred.ols[[j]])}#, col=subject.col[j])}
lines(dati$Days, pooled.pred, col="red", lwd=2)


#' fitting a pooled quadratic ols model
x=c(0,seq(0:9))
dati$Days2=dati$Days^2 
ols.pooled=lm(Reaction~Days+Days2, dati) 
pooled.pred=predict(ols.pooled, list(Days=x, Days2=x^2))
subject.col <- c("#FF0000", "#00FF00", "#0000FF", "#FFFF00", "#FF00FF", "#00FFFF", "#000000", "#FFFFFF", "#FFA500", "#800080", "#008000", "#000080", "#FFC0CB", "#800000", "#808080", "#FFFFF0", "#F0FFF0", "#F0FFFF")
plot(dati$Days, dati$Reaction, type = "n", xlab="Days", ylab="Reaction")#, col=subject.col[dati$Subject])
lines(x, pooled.pred)
summary(ols.pooled)


#' fitting a quadratic OLS model for each subject
#' obtaining BETA.LS[1]=intercept, BETA.LS[2]=slope and S2.LS=variance for each subject
S2.LS=BETA.LS=NULL
fit.ols=R2.LS=list()
days.cont=seq(min(xj), max(xj), by=0.1)
react.pred.ols=list()
for(j in 1:m) { #'togliendo l'intercetta le linee partono da 2 circa
  fit.ols[[j]]<-lm(Y[[j]]~xj+xj2, as.data.frame(X[[j]])) #regressione lineare per ogni scuola usando xj
  BETA.LS<-rbind(BETA.LS,c(fit.ols[[j]]$coef)) #la regressione lineare produce due coefficienti, uno relativo all'intera scuola quindi il punteggio medio della scuola (intercetta), uno relativo al coefficiente dello status socio-economico di quella scuola (slope)
  S2.LS<-c(S2.LS, summary(fit.ols[[j]])$sigma^2) #varianza per ogni singola scuola
  R2.LS[j]=summary(fit.ols[[j]])$r.squared
  react.pred.ols[[j]]=predict(fit.ols[[j]], list(xj=days.cont, xj2=days.cont^2))
} 

subject.col <- c("#FF0000", "#00FF00", "#0000FF", "#FFFF00", "#FF00FF", "#00FFFF", "#000000", "#FFFFFF", "#FFA500", "#800080", "#008000", "#000080", "#FFC0CB", "#800000", "#808080", "#FFFFF0", "#F0FFF0", "#F0FFFF")
plot(dati$Days, dati$Reaction, type = "n", xlab = "Days", ylab = "Reaction")#, col=subject.col[dati$Subject])
for(j in 1:m) {lines(days.cont, react.pred.ols[[j]])}#, col=subject.col[j])}
lines(x, pooled.pred, col="red", lwd=2)

#### bayesian approach ####
#'number of parameters
p=dim(X[[1]])[2]

#'setting priors

#'prior distribution of iSigma
#'distribution: Wishart
#'parameters: eta0, S0
Sigma=cov(BETA.LS)
iSigma=solve(Sigma)
eta0=p+2 #Hoff page 200
S0=cov(BETA.LS) #prior sum of squares matrix, Hoff page 200

#'prior distribution of theta
#'distribution: Multivariate Normal
#'parameters: mu0, Lambda
theta=apply(BETA.LS, 2, mean)
mu0=apply(BETA.LS, 2, mean) #average values of intercept, b1 and b2, Hoff page 200
L0=cov(BETA.LS) #Hoff page 200
iL0=solve(L0)

#'diffusion of the theta prior
#'this is a really diffuse prior
th1.prior=mu0[1]+c(-1.96,1.96)*sqrt(L0[1,1])
th2.prior=mu0[2]+c(-1.96,1.96)*sqrt(L0[2,2])
th3.prior=mu0[3]+c(-1.96,1.96)*sqrt(L0[3,3])

min.th.prior=th1.prior[1]+th2.prior[1]*days.cont+th3.prior[1]*(days.cont^2)
max.th.prior=th1.prior[2]+th2.prior[2]*days.cont+th3.prior[2]*(days.cont^2)

plot(dati$Days, dati$Reaction, col=subject.col[dati$Subject], ylim = c(-100, 600))
for(j in 1:m) {lines(days.cont, react.pred.ols[[j]], col=subject.col[j])}
lines(days.cont, min.th.prior)
lines(days.cont, max.th.prior)


#'prior distribution of s2
#'distribution: Inverse Gamma
#'parameters: nu0, s20
s2=mean(S2.LS)
s20=mean(S2.LS) #average of within-group sample variance, Hoff page 200
nu0=1


#'sampling distribution of BETA
#'distribution: Multivariate Normale
#'parameters: Ej, Vj
BETA=BETA.LS

#'starting values
THETA.b=S2.b=NULL
Sigma.ps=matrix(0,p,p)
SIGMA.PS=NULL
BETA.ps=BETA*0
BETA.pp=NULL

#'Gibbs sampler
S=10000
for (s in 1:S) {
  #'update theta
  Lm=solve(iL0+m*iSigma)
  mum=Lm%*%(iL0%*%mu0+iSigma%*%apply(BETA,2,sum))
  theta=t(rmvnorm(1,mum,Lm))
  
  #'update Sigma
  mtheta=matrix(theta,m,p,byrow=TRUE)
  iSigma=matrix(rWishart(1, eta0+m, solve(S0+t(BETA-mtheta)%*%(BETA-mtheta))), ncol=3) 
  
  #'update beta_j with Gibbs sampling
  for (j in 1:m) {
    Vj=solve( iSigma + t(X[[j]])%*%X[[j]]/s2 )
    Ej=Vj%*%( iSigma%*%theta + t(X[[j]])%*%Y[[j]]/s2 )
    BETA[j,]=rmvnorm(1,Ej,Vj) 
  }
  
  
  #'update s2
  RSS=0 #residual sum of squared, la somma dei residui^2
  for(j in 1:m) { RSS=RSS+sum( (Y[[j]]-X[[j]]%*%BETA[j,] )^2 ) } #residui^2=Reaction vero - Reaction previsto da OLS, tutto al quadrato
  s2=1/rgamma(1,(nu0+sum(N))/2, (nu0*s20+RSS)/2 )
  
  #'update beta_j with MH
  #dSigma<-det(Sigma) 
  #for(j in 1:m) 
  #{ 
  #  beta.p<-t(rmvnorm(1,BETA[j ,],.5*Sigma)) 
  #  
  #  lr <-sum( dnorm(Y[[j]],beta.p, s2, log = TRUE) - dnorm(Y[[j]], BETA[j,],s2, log = TRUE) + ldmvnorm( t(beta.p),Sigma)  - ldmvnorm(t(BETA[j,]),Sigma))
  #  if( log(runif(1)) < lr ) { BETA[ j ,]<- beta.p }
  #}
  
  
  ##store results
  if(s%%10==0) 
  { 
    cat(s,s2,"\n")
    S2.b<-c(S2.b,s2)
    THETA.b<-rbind(THETA.b,t(theta))
    Sigma.ps<-Sigma.ps+solve(iSigma) #posterior of Sigma  
    BETA.ps<-rbind(BETA.ps, BETA) #posterior of betas
    SIGMA.PS<-rbind(SIGMA.PS,c(solve(iSigma)))
    BETA.pp<-rbind(BETA.pp,rmvnorm(1,theta,solve(iSigma)) ) #this is the beta_s in page 201 of Hoff
  }
}


#'obtaining posterior distribution of betas for each subject
vect.subjects=as.numeric(names(table(dati$Subject)))
BETA.ps[BETA.ps==0]=NA
BETA.ps=BETA.ps[complete.cases(BETA.ps),]
BETA.ps=data.frame(BETA.ps, subject=rep(vect.subjects, S))

BETA.mean.subject=data.frame(b1=rep(NA, 18), b2=rep(NA, 18), b3=rep(NA, 18), subjects=vect.subjects)
for (elem in vect.subjects) {
  BETA.mean.subject[which(BETA.mean.subject$subjects==elem),]=colMeans(BETA.ps[which(BETA.ps$subject==elem),])  
}


#'obtaining posterior distribution of theta 
THETA.mean=apply(THETA.b, 2, mean)


#### bayesian prediction ####
#'BETA
bayes.pred=data.frame(matrix(NA, ncol = 18, nrow = 91)) 
rownames(bayes.pred)=days.cont
for (i in 1:m) {
  bayes.pred[i]=BETA.mean.subject[i,1]+BETA.mean.subject[i,2]*days.cont+BETA.mean.subject[i,3]*(days.cont^2)
}
#'THETA
bayes.theta.pred=THETA.mean[1]+THETA.mean[2]*days.cont+THETA.mean[3]*(days.cont^2)
subject.col <- c("#FF0000", "#00FF00", "#0000FF", "#FFFF00", "#FF00FF", "#00FFFF", "#000000", "#FFFFFF", "#FFA500", "#800080", "#008000", "#000080", "#FFC0CB", "#800000", "#808080", "#FFFFF0", "#F0FFF0", "#F0FFFF")
plot(dati$Days, dati$Reaction, col=subject.col[dati$Subject])
for(j in 1:m) {lines(days.cont, unlist(bayes.pred[j]), col=subject.col[j])}
lines(days.cont, bayes.theta.pred, lwd=4)

#### comparison between bayesian and frequentist model ####
#'comparing bayesian hierarchial model with the frequentist one
par(mfrow=c(1,2))

subject.col <- c("#FF0000", "#00FF00", "#0000FF", "#FFFF00", "#FF00FF", "#00FFFF", "#000000", "#FFFFFF", "#FFA500", "#800080", "#008000", "#000080", "#FFC0CB", "#800000", "#808080", "#FFFFF0", "#F0FFF0", "#F0FFFF")
plot(dati$Days, dati$Reaction, type="n", xlab = "Days", ylab = "Reaction", main = "Bayesian regressions")#, col=subject.col[dati$Subject])
for(j in 1:m) {lines(days.cont, unlist(bayes.pred[j]), col=subject.col[j], lwd=1)}
lines(days.cont, bayes.theta.pred, lwd=4)

plot(dati$Days, dati$Reaction, type="n", xlab="Days", ylab = "Reaction", main="OLS regressions")#, col=subject.col[dati$Subject])
for(j in 1:m) {lines(days.cont, react.pred.ols[[j]], col=subject.col[j], lwd=1)}
react.pred.pooled=predict(ols.pooled, list(Days=days.cont, Days2= days.cont^2))
lines(days.cont, react.pred.pooled, lwd=4)

##### priors and posteriors of theta and Sigma ####

## THETA ##
par(mfrow=c(1,3))
#'prior of theta
prior.theta=rmvnorm(1000, mu0, L0)
#'posterior of theta
THETA.b
#'plot of theta1
plot(density(prior.theta[,1]), xlim=range(prior.theta[,1]), ylim=c(0, 0.06), main="Prior and posterior distribution of theta1", xlab="theta1")
lines(density(THETA.b[,1]), col="red")
legend(175,0.05, legend = c("Prior", "Posterior"), col = c("black", "red"), lty = 1, cex = 1)
#'plot of theta2
plot(density(prior.theta[,2]), xlim=range(prior.theta[,2]), ylim=c(0, 0.15), main="Prior and posterior distribution of theta3", xlab="theta2")
lines(density(THETA.b[,2]), col="red")
legend(-30,0.12, legend = c("Prior", "Posterior"), col = c("black", "red"), lty = 1, cex = 1)
#'plot of theta3
plot(density(prior.theta[,3]), xlim=range(prior.theta[,3]), ylim=c(0, 1.2), main="Prior and posterior distribution of theta3", xlab = "theta3")
lines(density(THETA.b[,3]), col="red")
legend(-4,1, legend = c("Prior", "Posterior"), col = c("black", "red"), lty = 1, cex = 1)

## SIGMA ##
#'prior of Sigma
#'rWishart(rWishart(10000, eta0, solve(S0)) crea una lista di 10000 matrici 3x3
mat=1/rWishart(1000, eta0, solve(S0))
#'io sono interessato a creare una matrice dove le colonne sono V(beta1), v(beta2), V(beta3), Cov(beta1, beta2), Cov(beta2, beta3), Cov(beta1, beta3)
#'sono quindi interessato solamente a quelle 6 colonne, ma essendo la 3x3 estraendo solamente quei valori otteniamo 9 colonne
#'le eccessive 3 colonne sono quelle derivanti dall'angolo della matrice, cioè le covarianze, e infatti sono uguali alle altre colonne relative alla covarianza
##
#'l'ultima colonna è la varianza di beta3
#'la penultima colonna è la covarianza tra beta3 e beta2
#'la terzultima colonna è la covarianza tra beta3 e beta1
#'la quarta colonna è la covarianza tra beta1 e beta2
#'la prima colonna è la varianza di beta1
#'la quinta colonna è la varianza di beta2
prior.Sigma=matrix(mat, ncol = 9, byrow = TRUE)
#'posterior of Sigma
#'SIGMA.PS
#'plot of beta1 variance
par(mfrow=c(2,1))
plot(density(prior.Sigma[,1]), main = "Prior distributon of V(beta1)")
plot(density(SIGMA.PS[,1]), col="red", main = "Posterior distribution of V(beta1)")
#'plot of beta2 variance
par(mfrow=c(2,1))
plot(density(prior.Sigma[,5]), main = "Prior distributon of V(beta2)")
plot(density(SIGMA.PS[,5]), col="red", main = "Posterior distribution of V(beta2)")
#'plot of beta3 variance
par(mfrow=c(2,1))
plot(density(prior.Sigma[,9]), main = "Prior distributon of V(beta3)")
plot(density(SIGMA.PS[,9]), col="red", main = "Posterior distribution of V(beta3)")
#'plot of beta1 and beta 2 covariance
par(mfrow=c(2,1))
plot(density(prior.Sigma[,4]), main = "Prior distributon of Cov(beta1, beta2)")
plot(density(SIGMA.PS[,4]), col="red", main = "Posterior distribution of Cov(beta1, beta2)")
#'plot of beta1 and beta3 covariance
par(mfrow=c(2,1))
plot(density(prior.Sigma[,7]), main = "Prior distributon of Cov(beta1, beta3)")
plot(density(SIGMA.PS[,7]), col="red", main = "Posterior distribution of Cov(beta1, beta3)")
#'plot of beta2 and beta 3 covariance
par(mfrow=c(2,1))
plot(density(prior.Sigma[,8]), main = "Prior distributon of Cov(beta2, beta3)")
plot(density(SIGMA.PS[,8]), col="red", main = "Posterior distribution of Cov(beta2, beta3)")

## BETA ##
#'prior of beta
BETA.LS
#'posterior of beta
BETA.ps
par(mfrow=c(1,3))
#'plot of prior and posterior of beta1
plot(density(BETA.LS[,1]), main = "Prior and posterior distributon of beta1", ylim=c(0,0.02))
lines(density(BETA.ps[,1]), col="red")
legend(150,0.02, legend = c("Prior", "Posterior"), col = c("black", "red"), lty = 1, cex = 0.8)
#'plot of prior and posterior of beta2
plot(density(BETA.LS[,2]), main = "Prior and posterior distributon of beta2", ylim=c(0,0.08))
lines(density(BETA.ps[,2]), col="red")
legend(-30,0.06, legend = c("Prior", "Posterior"), col = c("black", "red"), lty = 1, cex = 0.8)
#'plot of prior and posterior of beta3
plot(density(BETA.LS[,3]), main = "Prior and posterior distributon of beta3", ylim=c(0,0.5))
lines(density(BETA.ps[,3]), col="red")
legend(-30,0.06, legend = c("Prior", "Posterior"), col = c("black", "red"), lty = 1, cex = 0.8)


#### BETA VS PREDICTIVE BETA #### 
#'non ha molto senso farlo
#'plotting prior and predictive of beta1
plot(density(BETA[,1]), main = "Prior and predictive distribution of beta1", ylim=c(0,0.04))
lines(density(BETA.pp[,1]), col="red")
legend(200,0.03, legend = c("Prior", "Predictive"), col = c("black", "red"), lty = 1, cex = 0.8)
#'plotting prior and predictive of beta2
plot(density(BETA.LS[,2]), main = "Prior and predictive distribution of beta2", ylim=c(0,0.08))
lines(density(BETA.pp[,2]), col="red")
legend(-30,0.06, legend = c("Prior", "Predictive"), col = c("black", "red"), lty = 1, cex = 0.8)
#'plotting prior and predictive of beta3
plot(density(BETA.LS[,3]), main = "Prior and predictive distribution of beta3", ylim=c(0,0.4), xlim=c(-8,8))
lines(density(BETA.pp[,3]), col="red")
legend(-30,0.06, legend = c("Prior", "Predictive"), col = c("black", "red"), lty = 1, cex = 0.8)


#### THETA VS PREDICTIVE BETA ####
#plotting posterior of theta1 and predictive of beta1
plot(density(THETA.b[,1]), main = "Posterior theta1 and predictive beta1", xlim=c(100,400))
lines(density(BETA.pp[,1]), col="red")
legend(150,0.04, legend = c("Posterior theta1", "Predictive beta1"), col = c("black", "red"), lty = 1, cex = 0.8)
#plotting posterior of theta2 and predictive of beta2
plot(density(THETA.b[,2]), main = "Posterior theta2 and predictive beta2", ylim=c(0,0.12), xlim=c(-50,50))
lines(density(BETA.pp[,2]), col="red")
legend(-30,0.06, legend = c("Posterior theta2", "Predictive beta2"), col = c("black", "red"), lty = 1, cex = 0.8)
#plotting posterior of theta3 and predictive of beta3
plot(density(THETA.b[,3]), main = "Posterior theta3 and predictive beta3", ylim=c(0,1.2), xlim=c(-6,6))
lines(density(BETA.pp[,3]), col="red")
legend(-4,0.6, legend = c("Posterior theta3", "Predictive beta3"), col = c("black", "red"), lty = 1, cex = 0.8)


#### MCMC DIAGNOSTICS ####
#THETA[,1]
par(mfrow=c(1,2))
plot(THETA.b[,1], type = "l", main="Theta1")
hist(THETA.b[,1], main="Theta1")
#THETA[,2]
par(mfrow=c(1,2))
plot(THETA.b[,2], type = "l", main="Theta2")
hist(THETA.b[,2], main="Theta2")
#THETA[,3]
par(mfrow=c(1,2), cex.main=0.8)
plot(THETA.b[,3], type = "l", main="Theta3")
hist(THETA.b[,3], main="Theta3")
#SIGMA[1,1]
par(mfrow=c(1,2), cex.main=0.8)
plot(SIGMA.PS[,1], type = "l", main="SIGMA[1,1]")
hist(SIGMA.PS[,1], main="SIGMA[1,1]")
#SIGMA[2,2]
par(mfrow=c(1,2), cex.main=0.8)
plot(SIGMA.PS[,5], type = "l", main="SIGMA[2,2]")
hist(SIGMA.PS[,5], main="SIGMA[2,2]")
#SIGMA[3,3]
par(mfrow=c(1,2), cex.main=0.8)
plot(SIGMA.PS[,9], type = "l", main="SIGMA[3,3]")
hist(SIGMA.PS[,9], main="SIGMA[3,3]")
#SIGMA[1,2]
par(mfrow=c(1,2), cex.main=0.8)
plot(SIGMA.PS[,4], type = "l", main="SIGMA[1,2]")
hist(SIGMA.PS[,4], main="SIGMA[1,2]")
#SIGMA[1,3]
par(mfrow=c(1,2), cex.main=0.8)
plot(SIGMA.PS[,7], type = "l", main="SIGMA[1,3]")
hist(SIGMA.PS[,7], main="SIGMA[1,3]")
#SIGMA[2,3]
par(mfrow=c(1,2), cex.main=0.8)
plot(SIGMA.PS[,8], type = "l", main="SIGMA[2,3]")
hist(SIGMA.PS[,8], main="SIGMA[2,3]")





