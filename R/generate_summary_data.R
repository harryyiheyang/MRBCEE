#' Generate Simulated Data for Mendelian Randomization Analysis
#'
#' This function generates simulated data for Mendelian Randomization (MR) analysis,
#' considering genetic effects, estimation errors, and horizontal pleiotropy.
#' It allows for different distributions of genetic effects and pleiotropy,
#' and accommodates both independent and correlated instrumental variables (IVs).
#'
#' @param theta An (px1) vector of causal effects.
#' @param m The number of instrumental variables (IVs).
#' @param Rbb An (pxp) correlation matrix of genetic effects.
#' @param Ruv An ((p+1)x(p+1)) correlation matrix of residuals in outcome and exposures; the outcome is the first one.
#' @param Rnn An ((p+1)x(p+1)) correlation matrix of sample overlap; the outcome is the first one.
#' @param Nxy An ((p+1)x1) vector of GWAS sample sizes; the outcome is the first one.
#' @param Hxy An ((p+1)x1) vector of heritabilities; the outcome is the first one.
#' @param LD An (mxm) correlation matrix of the IVs or "identify" indicating independent IVs.
#' @param zero.frac An (px1) vector with all entries in (0,1]; each entry is the probability of deltaj such that betaj=betaj'*deltaj.
#' @param pleiotropy.frac A number in [0,0.5) indicating the fraction of IVs affected by pleiotropy.
#' @param pleiotropy.var A number in [0,+∞) indicating the variance attributed to pleiotropy.
#' @param effect.dis Distribution of genetic effects: "normal" (default), "uniform", or "t" distribution (with degree of freedom 5).
#' @param pleiotropy.dis Distribution of pleiotropy effects: "normal" (default), "uniform", or "t" distribution (with degree of freedom 5).
#'
#' @return A list containing simulated GWAS effect sizes for exposures (bX), their standard errors (bXse),
#'         the GWAS effect size for the outcome (by), its standard error (byse), and the pleiotropy effects (pleiotropy).
#'
#' @importFrom MASS mvrnorm
#' @importFrom mvtnorm rmvt
#' @importFrom stats rt sd var runif rnorm
#'
#' @note This function requires the `MASS`, `mvtnorm`, `stats`, and `Matrix` packages.
#' @export

generate_summary_data=function(theta,m,Rbb,Ruv,Rnn,Nxy,Hxy,LD="identify",non.zero.frac,pleiotropy.frac=0,pleiotropy.var=0.5,effect.dis="normal",pleiotropy.dis="uniform"){

p=length(theta)
hy=Hxy[p+1];hx=Hxy[1:p];Nx=Nxy[p+1];Ny=Nxy[1:p];

################ generate genetic effect ##########################
if(effect.dis=="normal"){
b0=MASS::mvrnorm(m,rep(0,p),Rbb)
b0=b0/sqrt(m);b0=t(t(b0)*sqrt(hx/non.zero.frac))
cutoff=qnorm(1-non.zero.frac)
b00=MASS::mvrnorm(m,rep(0,p),Rbb)
for(i in 1:p){
b00[,i]=b00[,i]>cutoff[i]
}
b00=checkzero(b00,nonzero.frac=non.zero.frac*0.9)
### Ensure that there must be 0.9*non.zero.frac nonzero elements
b0=b00*b0
}
if(effect.dis=="uniform"){
b0=MASS::mvrnorm(m,rep(0,p),Rbb);b0=(pnorm(b0)-0.5)*sqrt(12);b0=b0/sqrt(m);b0=t(t(b0)*sqrt(hx/non.zero.frac))
cutoff=qnorm(1-non.zero.frac)
b00=MASS::mvrnorm(m,rep(0,p),Rbb)
for(i in 1:p){
b00[,i]=b00[,i]>cutoff[i]
}
b00=checkzero(b00,nonzero.frac=non.zero.frac*0.9)
### Ensure that there must be 0.9*non.zero.frac nonzero elements
b0=b00*b0
}
if(effect.dis=="t"){
b0=mvtnorm::rmvt(m,Rbb,df=5);b0=b0/sqrt(m)/sqrt(5/3);b0=t(t(b0)*sqrt(hx/non.zero.frac))
cutoff=qnorm(1-non.zero.frac)
b00=MASS::mvrnorm(m,rep(0,p),Rbb)
for(i in 1:p){
b00[,i]=b00[,i]>cutoff[i]
}
b00=checkzero(b00,nonzero.frac=non.zero.frac*0.9)
### Ensure that there must be 0.9*non.zero.frac nonzero elements
b0=b00*b0
}

################### Generate Estimation Error #######################
Svv=verrorvar(Rbb=Rbb,Suu=Ruv[-1,-1],Suv=Ruv[1,-1],hx=diag(hx),theta=theta,hy=hy,pleiotropy.var=pleiotropy.var)
Duv=c(sqrt(1-hx),sqrt(Svv))
Sigmauv=diag(Duv)%*%Ruv%*%diag(Duv)
Vxy=Sigmauv*Rnn

E=MASS::mvrnorm(m,rep(0,p+1),Vxy)
E=E%*%diag(1/sqrt(Nxy))

if(LD[1]=="identify"){
E1=E
bX=b0+E1[,-(p+1)]
by=b0%*%theta+E1[,p+1]
bX=as.matrix(bX)
by=as.vector(by)
}

if(LD[1]!="identify"){
C=matrixsqrt(LD)$w
E1=C%*%E
bX=LD%*%b0+E1[,-(p+1)]
by=LD%*%b0%*%theta+E1[,p+1]
bX=as.matrix(bX)
by=as.vector(by)
}

############################    generate pleiotropy    #####################

if(pleiotropy.frac!=0){
indpleio=sample(1:m,pleiotropy.frac*m,replace=F)
bu=b0[,1]*0

if(pleiotropy.dis=="normal"){
bu1=rnorm(length(indpleio),0,1)
bu[indpleio]=bu1
}
if(pleiotropy.dis=="uniform"){
bu1=runif(length(indpleio),0,1)
bu[indpleio]=bu1
}
if(pleiotropy.dis=="t"){
bu1=stats::rt(length(indpleio),5,ncp=0)
bu[indpleio]=bu1
}

s=as.vector(bu)
s1=as.vector(b0%*%theta)
ra=var(s)/var(s1)/pleiotropy.var
if(LD[1]=="identify"){by=by+s/sqrt(ra)}
if(LD[1]!="identify"){by=by+LD%*%s/sqrt(ra)}
}
by=c(by)
byse=rep(sd(E1[,p+1]),m)
bXse=matrix(1,m,p)
for(jj in 1:p){bXse[,jj]=sd(E1[,jj])}

A=list(bX=bX,bXse=bXse,by=by,byse=byse,pleiotropy=s/sqrt(ra))
return(A)
}
