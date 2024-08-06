vec=function(A){
return(as.vector(A))
}

soft=function(a,b){
c=abs(a)-b
c[c<0]=0
c=c*sign(a)
return(c)
}

bimin=function(mat){
min_element <- min(mat)
min_indices <- which(mat == min_element, arr.ind = TRUE)
if (nrow(min_indices) > 1) {
min_indices <- min_indices[nrow(min_indices), ]
}
return(min_indices)
}

#' @importFrom CppMatrix matrixInverse matrixMultiply matrixVectorMultiply matrixEigen
positiveadj=function(A,min.eps=0.001){
a=matrixEigen(A)
d=c(a$values)
d[d<min.eps]=min.eps
B=matrixMultiply(a$vectors,t(a$vectors)*d)
return(B)
}

dmcp=function(x,lam,a=3){
b=lam-abs(x)/a
b[abs(x/a)>lam]=0
d=0*b
d[abs(x/a)<=lam]=-1/a
return(A=data.frame(d1=b,d2=d))
}

mcp=function(x,lam,a=3){
b=abs(x)
z=soft(x,lam)/(1-1/a)
z[which(b>(a*lam))]=x[which(b>(a*lam))]
return(z)
}

mad=function(a){
b=1.483*median(abs(a-median(a)))
return(b)
}

trace=function(A){
a=sum(diag(A))
return(a)
}

bimin=function(mat){
min_element <- min(mat)
min_indices <- which(mat == min_element, arr.ind = TRUE)
if (nrow(min_indices) > 1) {
min_indices <- min_indices[nrow(min_indices), ]
}
return(min_indices)
}

#' @importFrom CppMatrix matrixInverse matrixMultiply matrixVectorMultiply matrixEigen
matrixsqrt=function(A){
fit=matrixEigen(A)
d=c(fit$value)
d1=d*0
d1[d>0]=1/d[d>0]
d=sqrt(d)
d1=sqrt(d1)
A=matrixMultiply(fit$vector,t(fit$vector)*d)
B=matrixMultiply(fit$vector,t(fit$vector)*d1)
C=list(w=A,wi=B,eigenfit=fit)
return(C)
}

IVweight=function(byse,bXse,Rxy){
bZse=cbind(bXse,byse)
p=dim(bZse)[2]
n=dim(bZse)[1]
RxyList=array(0,c(n,p,p))
for(i in 1:n){
s=bZse[i,]
RxyList[i,,]=t(t(Rxy)*s)*s
}
return(RxyList)
}

biasterm=function(RxyList,indvalid){
X=RxyList[1,,]*0
for(i in indvalid){
X=X+RxyList[i,,]
}
return(X)
}

standardized.residual=function(res,RxyList,theta,adjust=1){
n=length(res)
tilde.theta=c(theta,-1)
vars=res
p=length(theta)
for(i in 1:n){
Rxyi=RxyList[i,,]
adjusti=sqrt(c(rep(adjust,p),1))
Rxyi=t(t(Rxyi)*adjusti)*adjusti
vars[i]=max(0.05,sum(tilde.theta*c(Rxyi%*%tilde.theta)))
}
return(res/sqrt(vars))
}

checkzero=function(A,nonzero.frac){
p=dim(A)[2]
for(i in 1:p){
a=mean(A[,i]!=0)-nonzero.frac[i]
if(a<0){
b=abs(a)
c=sample(which(A[,i]==0),b*length(which(A[,i]==0)),replace=F)
A[c,i]=1
}
}
return(A)
}

verrorvar=function(Rbb,Suu,Suv,hx,theta,hy,pleiotropy.var=0){
p=dim(Rbb)[2]
Sigmabb=sqrt(hx)%*%Rbb%*%sqrt(hx)
Sigmauu=sqrt(diag(p)-hx)%*%Suu%*%sqrt(diag(p)-hx)
Sigmaxx=Sigmabb+Sigmauu
plevar=as.numeric(pleiotropy.var*t(theta)%*%Sigmabb%*%theta)
c=t(theta)%*%Sigmaxx%*%theta+plevar-t(theta)%*%Sigmabb%*%theta/hy*(1+pleiotropy.var);c=c[1,1]
b=-2*t(theta)%*%sqrt(diag(p)-hx)%*%Suv;b=b[1,1]
svv=(-b+sqrt(b^2-4*c))/2
return(svv^2)
}

hclust_cut <- function(corr_matrix, k) {
dist_matrix <- as.dist(1 - corr_matrix)
hc <- hclust(dist_matrix, method = "ward.D2")
cluster_result <- cutree(hc, k = k)
return(cluster_result)
}

block_subsampling <- function(cluster_result, rho) {
unique_clusters <- unique(cluster_result)
num_selected_clusters <- ceiling(length(unique_clusters) * rho)
selected_clusters <- sample(unique_clusters, num_selected_clusters)
selected_indices <- which(cluster_result %in% selected_clusters)
return(selected_indices)
}

colSD=function(A){
a=A[1,]
for(i in 1:ncol(A)){
a[i]=sd(A[,i])
}
return(a)
}

subsample_indices <- function(cluster.index, sampling.ratio) {
unique_clusters <- unique(cluster.index)
n_subsample <- round(length(unique_clusters) * sampling.ratio)
subsampled_clusters <- sample(unique_clusters, n_subsample)
subsample_indices <- which(cluster.index %in% subsampled_clusters)
return(subsample_indices)
}

generate_default_cluster <- function(n) {
max_category <- ceiling(n / 2)
cluster_index <- rep(1:max_category, each = 2, length.out = n)
return(cluster_index)
}

bimin=function(mat){
min_element <- min(mat)
min_indices <- which(mat == min_element, arr.ind = TRUE)
if (nrow(min_indices) > 1) {
min_indices <- min_indices[nrow(min_indices), ]
}
return(min_indices)
}

parametric.bootstrap=function(bX,bXse,Rxy,theta,LD,RC,var.inf=0){
n=nrow(bX)
p=ncol(bX)
E=matrix(0,n,p+1)
for(i in 1:n){
sei=c(bXse[i,],1)
Vi=t(t(Rxy)*sei)*sei
e=MASS::mvrnorm(n=1,mu=rep(0,p+1),Sigma=Vi)
E[i,]=e
}
hatX=bX+matrixMultiply(RC,E[,1:p])
haty=matrixVectorMultiply(bX,theta)+matrixVectorMultiply(RC,E[,p+1])
if(var.inf>0){
haty=haty+matrixVectorMultiply(LD,rnorm(n=n,mean=0,sd=sqrt(var.inf)))
}
return(list(hatX=hatX,haty=haty))
}

reliability.adj.uv=function(bx,bxse,Theta="identify",thres=0.7){
if(Theta[1]=="identify"){
total.var=mean(bx^2)
error.var=mean(bxse^2)
reliability=(total.var-error.var)/total.var
r=1
if(reliability<thres){
r=total.var/error.var*(1-thres)
}
r=sqrt(r)
}else{
r=1
total.var=mean(bx*matrixVectorMultiply(Theta,bx))
error.var=mean(bxse^2)
reliability=(total.var-error.var)/total.var
if(reliability<thres){
r=total.var/error.var*(1-thres)
}
r=sqrt(r)
}
return(r)
}

reliability.adj=function(bX,bXse,Theta="identify",thres=0.7){
if(Theta[1]=="identify"){
p=ncol(bX)
r=rep(1,p)
total.var=colMeans(bX^2)
error.var=colMeans(bXse^2)
reliability=(total.var-error.var)/total.var
ind=which(reliability<thres)
if(length(ind)>0){
r[ind]=total.var[ind]/error.var[ind]*(1-thres)
}
r=sqrt(r)
}else{
p=ncol(bX)
r=rep(1,p)
m=length(bX[,1])
total.var=diag(matrixListProduct(list(t(bX),Theta,bX)))/m
error.var=colMeans(bXse^2)
reliability=(total.var-error.var)/total.var
ind=which(reliability<thres)
if(length(ind)>0){
r[ind]=total.var[ind]/error.var[ind]*(1-thres)
}
r=sqrt(r)
}
return(r)
}
