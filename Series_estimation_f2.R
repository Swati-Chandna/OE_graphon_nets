###################SERIES ESTIMATION for f_2 example in Chandna & Maugis###########
#INPUT: adjacencies simulated from f_2 with n=100, m=150, MC=50
#(available as .mat and .RData files)

#OUTPUT: estimated graphon matrices and estimated latent node positions:
#(a) Oracle based - preds11nc_all_orac2 (n by n by MC)
#(b) Actual - preds11nc_all (n by n by MC)
#(c) xihatt_all  (n by MC vector of estimated node positions)
#############################

require(reshape2);
require(ggplot2);
require(grid);
require(loe);
require(vegan);
require(mgcv);
require(Matrix);
require(rgl);
require(GoFKernel)
require(fields)##for colorbar.plot
require(igraph)
library(R.matlab)#

Gc_all_I=readMat("Adj_f2_n100m150_MC50_I.mat")
Gc_all_II=readMat("Adj_f2_n100m150_MC50_II.mat")

A1=Gc_all_I$Gc.MA
A2=Gc_all_II$Gc.MA
dA1=dim(A1)
n=dA1[1]
m=150
MC=50 

A=array(NA, dim=c(n,n,m,MC))
A[,,1:m,1:25]=A1[,,1:m,];
A[,,1:m,26:50]=A2[,,1:m,];

rm(A1)
rm(A2)
rm(Gc_all_I)
rm(Gc_all_II)

load("f2_n100_xi_all.RData")


#MC=2###for a quick run

perm_xitr=matrix(0,n,MC);
Atr=array(NA,dim=c(n,n,m,MC))
mean_Atr=array(NA,dim=c(n,n,MC))

for (mc in 1:MC)
{
  perm_xitr[,mc]=order(xi_all[,mc])
  Atr[,,,mc]=A[perm_xitr[,mc],perm_xitr[,mc],,mc]
}

mean_Atr[,,]=apply(Atr, c(1, 2, 4), mean)



########################ESTIMATION

graph_smpl_dist_base <- function(G){#Vertex distances based on neighborhood similarity
  ## Input:
  # G : a graph sample; i.e., a list of adjacency matrices of simple graphs where vertices can be identified across realization.
  # 
  ## Output:
  # A n xn matrix dst, where n is the order of the entris of G. dst_{ij} is the distance between 
  # vertices i and j computed using G. The metric is that of Airoldi et. al. (2013) simplified
  # to address undirected graphs, and randomized for stability. Assuming G is an iid sample drawn
  # from a fixed graphon f and that each vertex i is associated with a latent x_i, this metric estimates
  # consistentli \int_{[0,1]} (f(x_i,x)-f(x_j,x))dx and therefore measures how the neighborhood of i
  # differs from that of j.
  #
  ## Bibliography:
  # E. M. Airoldi, T. B. Costa and S. H. Chan, (2013) "Stochastic blockmodel approximation of a graphon: Theory and consistent estimation", NIPS
  #
  #
  ## Allocating variables
  m <- length(G);
  N=m;
  n <- ncol(G[[1]]);
  Gg   <- array(unlist(G),c(n,n,m));
  # rijk <- array(0,c(n,n,n));
  #
  ## Core loop
  R <- sample.int(m,floor(m/2)); #Randomization, conrasting with Airoldi et. al. (2013)
  
  #view(R)
  
  rij <- (apply(Gg[,,R],c(1,2),mean)%*%apply(Gg[,,-R],c(1,2),mean))/(n-2)
  ## Metric computation
  rhohat <- sum(Gg)/(n*(n-1)*m)
  #view(rhohat)
  dij <- (outer(diag(rij),diag(rij),'+')-rij-t(rij))/(rhohat^2);#(N^2*(n-2)); #Only the first 4 terms of eq. (5) of Airoldi et. al. (2013)
  dij[which(dij<0)] <- 0;
  #
  ## Output
  return(list(dij=dij,rij=rij))
}
#
graph_smpl_dist <- function(G,r=20,gram=F,low_rank=0){#Further randomizes vertex distance
  ## Input:
  # G : a graph sample; i.e., a list of adjacency matrices of simple graphs where vertices can be identified across realization.
  # r : number of independent replications
  #
  ## Output:
  # A n xn matrix dst, where n is the order of the entris of G. dst_{ij} is the average distance
  # between vertices i and j computed using graph_smpl_dist_base in r independent replications. 
  #
  #
  ## Main loop (vectorized)
  if (!gram){
    dst <- replicate(r,graph_smpl_dist_base(G)$dij,simplify=F)
  } else {
    dst <- replicate(r,graph_smpl_dist_base(G)$rij,simplify=F)
  }
  #
  ## Allocating and averaging
  dst <- array(unlist(dst),c(ncol(G[[1]]),ncol(G[[1]]),length(G)));
  dst <- apply(dst,c(1,2),mean);
  ## Eigenvalue thresholding
  if (low_rank>0){
    edc <- eigen(dst);
    sbs <- if (gram) {head(1:nrow(dst),low_rank)} else {c(1,tail(1:nrow(dst),low_rank))};
    dst <- edc$vectors[,sbs]%*%diag(edc$values[sbs])%*%t(edc$vectors[,sbs]);
  }
  #
  ## Output
  return(dst)
}


graph_smpl_order <- function(G, type = 'i', ini = 'rand', ...) {
  # Ordering of vertices using metric or ordinal embedding method
  #
  # Input:
  # G : a graph sample; a list of adjacency matrices of simple graphs where vertices can be identified across realization.
  # ... : additional arguments to pass on to graph_smpl_dist
  #
  # Output:
  # A list containing:
  # - `perm`: A permutation of the vertices in G such that adjacent vertices in the ordering are close
  # - `dst`: The computed distance matrix
  
  ## Computing distance matrix dst if not provided
  if (!('dst' %in% names(list(...)))) {
    dst <- graph_smpl_dist(G, ...)
  }
  
  if (type == 'i') {
    # Calling Isomap function
    perm <- order(scores(isomap(dst, k = ncol(G[[1]]), ndim = 1)))
    res_SOE <- NULL  # No SOE result for type 'i'
  } else {
    # Calling SOE function
    comp_dst <- get.order(dst)
    res_SOE <- SOE(comp_dst, ncol(G[[1]]), p = 1, iniX = ini, rnd = nrow(comp_dst))
    perm <- order(res_SOE$X)  
  }
  
  # Output:
  return(list(perm = perm, dst = dst, res_SOE = res_SOE))
}

####################
xihatt=matrix(0,n,MC);
minstress=matrix(0,MC,1);
preds11nc_all_orac2=array(NA, dim=c(n,n,MC))
preds11nc_all=array(NA, dim=c(n,n,MC))

P=20
perm_pp<-matrix(NA,n,P);

xihatt_pp<-matrix(NA,n,P);
stressfn_pp<-matrix(NA,P,1);
dsthat<-array(NA,dim=c(n,n,MC));
dsthat_pp<-array(NA,dim=c(n,n,P));
Atr_xi=array(NA,dim=c(n,n,m,MC));


for (mc in 1:MC)
{
  Amc=array(0,dim=c(n,n,m));
  Amc=A[,,,mc];
  G=lapply(seq(dim(Amc)[3]), function(x) Amc[ , , x])
  {
    for (p in 1:P){
      reslt=graph_smpl_order(G,type='j') #
      res_SOE=reslt$res_SOE
      dsthat_pp[,,p]=reslt$dst
      xihatt_pp[,p]=res_SOE$X 
      perm_pp[,p]<-order(res_SOE$X)
      stressfn_pp[p]<-res_SOE$str
      print(stressfn_pp[p])
    }
    
    
    minstress[mc]=which.min(stressfn_pp);
    xihatt[,mc]=xihatt_pp[,minstress[mc]]
    perm<-perm_pp[,minstress[mc]];
    dsthat[,,mc]=dsthat_pp[,,minstress[mc]]
  }
}

Atr_estd=array(NA,dim=c(n,n,m,MC))
mean_Atr_estd=array(NA,dim=c(n,n,MC))
perm_xiestd=matrix(0,n,MC);

for (mc in 1:MC)
{
  perm_xiestd[,mc]=order(xihatt[,mc]) 
  Atr_estd[,,,mc]=A[perm_xiestd[,mc],perm_xiestd[,mc],,mc]
}

mean_Atr_estd[,,]=apply(Atr_estd, c(1, 2, 4), mean)




Gshatavg=matrix(NA,n,n);
for (mc in 1:MC)
{
  Gshatavg=mean_Atr_estd[,,mc]
  sxihatt=sort(xihatt[,mc])
  xi=xi_all[,mc]
  sxi=sort(xi)
  aa=seq(from=0.01,to=0.98,length.out=n);
  
  y1<-cbind(c(Gshatavg))
  v1<-rep(sxihatt,times=n)
  
  v2<-rep(sxihatt,each=n);
  
  xxx11nc <- bam(y1~s(v1,v2),family=gaussian(link = "identity"))
  
  #summary(xxx11nc)###check adequacy of default K
  
  loc1P<- data.frame(v1=rep(sxihatt,times=n),v2=rep(sxihatt,each=n))
  
  preds11nc<- matrix(predict.gam(xxx11nc,newdata=loc1P,type='response'),n)
  
  preds11nc_all[,,mc]=preds11nc;
  
  pxi=sort.list(xi);

  xihatt_os=pxi/(n+1); 
  sxi2=sort(xihatt_os);

  
  Atr_xi[,,,mc]=A[perm_xitr[,mc],perm_xitr[,mc],,mc]
  mean_Atr[,,]=apply(Atr_xi, c(1, 2, 4), mean)
  Gsavg=mean_Atr[,,mc]
  
  
  y1<-cbind(c(Gsavg))  
  v1_os<-rep(sxi2,times=n)
  v2_os<-rep(sxi2,each=n);
  
  xxx11nc_orac2 <- bam(y1~s(v1_os,v2_os),family=gaussian(link = "identity"))
  loc1P2<- data.frame(v1_os=rep(sxi,times=n),v2_os=rep(sxi,each=n))
  
  preds11nc_orac2<- matrix(predict.gam(xxx11nc_orac2,newdata=loc1P2,type='response'),n)
  preds11nc_all_orac2[,,mc]=preds11nc_orac2;
  
} 
xihatt_all=xihatt;
N=m

#######################visualize the true and estimated graphons (oracle vs actual)
mcc=1
par(mfrow = c(1, 3)) 
image.plot(0.84*fijtr[,,mcc],main = "True P")
image.plot(preds11nc_all_orac2[,,mcc], main = "Estimated P (oracle)")
image.plot(preds11nc_all[,,mcc],main = "Estimated P (proposed)")