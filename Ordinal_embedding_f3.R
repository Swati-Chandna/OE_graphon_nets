
######################################OE based latent node position estimation for a 2-blockmodel f2 as in [Chandna and Maugis, 2025*]
##############Code below proceeds by generation of network samples A_1,...A_m for a specified (n,m)
#############Output: the latent node position vector: xihatt which is subsequently used to visualize the unknown true network structure
require(reshape2);
require(ggplot2);
require(grid);
require(loe);
require(vegan);
require(mgcv);
require(Matrix);
require(rgl);
require(GoFKernel)
require(fields)
require(igraph)
library(R.matlab)



n=50
#n=100
#n=150

m=5
#m=150


MC=2 ##Monte Carlo replications


########################DATA GENERATION 
pref.matrix= cbind( c(.7, .1), c(.1, .3) )
###unequal block sizes: 
h1=ceiling(0.6*n);
h2=n-h1
block.sizes=c(h1,h2);
A=array(0,dim=c(n,n,m,MC));
uu=matrix(0,n,MC)
for (mc in 1:MC)
{
  uu[,mc]=runif(n);
  g=sample_sbm(n, pref.matrix, block.sizes, directed = FALSE, loops = FALSE)

  for (l in 1:m)
  {
    aa=as_adjacency_matrix(g,type="both",sparse = FALSE)
    A[,,l,mc]=aa[sort.list(uu[,mc]),sort.list(uu[,mc])]##randomize the true ordering

  }

}


###############################ESTIMATION:

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


graph_smpl_order <- function(G,type='i',ini='rand',...){# Ordering of vertices using metric or ordinal embedding method
  ## Input:
  # G : a graph sample; i.e., a list of adjacency matrices of simple graphs where vertices can be identified across realization.
  # ... : aditional arguments to pass on to graph_smpl_dist
  #
  # Output:
  # A permutaion perm of the vertices in G such that adjacent vertices in the ordering are close
  # according the the distance matrix computed with graph_smpl_dist. Uses the SOE function from the
  # 'loe' package.
  #
  #
  ## Computing distance matrix dst if not provided
  if (!('dst' %in% names(list(...)))){
    dst  <- graph_smpl_dist(G,...);
  };
  #
  if (type == 'i'){
    ## Calling Isomap function
    perm <- order(scores(isomap(dst,k=ncol(G[[1]]),ndim=1)));
  } else {
    ## Calling SOE function
    comp_dst <- get.order(dst)
    res_SOE <- (SOE(comp_dst,ncol(G[[1]]),p=1,iniX=ini,rnd=nrow(comp_dst)));
    
    
  }
  ## Output:
  return(res_SOE)
}




####################
xihatt=matrix(0,n,MC);####final n x 1 vector of latent node position estimates for each Monte Carlo replication
minstress=matrix(0,MC,1);###the corresponding minimized value of the SOE stress function


P=20#############to remove sensitivity of the SOE algorithm to its random initialization

xihatt_pp<-matrix(NA,n,P);
stressfn_pp<-matrix(NA,P,1);

for (mc in 1:MC)
{
  Amc=array(0,dim=c(n,n,m));
  Amc=A[,,,mc];
  G=lapply(seq(dim(Amc)[3]), function(x) Amc[ , , x])
  {
    for (p in 1:P){
      res_SOE<- graph_smpl_order(G,type='j')
      
      xihatt_pp[,p]=res_SOE$X  
    
      stressfn_pp[p]<-res_SOE$str
      print(stressfn_pp[p])
    }
    
    
    minstress[mc]=which.min(stressfn_pp);
    xihatt[,mc]=xihatt_pp[,minstress[mc]]###select the final embedding as the one which corresponds to the minimum stress function
  }
}

##################visualize one of the adjacencies as observed vs with the estimated ordering of nodes
mcc=2
ll=2 ##visualize one of the m networks
permcc=order(xihatt[,mcc])
par(mfrow = c(1, 2))  
image.plot(A[,,ll,mcc])
image.plot(A[permcc,permcc,ll,mcc])


