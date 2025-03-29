# ordinal-embedding-based-graphon-estimation

The code implements the graphon estimator as proposed in 'Ordinal embedding for network estimation via graphon' by S. Chandna and P.A. Maugis. 

Input: A collection of m binary networks (adjacencies), each of size n x n , observed on the same set of nodes. Output: Latent node position estimates which are used to lead to a graphon matrix estimate (of size n x n) either via series or local constant approximation. 

Example code illustrates the method for a 2-block stochastic blockmodel.
