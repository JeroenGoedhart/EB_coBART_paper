FiniteSum <- function(x) {
  sum(x[is.finite(x)])
} 

EstimateLeafNode <- function(Trees, ntree) {
  ids <- which(trees$var==-1) #check which rows correspond to leaf nodes
  samples <- trees$value[ids] #obtain samples of leaf nodes
  varhat <- (1/length(samples))*FiniteSum(samples^2) #maximum likelihood estimate of variance for known mean (equals 0)
  k_hat <- 3/(sqrt(varhat)*sqrt(ntree))
  return(c(varhat=varhat,k_hat=k_hat))
}

LikelihoodBin <- function(Ypred,Y){
  result <- apply(Ypred,1, function(x) Y*log(x)+(1-Y)*log(1-x))
  return(t(result))
}

getDepth <- function(tree) {
  getDepthRecurse <- function(tree, depth) {
    node <- list(
      depth = depth
    )
    if (tree$var[1] == -1) {
      node$n_nodes <- 1
      return(node)
    }
    
    headOfLeftBranch <- tree[-1,]
    left <- getDepthRecurse(headOfLeftBranch, depth + 1)
    n_nodes.left <- left$n_nodes
    left$n_nodes <- NULL
    
    headOfRightBranch <- tree[seq.int(2 + n_nodes.left, nrow(tree)),]
    right <- getDepthRecurse(headOfRightBranch, depth + 1)
    n_nodes.right <- right$n_nodes
    right$n_nodes <- NULL
    
    node$n_nodes <- 1 + n_nodes.left + n_nodes.right
    node$depth <- c(node$depth, left$depth, right$depth)
    return(node)
  }
  result <- getDepthRecurse(tree, 0)
  
  return(result$depth)
}

LikelihoodTreeStructure <- function(alpha,beta, Trees) {
  LogLike <- ifelse(Trees$var==-1,log(1-alpha*(1+Trees$depth)^(-beta)),log(alpha*(1+Trees$depth)^(-beta)))
  S <- FiniteSum(LogLike)
  return(-S)
}

MCMC_convergence <- function(Samples){
  if (!require(posterior, quietly = T)) {stop("Package posterior not installed")}
  if (class(Samples)[1] != "matrix") {stop("Samples should be specified as matrix with nsample rows and nchain columns")}
  if(nrow(Samples)<ncol(Samples)){print("Are you sure Samples is specified by nsample rows and nchain columns")}
  Rhat = rhat(Samples)
  return(Rhat)
}
