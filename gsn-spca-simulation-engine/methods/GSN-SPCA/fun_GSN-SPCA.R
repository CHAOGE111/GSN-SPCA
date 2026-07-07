library(reticulate)
use_python("D:/anaconda3/python.exe", required = TRUE)

# ========================================================================
#  ESPCA Function (Main Entry)
#  Performs Eigen Sparse Principal Component Analysis with group overlap
# ========================================================================
ESPCA = function(X, k=2, overlap.group, k.group=2, we=0.5, t = 0.1, niter=5, err=0.0001, Num.init=5){
  # X: input data matrix (n x p), n = samples, p = features
  n = nrow(X) 
  p = ncol(X)
  
  # Initialize matrices
  U = matrix(0, n, k) 
  D = matrix(0, k, k) 
  V = matrix(0, p, k)
  
  tX = X
  
  # First rank-1 decomposition
  out = rank1.ESPCA(tX, overlap.group, k.group, we, t, niter, err, Num.init)
  U[,1] = out$u; V[,1] = out$v; D[1,1] = out$d 
  
  if(k < 2) return (list(U=U, D=D, V=V))
  
  # Subsequent components
  for(i in 2:k){
    # Deflate previous component
    tX = tX - c(out$d) * out$u %*% t(out$v)
    UU = U %*% t(U)
    
    # Run refinement cycle
    out = cycleFun2(tX, UU, overlap.group, k.group, we, t, niter, err, Num.init)
    U[,i] = out$u; V[,i] = out$v; D[i,i] = out$d
  }
  return (list(U=U, D=D, V=V))
}


# ========================================================================
#  Rank-1 ESPCA (Base Solver)
#  Finds one sparse eigenvector pair (u,v) under group constraints
# ========================================================================
rank1.ESPCA = function(X, overlap.group, k.group, we, t, niter=1000, err=0.0001, Num.init=5){
  n = nrow(X)
  p = ncol(X)
  d.opt = -100
  
  # Try multiple random initializations
  for(ii in 1:Num.init){
    we1 = we
    set.seed(ii*100)
    
    v0 = matrix(rnorm(p,0,1), ncol=1); v0 = v0 / norm(v0, 'E')
    u0 = matrix(rnorm(n,0,1), ncol=1); u0 = u0 / norm(u0, 'E')
    
    # Iterative update of u and v
    for (i in 1:niter){
      u = u.project2(X %*% v0)
      v = overlap.group.penalty(t(X) %*% u, overlap.group, k.group, we1)
      
      # Decrease penalty weight
      we1 = max(0, we1 - t)
      
      # Stopping condition: convergence of u and v
      if ((norm(u - u0, 'E') <= err) & (norm(v - v0, 'E') <= err)) break
      u0 = u; v0 = v
    }
    
    # Evaluate decomposition strength
    d = t(u) %*% X %*% v
    if(d > d.opt){
      d.opt = d; u.opt = u; v.opt = v
    }
  }
  return (list(u=u.opt, v=v.opt, d=d.opt))
}


# ========================================================================
#  Cycle Function (Deflation Step)
#  Used for computing subsequent components after the first
# ========================================================================
cycleFun2 = function(X, UU, overlap.group, k.group, we, t, niter=1000, err=0.0001, Num.init=5){
  n = nrow(X)
  p = ncol(X)
  d.opt = -100
  
  for(ii in 1:Num.init){
    we1 = we
    set.seed(ii*100)
    
    v0 = matrix(rnorm(p,0,1), ncol=1); v0 = v0 / norm(v0,'E')
    u0 = matrix(rnorm(n,0,1), ncol=1); u0 = u0 / norm(u0,'E')
    
    # Iterative update of u and v
    for(i in 1:niter){
      # Project orthogonal to previous components
      u = (diag(n) - UU) %*% (X %*% v0)
      u = u.project2(u)
      
      # Apply group overlap penalty
      v = overlap.group.penalty(t(X) %*% u, overlap.group, k.group, we1)
      
      # Decrease penalty weight
      we1 = max(0, we1 - t)
      
      # Stopping condition
      if ((norm(u - u0,'E') <= err) & (norm(v - v0,'E') <= err)) break
      u0 = u; v0 = v
    }
    
    d = t(u) %*% X %*% v
    if(d > d.opt){
      d.opt = d; u.opt = u; v.opt = v
    }
  }
  return (list(u=u.opt, v=v.opt, d=d.opt))
}


# ========================================================================
#  Projection of vector u onto L2 unit sphere
# ========================================================================
u.project2 = function(z){  
  if(sum(z^2) == 0){return(rep(0, length(z)))}
  else return(z / sqrt(sum(z^2)))
}


# ========================================================================
#  Overlap Group Penalty Function
#  Applies group-sparse projection with optional weighting from external file
# ========================================================================
overlap.group.penalty = function(u, overlap.group, k0=1.0, we=0.5){
  group.num = length(overlap.group)
  group.norm = rep(0, group.num)
  
  # Compute group-wise norms (weighted)
  for(i in 1:group.num){
    g.set = overlap.group[[i]]
    g.set_clean <- g.set[!is.na(g.set)]
    w_i = 1 / sqrt(length(g.set_clean))
    group.norm[i] = norm(w_i * u[g.set_clean], "2")
  }
  
  # If external weights file exists, apply adjustment
  file_path <- "result_max_values.txt"
  if (file.exists(file_path)) {
    file_content <- readLines(file_path)
    if (length(file_content) == length(group.norm)) {
      for (j in 1:length(file_content)) {
        group.norm[j] <- group.norm[j] * as.numeric(file_content[j])
      }
    } else {
      stop("Row count of file does not match length of group.norm.")
    }
  } else {
    print("Weights file not found, skipping adjustment.")
  }
  
  # Determine number of groups to keep
  if(we > 0){
    k1 = ceiling(k0 * (1 + we))
  } else {
    k1 = k0
  }
  
  # Select top-k groups
  if(we != 0){
    ID1 = order(group.norm, decreasing=TRUE)[1:k1]
    ID  = sample(ID1, k0)
  } else {
    ID = order(group.norm, decreasing=TRUE)[1:k0]
  }
  
  # Collect selected features
  select.features = c()
  for(i in 1:length(ID)){
    temp = overlap.group[[ID[i]]]
    temp_clean <- temp[!is.na(temp)]
    select.features = c(select.features, temp_clean)
  }
  
  index = sort(unique(na.omit(select.features)))
  
  # Zero out unselected features
  x.opt = u; x.opt[-c(index)] = 0
  
  # Save selected feature indices for inspection
  write.csv(x.opt, file="x.opt.csv")
  opt <- read.csv("x.opt.csv", header=TRUE)
  a = opt[,1][opt[,2] != 0]
  write.csv(a, file="a.csv")
  
  # Normalize the output vector
  if(sum(abs(x.opt)) == 0) return(x.opt)
  else return(x.opt / norm(x.opt, "E"))
}
