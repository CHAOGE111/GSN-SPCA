library(reticulate)
# Python environment setting
use_python("D:/anaconda3/python.exe", required = TRUE)

ESPCA = function(X, k=2, overlap.group, k.group=2, we=0.5, t = 0.1, niter=1000, err=0.0001, Num.init=5){
  # --------------------------------------------------------------------------
  # [n,p] = dim(X), n is the number of samples and p is the number of features
  # --------------------------------------------------------------------------
  n = nrow(X); # number of samples
  p = ncol(X); # number of features
  U = matrix(0,n,k); D = matrix(0,k,k); V = matrix(0,p,k)
  tX = X
  out = rank1.ESPCA(tX, overlap.group, k.group, we, t, niter, err, Num.init)
  U[,1] = out$u; V[,1] = out$v; D[1,1] = out$d 
  if(k<2) return (list(U=U, D=D, V=V))
  
  # --------------------------------------------------------------------------
  # Iteratively extract components
  # --------------------------------------------------------------------------
  for(i in 2:k){
    tX = tX-c(out$d)*out$u%*%t(out$v) 
    UU = U%*%t(U)
    out = cycleFun2(tX, UU, overlap.group, k.group,we, t, niter, err, Num.init)
    U[,i] = out$u; V[,i] = out$v; D[i,i] = out$d
  }
  return (list(U=U, D=D, V=V))
}

# ----------------------------------------------------------------------------
rank1.ESPCA = function(X, overlap.group, k.group, we, t, niter=1000, err=0.0001, Num.init = 5){
  n = nrow(X); # number of samples
  p = ncol(X); # number of features
  d.opt = -100
  
  # Try multiple initializations
  for(ii in 1:Num.init){
    we1 = we
    print("rank")
    set.seed(ii*100)
    v0 = matrix(rnorm(p,0,1),ncol=1);v0 = v0/norm(v0,'E')
    u0 = matrix(rnorm(n,0,1),ncol=1);u0 = u0/norm(u0,'E')
    
    # Iterative algorithm to solve for u and v
    for (i in 1:niter){
      u = u.project2(X%*%v0)
      v = overlap.group.penalty(t(X)%*%u, overlap.group, k.group, we1)
      if(we1 > 0){
        we1 = we1 -t
      }else{
        we1 = 0
      }
      # Stopping condition
      if ((norm(u - u0,'E')<= err)&(norm(v - v0,'E')<= err)){break}
      else {
        u0 = u;v0 = v}
    }
    d =t(u)%*%X%*%v
    if(d>d.opt){
      d.opt = d
      u.opt = u
      v.opt = v
    }
  }
  return (list(u=u.opt, v=v.opt, d=d.opt))
}

# ----------------------------------------------------------------------------
cycleFun2 = function(X, UU, overlap.group, k.group,we, t, niter=1000, err=0.0001, Num.init = 5){
  n = nrow(X); # number of samples
  p = ncol(X); # number of features
  d.opt = -100
  
  for(ii in 1:Num.init){
    print("cyc")
    we1 = we
    set.seed(ii*100)
    v0 = matrix(rnorm(p,0,1),ncol=1);v0 = v0/norm(v0,'E')
    u0 = matrix(rnorm(n,0,1),ncol=1);u0 = u0/norm(u0,'E')
    
    # Iterative algorithm to solve for u and v
    for(i in 1:niter){
      u = (diag(n) - UU)%*%(X%*%v0); u = u.project2(u)
      v = overlap.group.penalty(t(X)%*%u, overlap.group, k.group, we1)
      if(we1 > 0){
        we1 = we1 -t
      }else{
        we1 = 0
      }
      # Stopping condition
      if ((norm(u - u0,'E')<= err)&(norm(v - v0,'E')<= err)){break}
      else {
        u0 = u;v0 = v}
    }
    d = t(u)%*%X%*%v
    if(d>d.opt){
      d.opt = d
      u.opt = u
      v.opt = v
    }
  }
  return (list(u=u.opt, v=v.opt, d=d.opt))
}

# ----------------------------------------------------------------------------
u.project2 = function(z){  
  u = z
  if(sum(u^2)==0){return(rep(0,length(u)))}
  else{
    u = u/sqrt(sum(u^2))
    return(u)} 
}

# ----------------------------------------------------------------------------
overlap.group.penalty = function(u, overlap.group, k0 = 1.0, we=0.5){
  # Greedy k-edge-sparse projection
  # k0 : number of variables to select
  # overlap.group : list of feature groups
  # --------------------------------------------------------------------------
  
  group.num = length(overlap.group)
  group.norm = rep(0,group.num)
  
  # Compute group norms
  for(i in 1:group.num){
    g.set = overlap.group[[i]]
    g.set_clean <- g.set[!is.na(g.set)]
    w_i = 1/sqrt(length(g.set_clean))
    group.norm[i] = norm(w_i*u[g.set_clean],"2")
  }
  
  # Read external weights from file if available
  file_path <- "result_max_values.txt"
  if (file.exists(file_path)) {
    file_content <- readLines(file_path)
    if (length(file_content) == length(group.norm)) {
      for (j in 1:length(file_content)) {
        group.norm[j] <- group.norm[j] * as.numeric(file_content[j])
      }
    } else {
      stop("The number of lines in the file does not match group.norm length.")
    }
  } else {
    print("The specified file does not exist. Skipping file reading step (normal for the first cycle).")
  }
  
  # Adjust number of groups selected
  if(we > 0){
    k1 = k0*(1+we)
    k1 = ceiling(k1)
  }else{
    k1 = k0
  }
  if(k1 < k0){
    k1 = k0
  }
  
  if(we != 0){
    ID1 = order(group.norm, decreasing = TRUE)[1:k1]
    ID = sample(ID1, k0)
  }else {
    ID = order(group.norm, decreasing = TRUE)[1:k0]
  }
  
  # Collect selected features
  select.features = c()
  for(i in 1:length(ID)){
    temp = overlap.group[[ID[i]]]
    temp_clean <- temp[!is.na(temp)]
    select.features = c(select.features, temp_clean)
  }
  
  # Final feature index set
  index= sort(unique(na.omit(select.features)))
  x.opt = u; x.opt[-c(index)] = 0
  
  # Save and reload selected features
  write.csv (x.opt, file ="x.opt.csv")
  opt <- read.csv("x.opt.csv", header=TRUE)
  a = opt[,1][opt[,2] != 0]
  write.csv (a, file ="a.csv")
  
  # Call external python function
  source_python("class1.py")
  
  if(sum(abs(x.opt))==0) return(x.opt)
  else return(x.opt/norm(x.opt,"E")) 
}
