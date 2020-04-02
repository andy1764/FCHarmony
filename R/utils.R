# convert covariance to Laplacian with a given threshold,
cov2lap <- function(cov, threshold = 0.5, gamma = 0.01) {
  corr <- cov2cor(as.matrix(nearPD(cov)$mat))
  corr[upper.tri(corr)] <- as.numeric(corr[upper.tri(corr)] >= threshold)
  corr[lower.tri(corr)] <- t(corr)[lower.tri(corr)]
  diag(corr) <- 0
  corr <- -corr
  # gamma addition ensures positive definite
  diag(corr) <- -apply(corr, 1, sum) + gamma
  return(corr)
}

# get matrix logarithm of SPD matrix
logm_eig <- function(x) {
  eig <- eigen(x, symmetric = TRUE)
  if (any(eig$values <= 0)) {stop("Input is not positive definite")}
  eig$vectors %*% diag(log(eig$values)) %*% t(eig$vectors)
}

# get local efficiency, based on code from brainGraph but fixed
local_eff <- function (g, weights = NULL, use.parallel = TRUE, A = NULL) {
  if (is.null(weights)) {
    if (is.null(A)) {
      A <- as_adj(g, names = FALSE, attr = "weight")
      weighted <- TRUE
    }
  } else {
    A <- as_adj(g, names = FALSE, sparse = FALSE)
    weighted <- NULL
  }
  eff <- rep(0, nrow(A))
  nodes <- which(rowSums((A > 0) + 0) > 1)
  X <- apply(A, 1, function(x) which(x > 0))
  if (length(nodes) > 0) {
    if (isTRUE(use.parallel)) {
      eff[nodes] <- foreach(i = nodes, .combine = "c") %dopar%
        {
          # originally used A[X[[i]], X[[i]]] which is a numeric value and
          # clearly incorrect, it should be getting the subgraph not including
          # the node of interest
          g.sub <- graph_from_adjacency_matrix(A[-X[[i]],
                                                 -X[[i]]], mode = "undirected", weighted = weighted)
          efficiency(g.sub, "global", weights = weights)
        }
    }
    else {
      for (i in nodes) {
        g.sub <- graph_from_adjacency_matrix(A[-X[[i]],
                                               -X[[i]]], mode = "undirected", weighted = weighted)
        eff[i] <- efficiency(g.sub, "global", weights = weights)
      }
    }
  }
  eff
}

# MatReg_QC_opt <- function(Y,X,method=c("CAP","CAP-C","CAP-C1"),max.itr=1000,tol=1e-4,trace=FALSE,score.return=TRUE,gamma0.mat=NULL,ninitial=NULL) {
#   n<-length(Y)
#   p<-ncol(Y[[1]])
#   Tvec<-rep(NA,n)
#
#   q<-ncol(X)
#
#   # Estimate covariance matrix for each subject
#   Sigma<-array(NA,c(p,p,n))
#   for(i in 1:n)
#   {
#     Tvec[i]<-nrow(Y[[i]])
#
#     Sigma[,,i]<-t(scale(Y[[i]],center=TRUE,scale=FALSE))%*%(scale(Y[[i]],center=TRUE,scale=FALSE))/nrow(Y[[i]])
#   }
#
#   # common PCA based method
#   if(method[1]=="CAP-C")
#   {
#     # find common eigenvectors and subject-specific eigenvalues
#     Ymat<-NULL
#     Group<-NULL
#     for(i in 1:n)
#     {
#       Tvec[i]<-nrow(Y[[i]])
#
#       Ymat<-rbind(Ymat,Y[[i]])
#
#       Group<-c(Group,rep(i,Tvec[i]))
#     }
#     re.FCPCA<-FCPCA(Data=Ymat,Group=Group)
#
#     # eigenvalues
#     lambda<-re.FCPCA$lambda
#     # common eigenvectors
#     phi<-re.FCPCA$loadings.common
#   }
#
#   # set initial values
#   if(is.null(gamma0.mat))
#   {
#     gamma0.mat<-matrix(NA,p,p+1+5)
#     for(j in 1:p)
#     {
#       gamma0.mat[,j]<-rep(0,p)
#       gamma0.mat[j,j]<-1
#     }
#     gamma0.mat[,p+1]<-rep(1,p)/sqrt(sum(rep(1,p)^2))
#
#     set.seed(500)
#     gamma.tmp<-matrix(rnorm(5*p,mean=0,sd=1),nrow=p)
#     gamma0.mat[,(p+2):(p+1+5)]<-apply(gamma.tmp,2,function(x){return(x/sqrt(sum(x^2)))})
#   }
#   if(is.null(ninitial))
#   {
#     ninitial<-min(ncol(gamma0.mat),10)
#   }else
#   {
#     if(ninitial>ncol(gamma0.mat))
#     {
#       ninitial<-ncol(gamma0.mat)
#     }
#   }
#   set.seed(500)
#   gamma0.mat<-matrix(gamma0.mat[,sort(sample(1:ncol(gamma0.mat),ninitial,replace=FALSE))],ncol=ninitial)
#
#   if(method[1]=="CAP-C1")
#   {
#     gamma0.mat<-apply(gamma0.mat,2,function(x){return(x/sqrt(sum(x^2)))})
#
#     re<-vector("list",ncol(gamma0.mat))
#     obj.func<-rep(NA,ncol(gamma0.mat))
#     for(kk in 1:ncol(gamma0.mat))
#     {
#       try(re[[kk]]<-MatReg_QC(Y,X,method=method[1],max.itr=max.itr,tol=tol,trace=trace,gamma0=gamma0.mat[,kk],score.return=score.return))
#
#       try(obj.func[kk]<-objfunc(Y,X,re[[kk]]$gamma,re[[kk]]$beta))
#     }
#
#     opt.idx<-which.min(obj.func)
#     re.opt<-re[[opt.idx]]
#
#     if(method[1]=="CAP-C")
#     {
#       dis<-rep(NA,1,p)
#       for(j in 1:p)
#       {
#         if(phi[1,j]<0)
#         {
#           phi.tmp<--phi[,j]
#         }else
#         {
#           phi.tmp<-phi[,j]
#         }
#         dis[j]<-sqrt(sum(re.opt$gamma-phi.tmp)^2)
#       }
#       re.opt$PC.idx<-which.min(dis)
#     }
#
#     return(re.opt)
#   }
#   if(method[1]=="CAP")
#   {
#     theta0.mat<-gamma0.mat
#
#     re<-vector("list",ncol(gamma0.mat))
#     obj.func<-rep(NA,ncol(gamma0.mat))
#
#     re.scale<-vector("list",ncol(gamma0.mat))
#     obj.func.scale<-rep(NA,ncol(gamma0.mat))
#
#     for(kk in 1:ncol(gamma0.mat))
#     {
#       try(re[[kk]]<-MatReg_QC(Y,X,method="CAP",max.itr=max.itr,tol=tol,trace=trace,gamma0=theta0.mat[,kk],score.return=score.return))
#
#       try(obj.func[kk]<-objfunc(Y,X,re[[kk]]$gamma,re[[kk]]$beta))
#
#       try(re.scale[[kk]]<-MatReg_QC_beta(Y,X,gamma=re[[kk]]$gamma/sqrt(sum((re[[kk]]$gamma)^2)),max.itr=max.itr,tol=tol,trace=trace,score.return=score.return))
#       try(obj.func.scale[kk]<-objfunc(Y,X,re.scale[[kk]]$gamma,re.scale[[kk]]$beta))
#     }
#
#     opt.idx<-which.min(obj.func)
#     opt.idx.scale<-which.min(obj.func.scale)
#
#     # re.opt<-list(unscale=re[[opt.idx]],scale=re.scale[[opt.idx]])
#     re.opt<-re.scale[[opt.idx]]
#
#     return(re.opt)
#   }
#   if(method[1]=="CAP-C")
#   {
#     optmat<-matrix(NA,p,p)
#     colnames(optmat)<-paste0("Dim",1:p)
#     rownames(optmat)<-paste0("BetaDim",1:p)
#     beta.est<-matrix(NA,q,p)
#     colnames(beta.est)<-paste0("Dim",1:p)
#     rownames(beta.est)<-colnames(X)
#     for(j in 1:p)
#     {
#       beta.tmp<-MatReg_QC_beta(Y,X,gamma=phi[,j])$beta
#       optmat[j,]<-(apply(lambda,2,function(x){return(sum(x*Tvec*exp(-X%*%beta.tmp)))})/apply(lambda,2,sum))*n/2+sum(Tvec*X%*%beta.tmp)/2
#
#       beta.est[,j]<-beta.tmp
#     }
#
#     min.idx<-apply(optmat,1,which.min)
#     sidx<-which(apply(cbind(min.idx,1:p),1,function(x){x[1]==x[2]})==TRUE)
#     if(length(sidx)>0)
#     {
#       svar<-rep(NA,length(sidx))
#       for(j in 1:length(sidx))
#       {
#         svar[j]<-sum(exp(X%*%beta.est[,sidx[j]]))
#       }
#       opt.idx<-sidx[which.max(svar)]
#       beta.opt<-beta.est[,opt.idx]
#       gamma.opt<-phi[,opt.idx]
#       if(gamma.opt[1]<0)
#       {
#         gamma.opt<--gamma.opt
#       }
#
#       if(score.return)
#       {
#         re.opt<-list(gamma=gamma.opt,beta=beta.opt,PC.index=opt.idx,score=c(lambda[,opt.idx]),minmix=TRUE)
#       }else
#       {
#         re.opt<-list(gamma=gamma.opt,beta=beta.opt,PC.index=opt.idx,minmix=TRUE)
#       }
#     }else
#     {
#       optvec<-rep(NA,p)
#       for(j in 1:p)
#       {
#         optvec[j]<-optmat[j,min.idx[j]]
#       }
#       opt.idx<-which.min(optvec)
#       beta.opt<-beta.est[,opt.idx]
#       gamma.opt<-phi[,opt.idx]
#       if(gamma.opt[1]<0)
#       {
#         gamma.opt<--gamma.opt
#       }
#
#       if(score.return)
#       {
#         re.opt<-list(gamma=gamma.opt,beta=beta.opt,PC.index=opt.idx,score=c(lambda[,opt.idx]),minmix=FALSE)
#       }else
#       {
#         re.opt<-list(gamma=gamma.opt,beta=beta.opt,PC.index=opt.idx,minmix=FALSE)
#       }
#     }
#
#     return(re.opt)
#   }
# }
#
# MatReg_QC_opt2 <- function(Y,X,Phi0=NULL,method=c("CAP","CAP-C","CAP-C1"),CAP.OC=FALSE,max.itr=1000,tol=1e-4,trace=FALSE,score.return=TRUE,gamma0.mat=NULL,ninitial=NULL) {
#   if(is.null(Phi0))
#   {
#     return(MatReg_QC_opt(Y,X,method=method,max.itr=max.itr,tol=tol,trace=trace,score.return=score.return,gamma0.mat=gamma0.mat,ninitial=ninitial))
#   }else
#   {
#     n<-length(Y)
#     q<-ncol(X)
#     p<-ncol(Y[[1]])
#
#     p0<-ncol(Phi0)
#     # estimate beta0
#     beta0<-rep(NA,p0)
#     for(j in 1:p0)
#     {
#       beta0[j]<-MatReg_QC_beta(Y,X,gamma=Phi0[,j],max.itr=max.itr,tol=tol,trace=FALSE)$beta[1]
#     }
#     Ytmp<-vector("list",length=n)
#     for(i in 1:n)
#     {
#       Y2tmp<-Y[[i]]-Y[[i]]%*%(Phi0%*%t(Phi0))
#
#       Y2tmp.svd<-svd(Y2tmp)
#
#       Ytmp[[i]]<-Y2tmp.svd$u%*%diag(c(Y2tmp.svd$d[1:(p-p0)],exp(beta0)))%*%t(Y2tmp.svd$v)
#
#       # Y2tmp0<-Y2tmp.svd$u%*%diag(c(rep(0,p-p0),exp(beta0)))%*%t(Y2tmp.svd$v)
#       # Ytmp[[i]]<-Y2tmp+Y2tmp0
#     }
#
#     if(method=="CAP")
#     {
#       if(CAP.OC==FALSE)
#       {
#         re.tmp<-MatReg_QC_opt(Ytmp,X,method=method,max.itr=max.itr,tol=tol,trace=FALSE,score.return=score.return,gamma0.mat=gamma0.mat,ninitial=ninitial)
#       }else
#       {
#         re.tmp<-MatReg_QC_RE(Ytmp,X,Phi0=Phi0,max.itr=max.itr,tol=tol,trace=FALSE,score.return=score.return)
#       }
#     }else
#     {
#       re.tmp<-MatReg_QC_opt(Ytmp,X,method=method,max.itr=max.itr,tol=tol,trace=FALSE,score.return=score.return,gamma0.mat=gamma0.mat,ninitial=ninitial)
#     }
#
#     re<-re.tmp
#     re$orthogonal<-c(t(re.tmp$gamma)%*%Phi0)
#
#     return(re)
#   }
# }
#
# MatReg_QC_RE <- function(Y,X,Phi0,max.itr=1000,tol=1e-4,trace=FALSE,score.return=TRUE) {
#   n<-length(Y)
#   p<-ncol(Y[[1]])
#   Tvec<-rep(NA,n)
#
#   q<-ncol(X)
#
#   # Estimate covariance matrix for each subject
#   Sigma<-array(NA,c(p,p,n))
#   for(i in 1:n)
#   {
#     Tvec[i]<-nrow(Y[[i]])
#
#     Sigma[,,i]<-t(scale(Y[[i]],center=TRUE,scale=FALSE))%*%(scale(Y[[i]],center=TRUE,scale=FALSE))/nrow(Y[[i]])
#   }
#
#   beta0<-rep(0,q)
#
#   if(trace)
#   {
#     gamma.trace<-NULL
#     beta.trace<-beta0
#   }
#
#   s<-0
#   diff<-100
#   while(s<=max.itr&diff>tol)
#   {
#     s<-s+1
#
#     # update gamma
#     A<-matrix(0,p,p)
#     for(i in 1:n)
#     {
#       A<-A+Tvec[i]*Sigma[,,i]/(exp(X[i,]%*%beta0)[1,1])
#     }
#     B<-apply(Sigma,c(1,2),mean)
#     B.inv<-solve(B)
#     P<-Phi0%*%ginv(t(Phi0)%*%B.inv%*%Phi0)%*%t(Phi0)%*%B.inv
#     B.inv.eigen<-eigen(B.inv)
#     B.inv.rt<-B.inv.eigen$vectors%*%diag(sqrt(B.inv.eigen$values))%*%solve(B.inv.eigen$vectors)
#     re.eigen<-eigen(B.inv.rt%*%(diag(rep(1,p))-P)%*%A%*%B.inv.rt)
#     re.eigen.vec<-Re(re.eigen$vectors)
#
#     x.tmp<-B.inv.rt%*%re.eigen.vec
#     # diag(t(x.tmp)%*%B%*%x.tmp)
#
#     # which.min(diag(t(x.tmp)%*%A%*%x.tmp))
#     # gamma.new<-x.tmp[,which.min(diag(t(x.tmp)%*%A%*%x.tmp))]
#
#     gamma.new.mat<-x.tmp
#
#     # update beta
#     # Q1<-t(X)%*%diag(Tvec)%*%X
#     beta.new.mat<-NULL
#     objfunc.tmp<-NULL
#     for(j in 1:ncol(gamma.new.mat))
#     {
#       gamma.new<-gamma.new.mat[,j]
#
#       Q1<-matrix(0,q,q)
#       Q2<-rep(0,q)
#       for(i in 1:n)
#       {
#         # Q1<-Q1+Tvec[i]*X[i,]%*%t(X[i,])
#
#         Q1<-Q1+(Tvec[i]*(t(gamma.new)%*%Sigma[,,i]%*%gamma.new)[1,1]*exp(-t(X[i,])%*%beta0)[1,1])*(X[i,]%*%t(X[i,]))
#
#         Q2<-Q2+Tvec[i]*(1-(t(gamma.new)%*%Sigma[,,i]%*%gamma.new)[1,1]*(exp(-t(X[i,])%*%beta0)[1,1]))*X[i,]
#       }
#       beta.new.mat<-cbind(beta.new.mat,beta0-ginv(Q1)%*%Q2)
#
#       objfunc.tmp<-c(objfunc.tmp,objfunc(Y,X,gamma.new,beta.new.mat[,j]))
#     }
#     beta.new<-beta.new.mat[,which.min(objfunc.tmp)]
#     gamma.new<-gamma.new.mat[,which.min(objfunc.tmp)]
#
#     if(trace)
#     {
#       gamma.trace<-cbind(gamma.trace,gamma.new)
#       beta.trace<-cbind(beta.trace,beta.new)
#     }
#
#     diff<-max(abs(beta.new-beta0))
#
#     beta0<-beta.new
#
#     # print(diff)
#   }
#
#   # scale gamma
#   gamma.new<-gamma.new/sqrt(sum(gamma.new^2))
#   if(gamma.new[1]<0)
#   {
#     gamma.new<--gamma.new
#   }
#   beta.new<-MatReg_QC_beta(Y,X,gamma=gamma.new,max.itr=max.itr,tol=tol,trace=trace,score.return=score.return)$beta
#
#   if(score.return)
#   {
#     score<-rep(NA,n)
#     for(i in 1:n)
#     {
#       score[i]<-t(gamma.new)%*%Sigma[,,i]%*%gamma.new
#     }
#   }
#
#   if(trace)
#   {
#     colnames(v.trace)<-NULL
#     colnames(beta.trace)<-NULL
#
#     if(score.return)
#     {
#       re<-list(gamma=c(gamma.new),beta=c(beta.new),convergence=(s<max.itr),score=score,gamma.trace=gamma.trace,beta.trace=beta.trace)
#     }else
#     {
#       re<-list(gamma=c(gamma.new),beta=c(beta.new),convergence=(s<max.itr),gamma.trace=gamma.trace,beta.trace=beta.trace)
#     }
#
#   }else
#   {
#     if(score.return)
#     {
#       re<-list(gamma=c(gamma.new),beta=c(beta.new),convergence=(s<max.itr),score=score)
#     }else
#     {
#       re<-list(gamma=c(gamma.new),beta=c(beta.new),convergence=(s<max.itr))
#     }
#   }
#
#   return(re)
# }
#
# MatReg_QC_beta <- function(Y,X,gamma,max.itr=1000,tol=1e-4,trace=FALSE,score.return=TRUE) {
#   n<-length(Y)
#   p<-ncol(Y[[1]])
#   Tvec<-rep(NA,n)
#
#   q<-ncol(X)
#
#   # Estimate covariance matrix for each subject
#   Sigma<-array(NA,c(p,p,n))
#   for(i in 1:n)
#   {
#     Tvec[i]<-nrow(Y[[i]])
#
#     Sigma[,,i]<-t(scale(Y[[i]],center=TRUE,scale=FALSE))%*%(scale(Y[[i]],center=TRUE,scale=FALSE))/nrow(Y[[i]])
#   }
#
#   if(score.return)
#   {
#     score<-rep(NA,n)
#     for(i in 1:n)
#     {
#       score[i]<-t(gamma)%*%Sigma[,,i]%*%gamma
#     }
#   }
#
#   beta0<-rep(0,q)
#
#   if(trace)
#   {
#     beta.trace<-beta0
#
#     obj<-objfunc(Y=Y,X=X,gamma=gamma,beta=beta0)
#   }
#
#   s<-0
#   diff<-100
#   while(s<=max.itr&diff>tol)
#   {
#     s<-s+1
#
#     # update beta
#     # Q1<-t(X)%*%diag(Tvec)%*%X
#     Q1<-matrix(0,q,q)
#     Q2<-rep(0,q)
#     for(i in 1:n)
#     {
#       # Q1<-Q1+Tvec[i]*X[i,]%*%t(X[i,])
#
#       Q1<-Q1+(Tvec[i]*(t(gamma)%*%Sigma[,,i]%*%gamma)[1,1]*exp(-t(X[i,])%*%beta0)[1,1])*(X[i,]%*%t(X[i,]))
#
#       Q2<-Q2+Tvec[i]*(1-(t(gamma)%*%Sigma[,,i]%*%gamma)[1,1]*(exp(-t(X[i,])%*%beta0)[1,1]))*X[i,]
#     }
#     # beta.new<-beta0-solve(Q1)%*%Q2
#     beta.new<-beta0-ginv(Q1)%*%Q2
#
#     if(trace)
#     {
#       beta.trace<-cbind(beta.trace,beta.new)
#
#       obj<-c(obj,objfunc(Y=Y,X=X,gamma=gamma,beta=beta.new))
#     }
#
#     diff<-max(abs(beta.new-beta0))
#
#     beta0<-beta.new
#
#     # print(diff)
#   }
#
#   if(trace)
#   {
#     if(score.return)
#     {
#       re<-list(gamma=c(gamma),beta=c(beta.new),convergence=(s<max.itr),score=score,beta.trace=beta.trace,obj=obj)
#     }else
#     {
#       re<-list(gamma=c(gamma),beta=c(beta.new),convergence=(s<max.itr),beta.trace=beta.trace,obj=obj)
#     }
#   }else
#   {
#     if(score.return)
#     {
#       re<-list(gamma=c(gamma),beta=c(beta.new),convergence=(s<max.itr),score=score)
#     }else
#     {
#       re<-list(gamma=c(gamma),beta=c(beta.new),convergence=(s<max.itr))
#     }
#   }
#
#   return(re)
# }
