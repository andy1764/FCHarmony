#' Common principal components using Flury's (1986) estimation
#'
#' Adapted from the FCPCA function in \emph{multigroup}
#'
#' @param dat an array of covariance matrices of dimension \eqn{(p \times p \times N)}
#' @param iter numeric, number of iterations
#'
#' @return \code{cpc_Flury} returns a list containing the following components:
#' \item{D}{eigenvalues for each input matrix as columns}
#' \item{CPC}{rotation matrix with eigenvectors in columns}
#' \item{ncomp}{number of components}
#' @export
#'
#' @examples

cpc_Flury <- function(dat, iter = 15) {
  N <- dim(dat)[3]
  P <- dim(dat)[2]
  L <- iter

  k <- N
  K <- N
  B <- diag(P)
  Bold <- diag(P)
  C <- dat
  T <- vector("list",k)
  d1 <- c(rep(0,k))
  d2 <- c(rep(0,k))
  for(l in 1:L){
    for(p in 1:(P-1)){
      for(e in (p+1):P){
        Q <- diag(2)
        #------ G step
        M <- diag(2)
        for(k in 1:K){
          H <- B[,c(p,e)]
          T[[k]] <- t(H)%*%C[,,k]%*%H
          d1[k] <- t(Q[,1])%*%T[[k]]%*%Q[,1]
          d2[k] <- t(Q[,2])%*%T[[k]]%*%Q[,2]
          M <-  M + (d1[k]-d2[k])/(d1[k]*d2[k])*T[[k]]
        }

        eig <- eigen(M)
        Q <- eig$vectors
        B[,c(p,e)] <- H%*%Q
      }
    }
  }
  W <- B

  # get eigenvalues
  lambda <- matrix(0, nrow=P, ncol=N)
  for(i in 1:N){
    lambda[,i] <- round(diag(t(W) %*% C[,,i] %*% W),3)
  }

  list(D = lambda,
       CPC = W,
       ncomp = P)
}
