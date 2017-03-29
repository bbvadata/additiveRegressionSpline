rk <- function(x,z) {
  ((z-0.5)^2 - 1/12)*((x-0.5)^2 - 1/12)/4 -
    ((abs(x-z)-0.5)^4-(abs(x-z)-0.5)^2/2 + 7/240) / 24
}

spl.X <- function(x,knots){
  q <- length(knots) + 2                # number of parameters
  n <- length(x)                        # number of observations
  X <- matrix(1, n, q)                  # initialized model matrix
  X[,2] <- x                            # set second column to x
  X[,3:q] <- outer(x, knots, FUN = rk)  # remaining to cubic spline
  X
}
prs.fit <- function(y, x, knots, lambda) {
  q = length(knots) + 2       # dimension of basis
  n = length(x)               # number of observations
  Xa = rbind(spl.X(x, knots), mat.sqrt(spl.S(knots))*sqrt(lambda))  # augmented model matrix
  y[(n+1):(n+q)] = 0          # augment the data vector
  lm(y ~ Xa - 1)              # fit and return penalized regression spline
}
spl.S <- function(knots) {
  q = length(knots) + 2
  S = matrix(0, q, q)                           # initialize matrix
  S[3:q, 3:q] = outer(knots, knots, FUN=rk)     # fill in non-zero part
  S
}
mat.sqrt <- function(S){
  d = eigen(S, symmetric=T)
  rS = d$vectors %*% diag(d$values^.5) %*% t(d$vectors)
  rS
}
GVC <- function(y,x,knots) {
  n <- length(x)
  # GCV (R)
  mod <- prs.fit(y, x, knots, lambda = 0.0001)
  trA <- sum(hatvalues(mod)[1:n])
  trA
  rss <- sum((y- fitted(mod)[1:n])^2)
  rss
  
  (length(x)*rss)/(length(x)-trA)^2
}
trunc <- function(x, ..., prec = 0) base::trunc(x * 10^prec, ...) / 10^prec

am.setup<-function(x,z,q=10)
  # Get X, S_1 and S_2 for a simple 2 term AM
{ # choose knots ...
  xk <- quantile(unique(x),1:(q-2)/(q-1))
  zk <- quantile(unique(z),1:(q-2)/(q-1))
  # get penalty matrices ...
  S <- list()
  S[[1]] <- S[[2]] <- matrix(0,2*q-1,2*q-1)
  S[[1]][2:q,2:q] <- spl.S(xk)[-1,-1]
  S[[2]][(q+1):(2*q-1),(q+1):(2*q-1)] <- spl.S(zk)[-1,-1]
  # get model matrix ...
  n<-length(x)
  X<-matrix(1,n,2*q-1)
  X[,2:q]<-spl.X(x,xk)[,-1]           # 1st smooth
  X[,(q+1):(2*q-1)]<-spl.X(z,zk)[,-1] # 2nd smooth
  list(X=X,S=S)
}
fit.am<-function(y,X,S,sp)
  # function to fit simple 2 term additive model
{ # get sqrt of total penalty matrix ...
  rS <- mat.sqrt(sp[1]*S[[1]]+sp[2]*S[[2]])
  q.tot <- ncol(X)                # number of params
  n <- nrow(X)                    # number of data
  X1 <- rbind(X,rS)               # augmented X
  y1 <- c(y,rep(0,q.tot))         # augment data
  b<-lm(y1~X1-1)                  # fit model
  trA<-sum(influence(b)$hat[1:n]) # tr(A)
  norm<-sum((y-fitted(b)[1:n])^2) # RSS
  list(model=b,gcv=norm*n/(n-trA)^2,sp=sp)
}


