Posdef <- function (n, ev = runif(n, 0, 10)) 
{
  Z <- matrix(ncol=n, rnorm(n^2))
  decomp <- qr(Z)
  Q <- qr.Q(decomp) 
  R <- qr.R(decomp)
  d <- diag(R)
  ph <- d / abs(d)
  O <- Q %*% diag(ph)
  Z <- t(O) %*% diag(ev) %*% O
  return(Z)
}

test = Posdef(4)
st = test/test[1,1]
Chol = chol(st)
ins = Chol[upper.tri(Chol,diag = T)][-1]
Mult = t(Chol)%*%Chol