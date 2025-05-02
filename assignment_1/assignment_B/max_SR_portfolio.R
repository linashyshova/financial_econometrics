max_SR_portfolio <- function(er,cov.mat){
    er <- as.vector(er)
    cov.mat <- as.matrix(cov.mat)
    N <- length(er)
    Dmat <- 2*cov.mat
    dvec <- rep.int(0, N)
    Amat <- cbind(er, diag(1,N))
    bvec <- c(1, rep(0,N))
    result <- quadprog::solve.QP(Dmat=Dmat,dvec=dvec,Amat=Amat,bvec=bvec,meq=1)
    w.t <- round(result$solution/sum(result$solution), 6)
    return(w.t)
}