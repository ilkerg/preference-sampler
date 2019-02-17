k1 <- 2
k2 <- 1
k3 <- 5
C12 <- 4
C13 <- 5
C23 <- 2

U <- function(q) -((k1+.1)*log(q[1]) + (k2+.1)*log(q[2]) + (k3+.1)*log(1-q[1]-q[2])
                   - C12*log(q[1]+q[2]) - C13*log(1-q[2]) - C23*log1p(1-q[1]))
K <- function(p) .5*sum(p^2)
E <- function(p, q) K(p) + U(q)

dU <- function(q) {
    dU_1 <- -(k1+.1)/q[1] + (k3+.1)/(1-q[1]-q[2]) + C12/(q[1]+q[2]) - C23/(1-q[1])
    dU_2 <- -(k2+.1)/q[2] + (k3+.1)/(1-q[1]-q[2]) + C12/(q[1]+q[2]) - C13/(1-q[2])
    c(dU_1, dU_2)
}

leapfrog <- function(p, q, grad, eps=.01, L=10) {

    for (i in 1:L) {
        p <- p - .5*eps*grad(q)
        q <- q + eps*p
        p <- p - .5*eps*grad(q)
    }

    list(p=p, q=q)
}

adaptive_lf <- function(p, q, grad, L=25) {
    eps = 1e-1 / max(abs(grad(q)))
    cat(sprintf("eps = %f\n", eps))
    leapfrog(p, q, grad, eps, L)
}

HMC <- function(q, grad) {
    p <- rnorm(length(q))

    lf <- adaptive_lf(p, q, grad, L=runif(1, min=10, max=15))
    if(any(lf$q<0))
        return(q)
    if(sum(lf$q)>1)
        return(q)

    E1 <- E(p, q)
    cat(sprintf("E1 = %f\n", E1))
    E2 <- E(lf$p, lf$q)
    cat(sprintf("E2 = %f\n", E2))

    alpha <- exp(E2-E1)
    cat(sprintf("alpha = %f\n", alpha))

    if(runif(1) < alpha) {
        cat("accept\n")
        lf$q
    }
    else {
        cat("reject\n")
        q
    }
}

cnt <- function() {
    f <- function(th1, th2) -((k1+.1)*log(th1) + (k2+.1)*log(th2) + (k3+.1)*log1p(-th1-th2)
                              - C12*log(th1 + th2) - C13*log1p(-th2) - C23 * log1p(-th1))

    x <- seq(0, 1, length.out=100)
    y <- seq(0, 1, length.out=100)
    z <- outer(x, y, FUN="f")
    z
}


test <- function(N, th0) {
    require(colorspace)

    th <- matrix(rep(c(th0, 1-sum(th0)), N), ncol=3, byrow=T)
    th_new <- t(apply(th[,1:2], 1, function(q) HMC(q, dU)))
    th_new <- cbind(th_new, 1-rowSums(th_new))

    z <- cnt()

    pal <- function(n) sequential_hcl(n, "Purples 3")

    png("1.png", width=1024, height=768)
    filled.contour(z, color=pal, nlevels=30,
                   plot.axis={
                       points(th, pch=19, col='#ff000055',
                              xlim=c(0,1), ylim=c(0,1))
                   })
    dev.off()

    png("2.png", width=1024, height=768)
    filled.contour(z, color=pal, nlevels=30,
                   plot.axis={
                       points(th_new, pch=19, col='#ff000055',
                              xlim=c(0,1), ylim=c(0,1))
                   })
    dev.off()
}

