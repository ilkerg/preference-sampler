#!/usr/bin/env Rscript

library(colorspace)

k1 <- 7
k2 <- 3
k3 <- 5
C12 <- 6
C13 <- 6
C23 <- 1


games <- matrix(c(1,2,
                  1,3,
                  1,2,
                  2,3,
                  2,1,
                  3,1), ncol=2, byrow=T)
winners <- c(1, 1, 2, 3, 2, 1)



# V = -log(posterior)
V <- function(th1, th2) {
    -(k1*log(th1) + k2*log(th2) + k3*log(1-th1-th2) - C12*log(th1+th2) - C13*log(1-th2) - C23*log(1-th1))
}


# posterior
p <- function(th1, th2) exp(-V(th1, th2))


x <- seq(0, 1, length.out=500)
y <- seq(0, 1, length.out=500)
z <- outer(x, y, V)



png("contour.png", width=1024, height=1024)
pal <- function(n) sequential_hcl(n, "Purples 3")
filled.contour(z, color.palette=pal, nlevels=100,
               xlab="th1", ylab="th2")
dev.off()

