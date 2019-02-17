#!/usr/bin/env Rscript

library(colorspace)

a1 <- 1.1
a2 <- 1.1
a3 <- 1.1

k1 <- 6
k2 <- 2
k3 <- 4
C12 <- 5
C13 <- 3
C23 <- 4

# V = -log(posterior)
V <- function(th1, th2) -((k1+a1-1)*log(th1) + (k2+a2-1)*log(th2) + (k3+a3-1)*log1p(-(th1+th2))
                          - C12*log(th1+th2) - C13*log(1-th2) - C23*log(1-th1))

# posterior
p <- function(th1, th2) exp(-V(th1, th2))


x <- seq(0, 1, length.out=1000)
y <- seq(0, 1, length.out=1000)
z <- outer(x, y, V)



png("contour.png", width=1024, height=1024)
pal <- function(n) sequential_hcl(n, "Purples 3")
filled.contour(log(z), color.palette=pal, nlevels=15,
               main="log V", xlab="th1", ylab="th2")
dev.off()

