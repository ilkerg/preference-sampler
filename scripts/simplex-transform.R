transform <- function(x) {
    y <- numeric(2)
    y[1] <- log(x[1]) - log(1-x[1]) + log(2)
    y[2] <- log(x[2]) - log(1-x[1]-x[2])
    y
}

inverse_transform <- function(y) {
    z <- numeric(2)
    x <- numeric(2)

    z[1] <- 1. / (1 + 2*exp(-y[1]))
    z[2] <- 1. / (1 + exp(-y[2]))

    x[1] <- z[1]
    x[2] <- (1-x[1])*z[2]

    x
}

disk_radius <- function(x0, radius=1) {
    y0 <- transform(x0)
    th <- seq(0, 2*pi, length.out=100)
    y <- cbind(y0[1] + radius*cos(th), y0[2] + radius*sin(th))
    t(apply(y, 1, inverse_transform))
}

unit_disk <- function(x0) disk_radius(x0, 1)

