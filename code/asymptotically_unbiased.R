par(mfrow=c(1,2))

beta <- 10

U = rnorm(500, mean=0, sd=10)
X <- data.frame(rnorm(n = 500, mean = 0, sd = 1))
Y <- beta * X + U
data <- data.frame(X, Y)
colnames(data) <- c("X", "Y")

windowsFonts(A = windowsFont("Times New Roman"))



plot(data[X <= 0, ], 
     ylim=c(-30,30) ,
     xlim = c(-3,3), 
     pch=19, 
     xlab="X",
     ylab="Y",cex.lab=1.2, cex.axis=1.2, cex.main=1.2, cex.sub=1.2,
     family="A",
     yaxt='n',
     bty = "n",
     axes=FALSE)
axis(1, col="white", tcl=0)
axis(2, col="white", tcl=0, las=1)
points(data[X > 0, ], col="orange", pch=19)
abline(lm(Y ~ X, data), lwd=2, col="orange")
abline(lm(Y ~ X, data[X <= 0, ]), lwd=2)


plot(data[Y <= 0, ],
     ylim=c(-30,30) ,
     xlim = c(-3,3), 
     pch=19, 
     xlab="X",
     ylab="Y",cex.lab=1.2, cex.axis=1.2, cex.main=1.2, cex.sub=1.2,
     family="A",
     yaxt='n',
     bty = "n",
     axes=FALSE)
axis(1, col="white", tcl=0)
axis(2, col="white", tcl=0, las=2)
points(data[Y > 0, ], col="orange", pch=19)
abline(lm(Y ~ X, data), lwd=2, col="orange")
abline(lm(Y ~ X, data[Y <= 0, ]), lwd=2)
