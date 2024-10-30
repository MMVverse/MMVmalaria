

X <- c(1,2,3,4,5,6,7)
Y <- c(10,4,1,1,6,10,20)


Y1 <- findMinMMV(X,Y,Method = "NoInterpolation")
Y2 <- findMinMMV(X,Y,Method = "CubicSpline")
Y3 <- findMinMMV(X,Y,Method = "Quadratic")

plot(X,Y,type="b",ylim=c(0,20))
points(Y1,pch=16,col="blue")
points(Y2,pch=16,col="green")
points(Y3,pch=16,col="red")
legend("topleft",legend=c("NoInterpolation","CubicSpline","Quadratic"),pch=rep(16,3),col=c("blue","green","red"))

