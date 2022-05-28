# Constant trend filtering (the 1d fused lasso)
set.seed(0)
n = 100
beta0 = rep(sample(1:10,5),each=n/5)
y = beta0 + rnorm(n,sd=0.8)
a = fusedlasso1d(y)
plot(a)

# Linear trend filtering
set.seed(0)
n = 100
beta0 = numeric(n)
beta0[1:20] = (0:19)*4/19+2
beta0[20:45] = (25:0)*3/25+3
beta0[45:80] = (0:35)*9/35+3
beta0[80:100] = (20:0)*4/20+8
y = beta0 + rnorm(n)
a = trendfilter(y,ord=1)
plot(a,df=c(10))

plot(y,xlab="Time",ylab="Trend filter estimate",main="Trend filtering example")
lines(a$fit[,30])