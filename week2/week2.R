# week 2

library(tissuesGeneExpression)
data(tissuesGeneExpression)


# problem 1
s <- svd(e)
m <- rowMeans(e)
cor(s$u[,1], m)

# problem 2
