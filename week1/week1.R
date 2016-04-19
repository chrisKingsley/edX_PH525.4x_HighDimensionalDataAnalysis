# week 1

library(tissuesGeneExpression)
data(tissuesGeneExpression)

# exercise 1
table(tissue)
sum(tissue=="hippocampus")

# exercise 2
sqrt(crossprod(e[,3] - e[,45]))

# exercise 3
sqrt(crossprod(e["210486_at",] - e["200805_at",]))

# exercise 4
nrow(e)^2

# exercise 5
d <- dist(t(e))
length(d)
