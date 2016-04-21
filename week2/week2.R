# week 2

library(tissuesGeneExpression)
data(tissuesGeneExpression)


# SVD problem 1
s <- svd(e)
m <- rowMeans(e)
cor(s$u[,1], m)


# SVD problem 4
y = e - rowMeans(e)
s = svd(y)
z = s$d * t(s$v)
sqrt(crossprod(e[,3]-e[,45])) - sqrt(crossprod(z[1:2,3]-z[1:2,45]))


# SVD problem 5
for(i in 1:10) {
    d <- sqrt(crossprod(e[,3]-e[,45]))
    dhat <- sqrt(crossprod(z[1:i,3]-z[1:i,45]))
    print(sprintf("i:%d d:%f frac:%f", i, d-dhat, (d-dhat)/d))
}


# SVD problem 6
distances1 = sqrt(apply(e[,-3]-e[,3],2,crossprod))
distances2 = sqrt(apply(z[1:2,-3]-z[1:2,3],2,crossprod))
cor(distances1, distances2, method="spearman")


# MDS problem 1
d = dist(t(e))
mds = cmdscale(d)
cor(z[1,], mds[,1])


# MDS problem 2
cor(z[2,], mds[,2])


# MDS problem 4-8
library(GSE5859Subset)
data(GSE5859Subset)
s = svd(geneExpression-rowMeans(geneExpression))
z = s$d * t(s$v)

rowGroupCors <- apply(z, 1, function(x) cor(x, sampleInfo$group))
which(rowGroupCors==max(rowGroupCors))
max(rowGroupCors)

which(rowGroupCors==sort(rowGroupCors, decreasing=T)[2])

month <- format(sampleInfo$date, "%m")
month <- as.numeric(month)
rowMonthCors <- apply(z, 1, function(x) cor(x, month))
max(rowMonthCors)
which(rowMonthCors==max(rowMonthCors))

good.idx <- which(geneAnnotation$CHR!="chrUn")
boxplot(s$u[good.idx, 5] ~ geneAnnotation$CHR[good.idx], las=3)
