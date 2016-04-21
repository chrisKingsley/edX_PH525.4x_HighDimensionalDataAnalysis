# week 3


# Hiearchichal Clustering Exercises #1
set.seed(1)
m <- 10000
n <- 24
x <- matrix(rnorm(m*n),m,n)
colnames(x) <- 1:n
distMat <- dist(t(x))
plot(hclust(distMat))
distMat


# Hiearchichal Clustering Exercises #2
set.seed(1)
numClusters <- rep(NA, 100)
for(i in 1:100) {
    m <- 10000
    n <- 24
    x <- matrix(rnorm(m*n),m,n)
    numClusters[i] <- length(unique(cutree(hclust(dist(t(x))), h=143)))
}
sd(numClusters)


# K-means Exercises #1
library(GSE5859Subset)
data(GSE5859Subset)

set.seed(10)
kClust <- kmeans(t(geneExpression), centers=5)
table(kClust$cluster, sampleInfo$group)
table(kClust$cluster, sampleInfo$date)
table(kClust$cluster, sampleInfo$ethnicity)


# Heat Maps Exercises #1
library(gplots)

rowMads <- apply(geneExpression, 1, mad, na.rm=T)
gene.idx <- order(rowMads, decreasing=T)[1:25]
heatmap.2(geneExpression[gene.idx,], scale="row",
          ColSideColors=as.character(sampleInfo$group),
          labCol=sampleInfo$date, labRow=geneAnnotation$CHR[gene.idx])


# Heat Maps Exercises #2
set.seed(17)
m <- nrow(geneExpression)
n <- ncol(geneExpression)
x <- matrix(rnorm(m*n),m,n)
g <- factor(sampleInfo$g)

rowtTests <- apply(x, 1, function(z) t.test(z ~ g)$p.val)
t.idx <- order(rowtTests)[1:50]
heatmap.2(x[t.idx,], ColSideColors=as.character(g),
          main="Feature Selection by T-test")

rowVars <- apply(x, 1, var)
var.idx <- order(rowVars, decreasing=T)[1:50]
heatmap.2(x[var.idx,], ColSideColors=as.character(g),
          main="Feature Selection by Variance")
