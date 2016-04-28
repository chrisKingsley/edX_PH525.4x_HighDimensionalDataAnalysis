# week 3
library(genefilter)


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


# Conditional Expectation Exercises #1
n = 10000
set.seed(1)
men = rnorm(n,176,7) #height in centimeters
women = rnorm(n,162,7) #height in centimeters
y = c(rep(0,n),rep(1,n))
x = round(c(men,women))
##mix it up
ind = sample(seq(along=y))
y = y[ind]
x = x[ind]

idx.176 <- which(x==176)
mean(y[idx.176])


# Conditional Expectation Exercises #2
prob.female <- sapply(160:178, function(height) mean(y[x==height]))
plot(160:178, prob.female, las=3)
abline(h=0.5, col="red", lty="dashed")
(160:178)[which(prob.female > 0.5)]


# Smoothing Exercises #1
set.seed(5)
N = 250
ind = sample(length(y),N)
Y = y[ind]
X = x[ind]

loess.model <- loess(Y ~ X)
predict(loess.model, 168)


# Smoothing Exercises #2
set.seed(5)
numPerms <- 1000
predicted.vals <- rep(NA, numPerms)
for(i in 1:numPerms) {
    N = 250
    ind = sample(length(y),N)
    Y = y[ind]
    X = x[ind]
    
    loess.model <- loess(Y ~ X)
    predicted.vals[i] <- predict(loess.model, 168)
}
sd(predicted.vals)


# kNN and Cross Validation Exercises #1
library(caret)

y <- factor(sampleInfo$group)
X <- t(geneExpression)
out <- which(geneAnnotation$CHR %in% c("chrX","chrY"))
X <- X[,-out]

set.seed(1)
folds <- createFolds(y, 10)
folds[[3]][2]


# kNN and Cross Validation Exercises #2
library(class)

train.idx <- unlist(folds[-2])
gene.t.tests <- colttests(X[train.idx, ], y[train.idx])
t.idx <- order(gene.t.tests$p.value)[1:8]

res <- knn(X[train.idx, t.idx], X[-train.idx, t.idx], y[train.idx], k=5)
sum(res!=y[-train.idx])


# kNN and Cross Validation Exercises #3
set.seed(1)
m <- 8
k <- 5
numErrors <- 0
for(i in 1:length(folds)) {
    train.idx <- unlist(folds[-i])
    gene.t.tests <- colttests(X[train.idx, ], y[train.idx])
    t.idx <- order(gene.t.tests$p.value)[1:m]
    
    res <- knn(train=X[train.idx, t.idx], test=X[-train.idx, t.idx],
               cl=y[train.idx], k=k)
    numErrors <- numErrors + sum(res!=y[-train.idx])
}
numErrors/length(y)


# kNN and Cross Validation Exercises #4
paramErrors <- expand.grid(k=seq(1,9,2), m=2^c(1:11))

for(i in 1:nrow(paramErrors)) {
    numErrors <- 0
    
    for(j in 1:length(folds)) {
        train.idx <- unlist(folds[-j])
        gene.t.tests <- colttests(X[train.idx, ], y[train.idx])
        t.idx <- order(gene.t.tests$p.value)[1:paramErrors$m[i]]
        
        res <- knn(train=X[train.idx, t.idx], test=X[-train.idx, t.idx],
                   cl=y[train.idx], k=paramErrors$k[i])
        numErrors <- numErrors + sum(res!=y[-train.idx])
    }
    
    paramErrors[i, "error"] <- numErrors/length(y)
}
paramErrors[paramErrors$error==min(paramErrors$error),]


# kNN and Cross Validation Exercises #5
paramErrors <- expand.grid(k=seq(1,9,2), m=2^c(1:11))

for(i in 1:nrow(paramErrors)) {
    numErrors <- 0
    gene.t.tests <- colttests(X, y)
    t.idx <- order(gene.t.tests$p.value)[1:paramErrors$m[i]]
    
    for(j in 1:length(folds)) {
        train.idx <- unlist(folds[-j])
        res <- knn(train=X[train.idx, t.idx], test=X[-train.idx, t.idx],
                   cl=y[train.idx], k=paramErrors$k[i])
        numErrors <- numErrors + sum(res!=y[-train.idx])
    }
    
    paramErrors[i, "error"] <- numErrors/length(y)
}
paramErrors[paramErrors$error==min(paramErrors$error),]


# kNN and Cross Validation Exercises #6
paramErrors <- expand.grid(k=seq(1,9,2), m=2^c(1:11))
y2 <- factor(as.numeric(format( sampleInfo$date, "%m")=="06"))

for(i in 1:nrow(paramErrors)) {
    numErrors <- 0
    
    for(j in 1:length(folds)) {
        train.idx <- unlist(folds[-j])
        gene.t.tests <- colttests(X[train.idx, ], y2[train.idx])
        t.idx <- order(gene.t.tests$p.value)[1:paramErrors$m[i]]
        
        res <- knn(train=X[train.idx, t.idx], test=X[-train.idx, t.idx],
                   cl=y2[train.idx], k=paramErrors$k[i])
        numErrors <- numErrors + sum(res!=y2[-train.idx])
    }
    
    paramErrors[i, "error"] <- numErrors/length(y2)
}
paramErrors[paramErrors$error==min(paramErrors$error),]


# differences between colttests and t.test - need to specify equal variances
# in t.test function
train.idx <- unlist(folds[-2])
coltTests <- apply(X, 2, function(z) t.test(z[train.idx] ~ y[train.idx],
                                            var.equal=T)$p.val)
t.idx <- order(coltTests)

gene.t.tests <- colttests(X[train.idx, ], y[train.idx])
t.idx2 <- order(gene.t.tests$p.value)

plot(t.idx, t.idx2)
sum(t.idx==t.idx2)
