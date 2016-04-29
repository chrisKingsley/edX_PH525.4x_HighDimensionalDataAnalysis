# week 4
library(genefilter)
library(dplyr)
library(Biobase)
library(qvalue)
library(GSE5859)
data(GSE5859)


# load admissions table
load('admissions.rda')


# Confounding Exercises #1
index <- which(admissions$Gender==0)
accepted <- sum(admissions$Number[index] * admissions$Percent[index]/100)
applied <- sum(admissions$Number[index])
accepted/applied


# Confounding Exercises #2
accTable <- admissions %>% 
    group_by(Gender) %>% 
    summarize(accepted=sum(Number*Percent/100),
              notAccept=sum(Number-Number*Percent/100))
chisq.test(accTable)


index <- admissions$Gender==1
men <- admissions[index,]
women <- admissions[!index,]
print( data.frame( major=admissions[1:6,1],men=men[,3], women=women[,3]) )


# Confounding Exercises #3/4
accTable2 <- admissions %>% 
    group_by(Major) %>% 
    summarize(totalApplicants=sum(Number),
              fracAdmitted=sum(Number*Percent/100)/sum(Number))
accTable2[accTable2$fracAdmitted==min(accTable2$fracAdmitted),]


# Confounding Exercises #5/6
cor(admissions$Number[admissions$Gender==1], accTable2$fracAdmitted)
cor(admissions$Number[admissions$Gender==0], accTable2$fracAdmitted)


# Confounding in Genomics Exercises #1/2
geneExpression <- exprs(e)
sampleInfo <- pData(e)
sampleInfo$year <- format(sampleInfo$date,"%y")
sampleInfo$month.year <- format(sampleInfo$date,"%m%y")

byYear <- sampleInfo %>% group_by(year) %>% 
    summarize(numEthnicities=length(unique(ethnicity)))
sum(byYear$numEthnicities > 1)

byMonthYear <- sampleInfo %>% group_by(month.year) %>% 
    summarize(numEthnicities=length(unique(ethnicity)))
sum(byMonthYear$numEthnicities > 1)/nrow(byMonthYear)


# Confounding in Genomics Exercises #3
sampleFactor <- rep(NA, nrow(sampleInfo))
sampleFactor[which(sampleInfo$ethnicity=="CEU" & sampleInfo$year=="03")] <- "03"
sampleFactor[which(sampleInfo$ethnicity=="CEU" & sampleInfo$year=="02")] <- "02"

tTests <- rowttests(geneExpression,  as.factor(sampleFactor))
qValues <- qvalue(tTests$p.value)
sum(qValues$qvalues < 0.05)
qValues$pi0


# Confounding in Genomics Exercises #4
sampleFactor <- rep(NA, nrow(sampleInfo))
sampleFactor[which(sampleInfo$ethnicity=="CEU" & sampleInfo$year=="03")] <- "03"
sampleFactor[which(sampleInfo$ethnicity=="CEU" & sampleInfo$year=="04")] <- "04"

tTests <- rowttests(geneExpression,  as.factor(sampleFactor))
qValues <- qvalue(tTests$p.value)
sum(qValues$qvalues < 0.05)


# Confounding in Genomics Exercises #5
sampleFactor <- rep(NA, nrow(sampleInfo))
sampleFactor[which(sampleInfo$ethnicity=="CEU")] <- "CEU"
sampleFactor[which(sampleInfo$ethnicity=="ASN")] <- "ASN"

tTests <- rowttests(geneExpression,  as.factor(sampleFactor))
qValues <- qvalue(tTests$p.value)
sum(qValues$qvalues < 0.05)


# Confounding in Genomics Exercises #6
sampleFactor <- rep(NA, nrow(sampleInfo))
sampleFactor[which(sampleInfo$ethnicity=="CEU" & sampleInfo$year=="05")] <- "CEU"
sampleFactor[which(sampleInfo$ethnicity=="ASN" & sampleInfo$year=="05")] <- "ASN"

tTests <- rowttests(geneExpression,  as.factor(sampleFactor))
qValues <- qvalue(tTests$p.value)
sum(qValues$qvalues < 0.05)


# Confounding in Genomics Exercises #7
set.seed(3)
sampleFactor <- rep(NA, nrow(sampleInfo))
sampleFactor[which(sampleInfo$ethnicity=="ASN" & sampleInfo$year=="05")] <- "ASN"
idx <- which(sampleInfo$ethnicity=="CEU" & sampleInfo$year=="02")
sampleFactor[sample(idx, 3)] <- "CEU"

tTests <- rowttests(geneExpression,  as.factor(sampleFactor))
qValues <- qvalue(tTests$p.value)
sum(qValues$qvalues < 0.05)

