library(rstudioapi)
    
setwd(dirname(getActiveDocumentContext()$path))

#'import training data
dati=read.csv("data_set_ALL_AML_train.csv", sep=",", quote = "", header = TRUE)
dati=dati[,-which(colnames(dati)=="Gene.Accession.Number")]
dati=dati[, -which(grepl("call", colnames(dati)))]

n=dati$Gene.Description
dati=dati[,-which(colnames(dati)=="Gene.Description")]

dati=apply(dati, 2, as.numeric)
dati=as.data.frame(t(dati))

colnames(dati)=n
dati=dati[,!duplicated(colnames(dati))]
dati=dati[order(row.names(dati)), ]

#'importing testing data
test.set=read.csv("data_set_ALL_AML_independent.csv", sep=",", quote = "", header = TRUE)
test.set=test.set[,-which(colnames(test.set)=="Gene.Accession.Number")]
test.set=test.set[, -which(grepl("call", colnames(test.set)))]

n.test=test.set$Gene.Description
test.set=test.set[,-which(colnames(test.set)=="Gene.Description")]

test.set=apply(test.set, 2, as.numeric)
test.set=as.data.frame(t(test.set))

colnames(test.set)=n.test
test.set=test.set[,!duplicated(colnames(test.set))]




#'importing response variable
#'training
y=read.csv("actual.csv", sep = ",", header = TRUE)
y.train=y[1:38,]
y.train=cbind(nomi=rownames(dati), y.train)
y.train=y.train[order(y.train$nomi), ]
#'testing
test.set=test.set[order(row.names(test.set)), ]
y.test=y[39:72,]
y.test=cbind(nomi=rownames(test.set), y.test)
y.test=y.test[order(y.test$nomi), ]



#'model
dati=cbind(dati, y=y.train$cancer)
which(sapply(dati, is.factor))
glm.fit=glm(y~., dati, family = "binomial")




summary(glm.fit)
predict(glm.fit, test.set)
sum(colnames(test.set)!=colnames(dati))
