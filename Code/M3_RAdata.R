install.packages(c("sparsebn"))
install.packages(c("caret"))
install.packages("tidyr")
install.packages("colortools")
install.packages("glmnet", repos = "http://cran.us.r-project.org")

library("glmnet")
library("colortools")
library("sparsebn")
library("caret")
library("tidyr")


MO3_data <- read.csv("data_MO3.csv")

#inverte variables/genes to columns
n <- MO3_data$Geneid
MO3_data <- as.data.frame(t(MO3_data[,-1]))
colnames(MO3_data) <- n
MO3_data$myfactor <- factor(row.names(MO3_data))
str(MO3_data) 

#vector w/ responders and non responders
MO3.response.vector <- as.vector(MO3_data[,25371])
MO3_data$myfactor <- NULL
MO3_data$RESPONDER <- NULL

#to find and delete variables with zero variance
apply(MO3_data, 2, var)
which(apply(MO3_data, 2, var) == 0)
MO3data_filtered <- MO3_data[ - as.numeric(which(apply(MO3_data, 2, var) == 0))]
str(MO3_data)

#Log transformation
MO3data_filtered <- log1p(MO3data_filtered)

#Standardization
#variables have been standardized with a mean value of zero
MO3_transformation <- preProcess(MO3data_filtered, method=c("center", "scale"))
MO3data_filtered <- predict(MO3_transformation, MO3data_filtered)
summary(MO3data_filtered)



                          #logistic regression: binomial model
MO3.response.vector = as.matrix(MO3.response.vector)
MO3data_filtered = as.matrix(MO3data_filtered)


#PARTIÇÃO (100x) - TESTE/TREINO

library(pROC)
degree = 1:100
auc_total = NULL


for (d in degree){
  index <- sample(nrow(MO3data_filtered), 46)
  datatreino <- MO3data_filtered[index, ]
  datateste <- MO3data_filtered[-index, ]
  vectortreino <- MO3.response.vector[index, ]
  vectorteste <- MO3.response.vector[-index, ]
  
  cvfit = cv.glmnet(datatreino, vectortreino, family = "binomial", type.measure = "class")
  
  modelteste <- glmnet(datatreino,
                       vectortreino,
                       alpha = 0.5,
                       lambda = cvfit$lambda.1se)
  #elastic_test <- coef(modelteste, s = 'lambda.1se')[,1] %>% {.[.!=0]} #to select jus variables with coef different from zero
  #genes.selected <- as.data.frame(names(elastic_test))
  #genes.selnames <- rbind(genes.selnames0, genes.selected)#two save gene names selected in a data.frame
  
  
  test <- predict(modelteste,datateste) #testar no conjunto teste
  testdf <- as.data.frame(test)
  vector.df <- as.data.frame(t(vectorteste))
  testdf$responder <- (vectorteste)
  
  #to find ROC
  treinoclasses<- data.frame(testdf$responder, stringsAsFactors=FALSE)
  treinoprediction <- data.frame(testdf$s0, stringsAsFactors=FALSE)
  obser.classes <- as.vector(t(treinoclasses)) #numeric vector check
  pred.probabilities <- as.vector(t(treinoprediction)) #numeric vector check
  
  res.roc <- roc(obser.classes, pred.probabilities)
  auc <- as.data.frame(res.roc$auc) #ou as.vector
  auc_total <- rbind(auc_total, auc)
  
}
#The user should save auc_total for each alpha used for further analysis