install.packages(c("sparsebn"))
install.packages(c("caret"))
install.packages("tidyr")
install.packages("tidyverse")
install.packages("colortools")
install.packages("glmnet", repos = "http://cran.us.r-project.org")
install.packages("readxl")
install.packages("xlsx")
install.packages("VennDiagram")


library("glmnet")
library("colortools")
library("sparsebn")
library("caret")
library("tidyr")
library("tidyverse")
library("readxl")
library("pROC")



baseline_data <- read_excel("data_baseline.xlsx")


# first remember the names
n <- baseline_data$Geneid
# transpose all but the first column (name)
baseline_data <- as.data.frame(t(baseline_data[,-1]))
colnames(baseline_data) <- n
baseline_data$myfactor <- factor(row.names(baseline_data))
str(baseline_data) # Check the column types


#vector w/ responders and non responders
response.vector <- as.vector(baseline_data[,25371])

baseline_data$myfactor <- NULL
baseline_data$RESPONDER <- NULL



#to find and delete variables with zero variance
baselinedata_filtered <- baseline_data[ - as.numeric(which(apply(baseline_data, 2, var) == 0))]
str(baseline_data)


#Log transformation
baselinedata_filtered <- log1p(baselinedata_filtered)

#Standardization
#variables have been standardized with a mean value of zero
transformation <- preProcess(baselinedata_filtered, method=c("center", "scale"))
baselinedata_filtered <- predict(transformation, baselinedata_filtered)
summary(baselinedata_filtered)



                            #LOGISTIC REGRESSION: binomial model

baselinedata_filtered = as.matrix(baselinedata_filtered)


#PARTITION (100x) - TEST/TRAIN

degree = 1:100
auc_total = NULL

for (d in degree){
  index <- sample(nrow(baselinedata_filtered), 43)
  datatreino <- baselinedata_filtered[index, ]
  datateste <- baselinedata_filtered[-index, ]
  vectortreino <- response.vector[index, ]
  vectorteste <- response.vector[-index, ]
  
  cvfit = cv.glmnet(datatreino, vectortreino, family = "binomial", type.measure = "class")
  
  modelteste <- glmnet(datatreino,
                       vectortreino,
                       alpha = 0,
                       lambda = cvfit$lambda.1se)

  
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


#LEAVE ONE OUT METHOD
#use alpha chosen with test/train method
degree = 1:63
response.vector = as.data.frame(response.vector)
loo <- NULL
loo_names <- NULL

for (d in degree){
  datan <- baselinedata_filtered[-d,]
  response.vectorn <- response.vector[-d,]
  cvfit = cv.glmnet(datan, response.vectorn, family = "binomial", type.measure = "class")
  
  model_loo <- glmnet(datan,
                      response.vectorn,
                      alpha = 0.4,
                      lambda = cvfit$lambda.1se)
  
  data_justn <- as.data.frame(baselinedata_filtered[d,])
  
  pred <- as.data.frame(predict(model_loo, t(data_justn)))
  loo <- rbind(loo, pred)
  
  loo_elastic <- coef(model_loo, s = 'lambda.1se')[,1] %>% {.[.!=0]}
  length(loo_elastic)
  loo_elasticnames <- as.data.frame(names(loo_elastic)) #ou as.vector
  loo_names <- rbind(loo_names, loo_elasticnames)
  
}

#To count the number of repeated genes in the 63 models
#with method Leave One Out 
loo_namest <- as.vector(t(loo_names))
loo_selectedgenes <- table(loo_namest)
loo_selectedgenes <- as.data.frame(loo_selectedgenes)
totalgenes_selected <- as.data.frame(loo_selectedgenes$Freq)
totalgenes_selected <- as.vector(t(totalgenes_selected)) #NR OF GENES ALWAYS SELECTED
table(totalgenes_selected)

#repeated.genes0.1 <- filter(loo_selectedgenes, Freq == 63)
repeated.genes0.4 <- filter(loo_selectedgenes, Freq == 63)


# Compute roc for leave one out 
library(pROC)
loo$responder <- (response.vector)
classes<- data.frame(loo$responder, stringsAsFactors=FALSE)
prediction <- data.frame(loo$s0, stringsAsFactors=FALSE)
observed.classes <- as.vector(t(classes)) #numeric vector check
prediction.probabilities <- as.vector(t(prediction)) #numeric vector check

res.roc <- roc(observed.classes, prediction.probabilities)
plot.roc(res.roc, print.auc = TRUE)


#fit the final model ALPHA = 0.1 and select the variables chosen
response.vector = t(response.vector)
baselinedata_filtered = as.matrix(baselinedata_filtered)

cvfit = cv.glmnet(baselinedata_filtered, response.vector, family = "binomial", type.measure = "class")
cvfit$lambda.1se

model <- glmnet(baselinedata_filtered,
                response.vector,
                alpha = 0.1,
                lambda = cvfit$lambda.1se)
elastic <- coef(model, s = 'lambda.1se')[,1] %>% {.[.!=0]} #to select jus variables with coef different from zero
length(elastic)
names(elastic)

#to select in the data the genes given by elastic net with a specific alpha
elastic.names <- as.vector(names(elastic))
baselinedata_filtered <- as.data.frame(baselinedata_filtered)
names.use <- names(baselinedata_filtered)[(names(baselinedata_filtered) %in% elastic.names)]
baselinedata_subset <- baselinedata_filtered[,names.use]

#split Responders from Non-Responders for Bayesian Network learning
baselinedata_subset$responder <- response.vector
R_BLdata<-baselinedata_subset[!(baselinedata_subset$responder=="0"),]
R_BLdata$responder <- NULL
NR_BLdata<-baselinedata_subset[!(baselinedata_subset$responder=="1"),]
NR_BLdata$responder <- NULL



                    #ANALYSE SEPARATELY RESPONDERS AND NON RESPONDERS

#For Responser data

names(R_BLdata)
RBLgene.data <- sparsebnData(R_BLdata,
                             type = "continuous")

#structure learning
RBLgene.learn <- estimate.dag(data = RBLgene.data, edge.threshold = 239) #edge.threshold to specify when the algorithm stops
RBLgene.learn


getPlotPackage()
plot(RBLgene.learn[[6]],
     vertex.size = 3,
     vertex.label = '',
     vertex.color = "yellowgreen",
     edge.color = gray(0),
     edge.arrow.size = 0.1)



RBLedges = as.matrix(RBLgene.learn[[6]][["edges"]])
RBLedges <- as.data.frame(RBLedges)

#eliminate variables without an edge
RBLedges <- RBLedges[ - as.numeric(which(apply(RBLedges, 2, var) == 0))]
RBLedges <- RBLedges[-which(apply(RBLedges[,-1], 1, function(x) all(x == 0)) == T),]
colSums(RBLedges != 0) #give me de genes with connection


#To get a matrix with genes edges to compare with STRING information
idx <- which(RBLedges == 1, arr.ind=TRUE)
RBLedges.compare <- cbind(rownames(RBLedges)[idx[,"row"]], colnames(RBLedges)[idx[,"col"]])
RBLedges.compare <- as.data.frame(RBLedges.compare)


#common genes interaction with STRING
STRINGdata <- read_excel("STRING_IDnames.xlsx")

RBLedges.compare$comp <- paste(RBLedges.compare$V1,RBLedges.compare$V2, sep="")
RBLedges.compare$compI <- paste(RBLedges.compare$V2,RBLedges.compare$V1, sep="")
STRINGdata$comp <- paste(STRINGdata$protein1,STRINGdata$protein2, sep="")

#count how many edges exist in String
summary(RBLedges.compare$comp %in% STRINGdata$comp)
commonR.S <- as.data.frame(RBLedges.compare$comp %in% STRINGdata$comp)
summary(RBLedges.compare$compI %in% STRINGdata$comp)
commonR.SI <- as.data.frame(RBLedges.compare$compI %in% STRINGdata$comp)



  #For Non Responser data

NRBLgene.data <- sparsebnData(NR_BLdata,
                              type = "continuous")

#structure learning
NRBLgene.learn <- estimate.dag(data = NRBLgene.data, edge.threshold = 239) #edge.threshold to specify when the algorithm stops
NRBLgene.learn


getPlotPackage()
plot(NRBLgene.learn[[6]],
     vertex.size = 3,
     vertex.label = '',
     vertex.color = "red",
     edge.color = gray(0),
     edge.arrow.size = 0.1)

NRBLedges <- as.matrix(NRBLgene.learn[[6]][["edges"]])
NRBLedges <- as.data.frame(NRBLedges)

#eliminate variables without an edge
NRBLedges <- NRBLedges[ - as.numeric(which(apply(NRBLedges, 2, var) == 0))]
NRBLedges <- NRBLedges[-which(apply(NRBLedges[,-1], 1, function(x) all(x == 0)) == T),]
colSums(NRBLedges != 0) #give me de genes with connection


#To get a matrix with genes edges to compare with STRING information
idx <- which(NRBLedges == 1, arr.ind=TRUE)
NRBLedges.compare <- cbind(rownames(NRBLedges)[idx[,"row"]], colnames(NRBLedges)[idx[,"col"]])
NRBLedges.compare <- as.data.frame(NRBLedges.compare)

#common genes interaction with STRING
NRBLedges.compare$comp <- paste(NRBLedges.compare$V1,NRBLedges.compare$V2, sep="")
NRBLedges.compare$compI <- paste(NRBLedges.compare$V2,NRBLedges.compare$V1, sep="")

#count how many common edges exist with/in String
summary(NRBLedges.compare$comp %in% STRINGdata$comp)
commonNR.S <- as.data.frame(NRBLedges.compare$comp %in% STRINGdata$comp)
summary(NRBLedges.compare$compI %in% STRINGdata$comp)
commonNR.SI <- as.data.frame(NRBLedges.compare$compI %in% STRINGdata$comp)



#COMPARE COMMON GENES BETWEEN R AND NR
summary(NRBLedges.compare$comp %in% RBLedges.compare$comp)
summary(NRBLedges.compare$compI %in% RBLedges.compare$comp)



#VENN DIAGRAM
library(VennDiagram)

venn.diagram(
  x = list(RBLedges.compare$compI, NRBLedges.compare$compI, STRINGdata$comp),
  category.names = c("R-BL" , "NR-BL " , "STRING"),
  filename = '#239_venn_diagramm.png',
  output = TRUE ,
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  col=c("#440154ff", '#21908dff', '#fde725ff'),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3)),
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.3,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  cat.col = c("#440154ff", '#21908dff', '#fde725ff'),
  #rotation = 3
)



