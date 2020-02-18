#TOP 100 GENES analysis 

library("glmnet")
library("colortools")
library("sparsebn")
library("caret")
library("tidyr")
library("tidyverse")
library("readxl")
#glmSparseNet:
library(dplyr)
library(ggplot2)
library(survival)
library(loose.rock)
library(futile.logger)


baseline_data <- read_excel("data_baseline.xlsx")


# first remember the names
n <- baseline_data$Geneid
# transpose all but the first column (name)
baseline_data <- as.data.frame(t(baseline_data[,-1]))
colnames(baseline_data) <- n
baseline_data$myfactor <- factor(row.names(baseline_data))
str(baseline_data) # Check the column types


#to find and delete variables with zero variance
data.filtered <- baseline_data[ - as.numeric(which(apply(baseline_data, 2, var) == 0))]


#SELECT TOP 100 GENES reported by Farutin et al.
topgenesdata <- data.filtered %>% select(ALPL, MGAM, PI3, CYP4F3, PGLYRP1, CA4, PROK2, ANXA3,
                                         ADM, MMP9, SLPI, BTNL8, CCNJL, ANPEP, MANSC1, KAZN, TRPM6, QPCT, MME,WLS,
                                         S100A12,TGM3, CAMP, MAK, WDFY3, MRVI1, KRT23, LGALSL, ACSL1, BASP1, AATK,
                                         S100P, IL1R2, LRRC4, TNFAIP6, DGAT2, FPR2, CMTM2, LTF, PYGL, CACNA1E, SLC26A8,
                                         DYSF, PLIN5, SLC22A4, AQP9, KREMEN1, HIST1H2BC, CDA, LRG1, LIN7A, IL1R1, FCAR,
                                         CEACAM3, PANX2, C5AR1, GLT1D1, TREM1, SIGLEC5, ROPN1L, TNFRSF10C, STEAP4,AOC3,
                                         VNN1, CEACAM4, VNN3, DSC2, HCAR3, NOV, AOC2, LUCAT1, SLC11A1, CLEC4D, CRISPLD2,
                                         LCN2, DOCK5, TCN1, C5AR2, SIPA1L2, PLIN4, HCAR2, PLB1, PADI4, PTGS2, HAL,
                                         REPS2, ARG1, NAMPT, FOLR3, CHI3L1, NLRP12, LSMEM1, MCEMP1, TNFRSF9, B3GNT8, 
                                         HRH2, "SEPT5-GP1BB", SLC8A1, FPR1, RESPONDER)

#vector w/ responders and non responders
topgenes.vector <- as.vector(topgenesdata$RESPONDER)

topgenesdata$RESPONDER <- NULL


#Log transformation
topgenesdata <- log1p(topgenesdata)


#Standardization
#variables have been standardized with a mean value of zero
transformation <- preProcess(topgenesdata, method=c("center", "scale"))
topgenesdata <- predict(transformation, topgenesdata)
#summary(topgenesdata)


class(topgenesdata) #check that is a data.frame
topgenes.vector = as.data.frame(topgenes.vector)
topgenes.vector = t(topgenes.vector)



#BINOMIAL REGRESSION 
topgenesdata = as.matrix(topgenesdata)

library(pROC)

topgenes.vector <- t(topgenes.vector)
topgenes.vector <- as.data.frame(topgenes.vector)

degree = 1:100
auc_total = NULL

for (d in degree){
  index <- sample(nrow(topgenesdata), 43)
  datatreino <- topgenesdata[index, ]
  datateste <- topgenesdata[-index, ]
  vectortreino <- topgenes.vector[index, ]
  vectorteste <- topgenes.vector[-index, ]
  
  cvfit = cv.glmnet(datatreino, vectortreino, family = "binomial")
  
  modelteste <- glmnet(datatreino,
                       vectortreino,
                       alpha = 0, #alpha=0 to the model use the 100 variables 
                       lambda = cvfit$lambda.1se)
  
  test <- predict(modelteste,datateste)#testar no conjunto teste
  
  testdf <- as.data.frame(test)
  vector.df <- as.data.frame(t(vectorteste))
  testdf$responder <- (vectorteste)
  
  #two find ROC
  treinoclasses<- data.frame(testdf$responder, stringsAsFactors=FALSE)
  treinoprediction <- data.frame(testdf$s0, stringsAsFactors=FALSE)
  obser.classes <- as.vector(t(treinoclasses)) #numeric vector check
  pred.probabilities <- as.vector(t(treinoprediction)) #numeric vector check
  
  res.roc <- roc(obser.classes, pred.probabilities)
  auc <- as.data.frame(res.roc$auc) #ou as.vector
  auc_total <- rbind(auc_total, auc)
  
}

#The user should save auc_total for each alpha used for further analysis




