library(jsonlite)
library(ggplot2)
library(plyr)
library(dplyr)
library(RColorBrewer)
library(pROC)

#Load Dataset 
data <- read.csv("C:/Users/Anny/Desktop/ASM/heart_data.csv")
#Pre-processing of Data
q2data <- data
summary(q2data)
sex <- factor(q2data$Sex)
summary(sex)
FastingBloodSugar<- factor(q2data$FastingBloodSugar)
summary(FastingBloodSugar)
ExerciseInduced<- factor(q2data$ExerciseInduced)
summary(ExerciseInduced)
Class<- factor(q2data$Class)
summary(Class)
Slope<- factor(q2data$Slope)
summary(Slope)
MajorVessels<- factor(q2data$MajorVessels)
summary(MajorVessels)
# Exclier processing
boxplot(q2data$Age)
boxplot(q2data$RestBloodPressure)
boxplot(q2data$SerumCholestoral)
boxplot(q2data$MaxHeartRate)
# build model
# model 1
q2data <-read.csv("C:/Users/Anny/Desktop/ASM/scale.csv")
head(q2data)
fit.model1 <- glm(Class~Age+Sex+RestBloodPressure+SerumCholestoral+FastingBloodSugar+MaxHeartRate+ExerciseInduced+Slope+MajorVessels,
                data=q2data,family = binomial())
summary(fit.model1)
# estimate model 1
pre <- predict(fit.model1,q2data)
modelroc <- roc(q2data$Class,pre)
plot(modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
     grid.col=c("green", "red"), max.auc.polygon=TRUE,
     auc.polygon.col="skyblue", print.thres=TRUE)
# model 2
fit.model2 <- glm(Class~Sex+RestBloodPressure+SerumCholestoral+FastingBloodSugar+MaxHeartRate+ExerciseInduced+Slope+MajorVessels,
                  data=q2data,family = binomial())
summary(fit.model2)
# estimate model 2
pre <- predict(fit.model2,q2data)
modelroc <- roc(q2data$Class,pre)
plot(modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
     grid.col=c("green", "red"), max.auc.polygon=TRUE,
     auc.polygon.col="skyblue", print.thres=TRUE)
