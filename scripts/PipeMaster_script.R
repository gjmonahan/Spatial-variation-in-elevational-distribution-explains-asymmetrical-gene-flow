###Running models pipemaster 
library(PipeMaster)
library(ggplot2)
library(caret)
library(doMC)
library(ggpubr)
library(vroom)
library(keras3)
library(matrixStats)
species="test_species"
setwd("/work/gjmonaha/Pipemaster/test")
#read pop assignemtspecies_pop.assign <- read.delim(paste(species, "_pop_list.txt", sep = ""), header = FALSE)
#draw models using PipeMaster
if (exists("m1_test")==F) {
  #No gene flow
  m1 <- main.menu(m1)
  #Bidirectional gene flow
  m2 <- main.menu(m2)
  #Unidirectional gene from the tropics to the subtropics
  m3 <- main.menu(m3)
  #Unidirectional gene from the subtropics to the tropics
  m4 <- main.menu(m4)
  
  m1 <- get.data.structure(m1,path.to.fasta = "/work/gjmonaha/Pipemaster/ipyrad_loci/test_species", pop.assign=vetula_pop.assign)
  m2 <- get.data.structure(m2,path.to.fasta = "/work/gjmonaha/Pipemaster/ipyrad_loci/test_species", pop.assign=vetula_pop.assign)
  m3 <- get.data.structure(m3,path.to.fasta = "/work/gjmonaha/Pipemaster/ipyrad_loci/test_species", pop.assign=vetula_pop.assign)
  m4 <- get.data.structure(m4,path.to.fasta = "/work/gjmonaha/Pipemaster/ipyrad_loci/test_species", pop.assign=vetula_pop.assign)
}

#Calculating summary stats for observed data (5,000 windows of 10kbp)
#Windows must be fasta alignments
if (file.exists("stats_obs_test.csv")==F) {
 stat <- obs.sumstat.ngs(model = m1, path.to.fasta = "/work/gjmonaha/Pipemaster/ipyrad_loci/test_species", 
                           pop.assign = test_species_pop.assign, moments = F)
}
 
setwd("/work/gjmonaha/Pipemaster/test")
if (file.exists("SIMS_m1.txt")==F) {
 
  ###m1###  No gene flow
  sim.msABC.sumstat(m1, path="./",
                    #mu.rates = list("rtnorm",1225,5e-10,5e-10,1e-11),
                    nsim.blocks = 50, use.alpha = F,
                    output.name = "m1",
                    append.sims = F, ncores=20, block.size = 200)
  
  ###m2###  Bidirectional gene flow
  sim.msABC.sumstat(m2, path="./",
                    #mu.rates = list("rtnorm",1225,5e-10,5e-10,1e-11),
                    nsim.blocks = 50, use.alpha = F,
                    output.name = "m2",
                    append.sims = F, ncores=20, block.size = 200)
  
  ###m3###  
  sim.msABC.sumstat(m3, path="./",
                    #mu.rates = list("rtnorm",1225,5e-10,5e-10,1e-11),
                    nsim.blocks = 50, use.alpha = F,
                    output.name = "m3",
                    append.sims = F, ncores=20, block.size = 200)
  
  ###m4###  
  sim.msABC.sumstat(m4, path="./",
                    #mu.rates = list("rtnorm",1225,5e-10,5e-10,1e-11),
                    nsim.blocks = 50, use.alpha = F,
                    output.name = "m4",
                    append.sims = F, ncores=20, block.size = 200)
  
  print("It's Done")
}


##############################################################################################################

##reading simulations
m1.sim <- read.delim("SIMS_m1.txt")
m2.sim <- read.delim("SIMS_m2.txt")
m3.sim <- read.delim("SIMS_m3.txt")
m4.sim <- read.delim("SIMS_m4.txt")
#m5.sim <- read.delim("SIMS_m5.txt")
###Removing some of the summary stats
m1.sim <- m1.sim[,grep("tajd|dvk|dvh|fwh|FayWuH", colnames(m1.sim), fixed = F, invert = T)]
m2.sim <- m2.sim[,grep("tajd|dvk|dvh|fwh|FayWuH", colnames(m2.sim), fixed = F, invert = T)]
m3.sim <- m3.sim[,grep("tajd|dvk|dvh|fwh|FayWuH", colnames(m3.sim), fixed = F, invert = T)]
m4.sim <- m4.sim[,grep("tajd|dvk|dvh|fwh|FayWuH", colnames(m4.sim), fixed = F, invert = T)]
#concatenating the data set

#assigning models to simulations
#whole data
data <- c(rep("m1",nrow(m1.sim)),
          rep("m2",nrow(m2.sim)),
          rep("m3",nrow(m3.sim)),
          rep("m4", nrow(m4.sim)))

data_pca <- c(rep("m1",1000),
                      rep("m2",1000),
                      rep("m3",1000),
                      rep("m4",1000))

#subset for caret/Keras
data_caret <- c(rep("m1",200000),
                        rep("m2",200000),
                        rep("m3",200000),
                        rep("m4",200000))

models <- rbind(
  m1.sim[, c(grep("sd.rate", colnames(m1.sim)) + 1):ncol(m1.sim)],
  m2.sim[, c(grep("sd.rate", colnames(m2.sim)) + 1):ncol(m2.sim)],
  m3.sim[, c(grep("sd.rate", colnames(m3.sim)) + 1):ncol(m3.sim)],
  m4.sim[, c(grep("sd.rate", colnames(m4.sim)) + 1):ncol(m4.sim)]
)

models_pca <- rbind(m1.sim[1:1000,c(grep("sd.rate", colnames(m1.sim))+1):ncol(m1.sim)],
                            m2.sim[1:1000,c(grep("sd.rate", colnames(m2.sim))+1):ncol(m2.sim)],
                            m3.sim[1:1000,c(grep("sd.rate", colnames(m3.sim))+1):ncol(m3.sim)],
                            m4.sim[1:1000,c(grep("sd.rate", colnames(m4.sim))+1):ncol(m4.sim)])

models_caret <- rbind(m1.sim[1:200000,c(grep("sd.rate", colnames(m1.sim))+1):ncol(m1.sim)],
                              m2.sim[1:200000,c(grep("sd.rate", colnames(m2.sim))+1):ncol(m2.sim)],
                              m3.sim[1:200000,c(grep("sd.rate", colnames(m3.sim))+1):ncol(m3.sim)],
                              m4.sim[1:200000,c(grep("sd.rate", colnames(m4.sim))+1):ncol(m4.sim)])


#making sure observed stats and simulations have the same columns
stat <- stat[, colnames(stat) %in% colnames(models_caret), drop = FALSE]

data_caret<-data_caret[complete.cases(models_caret)]
models_caret<-models_caret[complete.cases(models_caret),]

data_pca<-data_pca[complete.cases(models_pca)]
models_pca<-models_pca[complete.cases(models_pca),]


if (identical(colnames(stat), colnames(models_caret))==T) {
  print("all good - obs stats and simulations match")
}


######Running Keras####
###Testing models ####

dat <- data_caret
mdl <- models_caret

#assigning training and testing proportions
train = 0.75
test = 0.25

#calculating the mean of each variable to use as scaling point to normalize the data
mean_mdl <- colMeans(mdl)
#It might be a good idea to normalize the data. Keras have a function for that.
normalize = T
#Normalizing
if (normalize == T) {
  mdl <- apply(mdl,2,function(x){x/mean(x)})
  stat2 <- stat/mean_mdl
}
# first convert the data set to array
data <- as.numeric(as.factor(dat))-1
models <- as.matrix(mdl[,1:ncol(mdl)])

models <- cbind(models, data)
dimnames(models) <- NULL
set.seed(12345)

# Determine sample size
ind <- sample(2, nrow(models), replace=TRUE, prob=c(train, test))
models.training <- as.matrix(models[ind==1, 1:c(ncol(models) -1)])
dim(models.training)
models.test <- as.matrix(models[ind==2, 1:c(ncol(models) -1)])
dim(models.test)
models.trainingtarget <- as.matrix(models[ind==1, ncol(models)])
dim(models.trainingtarget)
models.testtarget <- as.matrix(models[ind==2, ncol(models)])
dim(models.testtarget)

##make sure that categories for the models start at 0 and not 1.
models.trainLabels <- as.matrix(to_categorical(models.trainingtarget))
models.testLabels <- as.matrix(to_categorical(models.testtarget))

#model
model <- keras_model_sequential()
# Add layers to the model
model %>% 
  layer_dense(units = 28, activation = 'relu', input_shape = c(ncol(models)-1)) %>% 
  layer_dense(units = 28, activation = 'relu') %>% 
  layer_dense(units = 28, activation = 'relu') %>%
  layer_dense(units = 4, activation = 'softmax')

#Define the optimizer
sgd <- optimizer_sgd(learning_rate = 0.01)
# Compile the model
model %>% compile(
  #loss = 'binary_crossentropy',
  loss = 'categorical_crossentropy',
  #optimizer = 'adam',
  optimizer = 'sgd',
  metrics = 'accuracy'
)
# Fit the model 
model %>% fit(
  models.training, 
  models.trainLabels, 
  epochs = 500, 
  batch_size = nrow(models.trainLabels)/600, 
  validation_split = 0.05
)


# Predict the classes for the test data
classes <- model %>% predict_classes(models.test)
classes <- model %>% predict(models.test)

classes <- model |>
  predict(models.test) |>
  apply(1, which.max)
# Confusion matrix
conf <- table(models.testtarget, classes)
write.csv(conf, "1_confusion_matrix_Keras_old.csv")
# Evaluate on test data and labels
score <- model %>% evaluate(models.test, models.testLabels)
write.csv(score, "1_accuracy_Keras_old.csv")

model$save("my_model.keras")
model$save_weights("my_model_weights.weights.h5")

# Predict the classes for the observed data
stat_nnet <- as.matrix(stat2)
dimnames(stat_nnet) <- NULL

classes_obs <- model %>%
  predict(stat_nnet) %>%
  apply(1, which.max)

classes_obs <- classes_obs +1
probs <- model %>% predict(stat_nnet)
colnames(probs) <- c("m1", "m2", "m3","m4")

write.csv(probs, "1_model_selection_table.csv")


######Estimating parameters########
#reading the classification probability of each model
probs <- read.csv("1_model_selection_table.csv")
probs <- probs[-1]

#Selecting model with the highest probability
best_model <- colnames(probs[match(max(t(probs)), probs)])

mdl2 <- models_caret[grep(best_model, data_caret),]

models_params <- get(paste(best_model, ".sim", sep = ""))
models_params<-models_params[complete.cases(models_params),]

models_params <- models_params[, 1:c(grep("mean.rate", colnames(models_params))-1)]

models_params <- models_params[,colMeans(models_params) != 0]

dat2 <- models_params

#mean is used to normalize the data
mean_mdl2 <- colMeans(mdl2)
mean_dat2 <- colMeans(dat2)

#Normalizing
if (normalize == T) {
  mdl2 <- apply(mdl2,2,function(x){x/mean(x)})
  dat2 <- apply(dat2,2,function(x){x/mean(x)})
  stat3 <- stat/mean_mdl2
}

# first convert the data set to array
data <- as.matrix(dat2[,1:ncol(dat2)])
models <- as.matrix(mdl2[,1:ncol(mdl2)])

models <- cbind(models, data)
dimnames(models) <- NULL
set.seed(12345)

#####
##Number of replicates. It is a good idea to run Keras multiple times and 
#calculate the mean of the stats and parameter estimation

#Parameters
#number of replicates
reps <- c(1:10)
#setting the proportion for train and test data
train = 0.75
test = 0.25
epochs = 500


##Running keras
parameters_keras <- function(rep) {
Res <- NULL
##here goes the loop for each variable 1:9
#i=8
for (i in 1:ncol(dat2)) {
# Determine sample size
ind <- sample(2, nrow(models), replace=TRUE, prob=c(train, test))
models.training <- as.matrix(models[ind==1, 1:c(ncol(models) - ncol(data))])
dim(models.training)
models.test <- as.matrix(models[ind==2, 1:c(ncol(models) - ncol(data))])
dim(models.test)
models.trainingtarget <- as.matrix(models[ind==1, c(ncol(mdl2)+i)])
dim(models.trainingtarget)
models.testtarget <- as.matrix(models[ind==2, c(ncol(mdl2)+i)])
dim(models.testtarget)

model <- keras_model_sequential() %>%
  layer_dense(units = 28, activation = "relu",
              input_shape = dim(models.training)[2]) %>%
  layer_dense(units = 28, activation = "relu") %>%
  layer_dense(units = 28, activation = "relu") %>%
  layer_dense(units = 1, activation = "relu")

model %>% compile(
  loss = "mae",
  optimizer = optimizer_rmsprop(),
  metrics = list("mean_absolute_error"))


model %>% fit(
  models.training,
  models.trainingtarget,
  epochs = epochs,
  batch_size = nrow(models.trainLabels)/600,
  validation_split = 0.05,
  verbose = 0)

mae %<-% (model %>% evaluate(models.test, models.testtarget, verbose = 0))

#mae <- mae*mean_dat2[i]

test_pred <- model %>% predict(models.test)
test_pred <- test_pred*mean_dat2[i]
models.testtarget2 <- models.testtarget*mean_dat2[i]



cor.coef <- cor.test(test_pred, models.testtarget)$estimate

if (dir.exists("Correlations_obs_vs_sim")==F) {
  dir.create("Correlations_obs_vs_sim")
}
pdf(paste("Correlations_obs_vs_sim/",colnames(dat2)[i],"_", rep,"_IBD.pdf",sep=""))
plot(test_pred ~ models.testtarget2, xlab = "true", ylab="estimated", main = colnames(dat2)[i])
lines(c(0:1000000),(0:1000000), col=2)
dev.off()


#predicting parameters one by one
predictions <- model %>% predict(stat_nnet)
predictions <- predictions*mean_dat2[i]

z <- c(predictions, cor.coef, mae[1])
if (is.na(z[2])) {
  z[1:3] <- NA
}
names(z)<- c(paste("estimate_", rep, sep = ""), paste("r2_", rep, sep = ""), paste("MAE_", rep, sep = ""))
Res <- rbind(Res, z)

}
rownames(Res) <- colnames(dat2)[1:ncol(dat2)]
print(paste("replicate ",rep, " is done!", sep = ""))
return(Res)
}


#####Running
system.time({
 parameters_keras_obs <- lapply(reps, parameters_keras)
})

system.time({
 parameters_keras <- lapply(reps, function(rep) {
  model %>% fit(
   models.training, 
    models.trainLabels, 
     epochs = epochs, 
      batch_size = nrow(models.trainLabels) / 600,  # explicitly named
      validation_split = 0.05, 
      verbose = 0
    )
  })
})

library(matrixStats)
parameters_keras_obs2 <- do.call(cbind.data.frame, parameters_keras_obs)

esti <- parameters_keras_obs2[,grep("estimate", colnames(parameters_keras_obs2))]
esti_numeric <- data.frame(lapply(esti, function(x) if(is.numeric(x)) x else as.numeric(as.character(x))))
esti_mean <- as.data.frame(rowMeans(esti_numeric, na.rm = TRUE))
esti_SD <- as.data.frame(apply(esti_numeric, 1, sd, na.rm = TRUE))


r2 <- parameters_keras_obs2[,grep("r2", colnames(parameters_keras_obs2))]
r2_numeric <- data.frame(lapply(r2, function(x) as.numeric(as.character(x))))
r2_mean <- as.data.frame(rowMeans(r2_numeric, na.rm = T))

MAE <- parameters_keras_obs2[,grep("MAE", colnames(parameters_keras_obs2))]
MAE_numeric <- data.frame(lapply(MAE, function(x) as.numeric(as.character(x))))
MAE_mean <- as.data.frame(rowMeans(MAE_numeric, na.rm = T))

mean_keras_replicates <- cbind(esti_mean, esti_SD, r2_mean, MAE_mean)

colnames(mean_keras_replicates) <- c(paste("mean of ", max(reps), " replicates", sep = ""), "SD", "R^2", "MAE")


write.table(round(mean_keras_replicates, 3), "Parameter_estimation_Keras_Results.csv",
            row.names = TRUE, sep = ",")


#############
#############
#############
#############
#PCA
#############
#############
#############



library(ggplot2)
library(ggpubr)
PCA <- prcomp(rbind(models_pca, stat), center = T, scale. = T, retx=T)
scores <- data.frame(PCA$x[,1:ncol(PCA$x)])
p1x2<-ggplot(scores, aes(x=PC1, y=PC2))+ theme(legend.position = "none") +  theme(legend.title=element_blank()) +
  geom_point(aes(colour=c(data_pca,"observed"), size=c(data_pca,"observed"), shape=c(data_pca,"observed")))+
  scale_size_manual(values=c(0.7,0.7,0.7,0.7,4))+
  scale_colour_manual(values = c("green", "blue", "red", "purple","black")) + scale_shape_manual(values = c(1,2,3,4,16))

p1x3<-ggplot(scores, aes(x=PC1, y=PC3))+ theme(legend.position = "none") + theme(legend.title=element_blank()) +
  geom_point(aes(colour=c(data_pca,"observed"), size=c(data_pca,"observed"), shape=c(data_pca,"observed")))+
  scale_size_manual(values=c(0.7,0.7,0.7,0.7,4))+
  scale_colour_manual(values = c("green", "blue", "red", "purple", "black")) + scale_shape_manual(values = c(1,2,3,4,16))

p1x4<-ggplot(scores, aes(x=PC1, y=PC4))+ theme(legend.position = "none") + theme(legend.title=element_blank()) +
  geom_point(aes(colour=c(data_pca,"observed"), size=c(data_pca,"observed"), shape=c(data_pca,"observed")))+
  scale_size_manual(values=c(0.7,0.7,0.7,0.7,4))+
  scale_colour_manual(values = c("green", "blue", "red", "purple","black")) + scale_shape_manual(values = c(1,2,3,4,16))

p2x3<-ggplot(scores, aes(x=PC2, y=PC3))+ theme(legend.position = "none") + theme(legend.title=element_blank()) +
  geom_point(aes(colour=c(data_pca,"observed"), size=c(data_pca,"observed"), shape=c(data_pca,"observed")))+
  scale_size_manual(values=c(0.7,0.7,0.7,0.7,4))+
  scale_colour_manual(values = c("green", "blue", "red", "purple", "black")) + scale_shape_manual(values = c(1,2,3,4,16))

p2x4<-ggplot(scores, aes(x=PC2, y=PC4))+ theme(legend.position = "none") + theme(legend.title=element_blank()) +
  geom_point(aes(colour=c(data_pca,"observed"), size=c(data_pca,"observed"), shape=c(data_pca,"observed")))+
  scale_size_manual(values=c(0.7,0.7,0.7,0.7,4))+
  scale_colour_manual(values = c("green", "blue", "red", "purple","black")) + scale_shape_manual(values = c(1,2,3,4,16))

p3x4<-ggplot(scores, aes(x=PC3, y=PC4))+ theme(legend.position = "none") + theme(legend.title=element_blank()) +
  geom_point(aes(colour=c(data_pca,"observed"), size=c(data_pca,"observed"), shape=c(data_pca,"observed")))+
  scale_size_manual(values=c(0.7,0.7,0.7,0.7,4))+
  scale_colour_manual(values = c("green", "blue", "red", "purple","black")) + scale_shape_manual(values = c(1,2,3,4,16))

PCA_Plot <- ggarrange(p1x2,p1x3, p1x4, p2x3, p2x4, p3x4, legend = "left", common.legend = T) + ggtitle(paste(species))
PCA_Plot
pdf(paste(species, "_PCA.pdf", sep = ""), width=12, height=6.69)
PCA_Plot
dev.off()

#Plotting variable contributions for summary stats
library("FactoMineR")
library("factoextra")

res.pca <- PCA(rbind(models_pca), graph = FALSE, scale.unit = TRUE, ncp=10)
eig.val <- get_eigenvalue(res.pca)
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))
pdf("1_PC_contribution.pdf")
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))
dev.off()

var <- get_pca_var(res.pca)

library("corrplot")
corrplot(var$cos2, is.corr=FALSE, tl.cex=0.5)
pdf("1_PCs_vs_variables_cor.pdf")
corrplot(var$cos2, is.corr=FALSE, tl.cex=0.5)
dev.off()

