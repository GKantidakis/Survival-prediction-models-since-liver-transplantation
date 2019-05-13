
install_packages <- c("foreign", "Hmisc", "plyr")
for (i in 1:length(install_packages)){
  if (!install_packages[i] %in% installed.packages()){
    install.packages(install_packages[i])
  }
}
library(foreign)
library(Hmisc)
library(plyr)

# read data file in sav format
unos <- read.spss("UNOS_file.sav", to.data.frame = TRUE)
# translate characters from upper to lower case
names(unos) <- casefold(names(unos))

response_variables <- c("patientsurvival", "death",
                        "gs_ffs", "gs_ffsstate")

# call created lists with donor, recipient and response variables
unos_selected <- unos[, c(donor_vars, patient_vars, response_variables)]



# Missing values - imputation technique
install_packages <- c("mice", "randomForestSRC", "Amelia",
                      "UpSetR", "VIM")
for (i in 1:length(install_packages)){
  if (!install_packages[i] %in% installed.packages()){
    install.packages(install_packages[i])
  }
}
library(mice)
library(randomForestSRC)
library(Amelia)
library(UpSetR)
library(VIM)


# missingness per row
which.max(rowSums(is.na(unos_selected)))

# sum of NA in the full data set is 2.5%
sum(is.na(unos_selected)) /
  (nrow(unos_selected)*ncol(unos_selected)) 

colSums(is.na(unos_selected)) # missing values per variable
colSums(is.na(unos_selected))[colSums(is.na(unos_selected))
                              > 10000 ]  # 15 variables with more than 10000 NA


# show which variables have missing values
vars_mis <- colnames(unos_selected)[sapply(colnames(unos_selected),
                                           function(x) sum(is.na(unos_selected[, x]))) > 0]


##################################################
# Data imputation
##################################################

# check Nelson-Aalen estimates of cumulative hazard as a replacement of patient survival
# for overall survival
H0_t <- nelsonaalen(data = unos_selected,
                    timevar = patientsurvival,
                    statusvar = death)
cor(cbind(H0_t, Survival = unos_selected$patientsurvival)) 
# 0.989 correlation

# for failure-free survival
H0_t2 <- nelsonaalen(data = unos_selected,
                     timevar = gs_ffs,
                     statusvar = gs_ffsstate)
cor(cbind(H0_t2, Survival = unos_selected$gs_ffs))
# 0.988 correlation
# correlation almost 1, so for these data it matters little
# whether we take Ho(t) or T as a predictor.

##############################################################
##############################################################
# Original missForest algorithm for data imputation using
# unsupervised splitting.
# We grow a forest to impute the data. To proceed split statistics
# are calculated. If a node splits on a variable with
# missing data, the variable's missing data is imputed by randomly
# drawing values from non-missing in-bag data. The purpose of this
# is to make it possible to assign cases to daughter nodes based
# on the split. We used all possible variable combinations as
# responses, and split by the rest of the variables using  
# multivariate composite splitting.
# Missing data for responses are imputed by prediction.
# The process is repeated using a new set of variables for responses 
# (mutually exclusive to the previous fit), until all variables
# have been imputed. The procedure is repeated until convergence
# of the algorithm. Maximum number of iterations was set to 5.

# This is the most accurate of all imputation procedures that the
# randomForestSCR offers, but also by far the most
# computationally expensive one. 

set.seed(12345)
time_rf_imput <- system.time(
  unos_complete <- impute.rfsrc(data = unos_selected,
                                ntree = 150, mf.q = 1,
                                max.iter = 5, do.trace = TRUE)
)

time_rf_imput / 3600 # elapsed time 8 hr 50 mins



################################################################################
# Application of the methods
################################################################################


#############################################################
# Cox proportional hazard models
#############################################################


install_packages <- c("survival", "survminer", "rms", "ggplot2",
                      "grid", "gridExtra", "pec", "forestmodel")

for (i in 1:length(install_packages)){
  if (!install_packages[i] %in% installed.packages()){
    install.packages(install_packages[i])
  }
}

library(survival)
library(survminer)
library(rms)
library(ggplot2)
library(grid)
library(gridExtra)
library(pec)
library(forestmodel)
source("https://bioconductor.org/biocLite.R")
biocLite("survcomp")
library(survcomp)

load("unos_failure_free.RData")
unos_failure_free$gs_ffs <- unos_failure_free$gs_ffs + 
  (1/365)

# Create 12 discrete time points
unos_failure_free$gs_ffs <- as.numeric(cut(unos_failure_free$gs_ffs,
                                           breaks = 12))

set.seed(12345) 
N <- nrow(unos_failure_free) 
index <- sample(1:N, round(N/3), replace = FALSE)
training_ffs <- unos_failure_free[-index, ] # create training set
test_ffs <- unos_failure_free[index, ] # create test set
# quantiles of the training set
quantile(training_ffs$gs_ffs) 



# create Cox model with all the variables

cox_ffs <- coxph(Surv(gs_ffs, gs_ffsstate) ~ ., data = training_ffs, method = "breslow")
res <- summary(cox_ffs)

res_nolist <- unlist(res$coefficients)
str(res_nolist)
# the factors sorted by the z-statistic
res_nolist[order(res_nolist[, 4], decreasing = T), ] 

# testing the proportional hazards assumption
test_ph <- cox.zph(cox_ffs)
round(test_ph$table[, 3], digits = 3)

# 14 variables violate the PH assumption
round(test_ph$table[, 3], digits = 3)[which(test_ph$table[, 3] < 0.05)]

#############################################################
# Cox backward model
#############################################################

cox_back <- selectCox(formula = Surv(gs_ffs, gs_ffsstate) ~.,
                      data = training_ffs, rule = "aic")
class (cox_back)
vars_back <- cox_back$In

#############################################################
# Cox Lasso model
#############################################################

response_pair <- c("gs_ffs", "gs_ffsstate")

# for the training set
X_train <- training_ffs[, !(colnames(training_ffs) %in% response_pair)]
X_train_ffs <- model.matrix(~., X_train)[, -1] # create martix X for glmnet
# dim(X_train_ffs) # 41529 x 129
Y_train_ffs <- training_ffs[, colnames(training_ffs) %in% response_pair] # create matrix Y
Y_survtrain_ffs <- Surv(Y_train_ffs$gs_ffs, Y_train_ffs$gs_ffsstate)

# for the test set
x_test <- test_ffs[, !(colnames(test_ffs) %in% response_pair)]
X_test_ffs <- model.matrix(~., x_test)[, -1]
Y_test_ffs <- test_ffs[, colnames(test_ffs) %in% response_pair]
# Surv function packages survival data into the form expected by glmnet
Y_survtest_ffs <- Surv(Y_test_ffs$gs_ffs, Y_test_ffs$gs_ffsstate)

# use these two datasets to fit Cox model in the selected features
training_ffs_extended <- as.data.frame(cbind(X_train_ffs, Y_train_ffs))
test_ffs_extended <- as.data.frame(cbind(X_test_ffs, Y_test_ffs))


seed <- 12345
set.seed(seed)
folds <- createFolds(y = training_ffs$gs_ffs, k = 5, list = FALSE) 

cv.fit <- cv.glmnet(x = X_train_ffs, y = Y_survtrain_ffs,
                    foldid = folds, parallel = TRUE,
                    family = "cox", grouped = TRUE,
                    maxit = 1000)

active_index <- as.numeric(unlist(predict(cv.fit, newx = X_test_ffs,
                                          s = "lambda.1se",
                                          type = "nonzero")))
names_coefs <- colnames(test_ffs_extended)[1:128]
vars_active <- names_coefs[active_index]
vars_active # 36 variables are active


# Figure in appendix: Survival and censoring functions for Data 


unos_ffs_km <- survfit(formula = Surv(gs_ffs, gs_ffsstate) ~ 1, data = training_ffs) # 13034 events
unos_ffs_cens <- survfit(formula = Surv(gs_ffs, gs_ffsstate==0) ~ 1, data = training_ffs)

oldpar <- par(no.readonly=TRUE) # save graphical parameters

layout(matrix(1:2, 1, 2),widths=c(10.25,9))
par(mar= c(5, 4, 4, 0.1) + 0.1)
plot(unos_ffs_km, mark.time = FALSE, conf.int = FALSE, lwd=2, xlim = c(0, 12),
     xlab = "Time (years since transplantation)", ylab = "Probability")
title(main="Survival")
par(mar = c(5, 0.1, 4, 1) + 0.1)
plot(unos_ffs_cens, mark.time = FALSE, conf.int = FALSE, lwd = 2, xlim = c(0, 12),
     xlab = "Time (years since transplantation)", ylab = "", axes = FALSE)
axis(1)
box()
title(main="Censoring")
par(oldpar) # reset graphical parameters


#############################################################
# Random survival forests
#############################################################



install_packages <- c("survival", "caret", "randomForestSRC", 
                      "parallel", "prodlim", "pec")

for (i in 1:length(install_packages)){
  if (!install_packages[i] %in% installed.packages()){
    install.packages(install_packages[i])
  }
}
library(survival)
library(caret)
library(randomForestSRC)
library(parallel)
library(prodlim)
library(pec)
options(rf.cores = 20, mc.cores = 20)

#############################################################
# 5-fold cross validation        
#############################################################

nodesize <- c(10, 20, 35, 50, 70, 85, 100, 120)
nsplit <- c(3, 4, 5, 6, 7)
mtry <- c(11, 17, 23, 29, 35, 41, 47)
combis <- expand.grid(nodesize, nsplit, mtry)

nfolds <- 5
set.seed(12345)
folds <- createFolds(training_ffs$gs_ffsstate, k = 5, list = TRUE)
cv_time <- matrix(0, nrow = nfolds, ncol = nrow(combis))
oob_error <- matrix(0, nrow = nfolds, ncol = nrow(combis))


for (i in 1:nfolds) {
  
  cat("Starting iteration i = ", i, "\n")
  indices <- folds[[i]]
  train_set <- training_ffs[-indices, ]
  validation_set <- training_ffs[indices, ]
  
  for (j in 1:nrow(combis)) {
    
    cat("Working on combination:", "nodesize: ",
        combis[j, 1],
        "nsplit: ", combis[j, 2],
        "and mtry: ", combis[j, 3], "\n")
    cv_time[i, j] <- { 
      system.time(
        fit_ffs <- rfsrc(Surv(gs_ffs, gs_ffsstate) ~ .,
                         splitrule = "logrank",
                         nodesize = combis[j, 1],
                         nsplit = combis[j, 2],
                         data = train_set, mtry = combis[j, 3],
                         ntree = 500, seed = -12345, 
                         ntime = 1000, sampsize = 1000))[3]
    }
    cat("Calculating the OOB validated prediction error ...",
        "\n")
    fit_val <- predict(fit_ffs, newdata = validation_set,
                       seed = -12345)
    oob_error[i, j] <- as.numeric(fit_val$err.rate[fit_val$ntree])
    cat("Error is: ", oob_error[i, j])
    
  }
}

df_logrank_ffs <- data.frame(Node_size = combis$Var1,
                             Nsplit = combis$Var2,
                             Mtry = combis$Var3,
                             Error = round(colMeans(oob_error), 
                                           4))

save.image("image_logrank_ffs_cv.Rdata")   
cv_time / 60 # time in minutes
mean(cv_time / 60)
df_logrank_ffs[which.min(df_logrank_ffs$Error), ]
sd(df_logrank_ffs$Error)

a <- list(
  title = 'Node size',
  autotick = FALSE,
  ticks = "outside",
  dtick = 10,
  ticklen = 2,
  tickwidth = 1,
  tickcolor = toRGB("blue")
)

b <-  list(
  title = 'Mtry',
  autotick = TRUE,
  ticks = "outside",
  tick0 = -1,
  dtick = 6,
  ticklen = 2,
  tickwidth = 1,
  range = c(10, 47),
  tickcolor = toRGB("blue")
)
c <- list(
  title = 'Nsplit',
  autotick = FALSE,
  ticks = "outside",
  dtick = 1,
  ticklen = 2,
  tickwidth = 1,
  tickcolor = toRGB("blue")
)

#############################################################
# Figure 5.7: Grid search on a 3D space for RSF
#############################################################

#create 3D scatterplot
library(plotly)
plot_ly(df_logrank_ffs, x = ~Node_size, y = ~Mtry, z = ~Nsplit,
        marker = list(color = ~Error,
                      colorscale = "Viridis", # default is "Reds",
                      showscale = TRUE)) %>%
  add_markers() %>%
  layout(scene = list(xaxis = a,
                      yaxis = b,
                      zaxis = c),
         annotations = list(
           x = 1.10,
           y = 1.10,
           z = 0.90,
           text = 'mean OOB Error',
           xref = 'paper',
           yref = 'paper',
           showarrow = FALSE
         ))


###################################
# Finding OOB training error   
###################################    

nodesize <- c(10, 20, 35, 50, 70, 85, 100, 120)
nsplit <- c(3, 4, 5, 6, 7)
mtry <- c(11, 17, 23, 29, 35, 41, 47)
combis <- expand.grid(nodesize, nsplit, mtry)
cv_time <- vector(mode = "numeric", length = nrow(combis))
oob_error <- vector(mode = "numeric", length = nrow(combis))

for (j in 1:nrow(combis)) {
  
  cat("Working on combination nr:", j, "with nodesize: ",
      combis[j, 1],
      "nsplit: ", combis[j, 2], "and mtry: ",
      combis[j, 3],  "\n")
  
  cv_time[j] <- { 
    system.time(
      fit_ffs <- rfsrc(Surv(gs_ffs, gs_ffsstate) ~ .,
                       splitrule = "logrank",
                       nodesize = combis[j, 1],
                       nsplit = combis[j, 2],
                       data = training_ffs,
                       mtry = combis[j, 3],
                       ntree = 500, seed = -12345,
                       sampsize = 5000,
                       ntime = 500, forest = FALSE))[3]  
  }
  
  cat("Calculating the OOB prediction error of combination  ...",
      "\n")
  oob_error[j] <- as.numeric(fit_ffs$err.rate[fit_ffs$ntree])
  cat("Error is: ", oob_error[j], "\n")
}


df_logrank_ffs <- data.frame(Node_size = combis$Var1,
                             Nsplit = combis$Var2,
                             Mtry = combis$Var3,
                             Error = round(oob_error, 4))


##############################################
# Finding OOB error per number of trees
##############################################  

fit_block <- rfsrc(Surv(gs_ffs, gs_ffsstate) ~ .,
                   splitrule = "logrank", nsplit = 5,
                   data = training_ffs, ntree = 500, 
                   split.depth = "all.trees",
                   var.used = "all.trees", seed = -12345,
                   mtry = 41, nodesize = 50,
                   ntime = 500, sampsize = 5000, block.size = 5, 
                   forest = TRUE
)

jpeg("plot_trees_ffs.jpg")
plot(fit_block)
dev.off()
cat("Plot of number of trees recorded!")

################################
# Fitting final model        
################################ 

cat("Fitting random forest... ")
# final fit ffs
fit <- rfsrc(Surv(gs_ffs, gs_ffsstate) ~ .,
             splitrule = "logrank", nsplit = 5,
             data = training_ffs, ntree = 300,
             split.depth = "all.trees",
             var.used = "all.trees", seed = -12345, 
             mtry = 41, nodesize = 50,
             ntime = 500, sampsize = 5000, forest = TRUE,
             importance = TRUE
)

fit$err.rate[fit$ntree] 

jpeg("plot_oob_mort.png", width = 600)
plot.survival.rfsrc(fit, plots.one.page = FALSE)
dev.off()


#######################################################
# number of terminal nodes for each tree in the forest
#######################################################
leaf_count <- fit$leaf.count 
jpeg("plot_nr_leaves.png")
ggplot(NULL, aes(x = leaf_count)) +
  geom_histogram(aes(y =..density..), 
                 bins = 20, 
                 fill=I("lightblue"), 
                 color = I("blue")) +
  geom_density(color = 2) +
  labs(x = "Number of leaves per tree") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

##################################################
# most frequent variables used for the forest
##################################################

vars_used <- tail(sort(fit$var.used), 10)
jpeg("plot_var_used.png", width = 750)
barplot(tail(sort(fit$var.used), 5), ylab = "Frequency")
dev.off()

# very important variables are used early
jpeg("plot_split_depth.png", width = 650)
hist(fit$split.depth, xlab = "Split depth per variable",
     col = "red", main = "", breaks = 20)
dev.off()



#####################################################
# Most important variables using VIMP on training set
#####################################################
vars2 <- fit$importance
head(sort(vars2, decreasing = TRUE), 15)
length(which(vars2 > 0.0005))

a <- head(sort(vars2, decreasing = TRUE), 20)
df <- data.frame(Variable = attr(a, "names")[1:10],
                 VIMP = round(a, 4)[1:10],
                 Variable = attr(a, "names")[11:20],
                 VIMP = round(a, 4)[11:20])
rownames(df) <- NULL


##############################################
# error on test set (1 - concordance)
##############################################

fit_test$err.rate[fit_test$ntree]

# always put the times in sorted order
probs_rf_ffs <- predictSurvProb(object = fit,
                                newdata = test_ffs,
                                times =
                                  sort(unique(test_ffs$gs_ffs)))

save(probs_rf_ffs, file = "probs_rf_ffs.Rdata")

metrics_rf_ffs <-  metrics_ffs_su(probs_rf_ffs, test_ffs)
pair_rfsrc <- as.data.frame(cbind(test_ffs$gs_ffs, 
                                  test_ffs$gs_ffsstate))
colnames(pair_rfsrc) <- c("time", "status")


library(doParallel)
if (!exists("cl")) {
  cl <- makeCluster(20)
  registerDoParallel(cl)
}

surv_f <- as.formula(Surv(gs_ffs, gs_ffsstate) ~ .)
pec_f <- as.formula(Hist(gs_ffs, gs_ffsstate) ~ 1)
## run cox/rfsrc models
## for illustration we use a small number of trees
cox_full <- coxph(surv_f, data = training_ffs,
                  x = TRUE, y = TRUE)




#######################################################################
# Survival neural networks
#######################################################################

##################################################
# Neural networks one hidden layer part 1
##################################################

# Hint: to download keras library for R you need to 
# have installed the Anaconda3

install_packages <- c("survival", "keras", "pec", 
                      "caret", "e1071", "doParallel")

for (i in 1:length(install_packages)){
  if (!install_packages[i] %in% installed.packages()){
    install.packages(install_packages[i])
  }
}

library(keras)
library(survival)
library(pec)
library(caret)
library(e1071)
library(doParallel)
cl <- makeCluster(3)
registerDoParallel(cl)
# parallel::stopCluster(cl)

# set this at the beginning of R session to get
# reproducible results (it is tensor-flow seed) 
# it disables GPU and enables parallelism only on CPU
# random number generator in tensor flow backend engine
use_session_with_seed(seed = 1235, disable_gpu = TRUE,
                      disable_parallel_cpu = TRUE,
                      quiet = FALSE)

# read the scaled training data for ffs survival
training_ffs_scaled <- read.table("training_ffs_scaled.txt")
# read the scaled test data for ffs survival
test_ffs_scaled <- read.table("test_ffs_scaled.txt")

# function that transforms the train data into long format
# for each patient create the interval (1 till time 12 years),
# survival and create an id number

data_train_creator <- function(data){
  
  N <- nrow(data)
  # assign survival times to 12 intervals, each is a 1-year period
  data$interval <- as.numeric(cut(data$gs_ffs,
                                  breaks = 12)) 
  data$survival <- as.numeric(cut(data$gs_ffs,
                                  breaks = 12)) 
  data$id <- 1:N
  n.times <- data$interval
  data_long <- data[rep(seq_len(N), times = n.times), ]
  
  # create the correct intervals
  for(i in unique(data_long$id)) {
    n_length <- length(data_long$interval[data_long$id == i])
    data_long$interval[data_long$id == i] <- 1:n_length
  }
  
  data_long$status <- vector(mode = "numeric",
                             length = nrow(data_long))
  
  # put indication 1 on status at the interval that patient dies
  for (i in 1:nrow(data_long)) {
    if (data_long$gs_ffsstate[i] != data_long$status[i] &&
        data_long$survival[i] == data_long$interval[i])
      data_long$status[i] <- 1
  }
  
  return(data_long)
  
}

# function that creates data in the right long format for
# test set for each patient the interval goes from 1 year
# till 12 years

data_test_creator <- function(data){
  
  N <- nrow(data)
  # assign survival times to 12 intervals
  data$interval <- max(as.numeric(cut(data$gs_ffs,
                                      breaks = 12)))
  # the true interval survival
  data$survival <- as.numeric(cut(data$gs_ffs,
                                  breaks = 12)) 
  data$id <- 70001:(70000 + N) # define the patient ids abstractly
  n.times <- data$interval
  data_long <- data[rep(seq_len(N), times = n.times), ]
  
  # create the correct intervals
  for(i in unique(data_long$id)) {
    n_length <- length(data_long$interval[data_long$id == i])
    data_long$interval[data_long$id == i] <- 1:n_length
    
  }
  
  data_long$status <- vector(mode = "numeric",
                             length = nrow(data_long))
  
  # put indication 1 on status at the intervals on
  # which a patient has died
  for (i in 1:nrow(data_long)) {
    if (data_long$gs_ffsstate[i] == 1 && 
        data_long$survival[i] <= data_long$interval[i]) 
      data_long$status[i] <- 1
  }
  
  return(data_long)
  
}


# function to calculate Brier and Integrated Brier scores for NN
brier_nnet <- function(pair_df, prob_matrix) {
  # calculating censoring distribution
  so <-  Surv(pair_df$true_surv , pair_df$status) 
  time <-  so[, 1]
  ot <-  order(time)
  cens <-  so[ot, 2 ]
  time <-  time[ot]
  N <-  nrow(so)
  hatcdist <- prodlim(Surv(time, cens) ~ 1, reverse = TRUE) 
  # censoring weights (reverse Kaplan-Meier)
  csurv <-  predict(hatcdist, times = time, type = "surv") 
  csurv [csurv == 0] <-  Inf # access infinite value to censored
  btime <-  time
  # matrix transpose of survival predictions
  survs <-  t(as.matrix(prob_matrix)) 
  btime.n <-  unique(time) # find the unique times
  bsc <- rep(0, length(btime.n))
  
  for (j in 1:length(btime.n)) {
    # indicator vectors for selecting relevant patients
    # the censored ones
    values1 <- as.integer(time <= btime.n[j] & cens == 1) 
    values2 <- as.integer(time > btime.n[j] )
    
    inb <- survs[j , ]
    inb[which(inb == 1)] <- inb[which(inb == 1)] 
    inb[which(inb == 0)] <- inb[which(inb == 0)]
    # Brier: bsc
    bsc[j] <- mean((0 - survs[j, ] )^2 * values1 * (1/csurv) +
                     (1 - survs [j, ])^2 * values2 * (1/csurv[j]))
  }
  bsc <- c(0, bsc) # Brier score from time 0
  
  # calculate integrated brier score
  timing <- 0:12
  indexa <- 2:length(timing)
  int_bs <- (diff(timing) %*% ((bsc[indexa - 1] + 
                                  bsc[indexa])/2))/diff(range(timing))
  
  return(list(Brier = bsc[-1], Int_brier = int_bs))
}


################################################################
# function that calculates weights, node size, cross entropy, 
# accuracy, sensitivity, specificity
# Brier score and integrated Brier score: uses function brier_nnet
# arguments are trained_model: aka the fit_keras,
# validation_set and real_status that is the event status
################################################################

measures_calculator <- function(trained_model, datanew, real_status) {
  
  df1 <- data.frame(hazard = predict_proba(trained_model,
                                           as.matrix(datanew[,
                                                             c(1:128, 131)]),
                                           batch_size = 500))
  df1$id <- datanew$id # ids of the patients
  df1$survival <- datanew$survival # survival time in years
  # create personalized groups per id
  groups <- split(df1, f = df1$id)
  
  true_surv <- unlist(lapply(groups, function(x) {
    surv_obj <- x$survival
    true_res <- surv_obj[1] 
    return(true_res)}
  ))
  group_probs <- lapply(groups, function(x) 
  {x <- cumprod(1 - x$hazard)})
  pred_mat <- do.call("rbind", group_probs)
  
  relative_probs <- vector(mode = "numeric",
                           length = length(group_probs))
  for (i in 1:length(group_probs)){ 
    temp <- group_probs[[i]]
    ind <- which(1:12 == true_surv[i])
    # find the relative survival probability of the
    # real survival interval
    relative_probs[i] <- temp[ind] 
  }
  N0 <- length(unique(df1$id)) # number of unique persons
  
  # in the data frame create random id numbers
  # to label the patients
  df2 <- data.frame(relative_probs, true_surv,
                    status = real_status,
                    id = (70001):(70000 + N0))
  df2$prediction <- 1 - round(df2$relative_probs, digits = 0)
  # auxiliary confusion table
  tabel <- table(df2$prediction, df2$status) 
  confusion_mat <- confusionMatrix(tabel, positive = "1")
  accuracy <- sum(df2$prediction == df2$status) / nrow(df2)
  sensitivity <- sum(df2$status == 1 & df2$prediction == 1) / 
    colSums(tabel)[2]
  specificity <- sum(df2$status == 0 & df2$prediction == 0) / 
    colSums(tabel)[1]
  precision <- as.numeric(confusion_mat$byClass[5])
  recall <- as.numeric(confusion_mat$byClass[6])
  f1score <- as.numeric(confusion_mat$byClass[7])
  
  brier_obj <- brier_nnet(df2, prob_matrix = pred_mat)
  int_brier <- as.numeric(brier_obj$Int_brier)
  brier_set <- brier_obj$Brier
  all_weights <- get_weights(fit_keras)
  nr_weights <- (length(all_weights[[1]]) + 
                   length(all_weights[[2]])
                 + length(all_weights[[3]]) + 
                   length(all_weights[[4]])) 
  
  return(list(weights = nr_weights,
              node_size = length(get_weights(fit_keras)[[2]]),
              cross_entropy = 
                result$metrics$loss[result$params$epochs],
              accuracy = accuracy,
              sensitivity = as.numeric(sensitivity),
              specificity = as.numeric(specificity),
              Precision = precision,
              Recall = recall,
              F1score = f1score,
              Integrated_brier = int_brier,
              Brier_scores = brier_set))
}

###########################################
# set up the cross-validation
###########################################

# most popular optimization algorithms used are the Stochastic 
# Gradient Descent (SGD), ADAM and RMSprop
# you need to tune certain parameters such as learning rate or 
# momentum
# things that can be tuned are node_size, lr, regularizer_l2,
# epochs, batch_size, decay

nfolds <- 5
set.seed(12345)
folds <- createFolds(training_ffs_scaled$gs_ffs,
                     k = 5, list = TRUE)

node_size <- seq(10, 100, by = 10) # grid of node sizes
dropout_rate <- c(0.2, 0.3, 0.4)
lr <- c(0.01, 0.1)
class_weights <- c(1, 2)
momentum <- c(0.8, 0.9)
combis <- expand.grid(node_size, dropout_rate,
                      lr, class_weights, momentum)

# initialize objects
cv_error <- matrix(0, nrow = nfolds, ncol = nrow(combis))
cv_accuracy <- matrix(0, nrow = nfolds, ncol = nrow(combis))
cv_specificity <- matrix(0, nrow = nfolds, ncol = nrow(combis))
cv_sensitivity <- matrix(0, nrow = nfolds, ncol = nrow(combis))
cv_precision <- matrix(0, nrow = nfolds, ncol = nrow(combis))
cv_recall <- matrix(0, nrow = nfolds, ncol = nrow(combis))
cv_f1score <- matrix(0, nrow = nfolds, ncol = nrow(combis))
cv_weights <- matrix(0, nrow = nfolds, ncol = nrow(combis))
cv_intbrier <- matrix(0, nrow = nfolds, ncol = nrow(combis))


# Learning rate controls how much to update the weight
# at the end of each batch and the momentum controls how
# much to let the previous update influence the current
# weight update. The number of neurons in a layer is an
# important parameter to tune. Generally the number of neurons
# in a layer controls the representational capacity of the 
# network, at least at that point in the topology.

for (i in 1:nfolds) {
  
  cat("Started iteration i = ", i, "\n")
  indices <- folds[[i]]
  cat("Creating the training set ...", "\n")
  # create the train set
  train_set <-  data_train_creator(
    training_ffs_scaled[-indices, ]) 
  cat("Creating the validation set ...", "\n")
  validation_set <- data_test_creator(data = 
                                        training_ffs_scaled[indices, ]) # create the validation set
  # real status for the validation set
  event_status <- training_ffs_scaled[indices, ]$gs_ffsstate 
  
  # create the matrices to be used for keras library
  # predictors: 128 variables + interval
  train_x <- as.matrix(train_set[, c(1:128, 131)]) 
  dimnames(train_x) <- NULL # the object must have empty dimnames
  train_y <- train_set$status
  validation_x <- as.matrix(validation_set[, c(1:128, 131)])
  # # the object must have empty dimnames
  dimnames(validation_x) <- NULL 
  validation_y <- validation_set$status
  
  for (j in 1:nrow(combis)) {
    
    cat("Testing combination number:", j, 
        "of repeat", i, " out of 5", "\n")
    cat("calculating for node size:", combis[j, 1], ", 
    dropout rate:", combis[j, 2], "\n",
        "and learning rate", combis[j, 3], "and weak class weight",
        combis[j, 4],
        "and momentum", combis[j, 5], "...", "\n")
    
    # start building the model
    fit_keras <- keras_model_sequential()
    # Add layers to the model
    # we create a densely connected ANN to the output
    fit_keras %>%
      layer_dense(units = combis[j, 1], activation = 'relu', 
                  input_shape = c(129)) %>% 
      layer_dropout(rate = combis[j, 2]) %>%
      layer_dense(units = 1, activation = 'sigmoid')
    # for binary class classification problem
    fit_keras %>% compile(
      loss = 'binary_crossentropy', 
      optimizer = optimizer_sgd(lr = combis[j, 3],
                                momentum = combis[j, 5])
    )
    
    result <- fit_keras %>% fit(
      train_x, 
      train_y, 
      epochs = 10, 
      batch_size = 500,
      validation_data = list(validation_x, validation_y),
      class_weight = list("0" = 1, "1" = combis[j, 4])
      #,callbacks = c(early_stopping)
    )
    
    # now that the model has run lets calculate the measures
    # the total weights are 129*node_size + bias_input 
    # + node_size + bias_node_size 
    values <- measures_calculator(trained_model = fit_keras,
                                  datanew = validation_set,
                                  real_status = event_status)
    
    cv_error[i, j] <- values$cross_entropy 
    cv_weights[i, j] <- values$weights
    cv_accuracy[i, j] <- values$accuracy
    cv_sensitivity[i, j] <- values$sensitivity
    cv_specificity[i, j] <- values$specificity
    cv_precision[i, j] <- values$Precision
    cv_recall[i, j] <- values$Recall
    cv_f1score[i, j] <- values$F1score
    cv_intbrier[i, j] <- values$Integrated_brier
    
  }
}

df_relu_ffs <- round(cbind(node_size = combis[, 1],
                           dropout_rate = combis[, 2],
                           learning_rate = combis[, 3],
                           momentum = combis[, 5],
                           weak_weight = combis[, 4],
                           weights = colMeans(cv_weights),
                           cross_entropy = colMeans(cv_error), 
                           accuracy = colMeans(cv_accuracy), 
                           sensitivity = 
                             colMeans(cv_sensitivity),
                           specificity =
                             colMeans(cv_specificity),
                           precision = colMeans(cv_precision),
                           recall = colMeans(cv_recall),
                           f1score = colMeans(cv_f1score),
                           integrated_brier = 
                             colMeans(cv_intbrier)), digits = 3) 

save(df_relu_ffs, file = "results_relu_ffs.Rdata")
save.image("image_relu_ffs.Rdata")

parallel::stopCluster(cl)

df_relu_ffs <- as.data.frame(df_relu_ffs)


#######################################################
# fitting the final models for survival neural networks 
#######################################################

# Hint: to download keras library you need to 
# have installed the Anaconda3 and TensorFlow for R

install_packages <- c("survival", "keras", "pec", 
                      "caret", "e1071", "gridExtra",
                      "grid", "ggpubr")

for (i in 1:length(install_packages)){
  if (!install_packages[i] %in% installed.packages()){
    install.packages(install_packages[i])
  }
}
library(keras)
library(survival)
library(pec)
library(caret)
library(e1071)
library(gridExtra)
library(grid)
library(ggpubr)

use_session_with_seed(seed = 1235, disable_gpu = TRUE,
                      disable_parallel_cpu = TRUE,
                      quiet = FALSE)


# read the scaled training data for ffs survival
training_ffs_scaled <- read.table("training_ffs_scaled.txt")
# read the scaled test data for ffs survival
test_ffs_scaled <- read.table("test_ffs_scaled.txt")

load("training_ffs_long.Rdata")
load("test_ffs_long.Rdata")

# function to calculate Brier and Integrated Brier scores for nnet
brier_nnet <- function(pair_df, prob_matrix) {
  # calculating censoring distribution
  so <-  Surv(pair_df$true_surv , pair_df$status) 
  time <-  so[, 1]
  ot <-  order(time)
  cens <-  so[ot, 2 ]
  time <-  time[ot]
  N <-  nrow(so)
  hatcdist <- prodlim(Surv(time, cens) ~ 1, reverse = TRUE) 
  # censoring weights (reverse Kaplan-Meier)
  csurv <-  predict(hatcdist, times = time, type = "surv") 
  csurv [csurv == 0] <-  Inf # access infinite value to censored
  btime <-  time
  # matrix transpose of survival predictions
  survs <-  t(as.matrix(prob_matrix)) 
  btime.n <-  unique(time) # find the unique times
  bsc <- rep(0, length(btime.n))
  
  for (j in 1:length(btime.n)) {
    # indicator vectors for selecting relevant patients
    # the censored ones
    values1 <- as.integer(time <= btime.n[j] & cens == 1) 
    values2 <- as.integer(time > btime.n[j] )
    
    inb <- survs[j , ]
    inb[which(inb == 1)] <- inb[which(inb == 1)] 
    inb[which(inb == 0)] <- inb[which(inb == 0)]
    # Brier: bsc
    bsc[j] <- mean((0 - survs[j, ] )^2 * values1 * (1/csurv) +
                     (1 - survs [j, ])^2 * values2 * (1/csurv[j]))
  }
  bsc <- c(0, bsc) # Brier score from time 0
  
  # calculate integrated brier score
  timing <- 0:12
  indexa <- 2:length(timing)
  int_bs <- (diff(timing) %*% ((bsc[indexa - 1] + 
                                  bsc[indexa])/2))/diff(range(timing))
  
  return(list(Brier = bsc[-1], Int_brier = int_bs))
}

# function that calculates all measures

measures_calculator <- function(trained_model,
                                datanew, real_status) {
  
  df1 <- data.frame(hazard = predict_proba(trained_model,
                                           as.matrix(datanew[, 
                                                             c(1:128, 131)]),
                                           batch_size = 500))
  df1$id <- datanew$id # ids of the patients
  df1$survival <- datanew$survival # survival time in years
  # create personalized groups per id
  groups <- split(df1, f = df1$id) 
  
  true_surv <- unlist(lapply(groups, function(x) {
    surv_obj <- x$survival
    true_res <- surv_obj[1] 
    return(true_res)}
  ))
  group_probs <- lapply(groups, function(x) {
    x <- cumprod(1 - x$hazard)})
  pred_mat <- do.call("rbind", group_probs)
  
  relative_probs <- vector(mode = "numeric",
                           length = length(group_probs))
  for (i in 1:length(group_probs)){ 
    temp <- group_probs[[i]]
    ind <- which(1:12 == true_surv[i])
    # find the relative survival probability of
    # the real survival interval
    relative_probs[i] <- temp[ind] 
  }
  N0 <- length(unique(df1$id)) # number of unique persons
  
  
  # in the data frame create random id numbers
  # to label the patients
  df2 <- data.frame(relative_probs, true_surv,
                    status = real_status, id = (70001):(70000 + N0))
  df2$prediction <- 1 - round(df2$relative_probs, digits = 0)
  # create possible classes
  classes <- c(0, 1)
  # auxiliary confusion table
  tabel <- table(factor(df2$prediction, levels = classes),
                 factor(df2$status, levels = classes)) 
  confusion_mat <- confusionMatrix(tabel, positive = "1")
  accuracy <- sum(df2$prediction == df2$status) / nrow(df2)
  sensitivity <- sum(df2$status == 1 & df2$prediction == 1) / 
    colSums(tabel)[2]
  specificity <- sum(df2$status == 0 & df2$prediction == 0)
  / colSums(tabel)[1]
  precision <- as.numeric(confusion_mat$byClass[5])
  recall <- as.numeric(confusion_mat$byClass[6])
  f1score <- as.numeric(confusion_mat$byClass[7])
  
  brier_obj <- brier_nnet(df2, prob_matrix = pred_mat)
  int_brier <- as.numeric(brier_obj$Int_brier)
  brier_set <- brier_obj$Brier
  all_weights <- get_weights(fit_keras)
  nr_weights <- (length(all_weights[[1]]) + 
                   length(all_weights[[2]]) 
                 + length(all_weights[[3]]) + 
                   length(all_weights[[4]])) 
  
  return(list(weights = nr_weights,
              node_size = length(get_weights(fit_keras)[[2]]),
              cross_entropy = 
                result$metrics$loss[result$params$epochs],
              accuracy = accuracy,
              sensitivity = as.numeric(sensitivity),
              specificity = as.numeric(specificity),
              Precision = precision,
              Recall = recall,
              F1score = f1score,
              Integrated_brier = int_brier,
              Brier_scores = brier_set))
}

fit_keras <- keras_model_sequential()
# Add layers to the model
# we create a densely connected ANN to the output
fit_keras %>%
  layer_dense(units = 75, activation = 'sigmoid', 
              input_shape = c(129)) %>% 
  layer_dropout(rate = 0.1) %>%
  layer_dense(units = 1, activation = 'sigmoid')
fit_keras %>% compile(
  loss = 'binary_crossentropy', # for binary class classification
  optimizer = optimizer_sgd(lr = 0.2, momentum = 0.9)
  #, metrics = c("accuracy")
)


event_status <- test_ffs_scaled$gs_ffsstate 

# create the matrices to be used for keras library
# predictors: 128 variables + interval
train_x <- as.matrix(training_ffs_long[, c(1:128, 131)]) 
dimnames(train_x) <- NULL # the object must have empty dimnames
train_y <- training_ffs_long$status
test_x <- as.matrix(test_ffs_long[, c(1:128, 131)])
dimnames(test_x) <- NULL # the object must have empty dimnames
test_y <- test_ffs_long$status

result <- fit_keras %>% fit(
  train_x, 
  train_y, 
  epochs = 10, 
  batch_size = 1000,
  validation_data = list(test_x, test_y),
  class_weight = list("0" = 1, "1" = 1)
)

# now that the model has run lets calculate the measures
values_final <- measures_calculator(trained_model = fit_keras,
                                    datanew = test_ffs_long,
                                    real_status = event_status)

df1 <- data.frame(hazard = predict_proba(fit_keras,
                                         as.matrix(test_ffs_long[, 
                                                                 c(1:128, 131)]),
                                         batch_size = 500))
df1$id <- test_ffs_long$id # ids of the patients
df1$survival <- test_ffs_long$survival # survival time in years
# create personalized groups per id
groups <- split(df1, f = df1$id) 

true_surv <- unlist(lapply(groups, function(x) {
  surv_obj <- x$survival
  true_res <- surv_obj[1]
  return(true_res)}
))
group_probs <- lapply(groups, function(x) {
  x <- cumprod(1 - x$hazard)})
pred_mat_nn1h <- do.call("rbind", group_probs)


#####################################
# mean survival probabilities
#####################################

plot1 <- rbind(df1, df2, df3, df4)

obj1 <- ggplot(plot1, aes(x = Time, y = Probs, color = Model)) + 
  geom_line(size = 0.9) + 
  scale_x_continuous(breaks = 0:12) + 
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) + 
  xlab("Time since transplantation in years") + 
  ylab("Survival probability") + 
  theme_classic() + 
  ggtitle(" ") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values=c("black", "red",
                              "green", "blue"), 
                     name = "NN 1 hidden")


##########################################################
# this is Garsons's connection weight method
# you can use it to get interpretation of the weights of
# the neural network. It returns a list with the variable
# names and the the relative importance measure
###########################################################

var_imp_keras <- function(model){
  
  list_weights <- get_weights(model)
  # list_weights[[1]] are the weights between input-hidden node
  # list_weights[[2]] are the weights between bias-hidden node
  # list_weights[[3]] are the weights between hidden node-output
  # list_weights[[4]] is the weight between bias-output
  # names input nodes , number input nodes
  names <- c(colnames(training_ffs_scaled)[1:128], "interval") 
  n_input <- length(list_weights[[1]]) / length(list_weights[[2]])
  # number hidden nodes , number output nodes
  n_nodes <- length(list_weights[[3]])
  n_outputs <- length(list_weights[[4]])
  
  # dimensions for weight matrices
  ncols1 <- (n_input + 1) 
  nrows1 <- n_nodes
  ncols2 <- (n_nodes + 1) 
  nrows2 <- n_outputs
  length1 <- ncols1 * nrows1 
  length2 <- ncols2 * nrows2
  
  # selecting weights
  mat1 <- rbind(list_weights[[2]], list_weights[[1]])
  mat1 <- t(mat1)
  weights1 <- mat1
  colnames(weights1) <- c("Bias", names)
  rownames(weights1) <- paste ("N" , seq (1:n_nodes) , sep = "")
  
  mat2 <- rbind(list_weights[[4]], list_weights[[3]])
  weights2 <- mat2
  rownames(weights2) <- c("Bias", paste ("N", seq (1:n_nodes),
                                         sep = ""))
  # calculating variable imp using connection weight method
  mega_mat <- matrix(NA, ncol = n_input, nrow = n_outputs)
  for (i in 1:n_input) {
    for (j in 1:n_outputs) {
      mega_mat[j, i] <- sum(weights1[, i + 1] * 
                              weights2[2:(n_nodes + 1), j ])
    }
  }
  colnames(mega_mat) <- names
  mega_mat_abs <- abs(mega_mat)
  totals <- rowSums (mega_mat_abs)
  mega_mat_rel <- as.data.frame(mega_mat_abs/ totals)
  rels <- as.vector(as.numeric(mega_mat_rel))
  
  return(list(Names = names, Rels = rels))
}

# find the variable importance of the last fit_keras model
var_imp <- var_imp_keras(model = fit_keras)

df2 <- data.frame(name = var_imp$Names,
                  variable_importance = var_imp$Rels)
df2 <- df2[order(df2$variable_importance, decreasing = TRUE), ]

rel_imp <- data.frame(Variable = df2$name[1:10],
                      Importance = df2$variable_importance[1:10],
                      Variable = df2$name[11:20],
                      Importance = df2$variable_importance[11:20])



###############################################
# final model with 2 hidden layers
###############################################
fit_keras <- keras_model_sequential()
# Add layers to the model
# we create a densely connected ANN to the output
fit_keras %>%
  layer_dense(units = 100, activation = 'tanh',
              input_shape = c(129)) %>% 
  layer_dropout(rate = 0.1) %>%
  layer_dense(units = 100, activation = 'tanh') %>%
  layer_dropout(rate = 0.1) %>%
  layer_dense(units = 1, activation = 'sigmoid')
fit_keras %>% compile(
  loss = 'binary_crossentropy', 
  optimizer = optimizer_sgd(lr = 0.2, momentum = 0.9)
)

event_status <- test_ffs_scaled$gs_ffsstate 

# create the matrices to be used for keras library
# predictors: 128 variables + interval
train_x <- as.matrix(training_ffs_long[, c(1:128, 131)]) 
dimnames(train_x) <- NULL # the object must have empty dimnames
train_y <- training_ffs_long$status
test_x <- as.matrix(test_ffs_long[, c(1:128, 131)])
dimnames(test_x) <- NULL # the object must have empty dimnames
test_y <- test_ffs_long$status

result <- fit_keras %>% fit(
  train_x, 
  train_y, 
  epochs = 10, 
  batch_size = 1000,
  validation_data = list(test_x, test_y),
  class_weight = list("0" = 1, "1" = 1)
)

# calculate the measures....
values_final <- measures_calculator(trained_model = fit_keras,
                                    datanew = test_ffs_long,
                                    real_status = event_status)

######################################################################
# Comparisons of the methods
######################################################################


##################################
# file "the_functions.R"     
##################################

install_packages <- c("survival", "pec", 
                      "caret", "e1071")

for (i in 1:length(install_packages)){
  if (!install_packages[i] %in% installed.packages()){
    install.packages(install_packages[i])
  }
}

# libraries to call
library(survival)
library(pec)
library(caret)
library(e1071)

# function that gets the probability matrix and returns
# the relative metrics

metrics_ffs <- function(prob_matrix, test_set, p = 0.5) {
  
  relative_probs <- vector(mode = "numeric", 
                           length = nrow(test_set))
  time <- sort(test_set$gs_ffs)
  for (i in 1:length(relative_probs)){ 
    temp <- prob_matrix[i, ]
    ind <- max(which(time == test_set$gs_ffs[i]))
    relative_probs[i] <- temp[ind]
  }
  
  # find expected status based on 0.5 probability cut-off
  expected_status <- ifelse(relative_probs > p, 0, 1)
  classes <- c(0, 1)
  # auxiliary confusion table
  tabel <- table(factor(test_set$gs_ffsstate, levels = classes),
                 factor(expected_status, levels = classes)) 
  confusion_mat <- confusionMatrix(tabel, positive = "1")
  confusion_mat <- confusionMatrix(tabel, positive = "1")
  accuracy <- as.numeric(confusion_mat$ffs[1])
  sensitivity <- as.numeric(confusion_mat$byClass[1])
  specificity <- as.numeric(confusion_mat$byClass[2])
  precision <- as.numeric(confusion_mat$byClass[5])
  recall <- as.numeric(confusion_mat$byClass[6])
  f1score <- as.numeric(confusion_mat$byClass[7])
  
  return(list(Accuracy = accuracy, Sensitivity = sensitivity,
              Specificity = specificity, Precision = precision,
              Recall = recall, F1score = f1score))
  
}

# this should be used in case you find the survival probabilities 
# at sort-unique time points
metrics_ffs_su <- function(prob_mat_su, test_set, p = 0.5) {
  
  relative_probs <- vector(mode = "numeric",
                           length = nrow(test_set))
  time <- sort(unique(test_set$gs_ffs))
  for (i in 1:length(relative_probs)){ 
    temp <- prob_mat_su[i, ]
    ind <- max(which(time == test_set$gs_ffs[i]))
    relative_probs[i] <- temp[ind] 
  }
  
  # find expected status based on 0.5 probability cut-off
  expected_status <- ifelse(relative_probs > p, 0, 1)
  classes <- c(0, 1)
  # auxiliary confusion table
  tabel <- table(factor(test_set$gs_ffsstate, levels = classes),
                 factor(expected_status, levels = classes)) 
  
  confusion_mat <- confusionMatrix(tabel, positive = "1")
  accuracy <- as.numeric(confusion_mat$ffs[1])
  sensitivity <- as.numeric(confusion_mat$byClass[1])
  specificity <- as.numeric(confusion_mat$byClass[2])
  precision <- as.numeric(confusion_mat$byClass[5])
  recall <- as.numeric(confusion_mat$byClass[6])
  f1score <- as.numeric(confusion_mat$byClass[7])
  
  return(list(Accuracy = accuracy, Sensitivity = sensitivity,
              Specificity = specificity, Precision = precision,
              Recall = recall, F1score = f1score))
}

# calculate Brier and Integrated Brier scores at all sorted times
brier_general <- function(pair_df, prob_matrix) {
  
  so <- Surv(pair_df$time, pair_df$status)
  time <- so[, 1]
  ot <- order(time)
  cens <- so[ot, 2]
  time <- time[ot] # ordered times
  N <- nrow(so)
  
  # find the reverse Kaplan-Meier distribution
  hatcdist <- prodlim(Surv(time, cens) ~ 1, reverse = TRUE) 
  csurv <- predict(hatcdist, times = time, type = "surv")
  csurv[csurv == 0] <- Inf
  btime <- time # duplicate time vector
  
  survs <- t(as.matrix(prob_matrix)) # put it as matrix  
  bsc <- rep(0, nrow(survs))
  
  for (j in 1:nrow(survs)) {
    values1 <- as.integer(time <= btime[j] & cens == 1)
    values2 <- as.integer(time > btime[j])
    inb <- survs[j,]
    
    inb[which(inb == 1)] <- inb[which(inb == 1)]
    inb[which(inb == 0)] <- inb[which(inb == 0)]
    bsc[j] <- mean((0 - survs[j, ])^2 * values1 * (1/csurv) + 
                     (1 - survs[j, ])^2 * values2 * (1/csurv[j]))
    
    # calculate the integrated brier score
  }
  
  pos <- sindex(jump.times = btime, eval.times = 1:12)
  brier <- c(0, bsc[pos])
  
  # calculate integrated brier score starting from time point 0
  timing <- 0:(length(brier) - 1)
  indexa <- 2:length(timing)
  int_bs <- (diff(timing) %*% ((brier[indexa - 1] + 
                                  brier[indexa])/2))/diff(range(timing))
  
  return(list(Brier = brier[-1], Int_brier = int_bs))
}

# this should be used when the probability matrix is found
# at sort unique time points
brier_general_su <- function(pair_df, prob_matrix_su) {
  
  so <- Surv(pair_df$time, pair_df$status)
  time <- so[, 1]
  ot <- order(time)
  cens <- so[ot, 2]
  time <- time[ot] # ordered times
  N <- nrow(so)
  
  # find the reverse Kaplan-Meier distribution
  hatcdist <- prodlim(Surv(time, cens) ~ 1, reverse = TRUE) 
  csurv <- predict(hatcdist, times = time, type = "surv")
  csurv[csurv == 0] <- Inf
  btime <- time # duplicate time vector
  
  survs <- t(as.matrix(prob_matrix_su)) # put it as matrix 
  
  bsc <- rep(0, nrow(survs))
  
  for (j in 1:nrow(survs)) {
    values1 <- as.integer(time <= btime[j] & cens == 1)
    values2 <- as.integer(time > btime[j])
    inb <- survs[j,]
    inb[which(inb == 1)] <- inb[which(inb == 1)]
    inb[which(inb == 0)] <- inb[which(inb == 0)]
    bsc[j] <- mean((0 - survs[j, ])^2 * values1 * (1/csurv) + 
                     (1 - survs[j, ])^2 * values2 * (1/csurv[j]))
    
    # calculate the integrated brier score
  }
  
  brier <- c(0, bsc)
  
  # calculate integrated brier score
  timing <- 0:(length(brier) - 1)
  indexa <- 2:length(timing)
  int_bs <- (diff(timing) %*% ((brier[indexa - 1] + 
                                  brier[indexa])/2))/diff(range(timing))
  
  return(list(Brier = brier[-1], Int_brier = int_bs))
}



install_packages <- c("glmnet", "keras", "randomForestSRC",
                      "pec", "survival", "caret")

for (i in 1:length(install_packages)){
  if (!install_packages[i] %in% installed.packages()){
    install.packages(install_packages[i])
  }
}
library(glmnet)
library(keras)
library(randomForestSRC)
library(pec)
library(survival)
library(caret)
source("https://bioconductor.org/biocLite.R")
biocLite("survcomp")
library(survcomp)

use_session_with_seed(seed = 1235, disable_gpu = TRUE,
                      disable_parallel_cpu = TRUE,
                      quiet = FALSE)

source("the_functions.R")
load("unos_failure_free.RData")
unos_failure_free$gs_ffs <- unos_failure_free$gs_ffs 
+ (1/365)

unos_failure_free$gs_ffs <- 
  as.numeric(cut(unos_failure_free$gs_ffs,
                 breaks = 12))

set.seed(12345) 
N <- nrow(unos_failure_free) 
index <- sample(1:N, round(N/3), replace = FALSE)
training_ffs <- unos_failure_free[-index, ] # create training set
test_ffs <- unos_failure_free[index, ] # create test set

# Cox full and Cox KM
cox_full <- coxph(Surv(gs_ffs, gs_ffsstate) ~., 
                  data = training_ffs, x = TRUE, y = TRUE)

cox_null <- coxph(Surv(gs_ffs, gs_ffsstate) ~ 1, 
                  data = training_ffs, x = TRUE, y = TRUE)

pi_cox_full <- predict(cox_full, newdata = test_ffs)

cindex_cox_full <-  concordance.index(x = pi_cox_full,
                                      surv.time =
                                        test_ffs$gs_ffs,
                                      surv.event = 
                                        test_ffs$gs_ffsstate)$c.index

# model cox backward
cox_back <- selectCox(formula = Surv(gs_ffs, gs_ffsstate) ~.,
                      data = training_ffs, rule = "aic")
vars_back <- cox_back$In # variables
f_back <- as.formula(paste("Surv(gs_ffs, gs_ffsstate) ~",
                           paste(vars_back, collapse = "+")))

# Cox backward in coxph object
cox_back <- coxph(formula = f_back, x = TRUE, y = TRUE,
                  data = training_ffs, method = "breslow")

pi_cox_back <- predict(cox_back, newdata = test_ffs)

cindex_cox_back <-  concordance.index(x = pi_cox_back,
                                      surv.time = 
                                        test_ffs$gs_ffs,
                                      surv.event = test_ffs$gs_ffsstate)$c.index 

# model Lasso
response_pair <- c("gs_ffs", "gs_ffsstate")

# for the training set
X_train <- training_ffs[, !(colnames(training_ffs) %in% response_pair)]
X_train_ffs <- model.matrix(~., X_train)[, -1] 
# dim(X_train_ffs) # 41529 x 129
Y_train_ffs <- training_ffs[, colnames(training_ffs) %in% response_pair] # create matrix Y
Y_survtrain_ffs <- Surv(Y_train_ffs$gs_ffs, 
                        Y_train_ffs$gs_ffsstate)

# for the test set
x_test <- test_ffs[, !(colnames(test_ffs) %in% response_pair)]
X_test_ffs <- model.matrix(~., x_test)[, -1]
Y_test_ffs <- test_ffs[, colnames(test_ffs) %in% response_pair]
# Surv function packages survival data into 
# the form expected by glmnet
Y_survtest_ffs <- Surv(Y_test_ffs$gs_ffs, 
                       Y_test_ffs$gs_ffsstate)

# use these two datasets to fit Cox model in the selected features
training_ffs_extended <- as.data.frame(cbind(X_train_ffs, 
                                             Y_train_ffs))
test_ffs_extended <- as.data.frame(cbind(X_test_ffs, Y_test_ffs))

seed <- 12345
set.seed(seed)
folds <- createFolds(y = training_ffs$gs_ffs,
                     k = 10, list = FALSE) 

cv.fit <- cv.glmnet(x = X_train_ffs, y = Y_survtrain_ffs,
                    foldid = folds, parallel = TRUE,
                    family = "cox", grouped = TRUE,
                    maxit = 1000)
cv.fit$lambda.1se

active_index <- as.numeric(unlist(predict(cv.fit,
                                          newx = X_test_ffs,
                                          s = "lambda.1se",
                                          type = "nonzero")))
names_coefs <- colnames(test_ffs_extended)[1:128]
vars_active <- names_coefs[active_index]  

f_lasso <- as.formula(paste("Surv(gs_ffs, gs_ffsstate) ~",
                            paste(vars_active, collapse = "+")))

cox_lasso <- coxph(formula = f_lasso, x = TRUE, y = TRUE,
                   data = training_ffs_extended)

pi_cox_lasso <- predict(cox_lasso, newdata = test_ffs_extended)

cindex_cox_lasso <-  concordance.index(x = pi_cox_lasso,
                                       surv.time = test_ffs$gs_ffs,
                                       surv.event = test_ffs$gs_ffsstate)$c.index

cox_pair <- as.data.frame(cbind(test_ffs_extended$gs_ffs,
                                test_ffs_extended$gs_ffsstate))
colnames(cox_pair) <- c("time", "status")

time_su <- sort(unique(test_ffs_extended$gs_ffs))

# find the probabilities of the models
cox_probs_null <- predictSurvProb(cox_null, test_ffs,
                                  times = time_su)

cox_probs_full <- predictSurvProb(cox_full, test_ffs,
                                  times = time_su)

brier_cox_full <-  brier_general_su(cox_pair, cox_probs_full)
metrics_cox_full <-  metrics_ffs_su(prob_mat_su = cox_probs_full,
                                    test_set = test_ffs,
                                    p = 0.5)

cox_probs_back <- predictSurvProb(cox_back, test_ffs,
                                  times = time_su)

brier_cox_back <-  brier_general_su(cox_pair, cox_probs_back)
metrics_cox_back <-  metrics_ffs_su(prob_mat_su = cox_probs_back,
                                    test_set = test_ffs,
                                    p = 0.5)

cox_probs_lasso <- predictSurvProb(cox_lasso, test_ffs_extended,
                                   times = time_su)

brier_cox_lasso <-  brier_general_su(cox_pair, cox_probs_lasso)
metrics_cox_lasso <-  metrics_ffs_su(prob_mat_su = 
                                       cox_probs_lasso,
                                     test_set = test_ffs_extended,
                                     p = 0.5)

# Random survival forests
model_rsf <- rfsrc(Surv(gs_ffs, gs_ffsstate) ~ .,
                   splitrule = "logrank", nsplit = 5,
                   data = training_ffs, ntree = 300,
                   seed = -12345, mtry = 41, nodesize = 50,
                   ntime = 500, sampsize = 5000,
                   forest = TRUE, importance = FALSE
)

fit_test <- predict(model_rsf,
                    newdata = test_ffs,
                    seed = -12345,
                    split.depth = "all.trees",
                    forest = FALSE)
cindex_rsf <-  1 - fit_test$err.rate[fit_test$ntree] # 0.6942

# always put the times in sorted order
probs_rf_ffs <- predictSurvProb(object = model_rsf,
                                newdata = test_ffs,
                                times = 
                                  sort(unique(test_ffs$gs_ffs)))

pair_rfsrc <- as.data.frame(cbind(test_ffs$gs_ffs,
                                  test_ffs$gs_ffsstate))
colnames(pair_rfsrc) <- c("time", "status")
brier_model_rsf <-  brier_general_su(pair_rfsrc, probs_rf_ffs)
metrics_model_rsf <- metrics_ffs_su(prob_mat_su = probs_rf_ffs,
                                    test_set = test_ffs,
                                    p = 0.5)                       

# Fit the neural networks

training_ffs_scaled <- read.table("training_ffs_scaled.txt")
test_ffs_scaled <- read.table("test_ffs_scaled.txt")
load("training_ffs_long.Rdata")
load("test_ffs_long.Rdata")

data_test_creator <- function(data){
  
  N <- nrow(data)
  # assign survival times to 12 intervals, each is a 1-year period
  data$interval <- max(as.numeric(cut(data$gs_ffs,
                                      breaks = 12))) 
  data$survival <- as.numeric(cut(data$gs_ffs,
                                  breaks = 12)) # the true interval survival 
  data$id <- 70001:(70000 + N) # define the patient ids abstractly
  n.times <- data$interval
  data_long <- data[rep(seq_len(N), times = n.times), ]
  
  # create the correct intervals
  for(i in unique(data_long$id)) {
    n_length <- length(data_long$interval[data_long$id == i])
    data_long$interval[data_long$id == i] <- 1:n_length   
  }
  
  data_long$status <- vector(mode = "numeric", 
                             length = nrow(data_long))
  
  # put indication 1 on status at the intervals
  # on which a patient has died
  for (i in 1:nrow(data_long)) {
    if (data_long$gs_ffsstate[i] == 1 && 
        data_long$survival[i] <= data_long$interval[i]) 
      data_long$status[i] <- 1
  }
  
  return(data_long)
}

measures_calculator <- function(trained_model,
                                datanew, real_status) {
  
  df1 <- data.frame(hazard = predict_proba(trained_model,
                                           as.matrix(datanew[,
                                              c(1:128, 131)]),
                                           batch_size = 500))
  df1$id <- datanew$id # ids of the patients
  df1$survival <- datanew$survival # survival time in years
  groups <- split(df1, f = df1$id) 
  
  true_surv <- unlist(lapply(groups, function(x) {
    surv_obj <- x$survival
    true_res <- surv_obj[1] 
    return(true_res)}
  ))
  group_probs <- lapply(groups, function(x) {
    x <- cumprod(1 - x$hazard)})
  pred_mat <- do.call("rbind", group_probs)
  
  relative_probs <- vector(mode = "numeric",
                           length = length(group_probs))
  for (i in 1:length(group_probs)){ 
    temp <- group_probs[[i]]
    ind <- which(1:12 == true_surv[i])
    relative_probs[i] <- temp[ind] 
  }
  N0 <- length(unique(df1$id)) # number of unique persons
  
  # in the data frame create random id numbers 
  # to label the patients
  df2 <- data.frame(relative_probs, true_surv,
                    status = real_status,
                    id = (70001):(70000 + N0))
  df2$prediction <- 1 - round(df2$relative_probs, digits = 0)
  # create possible classes
  classes <- c(0, 1)
  # auxiliary confusion table
  tabel <- table(factor(df2$prediction, levels = classes),
                 factor(df2$status, levels = classes)) 
  confusion_mat <- confusionMatrix(tabel, positive = "1")
  accuracy <- sum(df2$prediction == df2$status) / nrow(df2)
  sensitivity <- sum(df2$status == 1 & df2$prediction == 1) /
    colSums(tabel)[2]
  specificity <- sum(df2$status == 0 & df2$prediction == 0) / 
    colSums(tabel)[1]
  precision <- as.numeric(confusion_mat$byClass[5])
  recall <- as.numeric(confusion_mat$byClass[6])
  f1score <- as.numeric(confusion_mat$byClass[7])
  
  brier_obj <- brier_nnet(df2, prob_matrix = pred_mat)
  int_brier <- as.numeric(brier_obj$Int_brier)
  brier_set <- brier_obj$Brier
  all_weights <- get_weights(fit_keras)
  nr_weights <- (length(all_weights[[1]]) +
                   length(all_weights[[2]]) 
                 + length(all_weights[[3]]) 
                 + length(all_weights[[4]])) 
  
  return(list(weights = nr_weights, 
              node_size = length(get_weights(fit_keras)[[2]]),
              cross_entropy =
                result$metrics$loss[result$params$epochs],
              accuracy = accuracy,
              sensitivity = as.numeric(sensitivity), 
              specificity = as.numeric(specificity),
              Precision = precision, Recall = recall,
              F1score = f1score,
              Integrated_brier = int_brier,
              Brier_scores = brier_set))
}

# model with one hidden layer
set.seed(12345)

fit_keras <- keras_model_sequential()
# Add layers to the model

# we create a densely connected ANN to the output
fit_keras %>%
  layer_dense(units = 75, activation = 'sigmoid',
              input_shape = c(129)) %>% 
  layer_dropout(rate = 0.1) %>%
  layer_dense(units = 1, activation = 'sigmoid')
fit_keras %>% compile(
  loss = 'binary_crossentropy', # for binary class classification 
  optimizer = optimizer_sgd(lr = 0.2, momentum = 0.9)
)
event_status <- test_ffs_scaled$gs_ffsstate 

# create the matrices to be used for keras library
# predictors: 128 variables + interval
train_x <- as.matrix(training_ffs_long[, c(1:128, 131)]) 
dimnames(train_x) <- NULL # the object must have empty dimnames
train_y <- training_ffs_long$status
test_x <- as.matrix(test_ffs_long[, c(1:128, 131)])
dimnames(test_x) <- NULL # the object must have empty dimnames
test_y <- test_ffs_long$status

result <- fit_keras %>% fit(
  train_x, 
  train_y, 
  epochs = 10, 
  batch_size = 1000,
  validation_data = list(test_x, test_y),
  class_weight = list("0" = 1, "1" = 1)
  #,callbacks = c(early_stopping)
)

# now that the model has run lets calculate the measures
metrics_model_nn1h <- measures_calculator(trained_model =
                                            fit_keras,
                                          datanew = test_ffs_long,
                                          real_status = 
                                            event_status)

df1 <- data.frame(hazard = predict_proba(fit_keras,
                                         as.matrix(test_ffs_long[,
                                                  c(1:128, 131)]),
                                         batch_size = 500))
df1$id <- test_ffs_long$id # ids of the patients
groups <- split(df1, f = df1$id) 

group_probs <- lapply(groups, function(x) {
  x <- cumprod(1 - x$hazard)})
pred_mat_nn1h <- do.call("rbind", group_probs)

set.seed(12345)
# model with 2 hidden layers
# find the metrics for combination node_size = 100,
#dropout rate = 0.1, 
# learning rate 0.2, momentum 0.9, weak class weight 1
fit_keras2 <- keras_model_sequential()
# we create a densely connected ANN to the output
fit_keras2 %>%
  layer_dense(units = 100, activation = 'tanh', 
              input_shape = c(129)) %>% 
  layer_dropout(rate = 0.1) %>%
  layer_dense(units = 100, activation = 'tanh') %>%
  layer_dropout(rate = 0.1) %>%
  layer_dense(units = 1, activation = 'sigmoid')
fit_keras2 %>% compile(
  loss = 'binary_crossentropy', # for binary class classification 
  optimizer = optimizer_sgd(lr = 0.2, momentum = 0.9)
)

event_status <- test_ffs_scaled$gs_ffsstate 

# create the matrices to be used for keras library
# predictors: 128 variables + interval
train_x <- as.matrix(training_ffs_long[, c(1:128, 131)]) 
dimnames(train_x) <- NULL # the object must have empty dimnames
train_y <- training_ffs_long$status
test_x <- as.matrix(test_ffs_long[, c(1:128, 131)])
dimnames(test_x) <- NULL # the object must have empty dimnames
test_y <- test_ffs_long$status

result <- fit_keras2 %>% fit(
  train_x, 
  train_y, 
  epochs = 10, 
  batch_size = 1000,
  validation_data = list(test_x, test_y),
  class_weight = list("0" = 1, "1" = 1)
  #,callbacks = c(early_stopping)
)

metrics_model_nn2h <- measures_calculator(trained_model =
                                            fit_keras2,
                                          datanew = test_ffs_long,
                                          real_status = 
                                            event_status)

df2 <- data.frame(hazard = predict_proba(fit_keras2,
                                         as.matrix(test_ffs_long[,
                                                                 c(1:128, 131)]),
                                         batch_size = 1000))

df2$id <- test_ffs_long$id # ids of the patients
df2$survival <- test_ffs_long$survival # survival time in years
groups2 <- split(df2, f = df2$id) 
group_probs2 <- lapply(groups2, function(x) {
  x <- cumprod(1 - x$hazard)})
pred_mat_nn2h <- do.call("rbind", group_probs2)

#########################################################
# Global performance measures at distinct time points
#########################################################

obj_cox_full <- c(metrics_cox_full$Accuracy,
                  metrics_cox_full$Sensitivity,
                  metrics_cox_full$Specificity,
                  metrics_cox_full$Precision,
                  metrics_cox_full$F1score,
                  brier_cox_full$Int_brier,
                  cindex_cox_full)

obj_cox_back <- c(metrics_cox_back$Accuracy,
                  metrics_cox_back$Sensitivity,
                  metrics_cox_back$Specificity,
                  metrics_cox_back$Precision,
                  metrics_cox_back$F1score,
                  brier_cox_back$Int_brier,
                  cindex_cox_back)

obj_cox_lasso <- c(metrics_cox_lasso$Accuracy,
                   metrics_cox_lasso$Sensitivity,
                   metrics_cox_lasso$Specificity,
                   metrics_cox_lasso$Precision,
                   metrics_cox_lasso$F1score,
                   brier_cox_lasso$Int_brier,
                   cindex_cox_lasso)

obj_model_rsf <- c(metrics_model_rsf$Accuracy,
                   metrics_model_rsf$Sensitivity,
                   metrics_model_rsf$Specificity,
                   metrics_model_rsf$Precision,
                   metrics_model_rsf$F1score,
                   brier_model_rsf$Int_brier,
                   cindex_rsf)

obj_model_nn1h <- c(metrics_model_nn1h$accuracy,
                    metrics_model_nn1h$sensitivity,
                    metrics_model_nn1h$specificity,
                    metrics_model_nn1h$Precision,
                    metrics_model_nn1h$F1score,
                    metrics_model_nn1h$Integrated_brier)

obj_model_nn2h <- c(metrics_model_nn2h$accuracy,
                    metrics_model_nn2h$sensitivity,
                    metrics_model_nn2h$specificity,
                    metrics_model_nn2h$Precision,
                    metrics_model_nn2h$F1score,
                    metrics_model_nn2h$Integrated_brier)

df_final <- round(data.frame(rbind(obj_cox_full, obj_cox_back,
                                   obj_cox_lasso, obj_model_rsf,
                                   obj_model_nn1h, obj_model_nn2h)),
                  4)

colnames(df_final) <- c("Accuracy", "Sensitivity", "Specificity",
                        "Precision", "F1 score", "IBS", "C-index")
rownames(df_final) <- c("Cox all variables", "Cox backward",
                        "Cox LASSO", "RSF",
                        "Neural Network 1h", "Neural Network 2h")
df_final$`C-index`[5:6] <- c("-", "-")

#########################################
# Time-dependent Brier score plot
#########################################

line1 <- data.frame(times = time_su,
                    brier = brier_cox_full$Brier,
                    Model = "Cox all variables")
line2 <- data.frame(times = time_su,
                    brier = brier_cox_back$Brier,
                    Model = "Cox backward")
line3 <- data.frame(times = time_su,
                    brier = brier_cox_lasso$Brier,
                    Model = "Cox LASSO")
line4 <- data.frame(times = time_su,
                    brier = brier_model_rsf$Brier,
                    Model = "RSF")
line5 <- data.frame(times = time_su,
                    brier = metrics_model_nn1h$Brier_scores,
                    Model = "Neural Network 1h")
line6 <- data.frame(times = time_su,
                    brier = metrics_model_nn2h$Brier_scores,
                    Model = "Neural Network 2h")

df <- rbind(line1, line2, line3, line4, line5, line6)

a1 <- ggplot(df, aes(x = times, y = brier, color = Model)) + 
  geom_line(size = 0.9) + 
  scale_x_continuous(breaks = 0:12) + 
  xlab("Time since transplantation in years") + 
  ylab("Prediction error (Brier score)") + 
  ylim(c(0, 0.3)) + 
  theme_classic() + 
  ggtitle(" ") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values=c("black", "red", "purple",
                              "blue", "green", "darkgreen"))

############################################
# mean survival probabilities plot
############################################

mean_probs_cox_null <- colMeans(cox_probs_null)
mean_probs_cox_full <- colMeans(cox_probs_full)
mean_probs_cox_back <- colMeans(cox_probs_back)
mean_probs_cox_lasso <- colMeans(cox_probs_lasso)
mean_probs_rsf <- colMeans(probs_rf_ffs)
mean_probs_nn1h <- colMeans(pred_mat_nn1h)
mean_probs_nn2h <- colMeans(pred_mat_nn2h)

line0 <- data.frame(times = time_su, prob = mean_probs_cox_null,
                    Model = "Reference")
line7 <- data.frame(times = time_su, prob = mean_probs_cox_full,
                    Model = "Cox all variables")
line8 <- data.frame(times = time_su, prob = mean_probs_cox_back,
                    Model = "Cox backward")
line9 <- data.frame(times = time_su, prob = mean_probs_cox_lasso,
                    Model = "Cox LASSO")
line10 <- data.frame(times = time_su, prob = mean_probs_rsf,
                     Model = "RSF")
line11 <- data.frame(times = time_su, prob = mean_probs_nn1h,
                     Model = "Neural Network 1h")
line12 <- data.frame(times = time_su, prob = mean_probs_nn2h,
                     Model = "Neural Network 2h")

df2 <- rbind(line7, line8, line9, line10, line11, line12)

a2 <- ggplot(df2, aes(x = times, y = prob, color = Model)) + 
  geom_line(size = 0.9) + 
  scale_x_continuous(breaks = 0:12) + 
  xlab("Time since transplantation in years") + 
  ylab("Survival probability") + 
  ylim(c(0, 1)) + 
  theme_classic() + 
  ggtitle(" ") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values=c("black", "red", "purple",
                              "blue", "green", "darkgreen"))
grid.arrange(a1, a2)

#################################################
# predictions on new synthetic patients
#################################################

tumor <- factor(levels(test_ffs$rec_tumor))
age <- quantile(training_ffs$recipientage, c(0.10, 0.5, 0.90))

# create new hypothetic patient
combis <- with(training_ffs, expand.grid(
  rec_tumor = tumor, # increasing age to different
  recipientage = age
))

df_synth_normal$recipientage <- combis$recipientage
df_synth_normal$rec_tumor <- combis$rec_tumor

df_synth_expanded$recipientage <- combis$recipientage
df_synth_expanded$rec_tumorY <- as.numeric(combis$rec_tumor) - 1

synth_cox_full <- predictSurvProb(object = cox_full,
                                  newdata = df_synth_normal,
                                  times = time_su)
synth_cox_back <- predictSurvProb(object = cox_back,
                                  newdata = df_synth_normal,
                                  times = time_su)
synth_cox_lasso <- predictSurvProb(object = cox_lasso,
                                   newdata = df_synth_expanded,
                                   times = time_su)

time_val <- time_su

# predictions with rsf
model_rsf <- rfsrc(Surv(gs_ffs, gs_ffsstate) ~ .,
                   splitrule = "logrank", nsplit = 7,
                   data = training_ffs, ntree = 300,
                   seed = -12345, mtry = 23, nodesize = 50,
                   ntime = 500, sampsize = 5000,
                   forest = TRUE, importance = FALSE
)

synth_rsf_ffs <- predictSurvProb(object = model_rsf,
                                 newdata = df_synth_normal,
                                 times = time_su)

# create hypothetical patient for synthetic scaled

tumor_scaled <- c(0, 1)
age_scaled <- quantile(training_ffs_scaled$recipientage,
                       c(0.10, 0.5, 0.90))

# create new hypothetic patient
combis_scaled <- with(training_ffs_scaled, expand.grid(
  rec_tumorY = tumor_scaled, # increasing age to different
  recipientage = age_scaled
))

df_synth_scaled$rec_tumorY <- combis_scaled$rec_tumorY
df_synth_scaled$recipientage <- combis_scaled$recipientage
# create random patient survival and gs_ffsstate status
df_synth_scaled$gs_ffs <- c(2, 4, 5, 1, 9, 4)
df_synth_scaled$gs_ffsstate <- c(1, 0, 1, 0, 0 ,0)
rownames(df_synth_scaled) <- NULL
data_synth_long <- data_test_creator(df_synth_scaled)

df_new <- data.frame(hazard = predict_proba(fit_keras, 
                                            as.matrix(data_synth_long[, 
                                                                      c(1:128, 131)]),
                                            batch_size = 1000))
df_new$id <- rep(1:6, each = 12) # ids of the patients
groups_new <- split(df_new, f = df_new$id) 
group_probs_new <- lapply(groups_new, function(x) {
  x <- cumprod(1 - x$hazard)})
pred_mat_new <- do.call("rbind", group_probs_new)

# for neural network with 2 hidden layers
df_new2 <- data.frame(hazard = predict_proba(fit_keras2,    
                                             as.matrix(data_synth_long[,
                                                                       c(1:128, 131)]),
                                             batch_size = 500))
df_new2$id <- rep(1:6, each = 12) # ids of the patients
groups_new2 <- split(df_new2, f = df_new2$id) 
group_probs_new2 <- lapply(groups_new2, function(x) {
  x <- cumprod(1 - x$hazard)})
pred_mat_new2 <- do.call("rbind", group_probs_new2)



################################################################
# Survival curves for rec_immuno_maint_meds and rec_postx_los
################################################################

immuno_meds <- 
  factor(levels(training_ffs$rec_immuno_maint_meds))
los <- quantile(training_ffs$rec_postx_los,
                c(0.10, 0.5, 0.90))

# create new hypothetic patient
combis <- with(training_ffs, expand.grid(
  rec_immuno_maint_meds = immuno_meds, 
  rec_postx_los = los
))

df_synth2_normal$rec_postx_los <- combis$rec_postx_los
df_synth2_normal$rec_immuno_maint_meds <-
  combis$rec_immuno_maint_meds

df_synth2_expanded$rec_postx_los <- combis$rec_postx_los
df_synth2_expanded$rec_immuno_maint_medsY <- 
  as.numeric(combis$rec_immuno_maint_meds) - 1

synth_cox_full2 <- predictSurvProb(object = cox_full,
                                   newdata = df_synth2_normal,
                                   times = time_su)
synth_cox_back2 <- predictSurvProb(object = cox_back,
                                   newdata = df_synth2_normal,
                                   times = time_su)
synth_cox_lasso2 <- predictSurvProb(object = cox_lasso,
                                    newdata = df_synth2_expanded,
                                    times = time_su)
time_val <- time_su

synth_rsf_ffs2 <- predictSurvProb(object = model_rsf,
                                  newdata = df_synth2_normal,
                                  times = time_su)

# create hypothetical patient for synthetic scaled

immuno_meds_scaled <- c(0, 1)
los_scaled <- quantile(training_ffs_scaled$rec_postx_los,
                       c(0.10, 0.5, 0.90))

# create new hypothetic patient
combis_scaled2 <- with(training_ffs_scaled, expand.grid(
  rec_immuno_maint_medsY = immuno_meds_scaled, 
  rec_postx_los = los_scaled
))

df_synth2_scaled$rec_immuno_maint_medsY <-
  combis_scaled2$rec_immuno_maint_medsY
df_synth2_scaled$rec_postx_los <- combis_scaled2$rec_postx_los
# create random patient survival and gs_ffsstate status
df_synth2_scaled$gs_ffs <- c(2, 4, 5, 1, 9, 4)
df_synth2_scaled$gs_ffsstate <- c(1, 0, 1, 0, 1 ,0)
rownames(df_synth2_scaled) <- NULL
data_synth_long <- data_test_creator(df_synth2_scaled)

df_new <- data.frame(hazard = predict_proba(fit_keras,
                                            as.matrix(data_synth_long[,
                                                                      c(1:128, 131)]),
                                            batch_size = 500))
df_new$id <- rep(1:6, each = 12) # ids of the patients
groups_new <- split(df_new, f = df_new$id) 
group_probs_new <- lapply(groups_new, function(x) {
  x <- cumprod(1 - x$hazard)})
pred_mat_new <- do.call("rbind", group_probs_new)

df_new2 <- data.frame(hazard = predict_proba(fit_keras2,  
                                             as.matrix(data_synth_long[,
                                                                       c(1:128, 131)]),
                                             batch_size = 1000))
df_new2$id <- rep(1:6, each = 12) # ids of the patients
groups_new2 <- split(df_new2, f = df_new2$id) 
group_probs_new2 <- lapply(groups_new2, function(x) {
  x <- cumprod(1 - x$hazard)})
pred_mat_new2 <- do.call("rbind", group_probs_new2)

cases <- c("Los:  5 , Immuno meds:  N",
           "Los:  5 , Immuno meds:  Y",
           "Los:  10 , Immuno meds:  N",
           "Los:  10 , Immuno meds:  Y",
           "Los:  33 , Immuno meds:  N",
           "Los:  33 , Immuno meds:  Y")

synth_cox_full2 <-  cbind(1, synth_cox_full2)
synth_cox_back2 <-  cbind(1, synth_cox_back2)
synth_cox_lasso2 <-  cbind(1, synth_cox_lasso2)
synth_rsf_ffs2 <-  cbind(1, synth_rsf_ffs2)
pred_mat_new <- cbind(1, pred_mat_new)
pred_mat_new2 <- cbind(1, pred_mat_new2)


pat1 <- data.frame(Times = time_val,
                   Prob1 = synth_cox_full2[1, ], 
                   Prob2 = synth_cox_back2[1, ],
                   Prob3 = synth_cox_lasso2[1, ],
                   Prob4 = synth_rsf_ffs2[1, ],
                   Prob5 = pred_mat_new[1, ],
                   Prob6 = pred_mat_new2[1, ],
                   Patient = cases[1])
pat2 <- data.frame(Times = time_val, 
                   Prob1 = synth_cox_full2[2, ], 
                   Prob2 = synth_cox_back2[2, ],
                   Prob3 = synth_cox_lasso2[2, ],
                   Prob4 = synth_rsf_ffs2[2, ], 
                   Prob5 = pred_mat_new[2, ],
                   Prob6 = pred_mat_new2[2, ],
                   Patient = cases[2])
pat3 <- data.frame(Times = time_val,
                   Prob1 = synth_cox_full2[3, ], 
                   Prob2 = synth_cox_back2[3, ],
                   Prob3 = synth_cox_lasso2[3, ],
                   Prob4 = synth_rsf_ffs2[3, ],
                   Prob5 = pred_mat_new[3, ],
                   Prob6 = pred_mat_new2[3, ],
                   Patient = cases[3])
pat4 <- data.frame(Times = time_val,
                   Prob1 = synth_cox_full2[4, ], 
                   Prob2 = synth_cox_back2[4, ],
                   Prob3 = synth_cox_lasso2[4, ],
                   Prob4 = synth_rsf_ffs2[4, ],
                   Prob5 = pred_mat_new[4, ],
                   Prob6 = pred_mat_new2[4, ],
                   Patient = cases[4])
pat5 <- data.frame(Times = time_val,
                   Prob1 = synth_cox_full2[5, ], 
                   Prob2 = synth_cox_back2[5, ],
                   Prob3 = synth_cox_lasso2[5, ],
                   Prob4 = synth_rsf_ffs2[5, ],
                   Prob5 = pred_mat_new[5, ],
                   Prob6 = pred_mat_new2[5, ],
                   Patient = cases[5])
pat6 <- data.frame(Times = time_val,
                   Prob1 = synth_cox_full2[6, ], 
                   Prob2 = synth_cox_back2[6, ],
                   Prob3 = synth_cox_lasso2[6, ],
                   Prob4 = synth_rsf_ffs2[6, ],
                   Prob5 = pred_mat_new[6, ],
                   Prob6 = pred_mat_new2[6, ], 
                   Patient = cases[6])

df_all_pat <- rbind(pat1, pat2, pat3, pat4, pat5, pat6)

ggplot(df_all_pat, aes(x = Times)) + 
  geom_line(size = 0.8, aes(y  = Prob1, col = "Cox_full")) + 
  geom_line(size = 0.8, aes(y  = Prob2, col = "Cox_backward")) + 
  geom_line(size = 0.8, aes(y  = Prob3, col = "Cox_LASSO")) +
  geom_line(size = 0.8, aes(y  = Prob4, col = "RSF")) +
  geom_line(size = 0.8, aes(y  = Prob5, col = "NN_1h")) +
  geom_line(size = 0.8, aes(y  = Prob6, col = "NN_2h")) + 
  scale_x_continuous(breaks = 0:12) + 
  xlab("Time since transplantation in years") + 
  ylab("Survival probability") + 
  ylim(c(0, 1)) + 
  theme_classic() + 
  ggtitle(" ") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~Patient) +
  scale_colour_manual(name = " ",
                      values=c(Cox_full = "black",
                               Cox_backward = "red",
                               Cox_LASSO = "purple",
                               RSF = "blue",
                               NN_1h = "green",
                               NN_2h = "darkgreen")) + 
  theme(axis.title=element_text(face="bold.italic", 
                                size="12", color="brown"),
        legend.position="right")




