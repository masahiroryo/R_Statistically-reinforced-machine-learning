###########################################################################
#
#  Masahiro Ryo and Matthias C. Rillig. 2017. 
#  Statistically-reinforced machine learning for 
#  nonlinear patterns and variable interactions. Ecosphere.
#    Source: https://github.com/masahiroryo/R_Statistically-rainforced-machine-learning
#
###########################################################################
###  The original script for "Functions.R" is available in Hapfelmeier and Ulm (2013) 
###  http://dx.doi.org/10.1016/j.csda.2012.09.020
###  It was modified for applying their script with parallel computing in "Windows OS". 
###  [IMPORTANT] THOSE WHO USE THESE FUNCTIONS MUST CITE THEIR AWESOME WORK.
###########################################################################
#  
#  To analyze the dataset using Statistically-oriented Machine Learnings (Lines 242-320)
#   0. reading data, packages and functions
#   1. conditional inference tree
#   2. model-based tree
#   3. permutation-based random forest
#   4. partial dependence plot
#   5. LM with AIC best model selection
#
#   INPUT: "data_example.csv"
#   OUTPUT: analyses using the models listed above
#
###########################################################################
seed = 1   # number for set.seed for random number generation

###-------------------------------------------------------------------------
### 0. reading data, packages and functions
###-------------------------------------------------------------------------
# reading the example data from csv
  data = read.csv("data_example.csv", header=T) # see "Data_generator.R" 

# reading packages (to be automatically installed otherwise)
  package.list <- c("party", "caret","mlr", "randomForest","foreach","doSNOW")
  tmp.install <- which(lapply(package.list, require, character.only = TRUE)==FALSE)
  if(length(tmp.install)>0) install.packages(package.list[tmp.install])
  lapply(package.list, require, character.only = TRUE)

# reading functions: myvarimp and RF_Hapfelmeier
  source("Functions.R")


###-------------------------------------------------------------------------
### 1. Conditional inference tree  # package party
###-------------------------------------------------------------------------
# Fig.3
  set.seed(seed)
  model.ct = ctree(Y~.,data=data, controls=ctree_control(testtype = "Bonferroni")) # see ctree_control for setting
  plot(model.ct)

###-------------------------------------------------------------------------
### 2. Model-based tree  # package partykit 
###-------------------------------------------------------------------------
# Fig.4
  formula = as.formula(paste("Y~x9", paste(colnames(data)[-which(colnames(data) %in% c("Y","x9"))], collapse=" + "), sep=" | ")) # creating a formula
  set.seed(seed)
  lm_tree = mob(formula, data = data, na.action = na.omit, model=linearModel,
                control = mob_control(alpha = 0.05, bonferroni = T))
  plot(lm_tree)

# Fig.7
  formula.2 = as.formula(paste("Y~x9+x7", paste(colnames(data)[-which(colnames(data) %in% c("Y","x7","x9"))], collapse=" + "), sep=" | ")) # creating a formula
  set.seed(seed)
  lm_tree.2 = mob(formula.2, data = data, na.action = na.omit, model=linearModel,
                control = mob_control(alpha = 0.05, bonferroni = T))
  plot(lm_tree.2)

###-------------------------------------------------------------------------
### 3. Permutation-based random forest
###        note that the function RF_Hapfelmeier may take extraordinally long time
###        depending on nperm and ntree. We firstly recommend to try with 
###        nperm=100, ntree=30 to run quickly.
###        In our study, we used nperm=2000, ntree=100  
###
###        ncore determines the number of cores used for parallel computing
###-------------------------------------------------------------------------
# Fig.5
  formula = as.formula(paste("Y", paste(colnames(data)[2:length(colnames(data))], collapse=" + "), sep=" ~ ")) # creating a formula
  set.seed(seed)
  model.Hapfelmeier.RF = RF_permutation(formula, data=data, nperm=100, ntree=30, ncore=4)
  set.seed(seed)
  variable.importance = cbind(model.Hapfelmeier.RF$varimp, model.Hapfelmeier.RF$p.values.bonf) # variable importance and p-value for each predictor
  plot(variable.importance[,1], xlab="predictor_x", ylab="variable importance", pch = ifelse(variable.importance[,2]<0.05, 16, 1))

  ## Note that, in Hapfelmeier and Ulm (2013), they used p-value estimate for variable selection and
  ## they build the best model using only statistically-significant predictors 
  ## the best model can be confirmed with "model.Hapfelmeier.RF$forest.bonf"
  
###-------------------------------------------------------------------------
### 4. Partial dependence plot  # package mlr, party, randomForest
###-------------------------------------------------------------------------
  task = makeRegrTask(data = data, target = "Y")
  learner = makeLearner(cl="regr.cforest", predict.type = "response", par.vals = list(ntree=300))
  fit = train(learner, task)

# Fig.6(a,b,c): x7, x8, x9
  pd = generatePartialDependenceData(fit, task, c("x7","x8","x9"), interaction=F)
  plotPartialDependence(pd,geom = "line") # overall averages

# Fig.6(d,e,f):x1*x4*x9
  pd.int = generatePartialDependenceData(fit, task, c("x1","x4","x9"), interaction=T)
  
  plot(0, 0, type = 'n', xlim = range(pd.int$data[,4]), ylim = range(pd.int$data[,1]),
       xlab = "x9", ylab = 'partial dependence')
  col.id = 0
  for(i in unique(pd.int$data[,2])){
    for(j in unique(pd.int$data[,3])){
      col.id = col.id + 1
      lines(subset(pd.int$data[,4],pd.int$data[,2]==i&pd.int$data[,3]==j),
            subset(pd.int$data[,1],pd.int$data[,2]==i&pd.int$data[,3]==j),
            col =col.id,lwd=1.5)
    }
  }
  legend(15,4.5, legend = c("{0,low}","{0,high}","{0,moderate}","{1,low}","{1,high}","{1,moderate}"), lty = 1, cex=0.7, col = c(1:6), ncol = 1)



###-------------------------------------------------------------------------
### 5. LM with AIC best model selection
###-------------------------------------------------------------------------
  formula = as.formula("Y~poly(x7, 2)+x8+x1*x4*x9+x2+x3+x5+x6+x10+x11+x12+x13+x14+x15+x16")
  model.lm <-step(lm(formula,data=na.omit(data)))
  summary(model.lm)
