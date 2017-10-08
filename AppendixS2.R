###########################################################################
#
#  Masahiro Ryo and Matthias C. Rillig. 2017. 
#  Statistically-reinforced machine learning for 
#  nonlinear patterns and variable interactions. Ecosphere.
#
#  Appendix S2 R script used in this study.
#  please feel free to contact for any questions. [email: masahiroryo@gmail.com]
###########################################################################
#  1. Setting required analytical conditions (Lines 27-127)
#   1.1. To read/install packages
#   1.2. To set up functions for permutation-based random forest
#
#  2.To generate an artifical dataset (Lines 128-241)
#   2.1. To generate binary, categorical, and numeric variables (predictors)
#   2.2. To generate a response variable based on the predictors
#
#  3. To analyze the dataset using Statistically-oriented Machine Learnings (Lines 242-320)
#   3.1. conditional inference tree
#   3.2. model-based tree
#   3.3. permutation-based random forest
#   3.4. partial dependence plot
#   3.5. LM with AIC best model selection
###########################################################################


###################################################################
## 1. Setting required analytical conditions                     ##
###################################################################
###-------------------------------------------------------------------------
### 1.1. To read/install packages
###-------------------------------------------------------------------------
# if any of these have not been installed yet, they are done automatically
  package.list <- c("party", "caret","mlr", "randomForest","foreach","doSNOW")
  tmp.install <- which(lapply(package.list, require, character.only = TRUE)==FALSE)
  if(length(tmp.install)>0) install.packages(package.list[tmp.install])
  lapply(package.list, require, character.only = TRUE)

###-------------------------------------------------------------------------
### 1.2. To set up functions of permutation-based random forest
###-------------------------------------------------------------------------
###  The original script is available in Hapfelmeier and Ulm (2013) 
###  http://dx.doi.org/10.1016/j.csda.2012.09.020
###  and modified for applying their script with parallel computing in "Windows OS". 
###  [IMPORTANT] THOSE WHO USE THIS APPENDIX INCLUDING THIS SECTION 1.2 MUST
###              CITE THEIR AWESOME WORK IN ANY PUBLICATIONS.
###-------------------------------------------------------------------------
### function for calculating variable importance measure (modified by Masahiro Ryo)
###-------------------------------------------------------------------------
    myvarimp <- function(object, mincriterion = 0, conditional = FALSE, threshold = 0.2,
                         pre1.0_0 = conditional, varID) {
      response <- object@responses
      input <- object@data@get("input")
      xnames <- colnames(input)
      inp <- initVariableFrame(input, trafo = NULL)
      y <- object@responses@variables[[1]]
      if (length(response@variables) != 1) stop("cannot compute variable importance measure")
      CLASS <- all(response@is_nominal)
      ORDERED <- all(response@is_ordinal)
      if (CLASS) {
        error <- function(x, oob) mean((levels(y)[sapply(x, which.max)] != y)[oob])
      }else {
        if (ORDERED) {
          error <- function(x, oob) mean((sapply(x, which.max) != y)[oob])
        } else {
          error <- function(x, oob) mean((unlist(x) - y)[oob]^2)}
      }
      perror <- rep(0, length(object@ensemble))
      for (b in 1:length(object@ensemble)) {
        tree <- object@ensemble[[b]]
        w <- object@weights[[b]]
        w[w == 0] <- 1
        oob <- object@weights[[b]] == 0
        p <- .Call("R_predict", tree, inp, mincriterion, -1L, PACKAGE = "party")
        eoob <- error(p, oob)
        j <- varID
        p <- .Call("R_predict", tree, inp, mincriterion, as.integer(j), PACKAGE = "party")
        perror[(b - 1)] <- (error(p,oob) - eoob)
      }
      return(MeanDecreaseAccuracy = mean(perror))
    }
    environment(myvarimp) <- environment(varimp)

###-------------------------------------------------------------------------
### Permutation-based random forest function (modified by Masahiro Ryo)
###-------------------------------------------------------------------------
RF_Hapfelmeier <- function(formula, data, nperm = 400, ntree = 100, ncore=4)
  # formula: object of class "formula".
  # data: data frame containing the variables.
  # nperm: Number of permutation steps used for the permutation test.
  # ntree: Number of trees in the Random Forest.
  # ncore: Number of cores used for parallel computing
{
  library("foreach")  # to parallelize the main processes
  library("doSNOW")   # to activate foreach in windows OS
      
  x.names <- all.vars(formula)[-1]
  y.names <- all.vars(formula)[1]
  terms. <- terms.formula(formula)
  x.formula <- attr(terms., "term.labels")
  y.formula <- as.character(attr(terms., "variables"))[2]
  mtry <- ceiling(sqrt(length(x.formula)))
  dat <- subset(data, select = c(y.names, x.names))
  forest <- party::cforest(formula, data = dat, controls = cforest_unbiased(mtry = mtry, ntree = ntree))
  obs.varimp <- varimp(forest)
  perm.mat <- matrix(NA, ncol = length(x.names), nrow = nperm, dimnames = list(1:nperm, x.names))
      
  cl<-makeCluster(ncore) #change to your number of CPU cores
  registerDoSNOW(cl)

  for (j in x.names) {
    cat("\r", "Processing variable ", which(j == x.names), " of ", length(x.names)); flush.console()
    perm.dat <- dat
    perm.mat[, j] <- unlist(foreach(i = 1:nperm, .packages='party',.export="myvarimp") %dopar% {
      perm.dat[, j] <- sample(perm.dat[, j]);
      myvarimp(party::cforest(formula, data = perm.dat, controls = cforest_unbiased(mtry = mtry, ntree = ntree)), varID = which(x.names == j))
    })
  }
  stopCluster(cl)
      
  p.vals <- sapply(x.names, function(x) sum(perm.mat[, x] >= obs.varimp[which(x == x.names)]) / nperm)
  p.vals.bonf <- p.vals * length(p.vals)

  cat("\n", "\r"); flush.console()
  return(list("forest" = forest, "p.values" = p.vals,"p.values.bonf" = p.vals.bonf))
}

#################################################################
## 2. To generate an artifical dataset                         ##
#################################################################

###-------------------------------------------------------------------------
### 2.1. To generate binary, categorical, and numeric variables (predictors)
###-------------------------------------------------------------------------

## 2.1.1. Settings of data sizes 
  n = 300               # sample size (we tried 100, 200, 300)
  n.bi = 3              # number of binary predictors
  n.ca = 3              # number of categorical predictors
  n.nu = 10             # number of numeric predictors
  NA_ratio = 0.0       # proportion of NA for each predictors (we did not include NA)
  seed = 1              # number for set.seed for random number generation



## 2.1.2. binary data generation 
  for(i in c(1:n.bi)){
    set.seed(seed+2*i)  # changing random seed
    p = rnorm(1,mean=0.5,sd=0.1) # a probability for binomial dist.
    set.seed(seed+2*i)  
    eval(parse(text=paste("x",i, " = rbinom(",n,",","1,",p,")",sep=""))) # e.g. x1=rbinom(300,1,0.4)
  }


## 2.1.3. categorical data generation
  cat.type = list(c("high","moderate","low"),c("A","B","C","D"), 
                  c("c1","c2","c3","c4","c5"))  # representing 3, 4 & 5 classes
  set.seed(seed)
  ca.n = sample(c(1:length(cat.type)),n.ca,replace=T) # selecting ones from the category types
  for(i in c(1:n.ca)){
    c.num = length(cat.type[[ca.n[i]]])
    set.seed(seed+i)  
    weight = rnorm(c.num,mean=1/c.num,sd=0.1/c.num) # random weight generation for selecting each category
    weight = weight/sum(weight)                     # summed weight=1 (e.g., proportion = c(0.2, 0.5, 0.3))
    eval(parse(text=paste("x",i+n.bi,"=sample(cat.type[[ca.n[i]]],n,prob=weight, replace=T)",sep=""))) # taking sample randomly with the determined proportions
    eval(parse(text=paste("x",i+n.bi,"=as.factor(x",i+n.bi, ")",sep="")))
  }


## 2.1.4. numeric data generation
  prob.dist = c("runif","rnorm","rpois")                  # probability distributional types
  prob.para = c("min=0,max=10","mean=10,sd=3","lambda=5") # parameters for the types
  set.seed(seed)
  nu.type = sample(c(1:length(prob.dist)),n.nu,replace=T) # selecting ones from the prob.dist types
  
  for(i in c(1:n.nu)){
    eval(parse(text=paste("x",i+n.ca+n.bi, "=", prob.dist[nu.type[i]],    #e.g. x8=rpois(300,lambda=5)
                          "(n,",prob.para[nu.type[i]],")", sep="")))    
    eval(parse(text=paste("round(x",i+n.ca+n.bi, ", digits=2)", sep=""))) #e.g. 1.4310798712 -> 1.43   
  }
  




###-------------------------------------------------------------------------
### 2.2.To generate a response variable based on the predictors
###-------------------------------------------------------------------------

## 2.2.1. creating individual associations between predictors and a response
# unimodal peak (x7)
  y_x7 = numeric(n)
  y_x7 = 10*cos(pi*(x7-median(x7))/max(x7))

# linear with a threshold (x8)
  y_x8 = numeric(n)
  y_x8[which(x8<=median(x8))] = y_x8[which(x8<=median(x8))] - 0.6*(median(x8)-x8[which(x8<=median(x8))])    

# non-additive 3way interaction (x9|x1,x4)
  y_x9 = numeric(n)
  for(i in c(1:n)){ if((x1[i]==0)&&(x4[i]!="low")){y_x9[i] = y_x9[i] - 0.5*x9[i]} }



## 2.2.2. combining them, adding an error term, and replacing partial data by NA
# setting the response variable
  y = numeric(n)
  y = y + y_x7+ y_x8 + y_x9
  set.seed(seed)
  y = y + rnorm(n,mean=0,sd=3)

# accomodating all variables into a single data frame
  data = data.frame(matrix(NA, nrow=n,ncol=n.bi+n.ca+n.nu+1))
  data[,1] = y
  col = c("Y")
  for(i in c(1:(n.bi+n.ca+n.nu))){
    eval(parse(text=paste("data[,i+1]=x",i, sep="")))    
    eval(parse(text=paste("col = append(col,'x",i,"')", sep="")))     
  }
  colnames(data)=col # Y, x1, x2, ...

## NA replacement
  for(i in c(1:(n.bi+n.ca+n.nu))){
    set.seed(seed+i)
    row.NA = sample(c(1:n),NA_ratio*n)
    data[row.NA,i+1]= NA    
  }

## 2.2.3. Summary in plots (Fig.2)
  par(mfrow = c(2,3))
  hist(y)
  boxplot(y~x1)
  boxplot(y~x4)
  plot(y~x7, ylim=c(-6,16), ann=F)
  par(new=T); plot(y_x7~x7, pch=3,col=2, ylim=c(-6,16), ann=F)
  plot(y~x8, ylim=c(-6,16), ann=F)
  par(new=T); plot(y_x8~x8, pch=3,col=2, ylim=c(-6,16))
  plot(y~x9, ylim=c(-6,16), ann=F)
  par(new=T); plot(y_x9~x9, pch=3,col=2, ylim=c(-6,16))
  par(mfrow = c(1,1))

#################################################################################
## 3. To analyze the dataset using Statistically-oriented Machine Learnings    ##
#################################################################################
###-------------------------------------------------------------------------
### 3.1. Conditional inference tree  # package party
###-------------------------------------------------------------------------
# Fig.3
  set.seed(seed)
  model.ct = ctree(Y~.,data=data, controls=ctree_control(testtype = "Bonferroni")) # see ctree_control for setting
  plot(model.ct)

###-------------------------------------------------------------------------
### 3.2. Model-based tree  # package partykit 
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
### 3.3. Permutation-based random forest
###        see section 1.2 for function
###        note that the function RF_Hapfelmeier may take extraordinally long time
###        depending on nperm and ntree. We firstly recommend to try with 
###        nperm=100, ntree=30 to run quickly.  
###-------------------------------------------------------------------------
# Fig.5
  formula = as.formula(paste("Y", paste(colnames(data)[2:length(colnames(data))], collapse=" + "), sep=" ~ ")) # creating a formula
  set.seed(seed)
  model.Hapfelmeier.RF = RF_Hapfelmeier(formula, data=data,nperm=2000, ntree=100)
  set.seed(seed)
  variable.importance = cbind(varimp(model.Hapfelmeier.RF$forest), model.Hapfelmeier.RF$p.values.bonf) # variable importance and p-value for each predictor
  plot(variable.importance[,1], xlab="predictor_x", ylab="variable importance", pch = ifelse(variable.importance[,2]<0.05, 16, 1))


###-------------------------------------------------------------------------
### 3.4. Partial dependence plot  # package mlr, party, randomForest
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
### 3.5. LM with AIC best model selection
###-------------------------------------------------------------------------
  formula = as.formula("Y~poly(x7, 2)+x8+x1*x4*x9+x2+x3+x5+x6+x10+x11+x12+x13+x14+x15+x16")
  model.lm <-step(lm(formula,data=na.omit(data)))
  summary(model.lm)
