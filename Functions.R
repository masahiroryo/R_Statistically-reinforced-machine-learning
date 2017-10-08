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
# To define the two functions for permutation-based random forest
#   1. myvarimp
#   2. RF_Hapfelmeier
#
#   INPUT: none
#   OUTPUT: none (This script is called by the main script: "SML_application.R")
#
###########################################################################

  package.list = c("party", "caret","foreach","doSNOW")
  tmp.install = which(lapply(package.list, require, character.only = TRUE)==FALSE)
  if(length(tmp.install)>0) install.packages(package.list[tmp.install])
  lapply(package.list, require, character.only = TRUE)

###-------------------------------------------------------------------------
### 1. myvarimp
###    function for calculating variable importance measure (modified by Masahiro Ryo)
###-------------------------------------------------------------------------
  myvarimp = function(object, mincriterion = 0, conditional = FALSE, 
                      pre1.0_0 = conditional, varID) {
    response = object@responses
    input = object@data@get("input")
    xnames = colnames(input)
    inp = initVariableFrame(input, trafo = NULL)
    y = object@responses@variables[[1]]
    if (length(response@variables) != 1) stop("cannot compute variable importance measure")
    CLASS = all(response@is_nominal)
    ORDERED = all(response@is_ordinal)
    if (CLASS) {
      error = function(x, oob) mean((levels(y)[sapply(x, which.max)] != y)[oob])
    }else {
      if (ORDERED) {
        error = function(x, oob) mean((sapply(x, which.max) != y)[oob])
      } else {
        error = function(x, oob) mean((unlist(x) - y)[oob]^2)}
    }
    perror = rep(0, length(object@ensemble))
    for (b in 1:length(object@ensemble)) {
      tree = object@ensemble[[b]]
      w = object@weights[[b]]
      w[w == 0] = 1
      oob = object@weights[[b]] == 0
      p = .Call("R_predict", tree, inp, mincriterion, -1L, PACKAGE = "party")
      eoob = error(p, oob)
      j = varID
      p = .Call("R_predict", tree, inp, mincriterion, as.integer(j), PACKAGE = "party")
      perror[(b - 1)] = (error(p,oob) - eoob)
    }
    return(MeanDecreaseAccuracy = mean(perror))
  }
  environment(myvarimp) = environment(varimp)
  
###-------------------------------------------------------------------------
### 2. Permutation-based random forest function (modified by Masahiro Ryo)
###   
###-------------------------------------------------------------------------
  RF_permutation = function(formula, data, nperm = 500, ntree = 100, ncore=4, alpha=0.05)
    # formula: object of class "formula".
    # data: data frame containing the variables.
    # nperm: Number of permutation steps used for the permutation test.
    # ntree: Number of trees in the Random Forest.
    # ncore: Number of cores used for parallel computing
  {
    x.names = all.vars(formula)[-1]
    y.names = all.vars(formula)[1]
    terms. = terms.formula(formula)
    x.formula = attr(terms., "term.labels")
    y.formula = as.character(attr(terms., "variables"))[2]
    mtry = ceiling(sqrt(length(x.formula)))
    dat = subset(data, select = c(y.names, x.names))
    forest = party::cforest(formula, data = dat, controls = cforest_unbiased(mtry = mtry, ntree = ntree))
    obs.varimp = varimp(forest)
    perm.mat = matrix(NA, ncol = length(x.names), nrow = nperm, dimnames = list(1:nperm, x.names))
    
    cl=makeCluster(ncore) #change to your number of CPU cores
    registerDoSNOW(cl)
    
    for (j in x.names) {
      cat("\r", "Processing variable ", which(j == x.names), " of ", length(x.names)); flush.console()
      perm.dat = dat
      perm.mat[, j] = unlist(foreach(i = 1:nperm, .packages='party',.export="myvarimp") %dopar% {
        perm.dat[, j] = sample(perm.dat[, j]);
        myvarimp(party::cforest(formula, data = perm.dat, controls = cforest_unbiased(mtry = mtry, ntree = ntree)), varID = which(x.names == j))
      })
    }
    stopCluster(cl)
    
    p.vals = sapply(x.names, function(x) sum(perm.mat[, x] >= obs.varimp[which(x == x.names)]) / nperm)
    p.vals.bonf = p.vals * length(p.vals)
    
    if (any(p.vals < alpha)) {
      selection = names(p.vals)[which(p.vals < alpha)]
      mtry = ceiling(sqrt(length(selection)))
      forest = cforest(as.formula(paste(y.formula, ".", sep = " ~ ")), data = subset(dat, select = c(y.names, selection)),
                       controls = cforest_unbiased(mtry = mtry, ntree = ntree))
      p = p.vals[which(p.vals < alpha)]
    }
    if (any(p.vals.bonf < alpha)) {             
      selection.bonf = names(p.vals.bonf)[which(p.vals.bonf < alpha)]                         
      mtry = ceiling(sqrt(length(selection.bonf)))
      forest.bonf = party::cforest(as.formula(paste(y.formula, ".", sep = " ~ ")), data = subset(dat, select = c(y.names, selection.bonf)),
                                   controls = cforest_unbiased(mtry = mtry, ntree = ntree))
      p.bonf = p.vals.bonf[which(p.vals.bonf < alpha)]
      varimp.bonf = varimp(forest.bonf)
      accuracy.fitting = postResample(predict(forest.bonf), subset(dat, select = y.names))
      accuracy.validation = caret:::cforestStats(forest.bonf)
      varimp.R2 = accuracy.fitting[2]*varimp.bonf/sum(varimp.bonf)
      residual = subset(dat, select = y.names) - predict(forest.bonf)
    }
    if (!any(p.vals < alpha)) {
      selection = c(); forest = c(); p = c()
    }
    if (!any(p.vals.bonf < alpha)) {
      selection.bonf = c(); forest.bonf = c(); p.bonf = c(); varimp.bonf = c();
      accuracy.fitting = c(); accuracy.validation = c(); varimp.R2 = c(); residual = c()
    }
    Y = as.numeric(as.character(dat[,y.names]))
    oob.error = ifelse(length(selection) != 0, mean((Y - as.numeric(as.character(predict(forest, OOB = T))))^2), mean((Y - ifelse(all(Y %in% 0:1), round(mean(Y)), mean(Y)))^2))
    oob.error.bonf = ifelse(length(selection.bonf) != 0, mean((Y - as.numeric(as.character(predict(forest.bonf, OOB = T))))^2), mean((Y - ifelse(all(Y %in% 0:1), round(mean(Y)), mean(Y)))^2))
    cat("\n", "\r"); flush.console()
    return(list("selection" = selection, "forest" = forest, "oob.error" = oob.error, "p.values" = p.vals,
                "selection.bonf" = selection.bonf, "forest.bonf" = forest.bonf, "oob.error.bonf" = oob.error.bonf, "p.values.bonf" = p.vals.bonf,
                "varimp" =obs.varimp, "varimp.selection"=varimp.bonf, "fitting" = accuracy.fitting, "validation" = accuracy.validation, "varimp.R2"=varimp.R2, "residual" =residual))
  }