###########################################################################
#
#  Masahiro Ryo and Matthias C. Rillig. 2017. 
#  Statistically-reinforced machine learning for 
#  nonlinear patterns and variable interactions. Ecosphere.
#    Source: https://github.com/masahiroryo/R_Statistically-rainforced-machine-learning
#
###########################################################################
#
# To generate an artifical dataset
#   1. To generate binary, categorical, and numeric variables (predictors)
#   2. To generate a response variable based on the predictors
#
#   INPUT: some numbers specified in section 1.1
#   OUTPUT: "data_example.csv" (This data is used in 'SML_application.R')
#
###########################################################################


###-------------------------------------------------------------------------
### 1. To generate binary, categorical, and numeric variables (predictors)
###-------------------------------------------------------------------------

## 1.1. Settings of data sizes 
  n = 300               # sample size (the study tried 100, 200, 300)
  n.bi = 3              # number of binary predictors
  n.ca = 3              # number of categorical predictors
  n.nu = 10             # number of numeric predictors
  NA_ratio = 0.0       # proportion of NA for each predictors (we did not include NA)
  seed = 1              # number for set.seed for random number generation

## 1.2. binary data generation 
  for(i in c(1:n.bi)){
    set.seed(seed+2*i)  # changing random seed
    p = rnorm(1,mean=0.5,sd=0.1) # a probability for binomial dist.
    set.seed(seed+2*i)  
    eval(parse(text=paste("x",i, " = rbinom(",n,",","1,",p,")",sep=""))) # e.g. x1=rbinom(300,1,0.4)
  }


## 1.3. categorical data generation
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


## 1.4. numeric data generation
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
### 2.To generate a response variable based on the predictors
###-------------------------------------------------------------------------

## 2.1. creating individual associations between predictors and a response
# unimodal peak (x7)
  y_x7 = numeric(n)
  y_x7 = 10*cos(pi*(x7-median(x7))/max(x7))

# linear with a threshold (x8)
  y_x8 = numeric(n)
  y_x8[which(x8<=median(x8))] = y_x8[which(x8<=median(x8))] - 0.6*(median(x8)-x8[which(x8<=median(x8))])    

# non-additive 3way interaction (x9|x1,x4)
  y_x9 = numeric(n)
  for(i in c(1:n)){ if((x1[i]==0)&&(x4[i]!="low")){y_x9[i] = y_x9[i] - 0.5*x9[i]} }



## 2.2. combining them, adding an error term, and replacing partial data by NA
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

  
  write.csv(data,file="data_example.csv", row.names = F)
###-------------------------------------------------------------------------
  
