library(readr)
library(lme4)
library(mada)
library(ggplot2)
library(dplyr)

X <- requi
X<-requi%>%filter(group=="S",cutpoint==100)
X<-as.data.frame(X)
N <- length(X$tpos)

X$n1 <- X$tpos+X$fneg
X$n0 <- X$fpos+X$tneg
X$true1 <- X$tpos
X$true0 <- X$tneg

X$study <- 1:N


Y = reshape(X, direction = "long", varying = list( c("n1" , "n0") , c( "true1","true0" ) ) ,
            timevar = "sens" , times = c(1,0) , v.names = c("n","true") ) 

# Sort data by study to cluster the 2 records per study together 
Y = Y[order(Y$id),]  
Y$spec<- 1-Y$sens

MA_Y = glmer( formula = cbind(  true , n - true ) ~ 0 + sens + spec + (0+sens + spec|study),
              data = Y , family = binomial , nAGQ=1 , verbose=0 )

# More detail can be obtained by using the summary command 
ma_Y = summary(MA_Y)
labels(ma_Y) 

lsens = ma_Y$coeff[1,1]
lspec = ma_Y$coeff[2,1]



# Start by calculating derived and set parameters
seB     <- ma_Y$coefficients[2,2]
seA     <- ma_Y$coefficients[1,2]
r       <- ma_Y$vcov@x[2] / (seA*seB)
varA    <- ma_Y$varcor$study[1,1]
varB    <- ma_Y$varcor$study[2,2]
sAB     <- ma_Y$varcor$study[1,2]
covAB   <- ma_Y$vcov@x[2]
sepredA <- sqrt(varA + seA^2)
sepredB <- sqrt(varB + seB^2)
rpredAB <- (sAB + covAB) / (sepredA*sepredB)
level   <- 95
f       <- qf(0.95, df1=2, df2=N-2)
croot   <- sqrt(2*f)


# Empty data frames
conf_region <- (rep(NA, 361))
pred_region <- (rep(NA, 361))

# Confidence region
for (i in seq(0, 2*pi, length.out=361)){
  confB <- lspec + (seB*croot*cos(i))
  confA <- lsens + (seA*croot*cos(i + acos(r)))
  confsens <- plogis(confA)
  confspec <- plogis(confB)
  conf_i <- data.frame(X=1-confspec, Y=confsens)
  # Add results for most recent value to a data frame which contains the results of previous values
  conf_region<-rbind(conf_region, conf_i)
}
conf_region <- conf_region[2:362,]

# Predictive region
for (i in seq(0, 2*pi, length.out=361)){
  predB <- lspec + (sepredB*croot*cos(i))
  predA <- lsens + (sepredA*croot*cos(i + acos(rpredAB)))
  predsens <- plogis(predA)
  predspec <- plogis(predB)
  pred_i <- data.frame(X=1-predspec, Y=predsens)
  # Add results for most recent value to a data frame which contains the results of previous values
  pred_region<-rbind(pred_region, pred_i)
}
pred_region <- pred_region[2:362,]

# Sensitivity and specificity calculations for each study
X$sens <- X$tpos / (X$tpos + X$fneg)
X$spec <- X$tneg / (X$fpos + X$tneg)


study_level <- data.frame(ID=X$study, TP=X$tpos, FN=X$fneg, FP=X$fpos, TN=X$tneg, N=(X$tpos+X$fneg+X$fpos+X$tneg)
                          ,Sensitivity=X$sens, Specificity=X$spec)

# Mean point
Sens = plogis(lsens) 
Spec = plogis(lspec) 


mean_point=data.frame(x=Spec, y=Sens)


c<-ggplot(study_level,aes(Specificity,Sensitivity)) + geom_path(aes((1-X),Y,lty="95% prediction region"),data=pred_region)+
  geom_path(aes((1-X),Y,lty="95% confidence region"), data=conf_region)
c + geom_point() +xlim(1,0)+ylim(0,1)+geom_point(aes(x,y,col="Summary Estimate"),data=mean_point)