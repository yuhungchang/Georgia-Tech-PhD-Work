############################################################################
#Import data
############################################################################
dat = read.csv("../Data/export.csv")
des = read.csv("../Data/desnorm.csv")
factNames = names(des)

############################################################################
#Kriging
############################################################################
library(GPfit)
nRuns = 30
curdat = dat[dat[,"time"]==50,c(3,4)]
attach(curdat)
attach(des)
hfit = GP_fit(des[1:nRuns,1:5],h) #fitted kriging model for h
anglefit = GP_fit(des[1:nRuns,1:5],angle) # ... for angle

summary.GP(hfit)
summary.GP(anglefit) 

library(sensitivity)
n <- 1000
X1 <- data.frame(matrix(runif(5 * n), nrow = n))
X2 <- data.frame(matrix(runif(5 * n), nrow = n))
X3 <- data.frame(matrix(runif(5 * n), nrow = n))
colnames(X1) <- colnames(X2) <- colnames(X3) <- colnames(des)

hpred.fun <- function (X) {
  predict.GP(hfit, X)$Y_hat
}
hsa <- sobolowen(model = hpred.fun, X1, X2, X3, nboot = 100)
print(hsa)
plot(hsa)
title("First order and total sensitivity indices \nfor thickness response")

anglepred.fun <- function (X) {
  predict.GP(anglefit, X)$Y_hat
}
anglesa <- sobolowen(model = anglepred.fun, X1, X2, X3, nboot = 100)
print(anglesa)
plot(anglesa)
title("First order and total sensitivity indices \nfor angle response")



##  Sobol’ scheme (Sobol, 1993) to compute the indices given by the variance decomposition up to a specified order
h_secondorder <- sobol(model = hpred.fun, X1 = X1, X2 = X2, order = 2, nboot = 100)
print(h_secondorder)
library(plotrix)
plotCI(1:nrow(h_secondorder$S), 
       h_secondorder$S[,"original"], 
       ui = h_secondorder$S[,"max. c.i."], li = h_secondorder$S[,"min. c.i."],
       ylim = c(0,1), xaxt = "n", xlab = "", ylab = "Sobol' index")
axis(1, 1:nrow(h_secondorder$S), rownames(h_secondorder$S), las=2)
title("Sensitivity indices \nfor thickness (30-runs)")

angle_secondorder <- sobol(model = anglepred.fun, X1 = X1, X2 = X2, order = 2, nboot = 100)
print(angle_secondorder)
plotCI(1:nrow(angle_secondorder$S), 
       angle_secondorder$S[,"original"], 
       ui = angle_secondorder$S[,"max. c.i."], li = angle_secondorder$S[,"min. c.i."],
       ylim = c(0,1), xaxt = "n", xlab = "", ylab = "Sobol' index")
axis(1, 1:nrow(angle_secondorder$S), rownames(angle_secondorder$S), las=2)
title("Sensitivity indices \nfor angle (30-runs)")


