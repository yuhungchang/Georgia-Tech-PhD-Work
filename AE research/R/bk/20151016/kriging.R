library(GPfit)
library(lhs)
library(lattice)
library(gridExtra)

############################################################################
#Import data
############################################################################
dat = read.csv("../Data/export.csv");
des = read.csv("../Data/desnorm.csv");
factNames = names(des);

############################################################################
#Fitting the kriging model
############################################################################
attach(dat)
nRuns = 30;
curdat = dat[time==50,c(3,4)]; 
attach(curdat);
attach(des);
hfit = GP_fit(des,h); #fitted kriging model for h
anglefit = GP_fit(des,angle); # ... for angle

summary.GP(hfit)
summary.GP(anglefit)

############################################################################
#Plot factors one-at-a-time, fixing other four at midpoint
############################################################################
#Factor to vary
fact = 4
#Grid to plot
numPred = 100;
#Remaining factors to fix
fixInd = (1:5)[-fact]; 
prd = matrix(0,nrow=numPred,ncol=5);
prd[,fixInd] = rep(0.5,numPred); #fix other 4 factors at midpoint
prd[,fact] = seq(from=0,to=1,length.out=numPred);
hpred = predict.GP(hfit,prd);
anglepred = predict.GP(anglefit,prd);

#plot the results
plot(seq(from=0,to=1,length.out=numPred),hpred$Y_hat,xlab=factNames[[fact]],ylab="Film Thickness",type='l',ylim = c(min(h),max(h)))
points(des[1:nRuns,fact],h,pch=16,col=2)
plot(seq(from=0,to=1,length.out=numPred),anglepred$Y_hat,xlab=factNames[[fact]],ylab="Angle",type='l',ylim = c(min(angle),max(angle)))
points(des[1:nRuns,fact],angle,pch=16,col=2)

############################################################################
#Plot factors two-at-a-time, fixing other 3 at midpoint
############################################################################
#Factors to vary
facts = c(2,3);
#Grid size to plot
numPred = 100;
#Factors to be fixed
fixInds = (1:5)[-facts];
xvect = seq(from=0,to=1,length.out=numPred);
grd = expand.grid(x=xvect, y=xvect);
grd = as.matrix(grd);
prd = matrix(0,nrow=numPred^2,ncol=5);

#Values to fix
prd[,fixInds] = rep(0.5,numPred^2);
prd[,facts] = grd;

#Do the predictions and plot
hpred = predict.GP(hfit,prd);
anglepred = predict.GP(anglefit,prd);
hfitframe = data.frame(cbind(grd,hpred$Y_hat));
trellis.par.set("axis.line",list(col=NA,lty=1,lwd=1)) #MOD: Sets outer box to invisible
wireframe(V3~x*y,data = hfitframe, scales = list(arrows=FALSE),
          screen = list(z=120,x=-70,y=0), drape = T,
          zlab = list("Film Thickness", rot = 90), xlab = factNames[[facts[1]]], ylab = factNames[[facts[2]]])
# MOD:
# zlab,xlab,ylab - Setting the caption on the z-axis, x-axis and y-axis
# screen         - Rotates the plot around (x rotates the x-axis, y rotates the y-axis, z rotates the z-axis)

anglefitframe = data.frame(cbind(grd,anglepred$Y_hat));
trellis.par.set("axis.line",list(col=NA,lty=1,lwd=1)) #Sets outer box to invisible
wireframe(V3~x*y,data = anglefitframe, scales = list(arrows = FALSE), screen = list(z=60,x=-60), drape = T,zlab = expression(alpha), xlab = factNames[[facts[1]]], ylab = factNames[[facts[2]]])

############################################################################
#Compare 18-run results with 30-run results (5*p = 25 recommended)
############################################################################
nRuns1 = 18;
nRuns2 = 30;
hfit1 = GP_fit(des[1:nRuns1,],h[1:nRuns1]); # kriging model for h
anglefit1 = GP_fit(des[1:nRuns1,],angle[1:nRuns1]); # kriging model for angle
hfit2 = GP_fit(des[1:nRuns2,],h[1:nRuns2]);
anglefit2 = GP_fit(des[1:nRuns2,],angle[1:nRuns2]);

#Coefficient summary
summary(hfit1)
summary(hfit2) 

summary(anglefit1)
summary(anglefit2)

#Sensitivity analysis summary (see SA.R)

#Factors to vary
facts = c(4,5);
#Factors to be fixed
fixInds = (1:5)[-facts];
numPred = 100; #Grid size to plot
xvect = seq(from=0,to=1,length.out=numPred);
grd = expand.grid(x=xvect, y=xvect);
grd = as.matrix(grd);
prd = matrix(0,nrow=numPred^2,ncol=5);

#Values to fix
prd[,fixInds] = rep(0.5,numPred^2);
prd[,facts] = grd;

#Do the predictions and plot
hpred1 = predict.GP(hfit1,prd);
anglepred1 = predict.GP(anglefit1,prd);
hpred2 = predict.GP(hfit2,prd);
anglepred2 = predict.GP(anglefit2,prd);

hfitframe1 = data.frame(cbind(grd,hpred1$Y_hat,sqrt(hpred1$MSE)));
plot1p = wireframe(V3~x*y,data = hfitframe1, scales = list(arrows=FALSE),screen = list(z=120,x=-60), drape = T,
                  zlab = "Film Thickness", xlab = factNames[[facts[1]]], ylab = factNames[[facts[2]]],
                  zlim = c(min(hpred1$Y_hat,hpred2$Y_hat),max(hpred1$Y_hat,hpred2$Y_hat)),
                  main = "18 runs - Predicted surface")
plot1m = wireframe(V4~x*y,data = hfitframe1, scales = list(arrows=FALSE),screen = list(z=120,x=-60), drape = T,
                   zlab = "Film Thickness", xlab = factNames[[facts[1]]], ylab = factNames[[facts[2]]],
                   zlim = c(min(sqrt(hpred1$MSE),sqrt(hpred2$MSE)),max(sqrt(hpred1$MSE),sqrt(hpred2$MSE))),
                   main = "18 runs - Prediction error (1se)")
hfitframe2 = data.frame(cbind(grd,hpred2$Y_hat,sqrt(hpred2$MSE)));
plot2p = wireframe(V3~x*y,data = hfitframe2, scales = list(arrows=FALSE),screen = list(z=120,x=-60), drape = T,
                  zlab = "Film Thickness", xlab = factNames[[facts[1]]], ylab = factNames[[facts[2]]],
                  zlim = c(min(hpred1$Y_hat,hpred2$Y_hat),max(hpred1$Y_hat,hpred2$Y_hat)),
                  main = "30 runs - Predicted surface")
plot2m = wireframe(V4~x*y,data = hfitframe2, scales = list(arrows=FALSE),screen = list(z=120,x=-60), drape = T,
                   zlab = "Film Thickness", xlab = factNames[[facts[1]]], ylab = factNames[[facts[2]]],
                   zlim = c(min(sqrt(hpred1$MSE),sqrt(hpred2$MSE)),max(sqrt(hpred1$MSE),sqrt(hpred2$MSE))),
                   main = "30 runs - Prediction error (1se)")

grid.arrange(plot1p,plot1m,plot2p,plot2m,nrow=2,ncol=2,byrow = T)
#              main=textGrob("18-run vs. 30-run",gp=gpar(font=2,fontsize=25)))
