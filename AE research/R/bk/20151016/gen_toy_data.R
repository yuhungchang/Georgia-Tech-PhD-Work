library(randtoolbox)
##############################################################################
#Franke function scaled to [0,1]^2
##############################################################################
frankesc <- function(xx)
{
  #Scales function to [0,1]^2
  if ((xx[1]<0)|(xx[1]>1)|(xx[2]<0)|(xx[2]>1)){
    return (0)
  }
  else{
    
    lower1 = 0
    upper1 = 1
    lower2 = 0
    upper2 = 1
    
    x1 <- xx[1]*(upper1-lower1) + lower1
    x2 <- xx[2]*(upper2-lower2) + lower2
    
    term1 <- 0.75 * exp(-(9*x1-2)^2/4 - (9*x2-2)^2/4)
    term2 <- 0.75 * exp(-(9*x1+1)^2/49 - (9*x2+1)/10)
    term3 <- 0.5 * exp(-(9*x1-7)^2/4 - (9*x2-3)^2/4)
    term4 <- -0.2 * exp(-(9*x1-4)^2 - (9*x2-7)^2)
    
    y <- term1 + term2 + term3 + term4
    return(y)
  }
}
##############################################################################
#Toy flow: Franke 
#                  + periodic perturbation (decreases with time)
#                   + effect for 2 control factors
##############################################################################
flowtemp <- function(xx,t,c){
                3*sqrt(c[1])*frankesc(xx) + c[2]*(1-t^2)*(sin(1/(xx[1]+0.05)) + 2*sin(1/(xx[2]+0.05))) + 5*sqrt(c[1])*sqrt(c[2])*(3-(xx[1]-0.8)^2-(xx[2]-0.7)^2)
}

##############################################################################
#Plot contours of toy flow
##############################################################################
co <- c(0.25,0.25)
t <- 1
fn <- function(xx){flowtemp(xx,t,co)}
numPred = 40 #number of prediction in each dim
xvect = seq(from=0,to=1,length.out=numPred);
yvect = seq(from=0,to=1,length.out=numPred);
grd = expand.grid(x=xvect, y=yvect);
grd$truefn = apply(grd,1,fn);
z = matrix(grd$truefn,nrow=numPred)
filled.contour(xvect,yvect,matrix(grd$truefn,nrow=numPred))

##############################################################################
#Output toy flow data
##############################################################################
timeSteps = 20 #Number of time steps
timeVal = seq(from=0,to=1,length.out=timeSteps)
gridSize = 1000 #Number of grid points
grd = sobol(gridSize,2)
contSize = 30 #Number of control samples
cont = sobol(contSize,2)

for (i in 1:contSize){
  for (j in 1:timeSteps){
    t <- timeVal[j]
    co <- cont[i,]
    fn <- function(xx){flowtemp(xx,t,co)}
    resp <- apply(grd,1,fn);
    toSave <- cbind(grd,resp);
    write.csv(toSave,file=sprintf("flow_%i_$i.csv",i,j))
  }
}
