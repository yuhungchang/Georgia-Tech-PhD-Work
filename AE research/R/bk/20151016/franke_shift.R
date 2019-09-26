library(randtoolbox)

##############################################################################
# 2-D Test function: Franke (scaled to [0,1]^2)
##############################################################################
frankesh <- function(xx,shift)
{
  #Scales function to [0,1]^2
  if ((xx[1]<0)|(xx[1]>1)|(xx[2]<0)|(xx[2]>1)){
    return (0)
  }
  else{
    
    lower1 = shift
    upper1 = 1+shift
    lower2 = shift
    upper2 = 1+shift
    
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

# #Testing plot
# shift = 0.05
# numPred = 40 #number of prediction in each dim
# xvect = seq(from=0,to=1,length.out=numPred);
# yvect = seq(from=0,to=1,length.out=numPred);
# grd = expand.grid(x=xvect, y=yvect);
# grd$truefn = apply(grd,1,frankesh,shift);
# z = matrix(grd$truefn,nrow=numPred)
# filled.contour(xvect,yvect,z)

n = 50 #Number of different control runs
m = 3  #Number of observations per control setting
d = 2  #Number of control factors
x = sobol(n,2) #Control settings
dst1 = as.matrix(dist(x[,1]))  #Distance matrix for dim. 1
dst2 = as.matrix(dist(x[,2]))  #Distance matrix for dim. 1

y1 = apply(x,1,frankesh,0) #First response: no shift
y2 = apply(x,1,frankesh,0.05) #Second response: shifted 0.05
y3 = -apply(x,1,frankesh,-0.05) #Third response: shifted -0.05 and negated
Y = t(rbind(y1,y2,y3)) #Vectorized response: (y11, ..., y1m, ..., yn1, ..., ynm)
X = x
y <- c(t(Y))
