setwd("C:/Users/SimonLenovo/Documents/correlations")

for (ind in 4:9){
  load( paste0("params_resp",ind,"_time1.RData") )
  
  corrs = matrix(0,nrow=500,ncol=5)
  for (i in 1:500){
    corrs[i,] = params[[i]]$corr
  }
  mincut = -10
  corrs = pmax(corrs,mincut)
  
  #Histogram 
  hist(corrs[,1], breaks=20)
  hist(corrs[,2], breaks=20)
  hist(corrs[,3], breaks=20)
  hist(corrs[,4], breaks=20)
  hist(corrs[,5], breaks=20)
  print(colMeans(corrs))
}

#Plot over time
cut = 50
plot(1:cut,corrs[1:cut,1],type="l",col="red")
lines(1:cut,corrs[1:cut,2],type="l",col="blue")
lines(1:cut,corrs[1:cut,3],type="l",col="green")
lines(1:cut,corrs[1:cut,4],type="l",col="pink")
lines(1:cut,corrs[1:cut,5],type="l",col="black")
