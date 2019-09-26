library(MaxPro)
MaxPro.FUN <- function(range.ls, numIter = numIter, n){
  
  # Generating points
  p <- length(range.ls) #dimensions
  
  InitialDesign <- MaxProLHD(n, p)$Design #global search for LHD
  DOX <- MaxPro(InitialDesign, iteration = numIter) #local search from opt LHD
  D <- DOX #design points
  
  for (i in 2:numIter){
    print(i)
    InitialDesign <- MaxProLHD(n, p)$Design #global search for LHD
    DOX <- MaxPro(InitialDesign, iteration = numIter) #local search from opt LHD 
    if (DOX$measure < D$measure){
      D <- DOX
    } #design points
  }
  
  design.df <- D$Design
  colnames(design.df) <- names(range.ls)
  
  design_norm.df <- design.df
  for(name in colnames(design.df)) design_norm.df[,name] <- range.ls[[name]][1] + design.df[,name] * diff(range.ls[[name]])
  return(design_norm.df)
}

# Parameter Setting
range.ls <- list("diameter (mm)" = c(8.40, 14.5), "delta (mm)" = c(0.41, 0.77), 
                 "theta (degree)" = c(34.8, 64.6), "length (mm)" = c(7.07, 13.1))
numIter <- 100 #repeat for numIter to find best points
n <- 40
Design1 <- MaxPro.FUN(range.ls, numIter = numIter, n = n)

#save
write.csv(Design1, file="Dropbox/Aerospace injector/Design/Design_parameter4_n40.csv")

