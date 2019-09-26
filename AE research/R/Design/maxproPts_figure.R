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
range.ls <- list("L (mm)" = c(20, 100), "R_n (mm)" = c(2, 5), 
                 "theta (degree)" = c(45, 75), "delta  (mm)" = c(0.5, 2.0), "dL  (mm)" = c(1,4))
numIter <- 100 #repeat for numIter to find best points
n <- 30
Design1 <- MaxPro.FUN(range.ls, numIter = numIter, n = n)
pdf("../JASA paper/Figures/Projection_Two_dimension.pdf")
pairs(Design1, pch = 12, labels = c("L (mm)", expression(paste(R[n]," (mm)")),
                                   expression(paste(theta, " (degree)")), expression(paste(delta, " (mm)")),
                                   expression(paste(Delta,"L (mm)"))))
dev.off()

#save
write.csv(Design1, file="Dropbox/Aerospace injector/Design/MaxDesign_n30.csv")

