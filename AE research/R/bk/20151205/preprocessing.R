library(xlsx)
data_path <- "../Data/ResultsOfThinckness_0406.xlsx"
## normalize design matrix to [0,1]^5 ##
design.df <- read.xlsx(data_path, sheetIndex = 1, header = TRUE)
design.df <- design.df[,-c(1,ncol(design.df))]
colnames(design.df) <- c("L", "R_n", "theta", "delta", "dL")
range.ls <- list("L" = c(20,100), "R_n" = c(2,5), 
                 "theta" = c(45, 75), "delta" = c(0.5, 2), "dL" = c(1, 4))
design_norm.df <- design.df
for(name in colnames(design.df)) design_norm.df[,name] <- (design.df[,name] - range.ls[[name]][1])/diff(range.ls[[name]]) 
write.csv(design_norm.df, "../Data/desnorm.csv", row.names = FALSE)

## get output data ##
csvdata.df <- NULL
sheetnumber <- 30
for (i in 1:sheetnumber){
  data.df <- read.xlsx(data_path, sheetIndex = i + 1, header = TRUE)
  data.df <- data.df[ , !is.element(colnames(data.df), c("NA..1", 
                                                        "Donstream.Angle..2alpha..5..10..4.5.10.6",
                                                        "Donstream.Angle..2alpha..1.10..4.5.10.5"))]
  colnames(data.df) <- c("Time (ms)", "r_mn (mm)", "phi", "h (mm)", "2 alpha_inv(deg)", "Angle (2alpha) All",  "Angle (2alpha) exit max", "r_mn/R_n", "R_n")
  
  if(all(is.na(data.df[nrow(data.df),]))) data.df <- data.df[-nrow(data.df),]
  data.df <- cbind("Run" = rep(i, nrow(data.df)), data.df)
  if(i == 1L)  csvdata.df <- data.df 
  else csvdata.df <- rbind(csvdata.df, data.df)
}
csvdata.df[,"Time (ms)"] <- as.numeric(unlist(strsplit(as.character(csvdata.df[,"Time (ms)"]), " ms")))
csvdata.df <- csvdata.df[,c("Run", "Time (ms)", "h (mm)", "Angle (2alpha) All")]
colnames(csvdata.df) <- c("sim", "time", "h", "angle")
write.csv(csvdata.df, "../Data/export.csv", row.names = FALSE)
