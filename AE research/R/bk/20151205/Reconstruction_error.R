qqVec = c(25,50,75,100,200,300,400,500)
plot(qqVec,rErr,type='l',pch = 16, lwd = 2, col = "red",
     xlab="# of POD modes", ylab = "Average abs. recon. error",
     main = "POD reconstruction error")
