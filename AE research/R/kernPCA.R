library(kernlab)
library(bigmemory)
ind <- 8
setwd(paste0("/home/proj/jeffwu/isye-ae/bigmemory_bk/", ind))
FTtmp1 <- attach.big.matrix(paste0('FTtmp1.dsc'))
start.time <- proc.time()
kpc <- kpca(t(FTtmp1[,1:999]), kernel="rbfdot",
            kpar=list(sigma=2),features=30)
print(proc.time() - start.time)
save(kpc, file = "../KernPCA/KernPCA_case1.RData")