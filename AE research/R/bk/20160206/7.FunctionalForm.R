getResp <- function(ind){
  if (ind == 1){
    return("u")
  }
  else if (ind == 2){
    return("v")
  }
  else if (ind == 3){
    return("w")
  }
  else if (ind == 4){
    return("P")
  }
  else if (ind == 5){
    return("T")
  }
  else{
    return("rho")
  }
}

getMode <- function(ind){
  bool <- T;
  cs <- cumsum(numModes)
  nn <- length(numModes)
  ii <- 1 #This will be the response
  while ((bool) && (ii<=nn)){
    if (ind > cs[ii]){
      ii = ii + 1
    }else{
      bool = F
    }
  }
  
  #Mode number
  if (ii == 1){
    num <- ind
  }else{
    num <- ind-cs[ii-1] 
  }
  
  return(list(ii,num))
  
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Computing cut-off for POD modes
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
indvec = 4:9

#Load design
des = as.matrix(read.csv('/home/proj/jeffwu/isye-ae/desnorm.csv'));

#Set cutoff for POD modes at 95% energy
numModes = rep(-1,length(indvec)) #Number of modes to choose from
BTlist = vector('list',length(indvec))
cut = 0.95

for (i in 1:length(indvec)){
  #Load data
  ind = indvec[i] #Current flow response
  setwd(paste0('/home/proj/jeffwu/isye-ae/bigmemory_bk/', ind))
  load('PODdat.RData')
  BTlist[[i]] = BTsvg
  
  #Compute cut-off
  cs = cumsum(eigval)/sum(eigval)
  cs = pmax(cs-cut,0)
  cs[cs==0] = Inf
  numModes[i] = which.min(cs)
}

#Print significant correlations at time tt
tt <- 1
setwd(paste0('/home/proj/jeffwu/isye-ae/'))
load(paste0('kriging_T/allparams_t',tt,'.RData'))
idxnz <- params$idxnz
testmat <- params$Ticov
for (k in 1:nrow(idxnz)){
  tmp1 <- getMode(idxnz[k,1])
  tmp2 <- getMode(idxnz[k,2])
  if (tmp1[[1]] != tmp2[[1]]){
    print( paste0( getResp(tmp1[[1]])," - Mode ",tmp1[[2]],
                   " <=> ",
                   getResp(tmp2[[1]])," - Mode ",tmp2[[2]] ) )
  }
  print(paste0(idxnz[k,1],', ',idxnz[k,2],' - Covariance: ',testmat[idxnz[k,1],idxnz[k,2]]))
}

#Get coefficient for functional form
ii = 6
nnM = c(0,cumsum(numModes))
curind = (nnM[ii]+1):nnM[ii+1]
othind = setdiff(1:sum(numModes),curind)
Tcov = params$Tcov

T11 = Tcov[curind,curind]
T12 = Tcov[curind,othind]
T22 = Tcov[othind,othind]
T22inv = solve(T22)

wtMat <- T12%*%T22inv
tmp <- wtMat[15,]
sort(tmp[tmp!=0])
tmpind = which(abs(tmp)>1e-4)
tmpind
tmp[tmpind]

# Testing:
# tmp = params$Tcov %*% params$Ticov
# head(t(head(tmp)))
# sum(abs(tmp-diag(nrow(params$Tcov))))

# Correlations at time 1
[1] "P - Mode 1 <=> u - Mode 1"
[1] "v - Mode 1 <=> u - Mode 2"
[1] "rho - Mode 4 <=> u - Mode 6"
[1] "v - Mode 130 <=> u - Mode 32"
[1] "u - Mode 2 <=> v - Mode 1"
[1] "P - Mode 2 <=> v - Mode 1"
[1] "rho - Mode 15 <=> v - Mode 11"
[1] "u - Mode 32 <=> v - Mode 130"
[1] "P - Mode 1 <=> w - Mode 1"
[1] "T - Mode 1 <=> w - Mode 2"
[1] "rho - Mode 1 <=> w - Mode 2"
[1] "rho - Mode 4 <=> w - Mode 3"
[1] "T - Mode 3 <=> w - Mode 4"
[1] "T - Mode 19 <=> w - Mode 9"
[1] "u - Mode 1 <=> P - Mode 1"
[1] "w - Mode 1 <=> P - Mode 1"
[1] "v - Mode 1 <=> P - Mode 2"
[1] "w - Mode 2 <=> T - Mode 1"
[1] "rho - Mode 1 <=> T - Mode 1"
[1] "w - Mode 4 <=> T - Mode 3"
[1] "rho - Mode 3 <=> T - Mode 3"
[1] "rho - Mode 7 <=> T - Mode 11"
[1] "w - Mode 9 <=> T - Mode 19"
[1] "w - Mode 2 <=> rho - Mode 1"
[1] "T - Mode 1 <=> rho - Mode 1"
[1] "T - Mode 3 <=> rho - Mode 3"
[1] "u - Mode 6 <=> rho - Mode 4"
[1] "w - Mode 3 <=> rho - Mode 4"
[1] "T - Mode 11 <=> rho - Mode 7"
[1] "v - Mode 11 <=> rho - Mode 15"