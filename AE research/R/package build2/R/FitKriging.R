FitKriging <- function(num.case, des, 
                       ind_assumption.fg = TRUE, 
                       POD_filepath, output_path){
  
  
  cat("3. Fitting Kriging Model .....\n")
  if(ind_assumption.fg) cat("  Assume Response are independent (w/o T matrix): .....\n")
  else cat("  Do not assume Response are independent (w/ T matrix): .....\n")
  
  # create output_path folder
  if(!dir.exists(output_path)) dir.create(output_path)
  if(ind_assumption.fg & !dir.exists(paste0(output_path, "/ind.assumption"))) dir.create(paste0(output_path, "/ind.assumption"))
  
  des <- as.matrix(des)
  pars <- matrix(0, nrow = length(Snapshot.Index), ncol = ncol(des))
  
  cat("  3.1. Tuning correlation parameters .....\n")
  
  BTlist = vector("list", length(Response.ColumnIndex))
  for (i in 1:length(Response.ColumnIndex)){
    load(paste0(POD_filepath,"/", Response.ColumnIndex[i], "/PODdat.RData"))
    BTlist[[i]] = BTsvg
  }
  
  for (j in 1:length(Snapshot.Index)){
    snapshot <- Snapshot.Index[j]
    
    cat("     Tuning corr. param. for t = ", snapshot, "...\n")
    
    yMat = matrix(0, nrow = num.case, ncol = 0) #Matrix of POD coefficients (n x (numModes.numResponse))
    slc = seq(from = j, by = length(Snapshot.Index), length.out = num.case)
  
    numModes <- rep(0, length(Response.ColumnIndex))
    for (i in 1:length(Response.ColumnIndex)){
      #Find indices for each case
      yMat = cbind( yMat, BTlist[[i]][slc,] ) 
      numModes[i] <- ncol(BTlist[[i]])
    }
    
    #   load(paste0('yMat_',t,'.RData'))
    
    if(ind_assumption.fg){
      params <- lapply(data.frame(yMat), FUN = function(y.vt) GP_fit(des, y.vt, corr = list(type="exponential",power=2)))
      save(params, numModes, file = paste0(output_path, "/ind.assumption/allparams_t", snapshot, ".RData"))
    }else{
      #If j == 1, do a global search for correlation parameters
      if (j == 1){
        lower = -5
        upper = 5
        numIni = 100
        if(ncol(des) == 1){
          Mm = (runif(numIni))*(upper-lower)+lower
          Mm <- matrix(Mm, ncol = 1)
        }else{
          Mm = (maximinSLHD(1,numIni,ncol(des))$StandDesign)*(upper-lower)+lower
        }
        svList = vector('list',numIni)
        allVal = rep(-1,numIni)
        allPar = matrix(-1,nrow=numIni,ncol=ncol(des))
        for (i in 1:numIni){
          svList[[i]] <- lbfgs(Mm[i,], fn=GPdev.ind, control = list(xtol_rel = 1e-2), 
                               lower = rep(lower,ncol(des)), upper = rep(upper,ncol(des)), xMat=des, yMat=yMat, pw=2)
          allVal[i] = svList[[i]]$value
          allPar[i,] = svList[[i]]$par
        }
        pars[j,] = allPar[which.min(allVal),]
      }
      else{
        pars[j,] <- lbfgs(pars[(j-1),], fn=GPdev.ind, control = list(xtol_rel = 1e-2), 
                          lower = rep(lower,ncol(des)), upper = rep(upper,ncol(des)), xMat=des, yMat=yMat, pw=2)$par
      }
      
      corr <- pars[j,]
      print(corr)
      pm <- computeS(corr, des, yMat, 2);   #params[1] - mu, params[2] - Rinv, params[3] - data matrix for T
      vars <- apply(pm$Rsd, 2,function(xx){mean(xx^2)}) 
      
      cat("  3.2. Tuning T matrix .....\n")
      glhuge = huge(x = pm$Rsd, method = 'glasso', lambda.min.ratio=0.1)
      glsel = huge.select(glhuge)
      optind = which.min(glsel$ebic.score[-1])+1 #Find the best model which gives correlations
      invMat = cov2cor(as.matrix(glsel$icov[[optind]]))
      Tinv = cor2cov( (invMat + t(invMat))/2, 1/sqrt(vars*diag(glsel$icov[[optind]])) )
      Tcov = solve(Tinv)
      
      #Find which correlations are significant
      testmat = (Tinv)
      idx = which(testmat!=0,arr.ind=T)  
      offdiag = apply(idx,1,function(ind){
        if (ind[1]==ind[2]){
          return (0)
        }
        else{
          return (1)
        }
      })
      idxnz = idx[which(offdiag==1),]
      
      params = list(mu = pm$mu, Ticov = Tinv, Tcov = Tcov, corr = corr, 
                    yMat = yMat, numModes = numModes)
      
      save(params, file = paste0(output_path, "/allparams_t", snapshot, ".RData"))
    }
  }
}
