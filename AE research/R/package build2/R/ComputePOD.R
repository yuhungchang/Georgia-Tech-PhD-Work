ComputePOD <- function(num.case, energy.cutoff, POD_maxnum.cutoff, golden_case,
                       CommonGrid_filepath, output_path = getwd()){
  
  cat("2. Running ComputePOD .....\n")

  # read common grid
  load(paste0(CommonGrid_filepath, "/CommonCoord.RData"))
  Coordinate <- CommonCoord[[golden_case]]
  
  # create output_path folder
  if(!dir.exists(output_path)) dir.create(output_path)
  for(kk in 1:length(Response.ColumnIndex)){
    if(!dir.exists(paste0(output_path,"/", Response.ColumnIndex[kk]))) 
      dir.create(paste0(output_path,"/", Response.ColumnIndex[kk]))
  }
  
  TT <- length(Snapshot.Index)
  dimPart <- num.case/split.num * TT
  for(ind in Response.ColumnIndex){
    cat("  Response.ColumnIndex =", ind, "\n")
    cat("  2.1. Computing K matrix .....\n")
    K <- matrix(0, nrow = TT * num.case, ncol = TT * num.case)
    for (i in 1:split.num){
      for (j in i:split.num){
        #Loading libraries
        FTtmp1 <- attach.big.matrix(paste0(working_dir, "/bigmemory_bk/", ind, "/FTtmp", i,".dsc"))
        FTtmp2 <- attach.big.matrix(paste0(working_dir, "/bigmemory_bk/", ind, "/FTtmp", j,".dsc"))

        K[ ((i-1)*dimPart+1):(i*dimPart), ((j-1)*dimPart+1):(j*dimPart) ] <- bigInnerProd2(FTtmp1@address,FTtmp2@address)
        if(i != j) K[ ((j-1)*dimPart+1):(j*dimPart), ((i-1)*dimPart+1):(i*dimPart)] <- t(K[ ((i-1)*dimPart+1):(i*dimPart), ((j-1)*dimPart+1):(j*dimPart)])
      }
    }
    
    cat("  2.2. Eigendecomposition .....\n")
    if(ncol(K) < 5000){ # add it due to the unstable "eigs_sym" function for small size matrix (08/04/2016)
      eigdat <- eigen(K)
    }else {
      eigdat <- eigs_sym(K, min(POD_maxnum.cutoff, ncol(K)))
    }

    FTlist <- vector("list",split.num)
    tFTlist <- vector("list",split.num)
    for (i in 1:split.num){
      FTlist[[i]] <- attach.big.matrix(paste0(working_dir, "/bigmemory_bk/", ind, "/FTtmp",i,".dsc"))
      #if(ncol(FTlist[[i]]) > 10000){
      if(FALSE){
        tFTlist[[i]] <- big.t(FTlist[[i]], 
                              dir = paste0(working_dir, "/bigmemory_bk/", ind),
                              name = paste0("t.bigMat",i))
      }else{
        tFTlist[[i]] <- t(FTlist[[i]][])
        tFTlist[[i]] <- as.big.matrix(tFTlist[[i]], 
                                      backingfile = paste0("t.bigMat",i,".bck"), 
                                      backingpath = paste0(working_dir, "/bigmemory_bk/", ind), 
                                      descriptorfile = paste0("t.bigMat",i,".dsc"))
      }

    }

    eiglist <- vector("list",split.num)
    for (i in 1:split.num){
      eiglist[[i]] <- as.big.matrix(eigdat$vector[((i-1)*dimPart+1):(i*dimPart),], 
                                    backingfile = paste0("eiglist",i,".bck"), 
                                    backingpath = paste0(working_dir, "/bigmemory_bk/", ind), 
                                    descriptorfile = paste0("eiglist",i,".dsc"))
    }

    #POD modes
    PsiTlist <- vector("list",split.num)
    for (i in 1:split.num){
      PsiTlist[[i]] <- bigInnerProd2(tFTlist[[i]]@address,eiglist[[i]]@address)
    }
    PsiT <- matrix(0,nrow = nrow(FTlist[[1]]),ncol = ncol(PsiTlist[[1]])) # fixed bug on 08/05/2016
    for (i in 1:split.num){
      PsiT = PsiT + PsiTlist[[i]]
    }
    PsiT <- t(t(PsiT) / apply(PsiT, 2, FUN = function(x) sqrt(sum(x^2)))) # Normalize
    #POD coefficients
    PsiTbig <- as.big.matrix(PsiT, 
                             backingfile = "PsiT.bck", 
                             backingpath = paste0(working_dir, "/bigmemory_bk/", ind), 
                             descriptorfile = "PsiT.dsc")
    BT <- matrix(0, nrow = 0, ncol = ncol(PsiTlist[[1]])) # fixed bug on 08/05/2016
    for (i in 1:split.num){
      BT <- rbind( BT, bigInnerProd2(FTlist[[i]]@address,PsiTbig@address) )
    }
    
    #Cut POD by energy setting
    #cs = cumsum(eigdat$values)/sum(eigdat$values)
    cs = cumsum(eigdat$values)/sum(diag(K)) # edit on 07/20/2016 by chihli
    cs = pmax(cs - energy.cutoff,0)
    cs[cs==0] = Inf
    numModes = which.min(cs)
    cat("  2.3. Number of Modes with Energy: ", energy.cutoff, "\n")
    cat("       Response ", ind, ":", numModes, "\n")
    
    #Save POD information
    PsiTsvg <- PsiT[,1:numModes]
    BTsvg <- BT[,1:numModes]

    save(PsiTsvg, BTsvg, file = 
           paste0(output_path, "/", ind, "/PODdat.RData"))
    
    load(paste0(working_dir,"/bigmemory_bk/", ind, "/na.index.RData"))
    if(length(na.index) > 0) {
      POD_mode <- cbind(Coordinate[-na.index,], PsiTsvg)
    }else{
      POD_mode <- cbind(Coordinate, PsiTsvg)
    }
    write.table(POD_mode, file = paste0(output_path, "/", ind, "/POD_Modes.dat"), row.names = F)
    
    cat("  2.3. Cleaning some backup files .....\n")  
    for(i in 1:split.num){
      file.remove(paste0(working_dir, "/bigmemory_bk/", ind, "/eiglist",i,".bck"))
      file.remove(paste0(working_dir, "/bigmemory_bk/", ind, "/eiglist",i,".dsc"))
      file.remove(paste0(working_dir, "/bigmemory_bk/", ind, "/t.bigMat",i,".bck"))
      file.remove(paste0(working_dir, "/bigmemory_bk/", ind, "/t.bigMat",i,".dsc"))
    }
    file.remove(paste0(working_dir, "/bigmemory_bk/", ind, "/PsiT.bck"))
    file.remove(paste0(working_dir, "/bigmemory_bk/", ind, "/PsiT.dsc"))
  }
}
