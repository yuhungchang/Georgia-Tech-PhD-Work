CommonGrid <- function(num.case, Number.NN = 10, output_path = working_dir){
  
  if( num.case %% split.num != 0) stop("The setting of split.num was Wrong!")
  
  # create bigmemory_bk folder
  if(!dir.exists(paste0(working_dir, "/bigmemory_bk"))) dir.create(paste0(working_dir, "/bigmemory_bk"))
  for(kk in 1:length(Response.ColumnIndex)){
    if(!dir.exists(paste0(working_dir, "/bigmemory_bk/", Response.ColumnIndex[kk]))) 
      dir.create(paste0(working_dir, "/bigmemory_bk/", Response.ColumnIndex[kk]))
  }  
  
  cat("1. Running CommonGrid .....\n")
  
  ###################################################
  # Step 1: Find the densest grid among the 30 case #
  ###################################################
  cat("  1.1. Finding the densest grid: ")
  Coordinate = vector("list", num.case)
  Coordinate.index = vector("list", num.case)
  for (i in 1:num.case){
    Coordinate[[i]] <- read_rawdata(case = Case.Index[i], tm = 2, ColumnIndex = 1:2)
  }
  save(Coordinate, file = paste0(output_path, "/Coordinate.RData"))

  Grid.vt <- sapply(Coordinate, FUN = nrow) 
  gold_case <- which.max(Grid.vt)
  cat("Case: ", gold_case, "...\n")

  CommonCoord = vector("list",num.case)
  for (i in 1:num.case){
    CommonCoord[[i]] = Coordinate[[gold_case]]
  }
  save(CommonCoord, file = paste0(output_path, "/CommonCoord.RData"))
  
  ###################################################
  # Step 2: Fit common grid to rescaled grid and    #
  #         then Rescale back to original scale     #
  #         (interpolation)                         #
  ###################################################
  
  cat("  1.2. Interpolation to common grid: Case ")
  TT <- length(Snapshot.Index)
  

 #Set up clusters
  cl <- makeCluster(detectCores(),type="SOCK")
  registerDoSNOW(cl)

  for (i in 1:num.case){
    cat(i, ",")
    
    RealCoord.i <- Coordinate[[i]]
    CommonCoord.i <- CommonCoord[[i]]
    
    ptm <- proc.time()
    coords <- t(RealCoord.i[,1:2])
    
    out.Matrix <- foreach(k = 1:nrow(CommonCoord.i), .combine = rbind) %dopar% {
      sort.out <- sort.int(colSums((coords - as.numeric(CommonCoord.i[k,]))^2), 
                           decreasing = FALSE, index.return = TRUE)
      if(sort.out$x[1] < 10^(-16)){ #if exactly the same point, weight 1 o.w 0
        min.dist.ind <- sort.out$ix[1:Number.NN]
        output <- c(min.dist.ind, c(1, rep(0, Number.NN-1)))
      }else{
        min.dist <- sort.out$x[1:Number.NN]
        weight <- 1/min.dist
        weight <- weight/sum(weight)
        min.dist.ind <- sort.out$ix[1:Number.NN]
        output <- c(min.dist.ind, weight)
      }
      return(output)
    }
    min.dist.mx <- out.Matrix[,1:Number.NN]
    weight.mx  <- out.Matrix[,(Number.NN + 1) : (2 * Number.NN)]
    
    rownames(weight.mx) <- NULL
    rownames(min.dist.mx) <- NULL
    
    FT <- filebacked.big.matrix(nrow(min.dist.mx), length(Response.ColumnIndex) * length(Snapshot.Index),
                                backingfile = paste0("FT_tmp", i,".bck"), 
                                backingpath = paste0(working_dir, "/bigmemory_bk"), 
                                descriptorfile = paste0("FT_tmp", i,".dsc"))
    j <- 0
    for(tm in Snapshot.Index){
      #Read snapshot for case i at timestep tm
      j <- j + 1
      snapT <- read_rawdata(case = Case.Index[i], tm = tm, ColumnIndex = Response.ColumnIndex)
      if(any(Response.ColumnIndex == Rescale_by_1000.ColumnIndex)){
        snapT[, Response.ColumnIndex == Rescale_by_1000.ColumnIndex] <- snapT[, Response.ColumnIndex == Rescale_by_1000.ColumnIndex]/1000
      }
     
      if(length(Response.ColumnIndex) == 1) {
        interpolate.val <- rep(0, nrow(min.dist.mx))
        for(num in 1:ncol(min.dist.mx)){
          interpolate.val <- interpolate.val + snapT[min.dist.mx[,num]] * weight.mx[,num]
        }
        FT[, j] <- interpolate.val
      }else{
        interpolate.mx <- matrix(0, ncol = ncol(snapT), nrow = nrow(min.dist.mx))
        for(kk in 1:ncol(snapT)){
          interpolate.val <- rep(0, nrow(min.dist.mx))
          for(num in 1:ncol(min.dist.mx)){
            interpolate.val <- interpolate.val + snapT[min.dist.mx[,num], kk] * weight.mx[,num]
          }
          interpolate.mx[,kk] <- interpolate.val
        }
        FT[, (1+length(Response.ColumnIndex) * (j-1)):(length(Response.ColumnIndex)*j)] <- interpolate.mx
        
      }

      
    }   
    for(kk in 1:length(Response.ColumnIndex)){
      column.index <- kk + 0:(TT-1) * length(Response.ColumnIndex)
      FT.tmp <- as.big.matrix(FT[, column.index], backingfile = paste0("FT_",i,".bck"), backingpath = paste0(working_dir, "/bigmemory_bk/", Response.ColumnIndex[kk]), descriptorfile = paste0("FT_",i,".dsc"))
    }
    
    #file.remove(c(paste0(working_dir, "/bigmemory_bk/FT_tmp.bck"), paste0(working_dir, "/bigmemory_bk/FT_tmp.dsc")))
  }
  stopCluster(cl)
  
  grid.num <- nrow(CommonCoord[[1]]) 
  
  cat("\n  1.3. Splitting and Backing up the Big Matrix .....\n")  
  for(ind in Response.ColumnIndex){
    cat("       Response.ColumnIndex =", ind, "\n")
    FT <- filebacked.big.matrix(grid.num, TT * num.case,                      
                                backingfile = "FT_All.bck",  type='double',
                                backingpath = paste0(working_dir, "/bigmemory_bk/", ind), 
                                descriptorfile = "FT_All.dsc")
    
    for(i in 1:num.case){
      cat("       Making FT_All step",i,"\n")
      FT_part.desc <- attach.big.matrix(paste0(working_dir, "/bigmemory_bk/", ind, "/FT_",i,".dsc"))
      column.index <- (1 + (i-1) * TT) : (TT * i) 
      FT[,column.index] = FT_part.desc[]
    }

    
    ##### Compute mean at each grid
    #FTmean.vector <- bmcapply(FT, 1, mean, dir = paste0(working_dir, "/bigmemory_bk/", ind))
    FTmean.vector <- apply(FT[], 1, mean)
    na.index.fg <- is.na(FTmean.vector)
    na.index <- which(na.index.fg)
    
    if(any(na.index.fg)){
      cat("       Some NA values ....")
      print(which(na.index.fg))
    }
    save(na.index, file = paste0(working_dir, "/bigmemory_bk/", ind, "/na.index.RData"))
    FTmean.vector <- FTmean.vector[!na.index.fg]
    save(FTmean.vector, file = paste0(working_dir, "/bigmemory_bk/", ind, "/FTmean.vector.RData"))
    
    FT.centered <- filebacked.big.matrix(sum(!na.index.fg), TT * num.case,                      
                                         backingfile = "FT.centered.bck", 
                                         backingpath = paste0(working_dir, "/bigmemory_bk/", ind), 
                                         descriptorfile = "FT.centered.dsc")
    for(ii in 1:num.case){
      cat("       Centering FT.centered step", ii,"\n")
      col.index <- (1 + (ii - 1) * TT) : (ii * TT)
      FT.centered[,col.index] <- FT[!na.index.fg, col.index] - FTmean.vector
    }
    
    for(ii in 1:split.num){
      col.index <- (1 + (ii - 1) * (num.case/split.num) * TT) : (ii * (num.case/split.num) * TT)
      FTtmp <- filebacked.big.matrix(sum(!na.index.fg), length(col.index),                      
                                     backingfile = paste0("FTtmp",ii,".bck"),  type='double',
                                     backingpath = paste0(working_dir, "/bigmemory_bk/", ind), 
                                     descriptorfile = paste0("FTtmp",ii,".dsc"))
      
      for(kk in 1:(num.case/split.num)){
        cat("       Centering FTtmp", ii, "step", kk,"\n")
        col.tmp.index <- col.index[(1 + (kk - 1) * TT) : (kk * TT)]
        FTtmp[,(1 + (kk - 1) * TT) : (kk * TT)] <- FT.centered[,col.tmp.index]
      }
    }
  }
  
  ## clean some backup data
  cat("  1.4. Cleaning some backup files .....\n")  
  for(ind in Response.ColumnIndex){
    file.remove(paste0(working_dir, "/bigmemory_bk/", ind, "/FT.centered.bck"))
    file.remove(paste0(working_dir, "/bigmemory_bk/", ind, "/FT.centered.dsc"))
    file.remove(paste0(working_dir, "/bigmemory_bk/", ind, "/FT_All.bck"))
    file.remove(paste0(working_dir, "/bigmemory_bk/", ind, "/FT_All.dsc"))
    file.remove(paste0(working_dir, "/bigmemory_bk/", ind, "/FT_", 1:num.case,".bck"))    
    file.remove(paste0(working_dir, "/bigmemory_bk/", ind, "/FT_", 1:num.case,".dsc"))
  }
  
  return(gold_case)
}