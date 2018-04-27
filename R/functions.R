getNoise <- function(peaks, cwt2d, ridges){
  row_one <- row_one_del <- cwt2d[1,]
  del <- which(abs(row_one) < 10e-5)
  if (length(del)>0){
    row_one_del <- row_one[-del]
  }
  
  t <- 3*median(abs(row_one_del-median(row_one_del)))/0.67
  row_one[row_one > t] <- t
  row_one[row_one < -t] <- -t
  
  noises <- sapply(1:length(peaks),function(s){
    hf_win <- length(ridges$ridges_rows)
    win_s <- max(1, peaks[s] - hf_win)
    win_e <- min(ncol(cwt2d), peaks[s] + hf_win)
    return(as.numeric(quantile(abs(row_one[win_s:win_e]),0.9)))
  })
  return(noises)
}

localMaximum <- function (x, winSize = 5) {
  len <- length(x)
  rNum <- ceiling(len/winSize)
  
  ## Transform the vector as a matrix with column length equals winSize
  ##		and find the maximum position at each row.
  y <- matrix(c(x, rep(x[len], rNum * winSize - len)), nrow=winSize)
  y.maxInd <- apply(y, 2, which.max)
  ## Only keep the maximum value larger than the boundary values
  selInd <- which(apply(y, 2, function(x) max(x) > x[1] & max(x) > x[winSize]))
  
  ## keep the result
  localMax <- rep(0, len)
  localMax[(selInd-1) * winSize + y.maxInd[selInd]] <- 1
  
  ## Shift the vector with winSize/2 and do the same operation
  shift <- floor(winSize/2)
  rNum <- ceiling((len + shift)/winSize)	
  y <- matrix(c(rep(x[1], shift), x, rep(x[len], rNum * winSize - len - shift)), nrow=winSize)
  y.maxInd <- apply(y, 2, which.max)
  ## Only keep the maximum value larger than the boundary values
  selInd <- which(apply(y, 2, function(x) max(x) > x[1] & max(x) > x[winSize]))
  localMax[(selInd-1) * winSize + y.maxInd[selInd] - shift] <- 1
  
  ## Check whether there is some local maxima have in between distance less than winSize
  maxInd <- which(localMax > 0)
  selInd <- which(diff(maxInd) < winSize)
  if (length(selInd) > 0) {
    selMaxInd1 <- maxInd[selInd]
    selMaxInd2 <- maxInd[selInd + 1]
    temp <- x[selMaxInd1] - x[selMaxInd2]
    localMax[selMaxInd1[temp <= 0]] <- 0
    localMax[selMaxInd2[temp > 0]] <- 0
  }
  
  return(localMax)
}

peakDetection <- function(vec, scales=1:20, SNR.Th=5, amp.Th=0, ScaleRange=5){
  cwt2d <- cwtft(vec, scales=scales)$cwt2d
  ridges <- ridgesDetection(cwt2d, vec)
  if (length(ridges$ridges_rows)<1){return(NULL)}
  peaks <- peaksPosition(vec, ridges, cwt2d)
  ridge_lens <- sapply(ridges$ridges_rows, length)
  
  keep <- NULL
  for (p in unique(peaks)){
    this <- which(peaks==p)
    this <- this[which.max(ridge_lens[this])][1]
    keep <- c(keep, this)
  }
  peaks <- peaks[keep]
  ridges$ridges_rows <- ridges$ridges_rows[keep]
  ridges$ridges_cols <- ridges$ridges_cols[keep]
  
  signals <- getSignal(cwt2d, ridges, peaks)
  lens <- signals$ridge_lens + 1
  signals <- signals$signals
  peaks <- peaks+1
  noises <- getNoise(peaks, cwt2d, ridges)
  snr <- (signals+10^-5)/(noises+10^-5)
  refine <- snr>SNR.Th & lens>3 & vec[peaks]>amp.Th & scales[lens]>ScaleRange
  
  info <- cbind(peaks, scales[lens], snr)
  info <- info[refine,]
  if (length(info)==0){return(NULL)
  } else if (length(info)>3){
    info <- info[order(info[,1]),]
    peakIndex=info[,1]; peakScale=info[,2]; snr=info[,3]; signals=vec[info[,1]]
  } else {
    peakIndex=info[1]; peakScale=info[2]; snr=info[3]; signals=vec[info[1]]
  }
  return(list(peakIndex=peakIndex, peakScale=peakScale, snr=snr, signals=signals))
}

haar <- function(ms, scales=1) {
  extendNBase <-  function(x, nLevel=1, base=2, ...) {
    if (!is.matrix(x)) x <- matrix(x, ncol=1)
    nR <- nrow(x)
    if (is.null(nLevel)) {
      nR1 <- nextn(nR, base)
    } else {
      nR1 <- ceiling(nR / base^nLevel) * base^nLevel
    }
    if (nR != nR1) {
      x <- extendLength(x, addLength=nR1-nR, ...)
    }
    return(x)
  }
  extendLength <- function(x, addLength=NULL, method=c('reflection', 'open', 'circular'), direction=c('right', 'left', 'both')) {
    if (is.null(addLength)) stop('Please provide the length to be added!')
    if (!is.matrix(x)) x <- matrix(x, ncol=1)
    method <- match.arg(method)
    direction <- match.arg(direction)
    nR <- nrow(x)
    nR1 <- nR + addLength
    if (direction == 'both') {
      left <- right <- addLength
    } else if (direction == 'right') {
      left <- 0
      right <- addLength
    } else if (direction == 'left') {
      left <- addLength
      right <- 0
    }
    if (right > 0) {
      x <- switch(method,
                  reflection =rbind(x, x[nR:(2 * nR - nR1 + 1), , drop=FALSE]),
                  open = rbind(x, matrix(rep(x[nR,], addLength), ncol=ncol(x), byrow=TRUE)),
                  circular = rbind(x, x[1:(nR1 - nR),, drop=FALSE]))
    }
    if (left > 0) {
      x <- switch(method,
                  reflection =rbind(x[addLength:1, , drop=FALSE], x),
                  open = rbind(matrix(rep(x[1,], addLength), ncol=ncol(x), byrow=TRUE), x),
                  circular = rbind(x[(2 * nR - nR1 + 1):nR,, drop=FALSE], x))
    }
    if (ncol(x) == 1)  x <- as.vector(x)
    return(x)
  }
  
  psi_xval <- seq(0,1,length=1024)
  psi <- c(0,rep(1,511),rep(-1,511),0)
  oldLen <- length(ms)
  ms <- extendNBase(ms, nLevel=NULL, base=2)
  len <- length(ms)
  nbscales <- length(scales)
  wCoefs <- NULL
  
  psi_xval <- psi_xval - psi_xval[1]
  dxval <- psi_xval[2]
  xmax  <- psi_xval[length(psi_xval)]
  for (i in 1:length(scales)) {
    scale.i <- scales[i]
    f <- rep(0, len)
    j <- 1 + floor((0:(scale.i * xmax))/(scale.i * dxval))
    if (length(j) == 1)		j <- c(1, 1)
    lenWave <- length(j)
    f[1:lenWave] <- rev(psi[j]) - mean(psi[j])
    if (length(f) > len) stop(paste('scale', scale.i, 'is too large!'))
    wCoefs.i <- 1/sqrt(scale.i) * convolve(ms, f)
    wCoefs.i <- c(wCoefs.i[(len-floor(lenWave/2) + 1) : len], wCoefs.i[1:(len-floor(lenWave/2))])
    wCoefs <- cbind(wCoefs, wCoefs.i)
  }
  if (length(scales) == 1) wCoefs <- matrix(wCoefs, ncol=1)
  colnames(wCoefs) <- scales
  wCoefs <- wCoefs[1:oldLen,,drop=FALSE]
  return(wCoefs)
}

getPeakWidth <- function(x, peaks){
  wCoefs_haar <- haar(x, 1: min(length(x), max(peaks$peakScale)))
  peakIndex <- peaks$peakIndex
  peakScale <- peaks$peakScale
  
  peakBound <- lapply(seq_along(peaks$peakIndex), function(i){
    peakIndex.i <- peakIndex[i]
    peakScale.i <- min(length(x), peakScale[i])
    
    wCoefs_haar.i <- wCoefs_haar[,peakScale.i]
    wCoefs_haar.i.abs <- abs(wCoefs_haar.i)
    
    localmax <- t(localMax(t(-wCoefs_haar)))
    localmax <- as.numeric(localmax)
    localmax[peakIndex] <- 0
    tmprange <- (peakIndex.i-peakScale.i/2+1):(peakIndex.i+peakScale.i/2-1)
    localmax[tmprange[tmprange>0]] <- 0
    
    peakScale.i.3 <- 2*peakScale.i
    if(i==1){
      maxIndexL <- max(c((peakIndex.i-peakScale.i.3),1))
    }else{
      maxIndexL <- max(c((peakIndex.i-peakScale.i.3),peakIndex[i-1]))
    }
    
    if(i==length(peakIndex)){
      minIndexR <- min(c((peakIndex.i+peakScale.i.3),length(localmax)), length(x))
    } else{
      minIndexR <- min(c((peakIndex.i+peakScale.i.3),peakIndex[i+1]), length(x))
    }
    
    ignoreL <- 1:maxIndexL
    ignoreR <- minIndexR:length(localmax)
    localmax[c(ignoreL,ignoreR)] <- 0
    tmprange <- c(peakIndex.i,(peakIndex.i-(peakScale.i/2)):(peakIndex.i+(peakScale.i/2)))
    localmax[tmprange[tmprange>0]] <- 0
    bi <- which(localmax==1)
    
    biLeft <- bi[bi<peakIndex.i]
    useL <-  maxIndexL:peakIndex.i
    minIndexLeft <- useL[which(min(x[useL])==x[useL])][1]
    
    if(length(biLeft)==0){
      Lef <- minIndexLeft
    }else{
      minbaselineIndexLeft <- biLeft[which(min(x[biLeft])==x[biLeft])]
      if(minIndexLeft>=(peakIndex.i-peakScale.i/2+1)){
        Lef <- minbaselineIndexLeft
      }else{
        Lef <- max(c(minIndexLeft,minbaselineIndexLeft))
      }
    }
    
    biRight <- bi[bi>peakIndex.i]
    useR <- peakIndex.i:minIndexR
    minIndexRight <- useR[which(min(x[useR])==x[useR])][1]
    
    if(length(biRight)==0){
      Rig <- minIndexRight
    }else{
      minbaselineIndexRight <- biRight[which(min(x[biRight])==x[biRight])]
      if(minIndexRight <= (peakIndex.i+peakScale.i/2-1)){
        Rig <- minbaselineIndexRight
      }else{
        Rig <- min(c(minIndexRight,minbaselineIndexRight))
      }
    }
    cbind(Lef, Rig)
  })
  
  return(do.call(rbind, peakBound))
}

peakFit <- function(peak,pic,iter){
  # define sub-functions
  Gaussian <- function(x,position,width){
    exp(-((x-position)/(0.6005612*width))^2);
  }

  fitGaussian <- function(lambda,x,y){
    A <- matrix(0,length(x),round(length(lambda)/2))
    for(j in 1:(length(lambda)/2))
    {
      position <-lambda[2*j-1] ; width <- lambda[2*j];
      A[ ,j] <- Gaussian(x,lambda[2*j-1],lambda[2*j]);
    }
    lf=lsfit(A,y)
    PEAKHEIGHTS=abs(lf$coef[-1])
    e=y-A%*%PEAKHEIGHTS
    return(base::norm(as.matrix(e),'2'))
  }

  fitness <- function(lambda,x,y){
    -fitGaussian(lambda,x,y)
  }

  peakrange <- function(position,width,pos,wid){
    lambda <- cbind(position,width);
    NumPeaks <- nrow(lambda);
    lb <- as.vector(lambda);
    ub <- as.vector(lambda);
    for(i in 1:NumPeaks){
      ub[2*i-1] <- lambda[i,1]+pos*lambda[i,1];
      lb[2*i-1] <- lambda[i,1]-pos*lambda[i,1];
      ub[2*i] <- lambda[i,2]+wid*lambda[i,2];
      lb[2*i] <- lambda[i,2]-wid*lambda[i,2];
    }
    output <- list(ub=ub,lb=lb)
    return(output)
  }

  PEAKHEIGHTS <- function(fitresults,xa,y){
    # The positionwidth is the results of the initial estimation or GAs for position or width.
    NumPeaks <- length(fitresults)/2;
    lambda <- matrix(fitresults,nrow=NumPeaks,byrow=TRUE)
    A <- matrix(0,length(xa),round(length(lambda)/2))
    for(j in 1:(length(lambda)/2))
    {
      A[ ,j] <- Gaussian(xa,lambda[j,1],lambda[j,2]);
    }
    lf=lsfit(A,y)
    PEAKHEIGHTS=abs(lf$coef[-1])
    return(matrix(PEAKHEIGHTS))
  }

  # end of sub-functions

  hh <- 0.5*(pic[c(peak$peakIndex),2])
  ranges <- whichAsIRanges(pic[,2]>min(hh))

  if (length(peak$peakIndex)>0){
    mids <- round(diff(peak$peakIndex)/2) + peak$peakIndex[1:(length(peak$peakIndex)-1)]
    starts <- c(min(start(ranges)),mids)
    ends <- c(mids,max(end(ranges)))

    widths <- (ends-starts)*1.7
    heights <- peak$signals
    positions <- peak$peakIndex

    lambda <- as.vector(t(cbind(positions,widths)));

    lb <- peakrange(positions,widths,0.05,0.95)$lb
    ub <- peakrange(positions,widths,0.05,0.95)$ub

    GA <- ga(type = "real-valued", fitness = fitness,
             x = 1:nrow(pic), y = pic[,2], min = c(lb), max = c(ub),
             popSize=300, maxiter=iter)

    fitresults <- matrix(GA@solution, nrow=1);
    heights <- matrix(PEAKHEIGHTS(fitresults,1:nrow(pic),pic[,2]));
    positions <- sapply(1:(length(fitresults)/2),function(s){fitresults[2*s-1]})
    widths <- sapply(1:(length(fitresults)/2),function(s){fitresults[2*s]})

    fitpics <- lapply(1:length(positions),function(s){heights[s]*Gaussian(1:nrow(pic),positions[s],widths[s])})

  } else {
    stop()
  }
  peak$peakIndex <- round(positions)
  peak$width <- round(widths)
  peak$signals <- heights
  return(list(peaks=peak, fitpics=fitpics))
}

split.intensity <- function(pic, info) {
  pic[,2] <- smooth.spline(pic[,2])$y
  int <- pic[,2]
  rts <- pic[,1]
  pks <- info[,'rt']
  if (length(pks)>1){
    res <- sapply(1:(length(pks)-1), function(s){
      wh <- which(rts>pks[s] & rts<pks[s+1])
      lmin <- wh[which.min(pic[wh, 2])]
      if (pic[lmin,2] < 0.5*info[s,'intensity'] && pic[lmin,2] < 0.5*info[s+1,'intensity']) {
        lmin
      } else {
        NULL
      }
    })
    res <- unlist(res)
  } else {
    res <- NULL
  }
  return(res)
}

split.mz <- function(pic, info){
  rts <- pic[,1]
  mzs <- pic[,3]
  winSize <- round(min(diff(info[,'rt'])/mean(diff(rts))))
  CoV <- rollapply(mzs, winSize, function(x) sqrt(var(x))/mean(x))
  lmax <- which(localMaximum(CoV, winSize)==1)
  qua <- quantile(CoV, 0.9)
  lmax <- lmax[CoV[lmax]>qua]
  if (length(lmax)>0){
    res <- round(0.5*winSize + lmax)
  } else {
    res <- NULL
  }
  
  return(res)
}

getArea <- function(rt, intensity, lb, rb){
  integration <- function(x,yf){
    n <- length(x)
    integral <- 0.5*sum((x[2:n] - x[1:(n-1)]) * (yf[2:n] + yf[1:(n-1)]))
    return(integral)
  }
  Area <- integration(rt[lb:rb], intensity[lb:rb])
  return(Area)
}

MAPD <- function(pic, scales = 1:20, SNR.Th = 5, amp.Th = 0, PeakRange = 20, pval.Th = 0.005, FitIter = 0){
  mzs <- pic[,3]
  int <- pic[,2]
  rts <- pic[,1]
  ScaleRange <- round(PeakRange/mean(diff(rts))/4)
  peaks <- peakDetection(int, scales=scales, SNR.Th=SNR.Th, amp.Th=amp.Th, ScaleRange=ScaleRange)
  
  if (length(peaks$peakIndex) > 0){
    bounds <- getPeakWidth(int, peaks)
    ind <- peaks$peakIndex
    lef <- bounds[,'Lef']
    rig <- bounds[,'Rig']
    info <- cbind(rts[ind], int[ind], mzs[ind], rts[lef], rts[rig])
    colnames(info) <- c('rt', 'intensity', 'mz', 'left', 'right')
  } else {
    info <- NULL
    return(info)
  }
  
  if (!is.na(pval.Th)) {
    sp1 <- split.intensity(pic, info)
    sp2 <- split.mz(pic, info)
    sp <- sort(c(sp1, sp2))
    sec <- findInterval(info[,'rt'], sp)
    
    keep <- sapply(unique(sec), function(s){
      wh <- which(sec==s)
      if (length(wh)==1){
        wh
      } else {
        wh[which.max(peaks$signals[wh])]
      }
    })
    
    peaks <- lapply(peaks, function(p) p[keep])
    bounds <- getPeakWidth(int, peaks)
    lef <- bounds[,'Lef']
    rig <- bounds[,'Rig']
  }

  if (FitIter<=0) {
    bounds <- getPeakWidth(int, peaks)
    ind <- peaks$peakIndex
    areas <- sapply(seq_along(ind) , function(s){
      getArea(pic[,1], pic[,2], lef[s], rig[s])
    })
    info <- cbind(rts[ind], int[ind], areas, mzs[ind], rts[lef], rts[rig])
    colnames(info) <- c('rt', 'intensity', 'area', 'mz', 'left', 'right')
    return(info)
  } else {
    pf <- peakFit(peaks, pic, FitIter)
    ind <- pf$peaks$peakIndex
    areas <- sapply(seq_along(ind) , function(s){
      getArea(pic[,1], pf$fitpics[[s]], 1, nrow(pic))
    })
    info <- cbind(rts[ind], int[ind], areas, mzs[ind], rts[lef], rts[rig])
    colnames(info) <- c('rt', 'intensity', 'area', 'mz', 'left', 'right')
    return(list(info=info, fits=pf$fitpics))
  }
}