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
  lens <- signals$ridge_lens
  signals <- signals$signals
  peaks <- peaks+1
  noises <- getNoise(peaks, cwt2d, ridges)
  snr <- (signals+10^-5)/(noises+10^-5)
  refine <- snr>SNR.Th & lens>3 & vec[peaks]>amp.Th & scales[lens]>ScaleRange

  info <- cbind(peaks, scales[lens], snr)
  info <- info[refine,]
  info <- unique(info)
  if (length(info)==0){return(NULL)
  } else if (length(info)>3){
    info <- info[order(info[,1]),]
    peakIndex=info[,1]; peakScale=info[,2]; snr=info[,3]; signals=vec[info[,1]]
  } else {
    peakIndex=info[1]; peakScale=info[2]; snr=info[3]; signals=vec[info[1]]
  }
  return(list(peakIndex=peakIndex, peakScale=peakScale, snr=snr, signals=signals))
}

integration <- function(x,yf){
  n <- length(x)
  integral <- 0.5*sum((x[2:n] - x[1:(n-1)]) * (yf[2:n] + yf[1:(n-1)]))
  return(integral)
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

isIsolated <- function(index, peaks, int) {
  vec <- smooth.spline(int)$y
  loc <- peaks$peakIndex[index]
  npeaks <- length(peaks$peakIndex)
  peakHeight <- vec[loc]
  adjPeak <- {
    adj <- c(index-1, index+1)
    adj <- adj[adj>0&adj<=npeaks]
    peaks$peakIndex[adj]
  }
  adjPeakHeight <- vec[adjPeak]
  if (peakHeight > max(adjPeakHeight)) {
    return(TRUE)
  } else {
    lmin <- sapply(adjPeak, function(adj){
      min(vec[loc:adj])
    })
    if (max(lmin) < 0.5*peakHeight){
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
}

isSignificant <- function(peaks, int, mzs, pval.Th){
  vec <- smooth.spline(int)$y
  helfHeight <- vec[peaks$peakIndex]*0.5
  splitPoint <- sapply(1:(length(peaks$peakIndex)-1), function(s){
    peaks$peakIndex[s]+which.min(vec[peaks$peakIndex[s]:peaks$peakIndex[s+1]])
  })
  splitPoint <- unique(c(1, splitPoint, length(vec)))
  peakRanges <- lapply(1:length(peaks$peakIndex),function(s){
    splitPoint[s]+which(vec[splitPoint[s]:splitPoint[s+1]]>helfHeight[s])-1
  })
  peakMzs <- lapply(peakRanges, function(peakRange){
    mzs[peakRange]
  })
  res <- sapply(seq_along(peaks$peakIndex), function(i){
    adj <- {
      adj <- c(i-1, i+1)
      adj <- adj[adj>0&adj<=npeaks]
    }
    pvals <- sapply(adj, function(s){
      t.test(peakMzs[[s]], peakMzs[[i]])$p.value
    })
    max(pvals) < pval.Th
  })
  
  return(res)
}

WMPD <- function(pic, scales = 1:20, SNR.Th = 5, amp.Th = 0, PeakRange = 5, pval.Th = 0.05, FitIter = 0){
  int <- pic[,2]
  rts <- pic[,1]
  mzs <- pic[,3]
  ScaleRange <- round(PeakRange/mean(diff(rts))/4)
  peaks <- peakDetection(vec, scales=scales, SNR.Th=SNR.Th, amp.Th=amp.Th, ScaleRange=ScaleRange)

  if (length(peaks$peakIndex) > 1){
    rule1 <- sapply(seq_along(peaks$peakIndex), function(i) isIsolated(i, peaks, int))
    rule2 <- isSignificant(peaks, int, mzs, pval.Th)
    truePeak <- rule1|rule2
  } else {
    truePeak <- TRUE
  }
  
  return(res)
}