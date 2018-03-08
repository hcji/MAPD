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
  refine <- snr>SNR.Th & lens>3 & vec[peaks]>amp.Th & scales>ScaleRange

  info <- cbind(peaks, scales, snr)
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

WhittakerSmooth <- function(y,lambda){
  M <- length(y)
  E <- sparseMatrix(i=1:M,j=1:M,x=1)
  D <- Matrix::diff(E)
  C <- chol(E+lambda*Matrix::t(D)%*%D)
  z <- solve(C,solve(t(C),y))
  return(as.numeric(z))
}

integration <- function(x,yf){
  n <- length(x)
  integral <- 0.5*sum((x[2:n] - x[1:(n-1)]) * (yf[2:n] + yf[1:(n-1)]))
  return(integral)
}

plot.resolve <- function(pic, res){
  library(plotly)
  rts <- pic[,1]
  raw.vec <- pic[,2]
  fit.pics <- do.call(rbind, res$fitpics)
  sum.vec <- colSums(fit.pics)

  p <- plot_ly(x = rts, y = raw.vec, type = 'scatter', name = 'raw') %>%
    layout(xaxis = list(title = 'Retention Time (s)'),
           yaxis = list (title = 'Intensity'))

  for (i in 1:nrow(fit.pics)) {
    name <- paste('peak ', i)
    p <- add_trace(p, x = rts, y = fit.pics[i,], line = list(dash = 'dash'), mode='line', name = name)
  }
  p <- add_trace(p, x = rts, y = sum.vec, mode='line' , name = 'fitted')
  fitmodel <- colSums(do.call(rbind, res$fitpics))
  fiterror <- sqrt(sum((fitmodel-pic[,2])^2))/sum(pic[,2])*100
  R.square <- sum((mean(pic[,2])-fitmodel)^2)/sum((pic[,2]-mean(pic[,2]))^2)
  show(p)
  cat('fiterror: ', fiterror, '%', '\n')
  cat('R-square: ', R.square, '\n')
}

peak.fit <- function(peak,pic,iter){
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

WMPD <- function(pic, SNR.Th, amp.Th, pval, iter, min_width){
  library(IRanges)
  library(Matrix)
  library(GA)

  vec <- pic[,2]
  rts <- pic[,1]
  mzs <- pic[,3]
  ScaleRange <- round(min_width/mean(diff(rts))/4)
  peaks <- peak_detection(vec, SNR.Th, amp.Th, ScaleRange)

  if (length(peaks$peakIndex) > 1){
    # vec <- WhittakerSmooth(vec, 2)
    helf.height <- vec[peaks$peakIndex]*0.5
    split.point <- sapply(1:(length(peaks$peakIndex)-1), function(s){
      peaks$peakIndex[s]+which.min(vec[peaks$peakIndex[s]:peaks$peakIndex[s+1]])
    })
    split.point <- unique(c(1, split.point, length(vec)))

    peak.ranges <- lapply(1:length(peaks$peakIndex),function(s){
      split.point[s]+which(vec[split.point[s]:split.point[s+1]]>helf.height[s])-1
    })

    peak.mzs <- lapply(peak.ranges, function(peak.range){
      mzs[peak.range]
    })

    pvals <- sapply(1:(length(peak.mzs)-1),function(s){
      t.test(peak.mzs[[s]], peak.mzs[[s+1]])$p.value
    })

    TP <- rep(FALSE, length(peaks$peakIndex))
    for (i in 1:length(pvals)){
      if (pvals[i] < pval || vec[split.point[i+1]] < min(helf.height[i:(i+1)])){
        if (peaks$signals[i] > peaks$signals[i+1]){
          TP[i+1] <- TRUE
        } else {
          TP[i] <- TRUE
        }
      }
    }
    peaks$TF <- TP

    new.peaks <- list()
    new.peaks$peakIndex <- peaks$peakIndex[TP]
    new.peaks$peakScale <- peaks$peakScale[TP]
    new.peaks$snr <- peaks$snr[TP]
    new.peaks$signals <- peaks$signals[TP]
  } else {
    peaks$TF <- TRUE
    new.peaks <- peaks
  }

  res <- peak.fit(new.peaks, pic, iter)
  positions <- rts[new.peaks$peakIndex]
  areas <- sapply(res$fitpics, function(fp){
    integration(rts, fp)
  })

  heights <- sapply(res$fitpics, function(fp){
    max(fp)
  })
  res$pic <- pic
  res$peaks <- peaks
  res$new.peaks$positions <- positions
  res$new.peaks$areas <- areas
  res$new.peaks$heights <- heights
  return(res)
}