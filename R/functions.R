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

peak_detection <- function(vec, min_snr, level=0){
  cwt2d <- cwtft(vec)
  sca <- cwt2d$scales
  cwt2d <- cwt2d$cwt2d
  ridges <- ridgesDetection(cwt2d, vec)
  if (length(ridges$ridges_rows)<1){return(NULL)}
  peaks <- peaksPosition(vec, ridges, cwt2d)
  signals <- getSignal(cwt2d, ridges, peaks)
  lens <- signals$ridge_lens
  lens[lens<0] <- 0
  scales <- sca[1+lens]
  lens <- signals$ridge_lens
  signals <- signals$signals
  peaks <- peaks+1
  noises <- getNoise(peaks, cwt2d, ridges)
  snr <- (signals+10^-5)/(noises+10^-5)
  refine <- snr>min_snr & lens>3 & vec[peaks]>level

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

  if (length(peak$peakIndex)>1){
    mids <- round(diff(peak$peakIndex)/2) + peak$peakIndex[1:(length(peak$peakIndex)-1)]
    starts <- c(min(start(ranges)),mids)
    ends <- c(mids,max(end(ranges)))

    widths <- (ends-starts)*1.7
    heights <- peak$signals
    positions <- peak$peakIndex

    lambda <- as.vector(t(cbind(positions,widths)));

    lb <- peakrange(positions,widths,0.05,0.9)$lb
    ub <- peakrange(positions,widths,0.05,0.9)$ub

    GA <- ga(type = "real-valued", fitness = fitness,
             x = 1:nrow(pic), y = pic[,2], min = c(lb), max = c(ub),
             popSize=300, maxiter=iter)

    fitresults <- matrix(GA@solution, nrow=1);
    heights <- matrix(PEAKHEIGHTS(fitresults,1:nrow(pic),pic[,2]));
    positions <- sapply(1:(length(fitresults)/2),function(s){fitresults[2*s-1]})
    widths <- sapply(1:(length(fitresults)/2),function(s){fitresults[2*s]})

    fitpics <- lapply(1:length(positions),function(s){heights[s]*Gaussian(1:nrow(pic),positions[s],widths[s])})

  } else {
    starts <- min(start(ranges))
    ends <- max(end(ranges))

    widths <- (ends-starts+1)*1.7
    heights <- max(pic[,2])
    positions <- which.max(pic[,2])

    fit <- nls(pic[,2]~heights*Gaussian(1:nrow(pic),positions,widths),start=list(heights=heights,positions=positions,widths=widths),
               control=list(maxiter=500,tol=1e-3,minFactor=1/512,printEval=F,warnOnly=T))

    fitpics <- list(predict(fit,1:nrow(pic)))
    heights <- max(fitpics[[1]])
    positions <- which.max(fitpics[[1]])
    widths <- summary(fit)$parameters['widths','Estimate']
  }
  peak$peakIndex <- round(positions)
  peak$width <- round(widths)
  peak$signals <- heights
  return(list(peaks=peak, fitpics=fitpics))
}

WMPD <- function(pic, min_snr, level, pval, iter){
  library(IRanges)
  library(Matrix)
  library(GA)

  vec <- pic[,2]
  rts <- pic[,1]
  mzs <- pic[,3]
  peaks <- peak_detection(vec, min_snr, level)

  vec.smooth <- WhittakerSmooth(vec, 2)
  helf.height <- vec.smooth[peaks$peakIndex]*0.5
  split.point <- sapply(1:(length(peaks$peakIndex)-1), function(s){
    peaks$peakIndex[s]+which.min(vec.smooth[peaks$peakIndex[s]:peaks$peakIndex[s+1]])
  })
  split.point <- unique(c(1, split.point, length(vec)))

  peak.ranges <- lapply(1:length(peaks$peakIndex),function(s){
    split.point[s]+which(vec.smooth[split.point[s]:split.point[s+1]]>helf.height[s])-1
  })

  peak.mzs <- lapply(peak.ranges, function(peak.range){
    mzs[peak.range]
  })

  pvals <- sapply(1:(length(peak.mzs)-1),function(s){
    t.test(peak.mzs[[s]], peak.mzs[[s+1]])$p.value
  })

  TP <- rep(TRUE, length(peaks$peakIndex))
  for (i in 1:length(pvals)){
    if (pvals[i] > pval){
      if (peaks$signals[i] > peaks$signals[i+1]){
        TP[i+1] <- FALSE
      } else {
        TP[i] <- FALSE
      }
    }
  }
  peaks$TF <- TP

  new.peaks <- list()
  new.peaks$peakIndex <- peaks$peakIndex[TP]
  new.peaks$peakScale <- peaks$peakScale[TP]
  new.peaks$snr <- peaks$snr[TP]
  new.peaks$signals <- peaks$signals[TP]

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

runWMPD <- function(){
  library(shiny)
  library(plotly)
  ui <- fluidPage(
    titlePanel("Wavelet-based Peak Detection Assisted by Mass"),
    sidebarLayout(
      sidebarPanel(
        fileInput("file1", "Choose CSV File",
                  multiple = TRUE,
                  accept = c("text/csv",
                             "text/comma-separated-values,text/plain",
                             ".csv")),
        numericInput("snr",
                     "Minimum SNR of peaks:",
                     value = 4),
        numericInput("threshold",
                     "Threshold of peak height:",
                     value = 0),
        numericInput("pval",
                     "p-value of m/z difference of different peaks:",
                     value = 0.05,
                     step = 0.01),
        actionButton("goButton", "Go!")

      ),

      mainPanel(
        plotlyOutput("Plot1"),
        plotlyOutput("Plot2"),
        plotlyOutput("Plot3"),
        tableOutput('Table1')
      )
    )
  )

  server <- function(input, output) {
    pic <- reactive({
      req(input$file1)
      read.csv(input$file1$datapath)
    })

    output$Plot1 <- renderPlotly({
      pic <- pic()
      plot_ly(x=pic[,1], y=pic[,2], type = 'scatter', mode = 'lines', color = 'blue', showlegend= FALSE) %>%
        layout(xaxis = list(tick0 = 0, title = 'Retention Time (s)'),
               yaxis = list(title = 'Intensity'))
    })

    output$Plot2 <- renderPlotly({
      pic <- pic()
      plot_ly(x=pic[,1], y=pic[,3], color = pic[,2], type = 'scatter')%>%
        layout(xaxis = list(tick0 = 0, title = 'Retention Time (s)'),
               yaxis = list(title = 'M/Z'))
    })

    Peaks <- eventReactive(input$goButton, {
      withProgress(message = 'Detecting peaks', value = 0.5, {
        WMPD(pic(), input$snr, input$threshold, input$pval, 200)
      })
    })

    sum.vec <- reactive({
      fit.pics <- do.call(rbind, Peaks()$fitpics)
      colSums(fit.pics)
    })

    output$Plot3 <- renderPlotly({
      pic <- pic()
      withProgress(message = 'Creating plot', value = 0.1, {
        p <- plot_ly(x=pic[,1], y=pic[,2], type = 'scatter', mode = 'lines', name = 'raw', showlegend= TRUE) %>%
          layout(xaxis = list(tick0 = 0, title = 'Retention Time (s)'),
                 yaxis = list(title = 'Intensity'))
        incProgress(0.1)
        for (f in 1:length(Peaks()$fitpics)) {
          p <- add_trace(p, x = pic[,1], y = Peaks()$fitpics[[f]], line = list(dash = 'dash'), mode='line', name = paste('peak',f))
          incProgress(0.1)
        }
        setProgress(0.9)
        p <- add_trace(p, x = pic[,1], y = sum.vec(), mode='line', name = 'fitted profile')
        p
      })
    })

    output$Table1 <- renderTable({
      res.old <- do.call(cbind, Peaks()$peaks)
      res.old[,1] <- pic()[res.old[,1],1]
      colnames(res.old) <- c('original.rt', 'original.scale', 'original.snr', 'original.height', 'T/F')

      res.new <- matrix(NA, nrow(res.old), 3)
      res.new[res.old[,5]==1,] <- do.call(cbind, Peaks()$new.peaks)
      res.new[,1] <- pic()[res.new[,1],1]
      colnames(res.new) <- c('fitted.rt', 'fitted.area', 'fitted.height')

      as.data.frame(cbind(res.new, res.old))
    })
  }
  shinyApp(ui = ui, server = server)
}



