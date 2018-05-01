## simulated MM48 dataset

library(KPIC)
library(MAPD)
path <- 'E:/project/pymass/python/MM48_MSS.mzML'
pics <- getPIC(LoadData(path), mztol=0.02, width=1, level=100)
peaks <- NULL
for (i in seq_along(pics$pics)) {
  pic <- pics$pics[[i]]
  pic[,1] <- pics$scantime[pic[,1]]
  this <- MAPD(pic, SNR.Th = 3, amp.Th = 100, PeakRange = 3, Filter = TRUE, winSize = 10)
  peaks <- rbind(peaks, this)
}

ground_truth <- read.csv("E:/project/pymass/python/groud_truth.csv", row.names=1)
true_pos <- lapply(1:nrow(peaks), function(s){
  rt <- peaks[s, 'rt']
  mz <- peaks[s, 'mz']
  
  wh <- which(ground_truth[,'rt_min'] < rt & ground_truth[,'rt_max'] > rt &
                ground_truth[,'mz_min']-0.01 < mz & ground_truth[,'mz_max']+0.01 > mz)
  wh[1]
})

true_pos <- unique(true_pos)

recall <- length(true_pos)/nrow(ground_truth)
precision <- length(true_pos)/nrow(peaks)
F_score <- (2*recall*precision)/(recall+precision)

c(recall=recall, precision=precision, F_score=F_score)


## Arabidopsis thaliana dataset

library(KPIC)
path <- 'E:/dataset/Arabidopsis/Sol_Seed_Leaf_0_50_50_1-A,6_01_2554.mzxml'
pics <- getPIC(LoadData(path), mztol=0.05, width=10, level=100)
viewPICs(pics)

pic <- pics$pics[[111]][539:684, ]
pic[,1] <- pics$scantime[pic[,1]]
peak1 <- MAPD(pic, SNR.Th = 5, amp.Th = 100, PeakRange = 3, Filter = FALSE, winSize = 20)
peak2 <- MAPD(pic, SNR.Th = 5, amp.Th = 100, PeakRange = 3, Filter = TRUE, winSize = 20)

p5 <- plot_ly() %>%
  layout(xaxis = list(tick0 = 0, title = 'Retention Time (s)'),
         yaxis = list(title = 'Intensity'),
         legend = list(x = 0.7, y = 0.9))
p5 <- add_lines(p5, x=pic[,1], y=pic[,2], color = I('black'), name = 'ion trace')
p5 <- add_markers(p5, x=peak1[,'rt'], y=peak1[,'intensity'], color=I('blue'), name='initial detected peaks')
p5 <- add_markers(p5, x=peak2[,'rt'], y=peak2[,'intensity'], color=I('red'), name='refined peaks')
p5

p6 <- plot_ly(showlegend= FALSE) %>%
  layout(xaxis = list(tick0 = 0, title = 'Retention Time (s)'),
         yaxis = list(title = 'm/z'))
p6 <- add_markers(p6, x=pic[,1], y=pic[,3])
p6

