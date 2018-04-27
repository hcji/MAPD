\name{MAPD}
\alias{MAPD}

\title{
Peak detection combining wavelet transform and mass spectrometry
}
\description{
MAPD is a peak detection method, which combines wavelet transform and mass spectrometry information. The main purpose is to identify whether overlapped peaks exist and whether they are true positive. If the overlapped peaks are true peaks, MAPD can also resolve them and calculate their peak heights, widths and areas respectively.
}
\usage{
MAPD(pic, scales = 1:20, SNR.Th = 5, amp.Th = 0, PeakRange = 5, pval.Th = 0.05, FitIter = 0)
}

\arguments{
  \item{pic}{
A pure ion chromatogram with rt, intensity and mz values.
}
  \item{scales}{
Scales of CWT
}
  \item{SNR.Th}{
Minimum signal to noise ratio.
}
  \item{amp.Th}{
Minimum intensity of a peak.
}
  \item{PeakRange}{
The minimum scale range of the peak.
}
  \item{pval.Th}{
P-value of the t-test.
}
  \item{FitIter}{
Number of iterations of GA, if 0, no fitness.
}
}

\value{
Result of peak detection after removing false positive.
}

\examples{
pic <- read.csv(system.file('example_pic.csv', package = 'MAPD'), row.names=1)
peaks <- MAPD(pic)
}