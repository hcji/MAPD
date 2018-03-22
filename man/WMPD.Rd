\name{MAPD}
\alias{MAPD}

\title{
Peak detection combining wavelet transform and mass spectrometry
}
\description{
MAPD is a peak detection method, which combines wavelet transform and mass spectrometry information. The main purpose is to identify whether overlapped peaks exist and whether they are true positive. If the overlapped peaks are true peaks, MAPD can also resolve them and calculate their peak heights, widths and areas respectively.
}
\usage{
MAPD(pic, min_snr, level, pval, iter)
}

\arguments{
  \item{pic}{
A pure ion chromatogram with rt, intensity and mz values.
}
  \item{min_snr}{
Minimum signal to noise ratio.
}
  \item{level}{
Minimum intensity of a peak.
}
  \item{pval}{
P-value of the t-test.
}
  \item{iter}{
Number of iterations of GA.
}
}

\value{
\item{pic}{The input pic.}
\item{peaks}{Result of initial peak detection.}
\item{new.peaks}{Result of peak detection after removing false positive.}
\item{fitpic}{Fit curve of real peaks.}
}

\examples{
pic <- read.csv(system.file('example_pic.csv', package = 'MAPD'))
peaks <- MAPD(pic, 4, 0, 0.05, 100)
}
