# MAPD
  Mass Spectrometry Assisted Peak Detection Algorithm of Pure Ion Chromatograms from LC-MS
  
## Installation
	
	install.packages(c("devtools", "Rcpp", "RcppArmadillo", "Matrix", "GA", "zoo"))
	library(devtools);  
	install_github("hcji/MAPD")
	
## Usage

	pic <- read.csv(system.file('example_pic.csv', package = 'MAPD'), row.names=1)
	peaks <- MAPD(pic)
	
## Contact
For any questions, please contact: ji.hongchao@foxmail.com
	