library(KPIC)
path <- 'E:/dataset/Arabidopsis/Sol_Seed_Leaf_0_0_100_1-A,2_01_2431.mzxml'
pics <- getPIC(LoadData(path), level=500)
viewPICs(pics)
