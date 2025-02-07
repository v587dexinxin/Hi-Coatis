library('idr')
mu <- 2.6
sigma <- 1.3
rho <- 0.8
p <- 0.7


###K562_SF-----hg19------------
fname = paste('H:/work/niulongjian/K562/peaks/common_peak/K562_SF_common_peaks.txt' , sep = "")
outname = paste('H:/work/niulongjian/K562/peaks/IDR/K562_SF_IDR_Selected_0.05.narrowPeak' , sep = "")
data_0 = read.table(fname)
data <- data_0[,c('V8','V18')]
data_M = as.matrix(data)
out = est.IDR(data_M, mu, sigma, rho, p, eps=0.001,max.ite = 20)
data.selected = select.IDR(data_0,out$IDR,0.05)
write.table(data.selected$x , outname , row.names = FALSE,quote = FALSE,col.names = FALSE , sep = "\t")

###K562_SF-----hg38------------
fname = paste('H:/work/niulongjian/K562/hg38/K562_SF/peaks/common_peaks/K562_SF_common_peaks.txt' , sep = "")
outname = paste('H:/work/niulongjian/K562/hg38/K562_SF/peaks/IDR/K562_SF_IDR_Selected_0.05.narrowPeak' , sep = "")
data_0 = read.table(fname)
data <- data_0[,c('V8','V18')]
data_M = as.matrix(data)
out = est.IDR(data_M, mu, sigma, rho, p, eps=0.001,max.ite = 20)
data.selected = select.IDR(data_0,out$IDR,0.05)
write.table(data.selected$x , outname , row.names = FALSE,quote = FALSE,col.names = FALSE , sep = "\t")

