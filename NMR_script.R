library(Rnmr1D)
data_dir <- system.file("extra", package = "Rnmr1D")
RAWDIR <- file.path(data_dir, "CD_BBI_16P02")
CMDFILE <- file.path(data_dir, "NP_macro_cmd.txt")
SAMPLEFILE <- file.path(data_dir, "Samples.txt")

#get the samples
samples <- read.table(SAMPLEFILE, sep="\t", header=T,stringsAsFactors=FALSE)

#set up the commands
CMDTEXT <- readLines(CMDFILE)
CMDTEXT[grep("^#$", CMDTEXT, invert=TRUE)]

#do the processing
out <- Rnmr1D::doProcessing(RAWDIR, cmdfile=CMDFILE, samplefile=SAMPLEFILE, ncpu=2)
out$infos
out_df<-data.frame(out$infos)
sample_df<-data.frame(out$samples)
out_df$treament<-sample_df$Treatment
#Stacked plot
plotSpecMat(out$specMat, ppm_lim=c(0.5,5), K=0.33)


#overlayed plot
plotSpecMat(out$specMat, ppm_lim=c(0.5,5), K=0, pY=0.1)

#further processing
specMat.new <- Rnmr1D::doProcCmd(out, 
                                 c( "bucket aibin 10.2 10.5 0.2 3 0", "9.5 4.9", "4.8 0.5", "EOL" ), ncpu=2, debug=TRUE)

out$specMat <- specMat.new
#normalization
outMat <- Rnmr1D::getBucketsDataset(out, norm_meth='CSN')
outMat[, 1:10]

#bucket_clustering
options(warn=-1)
clustcor <- Rnmr1D::getClusters(outMat, method='corr', cval=0, dC=0.003, ncpu=2)
# hierarchical  bucket clustering
H_cluscor<-Rnmr1D::getClusters(outMat, method='hca', vcutusr=0.12)

#clustering results .
#clusters List of the ppm value corresponding to each cluster. the length of the list equal to number of clusters
#clustertab the associations matrix that gives for each cluster (column 2) the corresponding buckets (column 1)
clustcor$clustertab[1:20, ]
Clust_df<-data.frame(clustcor$clustertab)

#comparing the results of the two clustering methods
g1 <- ggplotCriterion(clustcor)
g1

g2 <- ggplotCriterion(H_cluscor)
g2

layout(matrix(1:2, 1, 2,byrow = TRUE))

hist(simplify2array(lapply(clustcor$clusters, length)), 
     breaks=20, main="CORR", xlab="size", col="darkcyan")
mtext("clusters size distribution", side = 3)

hist(simplify2array(lapply(H_cluscor$clusters, length)), 
     breaks=20, main="HCA", xlab="size", col="darkcyan")
mtext("clusters size distribution", side = 3)

g3 <- ggplotClusters(outMat,clustcor)
#ggplotPlotly(g3, width=820, height=400)
g3

g4 <- ggplotClusters(outMat,H_cluscor)
#ggplotPlotly(g4, width=820, height=400)
g4


#PCA
pca <- prcomp(outMat,retx=TRUE,scale=T, rank=2)
sd <- pca$sdev
eigenvalues <- sd^2
evnorm <- (100*eigenvalues/sum(eigenvalues))[1:10]

#Plot PCA Loadings (based on the HCA method)

g6 <- ggplotLoadings(pca$rotation, 1, 2, associations=H_cluscor$clustertab, EV=evnorm, main=sprintf("Loadings - Crit=%s",H_cluscor$vcrit), gcontour="ellipse" )
#ggplotPlotly(g6, width=820, height=650)
g6

#exporting results


out_mat_df<-data.frame(outMat)
out_mat_df$treatment<-sample_df$Treatment

write.csv(out_mat_df,"C:/Users/ccape/Downloads/NMR_data/Spec_data.csv")

write.csv(Clust_df,"C:/Users/ccape/Downloads/NMR_data/Cluster_data.csv")


#orginal code came from https://cran.r-project.org/web/packages/Rnmr1D/vignettes/Rnmr1D.html  
#I made changes to add in reporting and exporting useful information.


