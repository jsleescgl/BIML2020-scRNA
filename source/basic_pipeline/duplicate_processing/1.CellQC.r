# BIML 2020 Tutorial
# Part 1. Cell QC
# 10x genomics or DropSeq-based data QC

dirExp1='Q:\\BIML2020\\data\\2.Fundamental_pipeline\\exp1\\raw10x'
dirExp2='Q:\\BIML2020\\data\\2.Fundamental_pipeline\\exp2\\raw10x'

setwd(dirExp2)

library(scRNAseq)
library(SingleCellExperiment)
library(DropletUtils)
library(scater)

# Data load from cellranger output directory contains barcode/features/mtx files
sce <- read10xCounts(dirExp2)
rownames(sce) = uniquifyFeatureNames(rowData(sce)$ID, rowData(sce)$Symbol)
colnames(sce) = sce$Barcode

my.counts=counts(sce)

# plotting of UMI count barcode rank 
# barcode rank by total UMI counts
br.out <- barcodeRanks(my.counts)

# Making a plot.
plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
o <- order(br.out$rank)
lines(br.out$rank[o], br.out$fitted[o], col="red")

abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), 
       legend=c("knee", "inflection"))

# Filtering empty droplets using DropletUtils package (Genome Biol, 2019 - Aaron Lun)
set.seed(100)
e.out <- emptyDrops(my.counts)
e.out

is.cell <- e.out$FDR < 0.01
sum(is.cell, na.rm=TRUE) # result of this cmd = after filter out empty droplets

# barcode rank by total UMI counts
br.out_neRemoved <- barcodeRanks(counts(sce_ne))

# Making a plot.
plot(br.out_neRemoved$rank, br.out_neRemoved$total, log="xy", xlab="Rank", ylab="Total")
o <- order(br.out_neRemoved$rank)
lines(br.out_neRemoved$rank[o], br.out_neRemoved$fitted[o], col="red")

abline(h=metadata(br.out_neRemoved)$knee, col="dodgerblue", lty=2)
abline(h=metadata(br.out_neRemoved)$inflection, col="forestgreen", lty=2)
legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), 
       legend=c("knee", "inflection"))

# remove empty droplet
sce_ne=sce[,which(is.cell)]
sce_ne=sce_ne[rowSums(counts(sce_ne))>0,]


### Remove low quality cells based on total UMI counts and Mito-Genome mapped reads ratio per each cells
### Quality control on the cells
library(scater)
per.cell <- perCellQCMetrics(sce_ne, subsets=list(Mito=grep("mt-", rownames(sce_ne))))
summary(per.cell$sum)
colData(sce_ne) <- cbind(colData(sce_ne), per.cell)
sce_ne$log10Sum=log10(sce_ne$sum+1)

plotColData(sce_ne, x = "log10Sum", y="subsets_Mito_percent")

sce_ne <- runColDataPCA(sce_ne, variables=list("sum", "detected", "percent_top_50","percent_top_100",
                                               "percent_top_200", "percent_top_500", "subsets_Mito_sum",
                                               "subsets_Mito_detected","subsets_Mito_percent","log10Sum"))


# Visualize each QC metric
library(dplyr)
library(RColorBrewer)
qc.metric='subsets_Mito_percent'
fig <- tibble(x = reducedDim(sce_ne)[,1],
              y = reducedDim(sce_ne)[,2],
              qc.metric = sce_ne@colData[, qc.metric]) %>%
  ggplot(aes(x=x, y=y, colour=qc.metric)) +
  geom_point(size=0.1) +
  scale_colour_gradientn(colours = rev(brewer.pal(11, "RdBu"))) +
  # scale_colour_gradient(low = "grey", high = "red") +
  # scale_colour_gradientn(colours = rev(c("#300000", "red","#eeeeee")),
  #                        breaks=c(0,max(inferCNV_20q)),
  #                        labels=c(0,4)) +
  xlab("Component 1") +
  ylab("Component 2") +
  theme_bw() +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        # legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())

print(fig)

filtering <- c(rep(0,dim(sce_ne)[2]))
filtering[which(colData(sce_ne)[,'sum'] > 100 & colData(sce_ne)[,'subsets_Mito_percent'] < 10)] <- 1
sce_ne@colData$QC=filtering

fig <- tibble(x = reducedDim(sce_ne)[,1],
              y = reducedDim(sce_ne)[,2],
              QC = sce_ne@colData$QC) %>%
  ggplot(aes(x=x, y=y, colour=factor(sce_ne@colData$QC))) +
  geom_point(size=0.1) +
  # scale_colour_gradientn(colours = rev(brewer.pal(11, "RdBu"))) +
  # scale_colour_gradient(low = "grey", high = "red") +
  # scale_colour_gradientn(colours = rev(c("#300000", "red","#eeeeee")),
  #                        breaks=c(0,max(inferCNV_20q)),
  #                        labels=c(0,4)) +
  xlab("Component 1") +
  ylab("Component 2") +
  theme_bw() +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())

print(fig)

sce_filt=sce_ne[,sce_ne$QC==1]

# Gene filtering: remove genes which have all zero value across whole cells
sce_filt=sce_filt[rowSums(counts(sce_filt))>0,]

saveRDS(sce_filt, file='Q:\\BIML2020\\data\\mouse_lung_SCLC\\sce_filt_2.rds')
