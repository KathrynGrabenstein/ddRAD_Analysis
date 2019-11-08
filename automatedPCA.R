#!/usr/bin/env Rscript

#a file that is accessed by PCA_r to conduct principle component analyses

args = commandArgs()

dataset <- args[6]
outDir <- args[7]
pops <- args[8]
name <- args[9]
sample <- args[10]

library(devtools)
#install_github("zhengxwen/gdsfmt")
#install_github("zhengxwen/SNPRelate")
library(gdsfmt)
library(SNPRelate)

#input
vcf.fn <- paste(dataset)
snpgdsVCF2GDS(vcf.fn, paste0(outDir,"/",name,".gds"), method="biallelic.only")

#summary
sink(paste0(outDir,"/",name,"_summary.txt"), append=FALSE, split=FALSE)
snpgdsSummary(paste0(outDir,"/",name,".gds"))
sink()

#open the GDS file
genofile <- snpgdsOpen(paste0(outDir,"/",name,".gds"))

#get population information
pop_code <- scan(paste(pops), what=character())
table(pop_code)

#Run pca
pca <- snpgdsPCA(genofile, num.thread=2, autosome.only = FALSE)

#to calculate percent of variation that is accounted for by top PCA components:
pc.percent <- pca$varprop*100
sink(paste0(outDir,"/",name,"_percentvariation.txt"), append=FALSE, split=FALSE)
head(round(pc.percent, 2))
sink()

#get sample id
#sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
sample.id <- read.delim(sample, header = FALSE)


#make data frame
tab <- data.frame(sample.id = sample.id$V1,
pop = factor(pop_code)[match(pca$sample.id, sample.id)],
EV1 = pca$eigenvect[,1],    # the first eigenvector
EV2 = pca$eigenvect[,2],    # the second eigenvector
stringsAsFactors = FALSE)

#save files:
write.table(tab, file = paste0(outDir,"/",name,".txt"), sep = "\t")

tab <- read.table(paste0(outDir,"/",name,".txt"), header = TRUE, sep = "\t")


#draw it:
pdf(file = paste0(outDir,"/",name,"_pca_populations_notsubsetted.pdf"), useDingbats=FALSE)

plot(tab$EV2, tab$EV1, col=as.integer(tab$pop), xlab="eigenvector 2", ylab="eigenvector 1")
legend("topright", legend=levels(tab$pop), pch="o", col=1:nlevels(tab$pop))
text(tab$EV2, tab$EV1,, labels=rownames(sample.id),data=tab, cex=0.9, font=2)

dev.off()

pdf(file = paste0(outDir,"/",name,"_pca_individuals_notsubsetted.pdf"), useDingbats=FALSE)

plot(tab$EV2, tab$EV1, col=as.integer(tab$sample.id), xlab="eigenvector 2", ylab="eigenvector 1")
legend("topright", legend=levels(tab$sample.id), pch="o", col=1:nlevels(tab$sample.id))

dev.off()
