#
# This script explores the SNP data and removes batch effects
# Cleaned SNP and haplotypes files are written out for further analysis
# Several population genetics analyses are carried out
#
library(vcfR)
library(adegenet)
library(ggplot2)
library(scales)
library(gridExtra)
library(ggplotify)
library(cowplot)
#library(vegan)
#library(hierfstat)
library(FinePop2)
library(reshape2)
library(ggtree)
library(UpSetR)
library(ComplexUpset)

#
# 1) Load SNP data and metadata
#

indir=commandArgs(trailingOnly=TRUE)[1]
indir="stacks_m3M7_n0/singletonsHD+VS_D2_POOL7_p3r0.5/"
#indir="stacks_m3M7_n0/singletonsHD+VS_D2_p5r0.5/"
outdir=file.path(indir, "R_output")
dir.create(outdir)

# metadata
# M7 and M9 are mixed
# M3, M4, M8 are pure
meta <- read.table("geography_Mariana.txt", header=T, stringsAsFactors=F)
meta$Label[rev(order(meta$Lat))] <- c("A", "B", "C", "D", "E")
pmap <- read.table("popmaps/final_group_Mariana_splitpool.txt", stringsAsFactors = F,
                   col.names = c("ID", "Pop2"))
pmap$Pop <- substr(pmap$Pop2, 1, 2)
pmap <- merge(pmap, meta, by="Pop", sort=F)
apply(meta, 1, function(x) pmap$Pop2 <<- gsub(x["Pop"], x["Label"], pmap$Pop2))
pools <- read.table("sample_pools.txt", header=T)
pmap <- merge(pmap, pools, by="ID")

saveRDS(pmap, file=file.path(outdir, "popmap_metadata.RDS"))

# colours for each population
#pcols <- c(A="darkblue", B="steelblue3", C="#FF9900", D="red2", E="chartreuse4")
pcols <- c(A="#FF9900", B="dodgerblue", C="black", D="chartreuse3", E="red2")
#pcols.split <- c(A.N="darkblue", B.N="steelblue3", C.S="#FF9900", C.N="orange", D.S="red2", D.N="red2", E.S="chartreuse4")
#pcols.split <- c(A.P5="darkblue", B.P5="steelblue3", C.P7="#FF9900", C.P8="#FF9900", D.P6="red2", D.P7="red2", E.P7="chartreuse4")

# load SNP data
infile=file.path(indir, "populations.snps.vcf")
vcf <- read.vcfR(infile)
vcf.gen <- vcfR2genlight(vcf)
vcf.gen <- vcf.gen[na.exclude(match(pmap$ID, vcf.gen@ind.names)),]

# impute missing data (mean)
vcf.imp <- tab(vcf.gen, freq=T, NA.method="mean")

# locus presence/absence 
vcf.na <- t(tab(vcf.gen, freq=F, NA.method="asis"))
vcf.na[!is.na(vcf.na)] <- 1
vcf.na[is.na(vcf.na)] <- 0
vcf.na <- t(rowsum(vcf.na, gsub(":.*", "", rownames(vcf.na))))
vcf.na[vcf.na>1] <- 1

# purge popmap if individuals are missing
pmap <- subset(pmap, ID %in% vcf.gen@ind.names)

# load HWE data
source("HWE.R")
hwe <- load_sumstats(paste(indir, "populations.sumstats.tsv", sep="/"))
hwe$SNP.ID <- paste0(hwe$Locus, ":", hwe$Pos)

#################################
# check for batch effect
# pool7 contains all of E and half of C and D
# this pool causes lots of problems
# use DAPC to identify SNPs that split pool7 from other pools in C and D
#################################

# skip batch effect when input file is only pool7 or only other pools
# DAPC by pool
source("PCA_DAPC.R")

findbatch <- function(dd, jac=F){
	# DAPC by pool on all SNPs
	if(jac==T){
		pc.all <- PCA(dist.binary(dd, method=1), "a)", "Pool")
	} else {
		pc.all <- PCA(dd, "a)", "Pool")	
	}
	pc.opt <- which(cumsum(pc.all[[2]]$eig/sum(pc.all[[2]]$eig))*100>=80)[1]
	da.pool <- dapc(dd, pmap$Pool, scale=F, n.da=2, center=T, n.pca=pc.opt)
	gg.pool <- plot.DAPC(data.frame(da.pool$ind.coord, pmap), "b)", "Pool")
	
	# DAPC on sites C and D only, because these are the only ones with two pools
	pmap.CD <- subset(pmap, Label %in% c("C", "D"))
	dd.CD <- dd[as.character(pmap.CD$ID),]
	pc <- PCA(dd.CD)
	pc.opt <- which(cumsum(pc[[2]]$eig/sum(pc[[2]]$eig))*100>=80)[1]
	da.CD <- dapc(dd.CD, pmap.CD$Pool, scale=F, n.da=2, center=T, n.pca=pc.opt)
	da.CD.df <- data.frame(da.CD$ind.coord, pmap.CD)
	gg.CD <- plot.DAPC(da.CD.df, "c)", "Pool")
	
	# identify SNPs driving batch effect (6/8 vs 7 but NOT 6 vs 8!)
	#batch <- !with(subset(da.CD.df, !Pool==7), coef(summary(lm(Pool~LD1+LD2)))[-1,4]<=0.05)
	#batch=c(1,2)
	#batch.SNPs <- unique(unlist(lapply(c("LD1", "LD2"), function(x) names(which(da.CD$var.contr[,x]>=quantile(da.CD$var.contr[,x], 0.80))))[batch]))
	
	# better: C and D independently
	bi <- sapply(c("C", "D"), function(site) {
		pmap.site <- subset(pmap, Label %in% site)
		dd.site <- dd[as.character(pmap.site$ID),]
		pc <- PCA(dd.site)
		pc.opt <- which(cumsum(pc[[2]]$eig/sum(pc[[2]]$eig))*100>=80)[1]
		da.site <- dapc(dd.site, pmap.site$Pool, scale=F, n.da=2, center=T, n.pca=pc.opt)
		da.site.df <- data.frame(da.site$ind.coord, pmap.site)
	
		gg <- ggplot(da.site.df, aes(x=LD1, group=factor(Pool), fill=Label, colour=Label)) +
		geom_vline(xintercept=0, colour="lightgray", linetype="dotted") +
		geom_density(alpha=0.5) +
		scale_fill_manual(name="Site", values=pcols) +
		scale_colour_manual(name="Site", values=pcols) +
		theme_bw() +
		theme(legend.position="none",
			panel.grid = element_blank(),
			axis.title = element_text(size=14),
		        axis.text = element_text(size=12),
		        plot.tag.position = c(0.008,0.99),
		        plot.tag = element_text(size=16, face="bold"))

		bs <- names(which(da.site$var.contr[,1]>=quantile(da.site$var.contr[,1], 0.90)))
		return(list(gg, bs))
	})
	#gg <- arrangeGrob(gg.pool, gg.CD, arrangeGrob(bi[[1]], bi[[3]], ncol=1), hwe.bar(subset(hwe, SNP.ID %in% colnames(dd))), ncol=4)
	gg <- arrangeGrob(pc.all[[1]], gg.pool, gg.CD, arrangeGrob(bi[[1]], bi[[3]], ncol=1), ncol=4, widths=c(0.28,0.28,0.28,0.16))
	bs <- unique(unlist(c(bi[[2]], bi[[4]])))
	return(list(SNPs=bs, Loci=unique(gsub(":.*", "", bs)), ggPlot=gg))
}

if(length(unique(pmap$Pool))<4) {	
	# skip
	vcf.imp.clean <- vcf.imp
} else {	
	# identify batch effect
	print("Identifying batch effect...")
	# 1) presence/absence
	batch.PA <- findbatch(vcf.na, jac=T)
	# 2) genotypes after cleaning presence/absence batch
	#batch.SNPs <- findbatch(vcf.imp[, ! colnames(vcf.imp) %in% batch.PA$SNPs])	# remove SNPs
	batch.SNPs <- findbatch(vcf.imp[, ! gsub(":.*", "", colnames(vcf.imp)) %in% batch.PA$Loci])	# remove loci

	# remove from dataset
	vcf.imp.clean <- vcf.imp[, ! gsub(":.*", "", colnames(vcf.imp)) %in% batch.PA$Loci]
	#vcf.imp.clean <- vcf.imp.clean[, ! colnames(vcf.imp.clean) %in% batch.SNPs$SNPs]
	vcf.imp.clean <- vcf.imp.clean[, ! gsub(":.*", "", colnames(vcf.imp.clean)) %in% batch.SNPs$Loci]
	print(dim(vcf.imp.clean))
	writeLines(colnames(vcf.imp.clean), con=file.path(outdir, "SNPnames.clean.txt"))

	# re-run DAPC on pools to verify that batch effect is gone
	tmp <- findbatch(vcf.imp.clean)

	# plot supplemental figures
	pdf(file.path(outdir, "batcheffect.pdf"), width=18, height=18)
		grid.arrange(batch.PA$ggPlot, batch.SNPs$ggPlot, tmp$ggPlot, ncol=1)
	dev.off()
	# overlap between HWE and batch effect
	ups <- data.frame(SNP=colnames(vcf.imp), Batch.SNPs=0, Batch.loci=0, FIS=0)
	ups$Batch.SNPs[ups$SNP %in% batch.SNPs$SNPs] <- 1
	ups$Batch.loci[gsub(":.*", "", ups$SNP) %in% batch.PA$Loci] <- 1
	ups$FIS[ups$SNP %in% subset(hwe, Label %in% c("C", "D") & Sig=="Sig")$SNP.ID] <- 1

	pdf(file.path(outdir, "batcheffect_intersections.pdf"), width=5, height=3)
		print(UpSetR::upset(ups))
	dev.off()
}

###################################
# Outlier analyses
###################################

# DAPC
pc.clean <- PCA(vcf.imp.clean, "a)")
pc.opt <- which(cumsum(pc.clean[[2]]$eig/sum(pc.clean[[2]]$eig))*100>=80)[1]
#a.opt <- optim.a.score(dapc(vcf.imp.clean, pmap$Label, scale=F, n.da=2, center=T, n.pca=100), n.sim=25, plot=F, smart=F)$best
da.clean <- dapc(vcf.imp.clean, pmap$Label, scale=F, n.da=4, center=T, n.pca=pc.opt)
da.clean.df <- data.frame(da.clean$ind.coord, pmap)

source("PCA_DAPC.R")
gg.da.clean <- plot.DAPC(da.clean.df, "b)")
gg.da.clean2 <- plot.DAPC(da.clean.df, "", xa="LD3", ya="LD2")

# Figure for manuscript
pdf(file.path(outdir, "PCA+DAPC.pdf"), width=15, height=5)
	grid.arrange(pc.clean[[1]], gg.da.clean, gg.da.clean2, ncol=3, widths=c(0.345, 0.345, 0.31))
dev.off()

# correlation plot of LD1,2,3 with depth/latitude

gg.cor <- mapply(function(ld, var) {
	ct <- cor.test(da.clean.df[,ld], da.clean.df[,var], data=da.clean.df, method="spearman")
	ct2 <- paste0("rho==", round(ct$estimate, 3), "*\";\"~italic(P)==", round(ct$p.value, 3))
	df2 <- aggregate(da.clean.df, list(da.clean.df$Label), median)
	g <- ggplot(da.clean.df, aes_string(x=var, y=ld)) +
	geom_hline(yintercept=0, colour="lightgray", linetype="dotted") +
	geom_point(aes(fill=Label), shape=21, size=8, alpha=0.3, color="white") +
	#geom_boxplot(aes(colour=Label), position="identity", width=ifelse(var=="Depth", 100, 0.01)) +
	#geom_text(aes(label=Label), colour="white", size=5) +
	geom_smooth(method="lm", colour="black") +
	annotate("label", colour="white", fill="grey40", y = Inf, x = ifelse(ct$estimate<0 || var=="Depth", Inf, -Inf), label = ct2, parse=T, hjust=ifelse(ct$estimate<0 || var=="Depth", 1.1, -0.1), vjust=1.5, size=5) +
	#geom_point(data=df2, aes(fill=Group.1), shape=21, size=8, alpha=1, color="white") +	
	#geom_text(data=df2, aes(colour=Group.1, label=Group.1), colour="white") +	
	geom_label(data=df2, aes(fill=Group.1, label=Group.1), colour="white", size=5) +
	scale_fill_manual(name="Site", values=pcols) +
	scale_colour_manual(name="Site", values=pcols) +
	xlab(ifelse(var=="Lat", "Latitude", "Depth")) +
	theme_bw() +
	theme(legend.position="none",
		panel.grid = element_blank(),
		axis.title.y = element_text(size=14),
	        axis.text = element_text(size=12),
		axis.title.x = element_text(size=14),
	        plot.tag.position = c(0.008,0.99),
	        plot.tag = element_text(size=16, face="bold"),
		plot.margin = margin(5.5, 7.5, 5.5, 5.5, "pt"))
	#if(ld=="LD3") g <- g + theme(axis.title.x = element_text(size=12))
	return(g)
}, rep(c("LD1", "LD2", "LD3"), times=2), rep(c("Lat", "Depth"), each=3), SIMPLIFY=F)

pdf(file.path(outdir, "DAPC_correlations.pdf"), width=14, height=9)
	do.call(grid.arrange, c(gg.cor, ncol=3))
dev.off()

# old correlation plot
#ct1 <- cor.test(~LD1+Lat, data=da.clean.df, method="spearman")
#ct2 <- cor.test(~LD2+Depth, data=da.clean.df, method="spearman")
#ct <- c(LD1=paste0("rho==", round(ct1$estimate, 3), "*\";\"~italic(P)==", round(ct1$p.value, 3)),
#LD2=paste0("rho==", round(ct2$estimate, 3), "*\";\"~italic(P)==", round(ct2$p.value, 3)))

#da.clean.df$Lat.jitter <- da.clean.df$Lat + runif(nrow(da.clean.df), min=-0.02, max=0.02)
#da.clean.df$Depth.jitter <- da.clean.df$Depth + runif(nrow(da.clean.df), min=-100, max=100)

#gg.depth <- do.call(arrangeGrob, c(lapply(c("LD1", "LD2"), function(ld) {
#	ggplot(da.clean.df, aes_string(x=ifelse(ld=="LD1", "Lat", "Depth"), y=ld)) +
#	geom_hline(yintercept=0, colour="lightgray", linetype="dotted") +
##	geom_point(aes(fill=Label), shape=21, size=8, alpha=1, color="white") +
#	geom_text(aes(label=Label), colour="white", size=5) +
#	geom_smooth(method="lm", colour="black") +
#	annotate("text", x = Inf, y = Inf, label = ct[ld], parse=T, vjust=2, hjust=1.1, size=5) +
#	scale_fill_manual(name="Site", values=pcols) +
#	scale_colour_manual(name="Site", values=pcols) +
#	xlab(ifelse(ld=="LD1", "Latitude", "Depth")) +
#	theme_bw() +
#	theme(legend.position="none",
#		panel.grid = element_blank(),
#		axis.title = element_text(size=14),
#	        axis.text = element_text(size=12),
#	        plot.tag.position = c(0.008,0.99),
#	        plot.tag = element_text(size=16, face="bold"))
#
#}), ncol=1))

# rough summary figure
#pdf(file.path(outdir, "DAPC.pdf"), width=16, height=8)
#	grid.arrange(pc.clean[[1]], gg.da.clean, gg.depth, ncol=3)
#dev.off()

# sweep over numbers of PCs and do correlations with latitude/depth
dapc.sweep1 <- DAPCsweep(vcf.imp.clean, "DAPC_sweep_clean")

print("Outlier tests...")

# DAPC outliers
m <- matrix(1:nrow(da.clean$var.contr), nrow=1)
colnames(m) <- paste0("dummy", 1:ncol(m))
sz <- snpzip(m, da.clean, plot=F, method="ward.D2")
sz <- lapply(sz$FS, function(x) x$`List of selected alleles`)[1:3]
# convert to snp and locus names
szconvert <- function(snps, fname){
	sz.snps <- gsub(":", "_", colnames(vcf.imp.clean)[snps])
	sz.loci <- sort(unique(gsub("_.*", "", sz.snps)))
	writeLines(sz.loci, con=file.path(outdir, paste0(fname, "_locusnames.txt")))
	writeLines(sz.snps, con=file.path(outdir, paste0(fname, "_snpnames.txt")))
	return(sz.loci)
}
sz.loci1 <- szconvert(sz[[1]], "DAPC_snpzip_axis1")
sz.loci2 <- szconvert(sz[[2]], "DAPC_snpzip_axis2")
sz.loci3 <- szconvert(sz[[3]], "DAPC_snpzip_axis3")
dapc.outlier.SNPs <- unique(unlist(sz))
dapc.outlier.loci <- szconvert(dapc.outlier.SNPs, "DAPC_snpzip_threeaxes")
# 2 axes only
dapc.outlier.SNPs <- unique(unlist(sz[1:2]))
dapc.outlier.loci <- szconvert(dapc.outlier.SNPs, "DAPC_snpzip_twoaxes")

# test plot
pc.tmp <- PCA(vcf.imp.clean[,dapc.outlier.SNPs], "a)")
pc.opt <- which(cumsum(pc.tmp [[2]]$eig/sum(pc.tmp [[2]]$eig))*100>=80)[1]
#a.opt <- optim.a.score(dapc(vcf.imp.clean, pmap$Label, scale=F, n.da=2, center=T, n.pca=100), n.sim=25, plot=F, smart=F)$best
da.tmp <- dapc(vcf.imp.clean[,dapc.outlier.SNPs], pmap$Label, scale=F, n.da=2, center=T, n.pca=pc.opt)
gg.tmp <- plot.DAPC(data.frame(da.tmp$ind.coord, pmap), "a)")
pdf(file.path(outdir, "DAPC_outliers.pdf"), width=12, height=8)
	grid.arrange(pc.tmp[[1]], gg.tmp, ncol=2)
dev.off()


# OUTFLANK outlier analysis
# based on Fst between populations
#
library(OutFLANK)

# extract allele codes from VCF
# make genotype codes (0, 1, 2)
# missing data: 9
vcf.clean <- vcf[vcf@fix[,"ID"] %in% colnames(vcf.imp.clean),]
alleles <- as.numeric(unlist(lapply(strsplit(vcf.clean@gt[,-1], split = c("[|///:]")), function(x) x[1:2])))
odd = seq(1, length(alleles), by=2)
gtypes <- alleles[odd] + alleles[odd+1]
gtypes <- matrix(gtypes, ncol=ncol(vcf.clean@gt[,-1]))
SNPdata <- t(gtypes)
SNPdata[is.na(SNPdata)] <- 9

# run OUTFLANK
FstDataFrame <- MakeDiploidFSTMat(SNPdata, vcf.clean@fix[,"ID"], pmap$Pop[match(colnames(vcf.clean@gt[,-1]), pmap$ID)])
outflank <- OutFLANK(FstDataFrame, NumberOfSamples=5, qthreshold=0.1, Hmin = 0.1)
outflank.sig <- subset(outflank$results, OutlierFlag==TRUE)
writeLines(as.character(outflank.sig$LocusName), con=file.path(outdir, "SNPnames.OUTFLANK.txt"))
#vcf.imp.outfl <- vcf.imp[,colnames(vcf.imp) %in% outflank.sig$LocusName]

pdf(file.path(outdir, "OUTFLANK.pdf"), width=8, height=8)
	OutFLANKResultsPlotter(outflank, withOutliers = TRUE,
		NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =
		FALSE, RightZoomFraction = 0.05, titletext = NULL)
dev.off()

# PCAdapt outlier analysis
# based on PCA between individuals
library(pcadapt)
#pcad.in <- read.pcadapt(infile, type="vcf", type.out="matrix")
pcad.in <- read.pcadapt(vcf.imp.clean, type="lfmm")
#x <- pcadapt(pcad.in, K = 20)
#plot(x, option = "screeplot")
#dev.off()
tre <- try(x <- pcadapt(pcad.in, K = 6))
if("try-error" %in% class(tre)) x <- pcadapt(pcad.in, K = 2)
pcad.outliers <- which(p.adjust(x$pvalues,method="BH") <= 0.1)
length(pcad.outliers)
writeLines(colnames(vcf.imp.clean)[pcad.outliers], con=file.path(outdir, "SNPnames.PCAdapt.txt"))

# upset plot to show outlier intersections (and HWE in all populations)
ups <- data.frame(SNP=colnames(vcf.imp.clean), PCAdapt=0, OUTFLANK=0, DAPC=0, HWE=0)
ups$PCAdapt[ups$SNP %in% colnames(vcf.imp.clean)[pcad.outliers]] <- 1
ups$OUTFLANK[ups$SNP %in% as.character(outflank.sig$LocusName)] <- 1
ups$DAPC[ups$SNP %in% colnames(vcf.imp.clean)[dapc.outlier.SNPs]] <- 1
ups$HWE[ups$SNP %in% subset(hwe, Sig=="Sig")$SNP.ID] <- 1
pdf(file.path(outdir, "outliers_intersections.pdf"), width=5, height=3)
	UpSetR::upset(ups, mb.ratio = c(0.6, 0.4))
dev.off()

###
# summarise outliers from all approaches
all.outlier.loci <- unique(c(unique(gsub(":.*", "", colnames(vcf.imp.clean)[pcad.outliers])),
			unique(gsub(":.*", "", as.character(outflank.sig$LocusName))),
			dapc.outlier.loci))
all.outlier.SNPs <- unique(c(colnames(vcf.imp.clean)[pcad.outliers],
			as.character(outflank.sig$LocusName),
			colnames(vcf.imp.clean)[dapc.outlier.SNPs]))
vcf.imp.outl <- vcf.imp.clean[,colnames(vcf.imp.clean) %in% all.outlier.SNPs]
vcf.imp.nonoutl <- vcf.imp.clean[,! colnames(vcf.imp.clean) %in% all.outlier.SNPs]
writeLines(all.outlier.SNPs, con=file.path(outdir, "SNPnames.outliers.txt"))

#### FST dendrograms
library(hierfstat)
library(FinePop2)
source("read.GENEPOP.fixed.R")
fst_dendro <- function(dd) {
	# 1) pairwise FST 
	dd <- vcf.gen[,vcf.gen@loc.names %in% colnames(dd)]
	dd <- dd[,glSum(dd)>0]
	dd <- tab(dd, freq=F, NA.method="asis")

	# hierfstat
	dd.hf <- dd
	dd.hf[dd.hf==0] <- 11
	dd.hf[dd.hf==1] <- 12
	dd.hf[dd.hf==2] <- 22
	dd.hf <- data.frame(Pop=as.numeric(factor(pmap$Label)), dd.hf)
	#fstat.b <- basic.stats(dd)
	#fst <- genet.dist(dd.hf, diploid=TRUE, method="Nei87")
	#fst <- as.matrix(fst)
	#colnames(fst) <- sort(unique(pmap$Label))
	#rownames(fst) <- sort(unique(pmap$Label))
	
	# 2) population-specific FST
	# convert to GenePop for FinePop2
	dd[dd==0] <- "0101"
	dd[dd==1] <- "0102"
	dd[dd==2] <- "0202"
	dd[is.na(dd)] <- "0000"
	writeLines(c("GENEPOP conversion", paste(colnames(dd), collapse=",")), con=file.path(outdir, "tmp.convert.gen"))
	for(p in sort(unique(pmap$Label))){
		d <- dd[rownames(dd) %in% subset(pmap, Label==p)$ID,]
		d <- cbind(paste0(rownames(d), ","), d)
		write("POP",file=file.path(outdir, "tmp.convert.gen"),append=TRUE)
		write.table(d, file=file.path(outdir, "tmp.convert.gen"), append=TRUE, row.names=F, col.names=F, quote=F, sep="\t")
	}
	# read GenePop and calculate
	genepop <- read.GENEPOP.fixed(file.path(outdir, "tmp.convert.gen"))
	genepop$pop_names <- sort(unique(pmap$Label))
	popfst <- pop_specificFST(genepop, cov=T)

	# pairwise FST 	with finepop2
	fst <- pop_pairwiseFST(genepop)
	colnames(fst) <- sort(unique(pmap$Label))
	rownames(fst) <- sort(unique(pmap$Label))
	write(fst, file=file.path(outdir, "fst.matrix.txt"))

	# GLS test
	cov <- apply(popfst$cov, 2, function(x) x/sqrt(diag(popfst$cov)))
	cor <- apply(cov, 1, function(x) x/sqrt(diag(popfst$cov)))
	sink(file.path(outdir, "fst.GLS.txt"), append=T)
	print(GLS(fst~Depth, data=data.frame(fst=popfst$fst[,1], meta[match(rownames(popfst$fst), meta$Label),c("Depth", "Lat")]), omega=cor))
	print(GLS(fst~Lat, data=data.frame(fst=popfst$fst[,1], meta[match(rownames(popfst$fst), meta$Label),c("Depth", "Lat")]), omega=cor))
	sink()

	# FST dendrogram
	fstd <- as.dist(fst)
	fstd <- hclust(fstd, method="ward.D2")
	#fstd <- ape::nj(fst)
	# root at pop with lowest FST
	#fstd <- ape::root(fstd, outgroup=genepop$pop_names[which.min(popfst$fstd[,1])[1]], resolve.root=T)
	return(list(FST=fst, Dendro=fstd, PopFST=data.frame(Label=rownames(popfst$fst),popfst$fst)))
}
write("Clean SNPs:", file=file.path(outdir, "fst.GLS.txt"))
fst.clean <- fst_dendro(vcf.imp.clean)
write("Outlier SNPs:", file=file.path(outdir, "fst.GLS.txt"), append=T)
fst.outl <- fst_dendro(vcf.imp.outl)
write("Non-outlier SNPs:", file=file.path(outdir, "fst.GLS.txt"), append=T)
fst.nonoutl <- fst_dendro(vcf.imp.nonoutl)

# new FST plot
fst.df <- data.frame(fst.outl$PopFST, fst.nonoutl$PopFST)
gg.fst <- ggplot(fst.df, 
		aes(x=FST, y=FST.1, colour=Label)) +
		geom_abline(intercept=0,slope=1, colour="lightgray", linetype="dotted") +
		geom_segment(aes(x = FST-(2*SE), y=FST.1, xend = FST+(2*SE), yend = FST.1), size=1) +
		geom_segment(aes(y = FST.1-(2*SE.1), x=FST, yend = FST.1+(2*SE.1), xend = FST), size=1) +
		geom_point(aes(fill=Label), shape=21, size=8, alpha=1, color="white") +
		geom_text(aes(label=Label), colour="white", size=5) +
		scale_fill_manual(values=pcols) +
		scale_colour_manual(values=pcols) +
		xlab(expression(Outliers~local~italic(F)[ST])) + 
		ylab(expression(Non-outliers~local~italic(F)[ST])) + 
		labs(tag="b)") +
		theme_bw() +
		theme(legend.position="none",
			panel.grid = element_blank(),
			axis.title = element_text(size=12),
		        axis.text = element_text(size=10),
		        plot.caption = element_text(size=14, face="bold", hjust = 0.5),
			plot.tag.position = c(0.008,0.99),
		        plot.tag = element_text(size=14, face="bold"))

# PCA of FST matrix
fstpca <- rbind(data.frame(Pop=c("A", "B", "C", "D", "E"), Loci="non-outliers", prcomp(as.dist(fst.nonoutl$FST))$x[,1:2]), 
		data.frame(Pop=c("A", "B", "C", "D", "E"), Loci="outliers", prcomp(as.dist(fst.outl$FST))$x[,1:2]))
fstpca <- rbind(data.frame(Pop=c("C", "D", "E"), Loci="non-outliers", prcomp(as.dist(fst.nonoutl$FST))$x[,1:2]), 
		data.frame(Pop=c("C", "D", "E"), Loci="outliers", prcomp(as.dist(fst.outl$FST))$x[,1:2]))

#fstpca <- rbind(data.frame(Pop=c("A", "B", "C", "D", "E"), Loci="non-outliers", cmdscale(as.dist(fst.nonoutl$FST))), 
#		data.frame(Pop=c("A", "B", "C", "D", "E"), Loci="outliers", cmdscale(as.dist(fst.outl$FST))))
fstpca$Loci <- relevel(fstpca$Loci, "outliers")
gg.fstpca <- ggplot(fstpca,	aes(x=PC1, y=PC2, colour=Pop)) +
		geom_segment(aes(xend=0, yend=0)) +
		geom_vline(xintercept=0, colour="lightgray", linetype="dotted") +
		geom_hline(yintercept=0, colour="lightgray", linetype="dotted") +
		geom_point(aes(fill=Pop), shape=21, size=8, alpha=1, color="white") +
		geom_text(aes(label=Pop), colour="white", size=5) +
		scale_fill_manual(values=pcols) +
		scale_colour_manual(values=pcols) +
		facet_wrap(~Loci) +
		labs(tag="a)") +
		theme_bw() +
		theme(legend.position="none",
			panel.grid = element_blank(),
			axis.title = element_text(size=12),
		        axis.text = element_text(size=10),
		        strip.text.x = element_text(size=14, face="bold"),
		        plot.tag.position = c(0.008,0.99),
		        plot.tag = element_text(size=14, face="bold"))

pdf(file.path(outdir, "FST_PCA2.pdf"), width=12, height=4)
	grid.arrange(gg.fstpca, gg.fst, ncol=2, widths=c(0.62, 0.38))
dev.off()

# Mantel/correlation test
mantel(as.dist(fst.outl$FST), dist(meta[order(meta$Label), c("Lat")]), permutations=9999)
mantel(as.dist(fst.outl$FST), dist(meta[order(meta$Label), c("Depth")]), permutations=9999)
mantel(as.dist(fst.nonoutl$FST), dist(meta[order(meta$Label), c("Lat")]), permutations=9999)
mantel(as.dist(fst.nonoutl$FST), dist(meta[order(meta$Label), c("Depth")]), permutations=9999)

mantel.partial(as.dist(fst.outl$FST), dist(meta[order(meta$Label), c("Depth")]), dist(meta[order(meta$Label), c("Lat", "Long")]), permutations=9999)
mantel.partial(as.dist(fst.nonoutl$FST), dist(meta[order(meta$Label), c("Depth")]), dist(meta[order(meta$Label), c("Lat", "Long")]), permutations=9999)


for(a in c("outliers", "non-outliers")){
	for(b in c("PC1", "PC2")){
		for(d in c("Lat", "Depth")){
			print(c(a, b, d))
			print(paste(cor.test(subset(fstpca, Loci==a)[,b], meta[order(meta$Label), d], method="spearman")[c("estimate","p.value")]))
		}
	}
}



# old pairwise FST matrix
o <- hclust(as.dist(fst.nonoutl$FST))$order
pfst.df <- melt(fst.nonoutl$FST)
pfst.df$Var1 <- factor(pfst.df$Var1, levels=c("A", "B", "C", "D", "E")[o])
pfst.df$Var2 <- factor(pfst.df$Var2, levels=c("A", "B", "C", "D", "E")[o])
ggplot(pfst.df, aes(x=Var1, y=Var2, fill=value)) + 
	geom_tile() +  
	theme_classic()
dev.off()

# old FST plot
fst_plot <- function(fst, rv=F){
	# dendrogram
	tr <- ape::as.phylo(fst$Dendro)
	gg.tr <- ggtree(tr, branch.length="none") %<+% data.frame(Taxon=as.character(pmap$Label), pmap) + 
		geom_tippoint(aes(fill=Label), size=8, shape=21, colour="white") +
		geom_tiplab(aes(label=Label), angle=0, size=5, colour="white", hjust=0.6) +
		scale_fill_manual(values=pcols) +
		scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
		theme(legend.position="none",
			plot.tag.position = c(0.005,0.99),
			plot.tag = element_text(size=16, face="bold"))
	if(rv==T) gg.tr <- gg.tr + scale_x_reverse(expand = expansion(mult = c(0.2, 0)))
	# FST plot (ordered by tree labels)
	fst.df <- data.frame(fst$PopFST, SNPs="clean")
	fst.df$Label <- factor(fst.df$Label, levels=rev(get_taxa_name(gg.tr)))
	gg.fst <- ggplot(fst.df, 
		aes(x=FST, y=Label, colour=Label)) +
		geom_hline(yintercept=1:5, colour="lightgray", linetype="dotted") +
		geom_vline(xintercept=0, colour="lightgray", linetype="dotted") +
		geom_point(size=4) +
		geom_segment(aes(x = FST-(2*SE), y=Label, xend = FST+(2*SE), yend = Label), size=1) +
		#geom_text(aes(label=Label), colour="white") +
		scale_fill_manual(values=pcols) +
		scale_colour_manual(values=pcols) +
		scale_y_discrete(position=ifelse(rv==T, "right", "left")) +
		labs(caption=ifelse(rv==F, "outlier SNPs", "non-outlier SNPs")) +
		xlab(expression(local~italic(F)[ST])) + 
		theme_bw() +
		theme(legend.position="none",
			panel.grid = element_blank(),
			axis.line.y=element_blank(),
		        axis.text.y=element_blank(),
		        axis.title.y=element_blank(),
			axis.title.x = element_text(size=12),
		        axis.text.x = element_text(size=10),
			plot.margin = margin(l=0, r=0),
		        plot.caption = element_text(size=14, face="bold", hjust = 0.5),
		        plot.tag = element_text(size=14, face="bold"))
		return(list(gg.tr, gg.fst))
}
gg.fst.clean <- fst_plot(fst.clean)
gg.fst.outl <- fst_plot(fst.outl)
gg.fst.nonoutl <- fst_plot(fst.nonoutl, rv=T)
pg <- plot_grid(gg.fst.clean [[1]], gg.fst.clean [[2]], gg.fst.nonoutl[[2]], gg.fst.nonoutl[[1]], ncol = 4, rel_widths=c(0.15, 0.35, 0.35, 0.15), align = "h") +
	labs(tag="c)") + theme(plot.tag.position = c(0.02,0.99), plot.tag = element_text(size=16, face="bold"))


pdf(file.path(outdir, "FST_new3.pdf"), width=12, height=6)
	grid.arrange(gg.fst, pg, ncol=2)
dev.off()


gg.up <- as.ggplot(upset(ups[,-1], colnames(ups)[-1], height_ratio = 0.6, name="Intersections", stripes='white',min_degree=1,keep_empty_groups=TRUE,
	set_sizes=(upset_set_size(mapping=aes(fill="bars_color")) + 
		scale_fill_manual(values=c('bars_color'='black'), guide='none') + 
		ylab("Total SNPs") + 
		theme(panel.grid=element_blank(), axis.ticks.x=element_line(), axis.line.x=element_line())),
	base_annotations=list('Shared SNPs'=intersection_size(mapping=aes(fill="bars_color")) + 
		scale_fill_manual(values=c('bars_color'='black'), guide='none') + 
		theme(panel.grid=element_blank(), axis.ticks.y=element_line(), axis.line.y=element_line())),
	matrix=(intersection_matrix(geom=geom_point(size=5), segment=geom_segment(size=3), outline_color=list(active='black', inactive='black')) + scale_color_manual(values=c('TRUE'='black', 'FALSE'='white'), guide='none'))
   )) + labs(tag="b)") + theme(plot.tag.position = c(0.02,0.99), plot.tag = element_text(size=16, face="bold"))
pdf(file.path(outdir, "FST.pdf"), width=12, height=10)
	grid.arrange(gg.da.clean, gg.depth, gg.up, pg, ncol=2, heights=c(0.6, 0.4))
dev.off()

###################################
# write filtered files for external analyses
###################################

print("Writing filtered files...")

# write clean files for RADPainter
vcf.hap <- read.vcfR(file.path(indir, "populations.haps.vcf"))
rad <- readLines(file.path(outdir, "..", "populations.haps.radpainter"))
rad.loci <- unique(gsub(":.*", "", vcf.hap@fix[,"CHROM"]))
rad.clean <- c(rad[1], rad[-1][which(rad.loci %in% unique(gsub(":.*", "", colnames(vcf.imp.nonoutl))))])
writeLines(rad.clean, con=file.path(outdir, "..", "populations.haps.nonoutl.radpainter"))

# write outlier loci (union of all three outlier approaches)
rad.outl <- c(rad[1], rad[-1][which(rad.loci %in% all.outlier.loci)])
writeLines(rad.outl, con=file.path(outdir, "..", "populations.haps.outliers.radpainter"))

# write clean PLINK files for fastSTRUCTURE
plink.map <- read.table(file.path(outdir, "..", "populations.plink.map"), colClasses = "character")
plink.ped <- read.table(file.path(outdir, "..", "populations.plink.ped"), colClasses = "character")

plink.clean <- which(plink.map[,2] %in% gsub(":", "_", colnames(vcf.imp.nonoutl)))
plink.clean.map <- plink.map[plink.clean,]
plink.clean.ped <- plink.ped[,c(1:6, sort(c(plink.clean*2+5, plink.clean*2+6)))]
write.table(plink.clean.map, file.path(outdir, "..", "populations.plink.nonoutl.map"), sep="\t", row.names=F, col.names=F, quote=F)
write.table(plink.clean.ped, file.path(outdir, "..", "populations.plink.nonoutl.ped"), sep="\t", row.names=F, col.names=F, quote=F)

# write outlier PLINK files for fastSTRUCTURE
plink.outl <- which(plink.map[,2] %in% gsub(":", "_", all.outlier.SNPs))
plink.map <- plink.map[plink.outl,]
plink.ped <- plink.ped[,c(1:6, sort(c(plink.outl*2+5, plink.outl*2+6)))]
write.table(plink.map, file.path(outdir, "..", "populations.plink.outliers.map"), sep="\t", row.names=F, col.names=F, quote=F)
write.table(plink.ped, file.path(outdir, "..", "populations.plink.outliers.ped"), sep="\t", row.names=F, col.names=F, quote=F)

############ STOP
q()

########### infer clusters from outliers
source("structureplot.R")
pc.outl <- PCA(vcf.imp.outl, "a)")
pc.outl.opt <- which(cumsum(pc.outl[[2]]$eig/sum(pc.outl[[2]]$eig))*100>=80)[1]
fc <- lapply(4:6, function(x) {
	grp <- find.clusters(vcf.imp.outl, n.pca=100, n.clust=x)$grp
	da <- dapc(vcf.imp.outl, grp, scale=F, n.da=2, center=T, n.pca=pc.outl.opt)
	structureplot(da$posterior, "", c("purple", "gray80", "darkslategray4", "#FF9900", "red", "green"))
})
pdf(file.path(outdir, "outliers_structureplot.pdf"))
	do.call(grid.arrange, c(fc, ncol=1))
dev.off()




###### summarise catalog annotations
# RepeatMasker
repm <- read.table(file.path(outdir, "../..", "catalog.fa.gz.RepeatMasker.out.summary"), colClasses="character")
# TSA BLASTN
trans <- readLines(file.path(outdir, "../..", "catalog.fa.gz.GEZX01.1_blastn.loci"))

# only keep records for loci in the current dataset
repm.clean <- unique(subset(repm, V1 %in% gsub(":.*", "", colnames(vcf.imp.clean))))
repm.outl <- unique(subset(repm, V1 %in% gsub(":.*", "", colnames(vcf.imp.outl))))
trans.clean <- trans[trans %in% gsub(":.*", "", colnames(vcf.imp.clean))]
trans.outl <- trans[trans %in% gsub(":.*", "", colnames(vcf.imp.outl))]
total.loci.clean <- length(unique(gsub(":.*", "", colnames(vcf.imp.clean))))
total.loci.outl <- length(unique(gsub(":.*", "", colnames(vcf.imp.outl ))))

repm.df <- rbind(
	data.frame(Set="clean", melt(table(repm.clean[,2])/total.loci.clean*100)),
	data.frame(Set="outliers", melt(table(repm.outl [,2])/total.loci.outl*100)),
	data.frame(Set=c("clean", "outliers"), Var1="Transcripts", value=c(length(trans.clean)/total.loci.clean*100,length(trans.outl)/total.loci.outl*100))	
)

pdf(file.path(outdir, "annotations.pdf"), width=12, height=10)
ggplot(repm.df, aes(x=Set, y=value)) +
	geom_bar(stat="identity") +
	facet_wrap(~Var1, scales="free_y") +
	theme_bw()
dev.off()




sum(repm$V1 %in% gsub(":.*", "", colnames(vcf.imp)))/length(unique(gsub(":.*", "", colnames(vcf.imp))))
sum(repm$V1 %in% dapc.outlier.loci)/length(dapc.outlier.loci)

sum(dapc.outlier.loci %in% trans)/length(dapc.outlier.loci)
sum(trans %in% gsub(":.*", "", colnames(vcf.imp)))/length(unique(gsub(":.*", "", colnames(vcf.imp))))




# subset PLINK datafiles by snpzip loci (for STRUCTURE)
plink.map <- read.table(file.path(outdir, "..", "populations.plink.map"), colClasses = "character")
plink.ped <- read.table(file.path(outdir, "..", "populations.plink.ped"), colClasses = "character")
rad <- readLines(file.path(outdir, "..", "populations.haps.radpainter"))
rad.loci <- unique(gsub("_.*", "", plink.map[,2]))
rad <- c(rad[1], rad[-1][which(rad.loci %in% sz.lociboth)])
writeLines(rad, con=file.path(outdir, "..", "populations.haps.radpainter.snpzip"))


plink.map <- plink.map[sort(sz),]
plink.ped <- plink.ped[,c(1:6, sort(c(sz*2+5, sz*2+6)))]
rad <- c(rad[1], rad[])
write.table(plink.map, file.path(outdir, "..", "populations.plink.snpzip.map"), sep="\t", row.names=F, col.names=F, quote=F)
write.table(plink.ped, file.path(outdir, "..", "populations.plink.snpzip.ped"), sep="\t", row.names=F, col.names=F, quote=F)









# use the DAPC outliers to make DAPC-based cluster inference
fc <- find.clusters(vcf.imp.clean[,dapc.outliers], n.clust=5, n.pca=100)
da.fc <- dapc(vcf.imp.clean[,dapc.outliers], fc$grp, scale=F, n.da=4, center=T, n.pca=pc80)

source("structureplot.R")
gg.fc <- structureplot(da.fc$posterior, "b)", c("white", "gray80", "black", "darkslategray4", "#FF9900"))

pdf(file.path(outdir, "DAPC_final.pdf"), width=12, height=8)
	grid.arrange(gg.ld, gg.fc, ncol=2)
dev.off()


