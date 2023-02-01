library(ggplot2)
library(gridExtra)
library(grid)
library(reshape2)
library(ggtree)

#
# generate figures from all analyses
#
indir="stacks_m3M7_n0/singletonsHD+VS_D2_p5r0.5/"
outdir=file.path(indir, "R_output")
dir.create(outdir)

pmap <- readRDS(file=file.path(outdir, "popmap_metadata.RDS"))

# STRUCTUREplots from fastSTRUCTURE
source("structureplot.R")

# load fastSTRUCTURE results
# all Ks are loaded into the same object

pcols <- c(A="#FF9900", B="dodgerblue", C="black", D="chartreuse3", E="red2")
struc.cols <- c("purple", "gray80", "darkslategray4", "#FF9900", "red", "green")
struc.cols <- pcols[c(1,2,4,5,3)]

compoundstructure <- function(dd,tt=""){
	struc <- lapply(list.files(file.path(indir, "fastStructure"), pattern=paste0(dd, ".*meanQ"), full.names=T), function(x) {
		x <- read.table(x)
		colnames(x) <- 1:ncol(x)
		rownames(x) <- pmap$ID
		return(x) 
	})
	# only last plot must have x labels
	gg.1 <- lapply(struc[1:3], structureplot, "", struc.cols)
	gg.2 <- structureplot(struc[[4]], "", struc.cols, T)
	arrangeGrob(gg.1[[1]], gg.1[[2]], gg.1[[3]], gg.2, ncol=1, top=textGrob(tt,gp=gpar(fontsize=16)))
}
#struc.gg <- compoundstructure("clean.ld.logistic")
#struc.gg2 <- compoundstructure("outliers.logistic")

pdf(file.path(outdir, "fastSTRUCTUREplots.pdf"), width=12, height=10)
grid.arrange(compoundstructure("nonoutl.logistic", "Non-outliers"),
	     compoundstructure("nonoutl.ld.logistic", "Non-outliers in linkage equilibrium"),
	     compoundstructure("outliers.logistic", "Outliers"), 
	     compoundstructure("outliers.ld.logistic", "Outliers in linkage equilibrium"), ncol=2)
grid.arrange(compoundstructure("nonoutl.simple", "Non-outliers"),
	     compoundstructure("nonoutl.ld.simple", "Non-outliers in linkage equilibrium"),
	     compoundstructure("outliers.simple", "Outliers"), 
	     compoundstructure("outliers.ld.simple", "Outliers in linkage equilibrium"), ncol=2)
dev.off()

# fineRADSTRUCTURE
source("PCA_DAPC.R")

# coancestry matrix
coa <- lapply(list.files(file.path(indir, "RADpainter"), pattern="*chunks.out$", full.names=T), read.table, header=T, check.names=F)
coa <- lapply(coa, function(x) as.dist(x[,-1]))
pc.rad <- lapply(coa, function(x) PCA(x)[[1]])
pc.rad <- list(PCA(coa[[1]])[[1]] + ggtitle("Non-outliers") + theme(plot.title = element_text(size=18, hjust=0.5)),
		PCA(coa[[2]])[[1]] + ggtitle("Outliers") + theme(plot.title = element_text(size=18, hjust=0.5)))

# dendrograms
#tr <- read.tree(file.path(indir, "RADpainter", "populations.haps_chunks.out.tree.nw"))

tr <- lapply(list.files(file.path(indir, "RADpainter"), pattern="*.out.tree.nw", full.names=T), read.tree)

gg.tr <- lapply(tr, function(x) ggtree(x, branch.length="none", layout="circular") %<+% data.frame(Taxon=as.character(pmap$ID), pmap) + 
	geom_tippoint(aes(fill=Label), size=8, shape=21, colour="white") +
	geom_tiplab(aes(label=Label), angle=0, size=4, colour="white", hjust=0.6) +
	scale_fill_manual(values=pcols) +
	#layout_dendrogram() +
	#labs(tag="d)") +
	theme(legend.position="none",
		plot.tag.position = c(0.005,0.99),
		plot.tag = element_text(size=16, face="bold")))


pdf(file.path(outdir, "radpainter.pdf"), width=12, height=12)
	#do.call(grid.arrange, c(list(structureplot(qmat[[1]], "c)", pcols)), pc.rad, gg.tr, ncol=2))
	do.call(grid.arrange, c(pc.rad, gg.tr, ncol=2))
dev.off()



# barplot
qmat <- lapply(list.files(file.path(indir, "RADpainter"), pattern="*.admix.Q$", full.names=T), function(x) {
	qmat <- read.csv(x, header=T)
	rownames(qmat) <- qmat[,1]
	qmat <- qmat[,-1]
	colnames(qmat) <- 1:ncol(qmat)
	qmat
})

source("structureplot.R")

pdf(file.path(outdir, "radpainter.pdf"), width=12, height=14)
	#do.call(grid.arrange, c(list(structureplot(qmat[[1]], "c)", pcols)), pc.rad, gg.tr, ncol=2))
	do.call(grid.arrange, c(pc.rad, gg.tr, list(structureplot(qmat[[2]], "", pcols)),  ncol=2))
dev.off()


pdf(file.path(outdir, "multipanel.pdf"), width=12, height=8)
	grid.arrange(pc.all[[1]], pc.rad[[1]], 
			structureplot(qmat, "c)", c("white", "gray80", "black", "darkslategray4", "#FF9900")), 
			gg.tr,
			layout_matrix=matrix(c(1,2,3,3,4,4), ncol=2, byrow=T), heights=c(0.5,0.25,0.25))
dev.off()

###### summarise locus annotations
# RepeatMasker hits
repm <- read.table(file.path(outdir, "../..", "catalog.fa.gz.RepeatMasker.out.summary"), colClasses="character")
# TSA BLASTN hits
trans <- readLines(file.path(outdir, "../..", "catalog.fa.gz.GEZX01.1_blastn.loci"))
# Uniprot hits
uniprot <- readLines(file.path(outdir, "../..", "catalog.fa.gz.uniprot_blastx.loci"))

# locus names
loci.singletons <- unique(gsub("_.*", "", unlist(read.table(file.path(indir, "populations.plink.map"), colClasses="character")[2])))
loci.clean <- unique(gsub(":.*", "", readLines(file.path(outdir, "SNPnames.clean.txt"))))
loci.outl <- unique(gsub(":.*", "", readLines(file.path(outdir, "SNPnames.outliers.txt"))))
loci.paralogs <- unique(gsub("_.*", "", unlist(read.table(file.path(gsub("singletons", "paralogs", indir), "populations.plink.map"), colClasses="character")[2])))
loci.nooutl <- loci.clean[! loci.clean %in% loci.outl]
# alternative paralogue set (without populations filters)
loci.paralogs <- readLines(file.path(indir, "..", "combined_paralogsHD+VS_D2_whitelist.txt"))

# only keep records for loci in the current dataset
repm.keep <- lapply(list(loci.singletons, loci.clean, loci.outl, loci.nooutl, loci.paralogs), function(x) unique(subset(repm, V1 %in% x)))
trans.keep <- lapply(list(loci.singletons, loci.clean, loci.outl, loci.nooutl, loci.paralogs), function(x) trans[trans %in% x])
uniprot.keep <- lapply(list(loci.singletons, loci.clean, loci.outl, loci.nooutl, loci.paralogs), function(x) uniprot[uniprot%in% x])

# total numbers of loci
total.loci.singletons <- length(loci.singletons)
total.loci.clean <- length(loci.clean)
total.loci.outl <- length(loci.outl)
total.loci.paralogs <- length(loci.paralogs)
total.loci.nooutl <- length(loci.nooutl)

# transcripts %
df <- data.frame(Set=c("singletons", "clean", "outliers", "non-outliers", "paralogs"), 
		 Var1="Transcripts", 
		 value=sapply(trans.keep, length)/c(total.loci.singletons, total.loci.clean, total.loci.outl, total.loci.nooutl, total.loci.paralogs),
	 	 Number.Var1="dummy",
		 Number.Freq=sapply(trans.keep, length))

# uniprot %
df <- rbind(df, data.frame(Set=c("singletons", "clean", "outliers", "non-outliers", "paralogs"), 
		 Var1="UniProt", 
		 value=sapply(uniprot.keep, length)/c(total.loci.singletons, total.loci.clean, total.loci.outl, total.loci.nooutl, total.loci.paralogs),
	 	 Number.Var1="dummy",
		 Number.Freq=sapply(uniprot.keep, length)))

# repeats %
df <- rbind(df,
	data.frame(Set="singletons", melt(table(repm.keep[[1]][,2])/total.loci.singletons), Number=table(repm.keep[[1]][,2])),
	data.frame(Set="clean", melt(table(repm.keep[[2]][,2])/total.loci.clean), Number=table(repm.keep[[2]][,2])),
	data.frame(Set="outliers", melt(table(repm.keep[[3]][,2])/total.loci.outl), Number=table(repm.keep[[3]][,2])),
	data.frame(Set="non-outliers", melt(table(repm.keep[[4]][,2])/total.loci.nooutl), Number=table(repm.keep[[4]][,2])),
	data.frame(Set="paralogs", melt(table(repm.keep[[5]][,2])/total.loci.paralogs), Number=table(repm.keep[[5]][,2]))
)

gg <- ggplot(subset(df, Set %in% c("outliers", "non-outliers")), aes(x=Set, y=value)) +
	geom_bar(stat="identity") +
	geom_label(aes(label=Number.Freq)) +
	facet_wrap(~Var1, scales="free_x", nrow=2) +
	coord_flip() +
	theme_bw()

cols2 <- c(`outliers`="#FF9900", `non-outliers`="darkgray", `paralogs`="firebrick")
df$Var1 <- factor(df$Var1, levels=rev(levels(df$Var1)))
gg <- ggplot(subset(df, Set %in% c("outliers", "non-outliers", "paralogs")), aes(x=Var1, y=value, fill=Set)) +
	geom_bar(stat="identity", position=position_dodge(width=0.6), width=0.3) +
	geom_text(aes(label=Number.Freq, y=-0.02, colour=Set), position=position_dodge(width=0.6), size=2, fontface="bold", hjust=1) +
	scale_y_continuous(labels = scales::percent, breaks=seq(0, 0.6, 0.2), limits=c(-0.07, 0.6)) +
	scale_fill_manual(name="Loci", values=cols2) +
	scale_colour_manual(name="Loci", values=cols2) +
	xlab("") +
	ylab("% of loci") +
	labs(tag="a)") +
	coord_flip() +
	theme_classic() +
	theme(legend.position="top",
		axis.title.y = element_blank())


pdf(file.path(outdir, "annotations.pdf"), width=12, height=6)
	gg
dev.off()

# proportion tests
# outliers vs. non-outliers
prop.test(sapply(repm.keep[3:4], nrow), c(total.loci.outl, total.loci.nooutl))
prop.test(sapply(trans.keep[3:4], length), c(total.loci.outl, total.loci.nooutl))
prop.test(sapply(uniprot.keep[3:4], length), c(total.loci.outl, total.loci.nooutl))

# paralogs vs. clean loci
prop.test(sapply(repm.keep[c(2,5)], nrow), c(total.loci.clean, total.loci.paralogs))
prop.test(sapply(trans.keep[c(2,5)], length), c(total.loci.clean, total.loci.paralogs))
prop.test(sapply(uniprot.keep[c(2,5)], length), c(total.loci.clean, total.loci.paralogs))

# 3-way comparison
prop.test(sapply(trans.keep[3:5], length), c(total.loci.outl, total.loci.nooutl, total.loci.paralogs))

#prop.test(length(trans.keep[[3]]), total.loci.outl, 	p=length(trans.keep[[2]])/total.loci.clean)
#prop.test(length(trans.keep[[4]]), total.loci.paralogs, p=length(trans.keep[[2]])/total.loci.clean)
#prop.test(length(uniprot.keep[[3]]), total.loci.outl, 	p=length(uniprot.keep[[2]])/total.loci.clean)
#prop.test(length(uniprot.keep[[4]]), total.loci.paralogs, p=length(uniprot.keep[[2]])/total.loci.clean)

# Geneontology enrichment
library(clusterProfiler)
blastx <- read.table(file.path(outdir, "../..", "catalog.fa.gz.uniprot_blastx"), stringsAsFactors=F)
# strict filtering
#blastx <- subset(blastx, V11<1e-10 & V3>=50)
blastx <- blastx[,1:2]
blastx$V1 <- as.character(blastx$V1)
blastx$Uniprot <- gsub("\\|", "", regmatches(blastx$V2,regexpr("\\|.*\\|",blastx$V2)))
ugo <- read.table("blast/uniprot-taxonomy__Arthropoda6656_GO.tsv.gz", sep="\t", header=F, colClasses="character", quote = "")
names(ugo) <- c("Uniprot", "Entry")
ugo$Entry <- gsub("^ ", "", ugo$Entry)
ugo <- merge(blastx, ugo, by="Uniprot")
write.table(ugo, file.path(outdir, "GO_uniprot_database.txt"), col.names=T, row.names=F, sep="\t", quote=F)

# background
bg <- unique(subset(ugo, V1 %in% loci.singletons & Entry !="")[,c("Entry", "V1")])
bg.para <- unique(subset(ugo, Entry !="")[,c("Entry", "V1")])

# enrich
enr.outl <- enricher(loci.outl, TERM2GENE=bg, minGSSize=1)
enr.nooutl <- enricher(loci.nooutl, TERM2GENE=bg, minGSSize=1)
enr.para <- enricher(loci.paralogs, TERM2GENE=bg.para, minGSSize=1)
enr.outl@result$UniProt <- sapply(enr.outl@result$geneID, function(x) paste(sapply(unlist(strsplit(x, "/")), function(y) paste(unique(subset(ugo, V1 == y & Entry != "")$Uniprot), collapse=",")), collapse=";"))
enr.nooutl@result$UniProt <- sapply(enr.nooutl@result$geneID, function(x) paste(sapply(unlist(strsplit(x, "/")), function(y) paste(unique(subset(ugo, V1 == y & Entry != "")$Uniprot), collapse=",")), collapse=";"))
enr.para@result$UniProt <- sapply(enr.para@result$geneID, function(x) paste(sapply(unlist(strsplit(x, "/")), function(y) paste(unique(subset(ugo, V1 == y & Entry != "")$Uniprot), collapse=",")), collapse=";"))
write.table(enr.outl@result, file.path(outdir, "GO_enrich_outliers.txt"), col.names=T, row.names=F, sep="\t", quote=F)
write.table(enr.nooutl@result, file.path(outdir, "GO_enrich_nooutliers.txt"), col.names=T, row.names=F, sep="\t", quote=F)
write.table(enr.para@result, file.path(outdir, "GO_enrich_paralogs.txt"), col.names=T, row.names=F, sep="\t", quote=F)



# examine overlaps of GO terms among the three sets
go.ups <- data.frame(term=unique(c(enr.outl@result$ID,enr.nooutl@result$ID,enr.para@result$ID)), outliers=0, `non-outliers`=0, paralogs=0)
go.ups$outliers[go.ups$term %in% enr.outl@result$ID] <- 1
go.ups$non.outliers[go.ups$term %in% enr.nooutl@result$ID] <- 1
go.ups$paralogs[go.ups$term %in% enr.para@result$ID] <- 1
library(UpSetR)
pdf(file.path(outdir, "annotations_upset2.pdf"), width=6, height=3)
	upset(go.ups)
dev.off()

#

# get gene ratios for top-6 GO terms
terms <- unique(c(head(enr.outl@result[rev(order(enr.outl@result$Count)),], n=6)$ID,
head(enr.nooutl@result[rev(order(enr.nooutl@result$Count)),], n=6)$ID,
head(enr.para@result[rev(order(enr.para@result$Count)),], n=6)$ID))

gr <- rbind(data.frame(Loci="outliers", Total=total.loci.outl, subset(enr.outl@result, ID %in% terms)),
data.frame(Loci="non-outliers", Total=total.loci.nooutl, subset(enr.nooutl@result, ID %in% terms)),
data.frame(Loci="paralogs", Total=total.loci.paralogs, subset(enr.para@result, ID %in% terms)))

gr$Ratio <- sapply(gr$GeneRatio, function(x) eval(parse(text=x)))
gr$Pct <- gr$Count/gr$Total
gr$ID <- gsub("\\[GO.*\\]", "", gr$ID)
gr$ID <- factor(gr$ID, levels=unique(gr$ID[order(gr$Ratio)]))

# proportion test
pr <- sapply(split(gr, gr$ID), function(x) prop.test(x$Count, as.numeric(gsub(".*/", "", x$GeneRatio)))$p.value)

gg.gr <- ggplot(gr, aes(x=ID, y=Ratio, fill=Loci)) +
	geom_bar(stat="identity", position=position_dodge(width=0.6), width=0.3) +
	geom_text(aes(label=Count, y=-0.02, colour=Loci), position=position_dodge(width=0.6), size=2, fontface="bold", hjust=1) +
	scale_y_continuous(labels = scales::percent, breaks=seq(0, 0.5, 0.1), limits=c(-0.04, 0.55)) +
	scale_fill_manual(name="Loci", values=cols2) +
	scale_colour_manual(name="Loci", values=cols2) +
	xlab("") +
	ylab("% of loci") +
	labs(tag="b)") +
	coord_flip() +
	theme_classic() +
	theme(legend.position="top",
		axis.title.y=element_blank())
# combined figure
cmbdd <- subset(df, Set %in% c("outliers", "non-outliers", "paralogs"))
cmbdd$Var1 <- as.character(cmbdd$Var1)
cmbdd[! cmbdd$Var1 %in% c("Transcripts", "UniProt"),"Var1"] <- "Repetitive elements"
cmbdd <- aggregate(cmbdd[,c("value", "Number.Freq")], list(cmbdd$Set, cmbdd$Var1), sum)
names(cmbdd) <- c("Loci", "Description", "Ratio", "Count")
cmbdd <- rbind(cmbdd, gr[,c("Loci", "Description", "Ratio", "Count")])
cmbdd$Description <- factor(cmbdd$Description, levels=c(unique(gr$Description[order(gr$Ratio)]), "Repetitive elements", "UniProt", "Transcripts"))

gg.gr <- ggplot(cmbdd, aes(x=Description, y=Ratio, fill=Loci)) +
	geom_bar(stat="identity", position=position_dodge(width=0.9), width=0.5) +
	geom_text(aes(label=Count, y=-0.02, colour=Loci), position=position_dodge(width=0.9), size=3, fontface="bold", hjust=1) +
	scale_y_continuous(labels = scales::label_percent(accuracy = 1L), breaks=seq(0, 0.6, 0.1), limits=c(-0.06, 0.60)) +
	scale_fill_manual(name="Loci", values=cols2) +
	scale_colour_manual(name="Loci", values=cols2) +
	xlab("") +
	ylab("% of loci") +
	#labs(tag="b)") +
	coord_flip() +
	theme_classic() +
	theme(legend.position="right",
		axis.title.y=element_blank(),
		axis.text.y = element_text(size=14, face="bold"),
		axis.title = element_text(size=14),
		axis.text.x = element_text(size=12),
		legend.title = element_text(size=14),
		legend.text = element_text(size=12))
pdf(file.path(outdir, "annotations4.pdf"), width=12, height=7)
	gg.gr
dev.off()


pdf(file.path(outdir, "annotations2.pdf"), width=12, height=6)
	grid.arrange(gg, gg.gr, ncol=2, widths=c(0.4, 0.6))
dev.off()
