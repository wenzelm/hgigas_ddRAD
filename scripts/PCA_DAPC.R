source("glPCAFast.R")
PCA <- function(dd, tt="", ll="Label"){
	if(class(dd)=="genlight"){
		pca <- glPcaFast(dd, nf=2)
		df <- data.frame(pca$scores, pmap)
	} else if(class(dd)=="dist"){
		pca <- prcomp(dd, center=T, scale=F)
		names(pca)[1] <- "eig"
		#pca <- cmdscale(dd)
		df <- data.frame(pca$x, subset(pmap, ID %in% attr(dd, "Labels")))
		#df <- data.frame(pca, subset(pmap, ID %in% attr(dd, "Labels")))
		colnames(df)[1:2] <- c("PC1", "PC2")	
	} else {
		pca <- dudi.pca(dd, center=T, scale=F, scannf=F, nf=2)
		df <- data.frame(pca$li, subset(pmap, ID %in% rownames(dd)))
		colnames(df)[1:2] <- c("PC1", "PC2")	
	}

	hull.data <- do.call(rbind, lapply(split(df, df$Pop), function(x) x[chull(x[,c("PC1", "PC2")]),]))
	# % variance
	pct.var <- round(pca$eig/sum(pca$eig)*100, 2)
	gg <- ggplot(df, aes(x=PC1, y=PC2, fill=Label)) +
		geom_vline(xintercept=0, colour="lightgray", linetype="dotted") +
		geom_hline(yintercept=0, colour="lightgray", linetype="dotted") +
		geom_polygon(data=hull.data, aes(colour=Label), alpha=0.1, linetype="dashed") +
		geom_point(shape=21, size=8, alpha=1, color="white") +
		geom_text(aes_string(label=ll), colour="white", size=5) +
		scale_fill_manual(name="Site", values=pcols) +
		scale_colour_manual(name="Site", values=pcols) +
		#ggtitle(paste(ncol(dd), "SNPs")) +
		labs(tag=tt) +
		xlab(paste0("PC1 (", pct.var[1], "%)")) + ylab(paste0("PC2 (", pct.var[2], "%)")) +
		theme_bw() +
		theme(legend.position=ifelse(ll=="Label", "none", "top"),
			panel.grid = element_blank(),
			axis.title = element_text(size=14),
		        axis.text = element_text(size=12),
		        plot.tag.position = c(0.008,0.99),
		        plot.tag = element_text(size=16, face="bold"))
	return(list(gg, pca))
}

plot.DAPC <- function(df, tt="", ll="Label", xa="LD1", ya="LD2"){
	hull.data <- do.call(rbind, lapply(split(df, df$Label), function(x) x[chull(x[,c(xa, ya)]),]))
	g <- ggplot(df, aes_string(x=xa, y=ya, fill="Label")) +
		geom_vline(xintercept=0, colour="lightgray", linetype="dotted") +
		geom_hline(yintercept=0, colour="lightgray", linetype="dotted") +
		geom_polygon(data=hull.data, aes(colour=Label), alpha=0.1, linetype="dashed") +
		geom_point(shape=21, size=8, alpha=1, color="white") +
		geom_text(aes_string(label=ll), colour="white", size=5) +
		scale_fill_manual(name="Site", values=pcols) +
		scale_colour_manual(name="Site", values=pcols) +
		xlab(xa) + ylab(ifelse(xa=="LD3", "", ya)) + labs(tag=tt) +
		theme_bw() +
		theme(legend.position=ifelse(ll=="Label", "none", "top"),
			panel.grid = element_blank(),
			axis.title = element_text(size=14),
		        axis.text = element_text(size=12),
		        plot.tag.position = c(0.008,0.99),
		        plot.tag = element_text(size=16, face="bold"),
			plot.margin = margin(5.5, 1, 5.5, 5.5, "pt"))
	if(xa=="LD3") g <- g + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), plot.margin = margin(5.5, 5.5, 5.5, 1, "pt"))
	return(g)
}

# DAPC.sweep
# sweep across DAs
# 2) sweep across numbers of PCs (15 to 46)
library(vegan)
DAPCsweep <- function(dd, ll){
	sw <- lapply(15:46, function(a) dapc(dd, pmap$Label, scale=F, n.da=2, center=T, n.pca=a))
	sw.df <- lapply(sw, function(x) data.frame(x$ind.coord, pmap))
	names(sw.df) <- as.character(15:46)
	pdf(file.path(outdir, paste0(ll,".pdf")), width=24, height=32)
		do.call(grid.arrange, c(lapply(15:46, function(a) plot.DAPC(sw.df[[as.character(a)]], a)), ncol=4, nrow=8))
	dev.off()
	# correlation tests sweep
	# partial Mantel tests
	ct <- do.call(rbind, lapply(sw, function(a){
		df <- data.frame(a$ind.coord, pmap)
		ct1 <- cor.test(df$LD1, df$Lat, method="spearman")
		ct2 <- cor.test(df$LD1, df$Depth, method="spearman")
		ct3 <- cor.test(df$LD2, df$Lat, method="spearman")
		ct4 <- cor.test(df$LD2, df$Depth, method="spearman")
		mt1 <- mantel.partial(dist(df$LD1), dist(df$Depth), dist(df[,c("Long", "Lat")]), permutations=9999)
		mt2 <- mantel.partial(dist(df$LD2), dist(df$Depth), dist(df[,c("Long", "Lat")]), permutations=9999)
		c(
		round(ct1$estimate, 3),
		round(ct1$p.value, 3),
		round(ct2$estimate, 3),
		round(ct2$p.value, 3),
		round(ct3$estimate, 3),
		round(ct3$p.value, 3),
		round(ct4$estimate, 3),
		round(ct4$p.value, 3),
		round(mt1$statistic, 3),
		round(mt1$signif, 3),
		round(mt2$statistic, 3),
		round(mt2$signif, 3))
	}))
	rownames(ct) <- 15:46
	colnames(ct) <- c("LD1.Lat.r", "LD1.Lat.P",
		"LD1.Depth.r", "LD1.Depth.P",
		"LD2.Lat.r", "LD2.Lat.P",
		"LD2.Depth.r", "LD2.Depth.P",
		"LD1.Mantel.r", "LD1.Mantel.P",
		"LD2.Mantel.r", "LD2.Mantel.P")
	write.table(ct, file.path(outdir, paste0(ll,".cortests.txt")), sep="\t", col.names=NA, quote=F)
	return(sw)
}

