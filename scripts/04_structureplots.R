library(ggplot2)
library(gridExtra)
library(reshape2)

indir=commandArgs(trailingOnly=TRUE)[1]
#indir="ustacks_m3M10_nogaps_HQ_allowN/populations_p5r0.2_QCsnps"
#indir="sweep_m3M5/singletonsHD+VS_p5r0.2"
#indir="sweep_m3M4_n0/singletonsHD+VS_p5r0.5/"
outdir=file.path(indir, "R_output")
dir.create(outdir)

# load DAPC results
da <- readRDS(list.files(path=outdir, pattern="DAPC_eigenval.*RDS", full.names=T))
da.clusters <- readRDS(list.files(path=outdir, pattern="DAPC.*findclusters.*RDS", full.names=T))

# load fastSTRUCTURE results
struc <- lapply(list.files(file.path(indir, "fastStructure"), pattern="*logistic.*meanQ", full.names=T), function(x) {
	x <- read.table(x)
	colnames(x) <- 1:ncol(x)
	rownames(x) <- sort(as.numeric(rownames(da$ind.coord)))
	return(x) 
})

# metadata

pmap <- readRDS(file.path(outdir, "popmap_metadata.RDS"))

# custom compoplot
structureplot <- function(memberprobs){
  # synchronise colours
  df <- cbind(ID=as.numeric(rownames(memberprobs)), memberprobs)
  df <- merge(pmap[,c("ID", "Label")], df, by="ID")
  df <- aggregate(df[,-c(1:2)], list(df$Label), median)[,-1]
  o <- unique(unlist(sapply(1:ncol(df), function(x){ which.max(df[c(1,3,5)[x],])  })))
  o <- c(o, setdiff(1:ncol(df), o))
  
  memberprobs <- memberprobs[,o]
  colnames(memberprobs) <- 1:ncol(memberprobs)
  df <- melt(as.matrix(memberprobs))
  df <- merge(pmap[,c("ID", "Label")], df, by.x="ID", by.y="Var1")
  df <- df[order(df$Label, df$ID, df$Var2),]
  df$X <- rep(1:(nrow(df)/max(df$Var2)), each=max(df$Var2))
  lbl <- data.frame(X=tapply(df$X, df$Label, mean))
  lbl$Label <- rownames(lbl)
  pcols <- c(A="darkblue", B="steelblue3", C="#FF9900", D="red2", E="chartreuse4")
  pcols.cl <- c(pcols[1], "gray88", pcols[5], pcols[3])
  names(pcols.cl) <- 1:4
  gg <- ggplot() + 
    geom_bar(data=df, aes(x=X, y=value, fill=factor(Var2)), position="stack", stat="identity", width=1) + 
    #geom_vline(xintercept=df$X[which(!duplicated(df$Label))[-1]]-0.5, size=1) +
    geom_segment(aes(x = df$X[which(!duplicated(df$Label))[-1]]-0.5, y = -0.1, xend = df$X[which(!duplicated(df$Label))[-1]]-0.5, yend = 1), size=1) +
    #geom_text(aes(label=Label, colour=factor(Label)), y=-0.05) +
    geom_text(data=lbl, aes(x=X, label=Label, colour=factor(Label), y=-0.1), fontface="bold") +
    #scale_fill_manual(values=c(`1`="black", `2`="gray88", `3`="darkmagenta", `4`="green")) +
    scale_fill_manual(values=pcols.cl) +
    scale_colour_manual(values=pcols) +
    scale_y_continuous(limits=c(-0.15,1.02)) +
    scale_x_discrete(expand=c(0,0)) +
    ylab(paste0("Q (K=", ncol(memberprobs), ")")) +
    theme_classic() +
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.line.x = element_blank(),
          #panel.grid = element_line(colour = "black", size=1),
          #panel.grid.major = element_blank(), 
          #panel.grid.minor = element_blank(),
          legend.position="none")
  return(gg)
}

pdf(file.path(outdir, "STRUCTUREplots_logistic.PDF"), width=8, height=4)

grid.arrange(structureplot(da.clusters[[1]]$posterior), structureplot(struc[[1]]),
             structureplot(da.clusters[[2]]$posterior), structureplot(struc[[2]]),
             structureplot(da.clusters[[3]]$posterior), structureplot(struc[[3]]), ncol=2)


dev.off()
