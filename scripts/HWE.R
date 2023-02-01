# load STACKS sumstats data
# this contains HWE statistics for each SNP in each population

load_sumstats <- function(f){
  dd <- read.table(f, sep="\t", header=F, stringsAsFactors = F)
  colnames(dd)[c(1,4,5,17,20)] <- c("Locus", "Pos", "Pop.ID", "Fis", "HWE.P.value")
  dd$Plog <- -log10(dd$HWE.P.value)
  dd$HWE <- cut(dd$Fis, c(-1, 0, 1), labels=c("-", "+"))
  dd$Sig=cut(dd$HWE.P.value, c(0,0.05, 1), labels = c("Sig", "N.S."))
  dd$Type[dd$Pop.ID %in% c("M3", "M4", "M8", "M3.P5", "M4.P5", "M8.P7")] <- "pure"
  dd$Type[dd$Pop.ID %in% c("M7", "M9", "M7.P6", "M7.P7", "M9.P7", "M9.P8")] <- "mixed" 
  if(length(grep("\\.", dd$Pop.ID))>0) {
    for(a in 1:nrow(meta)){
      x <- grep(meta[a,"Pop"], dd$Pop.ID)
      dd$Pop.ID[x] <- gsub(meta[a,"Pop"], meta[a,"Label"], dd$Pop.ID[x])
      pmap[pmap[,2]==meta[a,"Pop"], "Label"] <- meta[a,"Label"]
      dd$Label=dd$Pop.ID
    }  
  } else{
    dd$Label=meta$Label[match(dd$Pop.ID, meta$Pop)]
  }
  return(dd)
}

# HWE barplot by population
hwe.bar <- function(d, blank=F){
  bb <- aggregate(Locus~Label+HWE+Sig, d, length)  
  bb[bb$HWE=="-", "Locus"] <- -(bb[bb$HWE=="-", "Locus"])
  bb$V2 <- bb$Locus/table(d$Pop.ID)[1]
  bb <- bb[bb$Sig=="Sig",]
  if(length(unique(bb$Label))<=5) pc <- pcols
  if(length(unique(bb$Label))>5) pc <- pcols.split
  #print(bb)
  gg <- ggplot(bb, aes(x=Label, y=V2, fill=Label, alpha=HWE)) + 
    geom_bar(stat="identity") +
    geom_hline(yintercept=0) +
    geom_label(aes(x=Label, y=V2, label=Label, color=Label), data=bb[bb$HWE=="+",], 
               vjust=0.5, size=3, fill="white", label.r = unit(0, "lines")) +
    scale_y_continuous(labels=percent, limits=c(-0.01, 0.20), breaks=seq(-0.01, 0.1, 0.01)) +
    scale_x_discrete(breaks=NULL) +
    scale_fill_manual(values=pc) +
    scale_color_manual(values=pc) +
    scale_alpha_manual(values=c(0.3, 1)) +
    ylab("Proportion of loci") +
    theme_classic() +
    theme(legend.position = "none",
          axis.title.x = element_blank())
  if(blank) gg <- gg + theme(axis.title.y=element_blank(),
                             axis.text.y = element_blank())
  return(gg)
}
