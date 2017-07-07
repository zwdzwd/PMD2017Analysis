##########################
## sequence context plot
##########################

context <- read.table('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/cgcontext.hg19/cgcontext.merged.bed', header=F, stringsAsFactors=F, sep='\t', col.names=c('chrm','beg','end','orph35','gc35','context1','context2'))
save(context, file='/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/cgcontext.hg19/context.rda')

contextHMD <- read.table('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/cgcontext.hg19/cgcontext.merged.comHMD.bed', header=F, stringsAsFactors=F, sep='\t', col.names=c('chrm','beg','end','orph35','gc35','context1','context2'))
contextPMD <- read.table('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/cgcontext.hg19/cgcontext.merged.comPMD.bed', header=F, stringsAsFactors=F, sep='\t', col.names=c('chrm','beg','end','orph35','gc35','context1','context2'))
save(contextHMD, file='/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/cgcontext.hg19/context.comHMD.rda')
save(contextPMD, file='/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/cgcontext.hg19/context.comPMD.rda')

pdf('~/gallery/2016_12_27_comPMD_gc_dist.pdf', width=3.5, height=3)
par(mar=c(5,5,3,2))
plot(density(contextPMD$gc35, bw=0.05), main='common PMD', xlab='+/-35bp GC content', ylab='Density', cex.lab=1.3, cex.axis=1.3)
dev.off()

pdf('~/gallery/2016_12_27_comHMD_gc_dist.pdf', width=3.5, height=3)
par(mar=c(5,5,3,2))
plot(density(contextHMD$gc35, bw=0.05), main='common HMD', xlab='+/-35bp GC content', ylab='Density', cex.lab=1.3, cex.axis=1.3)
dev.off()

pdf('~/gallery/2016_12_27_all_gc_dist.pdf', width=3.5, height=3)
par(mar=c(5,5,3,2))
plot(density(context$gc35, bw=0.05), main='all', xlab='+/-35bp GC content', ylab='Density', cex.lab=1.3, cex.axis=1.3)
dev.off()

contextPMD$orph35class <- cut(contextPMD$orph35, breaks=c(0,1,2,3,4,Inf), right=F)
contextPMD$gc35class <- cut(contextPMD$gc35, breaks=c(0,0.35,0.5,0.65,1), include.lowest=T)
pdf('~/gallery/2016_12_27_orphan_pmd.pdf', width=7, height=5)
qplot(orph35class, data=contextPMD, fill=gc35class, geom='bar') + scale_fill_grey(start=0.85,end=0.1) + xlab('#Flanking CpGs') + ylab('Count') + scale_x_discrete(labels=c('0','1','2','3','4+')) + theme(text=element_text(size=16)) + guides(fill=guide_legend(title='GC content \n+/-35bp')) + ggtitle('common PMD')
dev.off()

contextHMD$orph35class <- cut(contextHMD$orph35, breaks=c(0,1,2,3,4,Inf), right=F)
contextHMD$gc35class <- cut(contextHMD$gc35, breaks=c(0,0.35,0.5,0.65,1), include.lowest=T)
pdf('~/gallery/2016_12_27_orphan_hmd.pdf', width=7, height=5)
qplot(orph35class, data=contextHMD, fill=gc35class, geom='bar') + scale_fill_grey(start=0.85,end=0.1) + xlab('#Flanking CpGs') + ylab('Count') + scale_x_discrete(labels=c('0','1','2','3','4+')) + theme(text=element_text(size=16)) + guides(fill=guide_legend(title='GC content \n+/-35bp')) + ggtitle('common HMD')
dev.off()

base1 <- ifelse(substr(contextPMD$context1,1,1) %in% c('A','T'), 'W', 'S')
base4 <- ifelse(substr(contextPMD$context1,4,4) %in% c('A','T'), 'W', 'S')
a <- paste0(base1, 'CG', base4)
a[base1 != base4] <- 'SCGW'
contextPMD$contextclass <- as.factor(a)
pdf('~/gallery/2016_12_27_orphan_pmd_context.pdf', width=7, height=5)
qplot(orph35class, data=contextPMD, fill=contextclass, geom='bar') + scale_fill_grey(start=0.85,end=0.1) + xlab('#Flanking CpGs') + ylab('Count') + scale_x_discrete(labels=c('0','1','2','3','4+')) + theme(text=element_text(size=16)) + guides(fill=guide_legend(title='Context')) + ggtitle('common PMD')
dev.off()

base1 <- ifelse(substr(contextHMD$context1,1,1) %in% c('A','T'), 'W', 'S')
base4 <- ifelse(substr(contextHMD$context1,4,4) %in% c('A','T'), 'W', 'S')
a <- paste0(base1, 'CG', base4)
a[base1 != base4] <- 'SCGW'
contextHMD$contextclass <- as.factor(a)
pdf('~/gallery/2016_12_27_orphan_hmd_context.pdf', width=7, height=5)
qplot(orph35class, data=contextHMD, fill=contextclass, geom='bar') + scale_fill_grey(start=0.85,end=0.1) + xlab('#Flanking CpGs') + ylab('Count') + scale_x_discrete(labels=c('0','1','2','3','4+')) + theme(text=element_text(size=16)) + guides(fill=guide_legend(title='Context')) + ggtitle('common HMD')
dev.off()

###########################
## sequence context count
###########################
setwd('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS')
load('allprone/bin_mean_100kb.rda')
alltcga.samples <- colnames(betas)[grep('TCGA_', colnames(betas))]
tumorsamples <- alltcga.samples[!grepl('_N', alltcga.samples)]
normalsamples <- alltcga.samples[grepl('_N', alltcga.samples)]
sd.binmeanT <- na.omit(apply(betas[,tumorsamples], 1, sd, na.rm=T))
sd.binmeanN <- na.omit(apply(betas[,normalsamples], 1, sd, na.rm=T))

coords <- strsplit(names(sd.binmeanT),'[:-]')
write.table(data.frame(chrm=sapply(coords, function(x) x[1]), beg=sapply(coords, function(x) x[2]), end=sapply(coords, function(x) as.integer(x[3])-1), sd=sd.binmeanT), file='/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/sd.100kbmeanT.bed', quote=F, sep='\t', col.names=F, row.names=F)
## sortbed /secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/sd.100kbmeanT.bed >1
## mv 1 /secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/sd.100kbmeanT.bed

chrms <- c(paste0('chr',1:22),'chrX','chrY','chrM')
chrm <- sapply(strsplit(names(sd.binmeanT),"[:-]"), function(x) x[1])
beg <- sapply(strsplit(names(sd.binmeanT),"[:-]"), function(x) as.integer(x[2]))
end <- sapply(strsplit(names(sd.binmeanT),"[:-]"), function(x) as.integer(x[3]))
gr <- GRanges(seqnames=chrm, ranges=IRanges(beg, end), seqinfo=Seqinfo(chrms))
gr$sd.binmeanT <- sd.binmeanT
gr <- sort(gr)
ctxt <- gr
save(ctxt, file='/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/cgcontext.hg19/ctxt_w_sdbinmeanT.rda')
gr$PMD <- TRUE
gr$PMD[gr$sd.binmeanT < 0.125] <- FALSE
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/cgcontext.hg19/ctxt.rda')
gr <- mergeOv(ctxt, gr)
cnt <- table(interaction(gr$ctxt, gr$orph35class, gr$PMD))

## > cnt / sum(cnt)
## SCGS.[0,1).FALSE   SCGW.[0,1).FALSE   WCGW.[0,1).FALSE
##   0.03031861         0.07276171         0.04617708
## SCGS.[1,2).FALSE   SCGW.[1,2).FALSE   WCGW.[1,2).FALSE
##   0.03455364         0.06305955         0.02871867
## SCGS.[2,3).FALSE   SCGW.[2,3).FALSE   WCGW.[2,3).FALSE
##   0.02769898         0.04078148         0.01410548
## SCGS.[3,Inf).FALSE SCGW.[3,Inf).FALSE WCGW.[3,Inf).FALSE
##   0.06464961         0.05608584         0.01429013
## SCGS.[0,1).TRUE    SCGW.[0,1).TRUE    WCGW.[0,1).TRUE  
##   0.03674141         0.10137136         0.07725688         
## SCGS.[1,2).TRUE    SCGW.[1,2).TRUE    WCGW.[1,2).TRUE
##   0.03143857         0.06467612         0.03639626 
## SCGS.[2,3).TRUE    SCGW.[2,3).TRUE    WCGW.[2,3).TRUE   
##   0.02105720         0.03337730         0.01350964 
## SCGS.[3,Inf).TRUE  SCGW.[3,Inf).TRUE  WCGW.[3,Inf).TRUE
##   0.04111493         0.03917702         0.01068252

> sum(!gr$PMD) / length(gr)
[1] 0.4932008
> 2066416 / 13555530
## [1] 0.1524408
> 1235114 / 13191808
## [1] 0.09362735

> sum(ctxt$sd.binmeanT > 0.15) / length(ctxt)
## [1] 0.4862212
## where the 49% 0.15 comes from

##################################
## preprocess all betas of cpgs
##################################

library(GenomicRanges)
library(parallel)

setwd('~/projects/laird-secondary/2016_12_26_TCGA_WGBS')
chrms <- c(paste0('chr',1:22),'chrX','chrY','chrM')
allfns <- list.files('bed/','.bed.gz')
mclapply(allfns, function(fn) {
  gc()
  a <- read.table(gzfile(paste0('bed/', fn)), header=F, stringsAsFactors=F, sep='\t')
  a <- a[a$V1 %in% chrms,]
  gr <- GRanges(seqnames=a$V1, ranges=IRanges(a$V2+1,a$V3), seqinfo=Seqinfo(chrms))
  gr$betas <- a$V4
  gr <- sort(gr)
  save(gr, file=paste0('rda/', gsub('.bed.gz','',fn), '.rda'))
}, mc.cores=30)

##################
## germ line RRBS
##################
setwd('~/projects/laird-secondary/2017_01_25_early_dev_germline/2014-Nature-Guo-human-embryo/RRBS')
chrms <- c(paste0('chr',1:22),'chrX','chrY','chrM')
allfns <- list.files('bed/','.bed.gz')
mclapply(allfns, function(fn) {
  gc()
  a <- read.table(gzfile(paste0('bed/', fn)), header=F, stringsAsFactors=F, sep='\t')
  a <- a[a$V1 %in% chrms,]
  gr <- GRanges(seqnames=a$V1, ranges=IRanges(a$V2+1,a$V3), seqinfo=Seqinfo(chrms))
  gr$betas <- a$V4
  gr <- sort(gr)
  save(gr, file=paste0('rda/', gsub('.bed.gz','',fn), '.rda'))
}, mc.cores=30)

load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/annotations/pronecpg.rda')
setwd('~/projects/laird-secondary/2017_01_25_early_dev_germline/2014-Nature-Guo-human-embryo/RRBS')
samples <- list.files('rda')
pmd2betas <- lapply(samples, function(x) {
  load(paste0('rda/',x))
  gr <- mergeOv(gr, pronecpg)
  split(gr$betas, gr$commonPMD)
})
names(pmd2betas) <- samples
save(pmd2betas, file='pmd2betas.rda')
betas <- do.call(rbind, lapply(names(pmd2betas), function(sn) {s <- pmd2betas[[sn]]; s <- wzbind.list(s); s$sample <- sn; s}))
betas$sample <- sub('RRBS_','',betas$sample)
betas$sample <- factor(betas$sample, levels=c(
                                       'Sperm1','Sperm2','Sperm3','Sperm4',
                                       '1st_PB1','1st_PB2','1st_PB3',
                                       '2nd_PB1','2nd_PB2',
                                       'MII_Oocyte1','MII_Oocyte2',
                                       'Zygote1','Zygote2',
                                       '2-cell1','2-cell2',
                                       '4-cell1','4-cell2',
                                       '8-cell1','8-cell2','8-cell3',
                                       'Morula1','Morula2','Morula3',
                                       'ICM1','ICM2','ICM3',
                                       'TE1','TE2','TE3',
                                       'Postimplantation_embryo1','Postimplantation_embryo2','Postimplantation_embryo3'
  ))
pdf('~/gallery/2017_01_26_preimplantation.pdf', width=14, height=4.5)
ggplot(betas[betas$cat %in% c('comHMD','comPMD'), ], aes(sample, x, fill=cat)) + geom_violin(scale='area', draw_quantiles=0.5, bw=0.1, alpha=0.8, size=0.3) + theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) + scale_fill_manual(values=c('#CC0000','#005AA0')) + ylab('Beta Value')
dev.off()

#####################
## cpg context
#####################
library(GenomicRanges)
setwd('~/projects/laird-secondary/2016_12_26_TCGA_WGBS')
chrms <- c(paste0('chr',1:22),'chrX','chrY','chrM')
a <- read.table(gzfile('cgcontext.hg19/cgcontext.merged.contextclass.bed'), header=F, stringsAsFactors=F, sep='\t')
a <- a[a$V1 %in% chrms,]
a <- a[!grepl('N',a$V6),]
gr <- GRanges(seqnames=a$V1, ranges=IRanges(a$V2+1,a$V3), seqinfo=Seqinfo(chrms))
## gr$orph35 <- a$V4
## gr$gc35 <- a$V5
## gr$ctxt1 <- as.factor(a$V6)
## gr$ctxt2 <- as.factor(a$V7)
gr$ctxt <- as.factor(a$V8)
gr$orph35class <- cut(a$V4, breaks=c(0,1,2,3,Inf), right=F)
gr$gc35class <- cut(a$V5, breaks=c(0,0.35,0.5,0.65,1), include.lowest=T)
gr <- sort(gr)
ctxt <- gr
save(ctxt, file='/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/cgcontext.hg19/ctxt.rda')

## with extra flanking context
a <- read.table(gzfile('cgcontext.hg19/cgcontext.merged.contextclass.bed'), header=F, stringsAsFactors=F, sep='\t')
a <- a[a$V1 %in% chrms,]
a <- a[!grepl('N',a$V6),]
gr <- GRanges(seqnames=a$V1, ranges=IRanges(a$V2+1,a$V3), seqinfo=Seqinfo(chrms))
## gr$orph35 <- a$V4
## gr$gc35 <- a$V5
## gr$ctxt1 <- as.factor(a$V6)
## gr$ctxt2 <- as.factor(a$V7)
gr$ctxt <- as.factor(a$V8)
gr$orph35class <- cut(a$V4, breaks=c(0,1,2,3,Inf), right=F)
gr$gc35class <- cut(a$V5, breaks=c(0,0.35,0.5,0.65,1), include.lowest=T)
gr$ctxtextra <- as.factor(a$V6)
gr <- sort(gr)
ctxt <- gr
save(ctxt, file='/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/cgcontext.hg19/ctxtextra.rda')

## ulra-prone/solo-WCGW cpg
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/archived_Huy/2015_10_12_Huy_segments/commonPMDanno.rda')
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/cgcontext.hg19/ctxt.rda')
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/annotations/hg19.cpg.cgi.rda')
gr <- mergeOv(mergeOv(ctxt, common.pmd), cgi)
gr <- gr[gr$cgi == 'NonCGI' & gr$ctxt=='WCGW' & gr$orph35class=='[0,1)']
pronecpg <- gr
save(pronecpg, file='/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/annotations/pronecpg.rda')

#####################
## common PMD
#####################
library(GenomicRanges)
setwd('~/projects/laird-secondary/2016_12_26_TCGA_WGBS')
chrms <- c(paste0('chr',1:22),'chrX','chrY','chrM')
a <- read.table(gzfile('archived_Huy/2015_10_12_Huy_segments/commonPMDanno.bed'), header=F, stringsAsFactors=F, sep='\t')
a <- a[a$V1 %in% chrms,]
gr <- GRanges(seqnames=a$V1, ranges=IRanges(a$V2+1,a$V3), seqinfo=Seqinfo(chrms))
gr$commonPMD <- as.factor(a$V4)
gr <- sort(gr)
common.pmd <- gr
save(common.pmd, file='~/projects/laird-secondary/2016_12_26_TCGA_WGBS/archived_Huy/2015_10_12_Huy_segments/commonPMDanno.rda')

load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/archived_Huy/2015_10_12_Huy_segments/commonPMDanno.rda')
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/cgcontext.hg19/ctxt.rda')
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/annotations/hg19.cpg.cgi.rda')

## sample <- 'LUAD_6840'
## TCGAsamples <- list.files('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/rda/', 'TCGA.*.rda')
## for (fn in TCGAsmaples) {
##   sample <- gsub('.rda','',gsub('TCGA_','',fn))
##   load(paste0('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/rda/', fn))
##   cancer <- gr
##   cancer <- mergeOv(mergeOv(mergeOv(cancer, ctxt), common.pmd), cgi)
##   cancer <- cancer[cancer$cgi == 'NonCGI' & cancer$commonPMD!='Neither']
##   df <- as.data.frame(mcols(cancer))
##   pdf(sprintf('~/gallery/2016_12_28_context_%s.pdf', sample), width=7, height=2.5)
##   print(ggplot(aes(orph35class, betas, alpha=ctxt, fill=commonPMD, color='white'), data=) + ylab('Beta Value') + xlab('# Flanking CpGs') + annotate("rect", fill = "grey", alpha = 0.3, xmin = c(1,3)-0.5, xmax = c(1,3)+0.5, ymin = -Inf, ymax = Inf) + geom_violin(draw_quantiles=c(0.5), bw=0.05, scale='width')  + scale_color_manual(values=c('white')) + scale_fill_manual(values=c('#CC0000','#005AA0')) + scale_alpha_manual(values=c(0.4,0.6,1)) + scale_x_discrete(labels=c('0','1','2','3+')) + theme(text=element_text(size=18, face='bold'), axis.text.x = element_text(size=18), axis.text.y = element_text(size=17)) + ggtitle(sample))
##   dev.off()
## }

## just COAD
for (fn in c('methbase_Berman_2012_Human_ColonCancer.rda','TCGA_COAD_A00R.rda', 'TCGA_COAD_3158.rda', 'TCGA_COAD_N3158.rda', 'methbase_Berman_2012_Human_ColonNormal.rda')) {
  sample <- gsub('.rda','',gsub('TCGA_','',fn))
  load(paste0('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/rda/', fn))
  cancer <- gr
  cancer <- mergeOv(mergeOv(mergeOv(cancer, ctxt), common.pmd), cgi)
  cancer <- cancer[cancer$cgi == 'NonCGI' & cancer$commonPMD!='Neither']
  df <- as.data.frame(mcols(cancer))
  df$order <- interaction(df$ctxt, -as.integer(df$orph35class), df$commonPMD)
  order2mean <- tapply(df$betas, df$order, mean)
  df$mbetas <- order2mean[df$order]
  df$mbetas[df$mbetas > 0.8] <- 0.8
  df$mbetas[df$mbetas < 0.4] <- 0.4
  pdf(sprintf('~/gallery/2017_01_06_WGBS_violin_%s.pdf', sample), width=10, height=3.4)
  ## print(ggplot(aes(order, betas), data=df) + geom_violin(aes(fill=mbetas), draw_quantiles=c(0.5), bw=0.05, size=0.8) + theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)))
  print(ggplot(aes(order, betas), data=df) + geom_violin(aes(fill=mbetas), draw_quantiles=c(0.5), bw=0.05, size=0.8) + theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) + scale_fill_gradient2(high='#CC0000',low='#005AA0', mid='white', midpoint=0.6, limits=c(0.4,0.8)))
  dev.off()
}

## all WGBS samples
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/archived_Huy/2015_10_12_Huy_segments/commonPMDanno.rda')
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/cgcontext.hg19/ctxt.rda')
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/annotations/hg19.cpg.cgi.rda')
allsamples <- list.files('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/rda/', '.rda')
ctxt1 <- mergeOv(mergeOv(ctxt, common.pmd), cgi)
for (fn in allsamples) {
  cat(fn,'\n')
  sample <- gsub('.rda','',fn)
  load(paste0('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/rda/', fn))
  cancer <- gr
  cancer <- mergeOv(cancer, ctxt1)
  cancer <- cancer[cancer$cgi == 'NonCGI' & cancer$commonPMD!='Neither']
  df <- as.data.frame(mcols(cancer))
  df$orph35class <- factor(df$orph35class, rev(levels(df$orph35class)))
  df$order <- interaction(df$ctxt, df$orph35class, df$commonPMD)
  order2mean <- tapply(df$betas, df$order, mean)
  df$mbetas <- order2mean[df$order]
  df$mbetas[df$mbetas > 0.8] <- 0.8
  df$mbetas[df$mbetas < 0.4] <- 0.4
  png(sprintf('~/gallery/2017_02_27_WGBS_violin_all/%s.png', sample), width=1200, height=150)
  par(mar=c(5,5,1,1))
  print(ggplot(aes(order, betas), data=df) + geom_violin(aes(fill=mbetas), draw_quantiles=c(0.5), bw=0.05, size=0.8) + theme(axis.text.x = element_text(size=0,angle=90, vjust=1, hjust=1)) + scale_fill_gradient2(high='#CC0000',low='#005AA0', mid='white', midpoint=0.6, limits=c(0.4,0.8)) + xlab('') + ggtitle(sample) + theme(plot.title = element_text(size=20)))
  dev.off()
}

## just mean heatmap
allcat.mean <- mclapply(allsamples, function(fn) {
  cat(fn,'\n')
  sample <- gsub('.rda','',fn)
  load(paste0('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/rda/', fn))
  cancer <- gr
  cancer <- mergeOv(cancer, ctxt1)
  cancer <- cancer[cancer$cgi == 'NonCGI' & cancer$commonPMD!='Neither']
  df <- as.data.frame(mcols(cancer))
  df$orph35class <- factor(df$orph35class, rev(levels(df$orph35class)))
  df$order <- droplevels(interaction(df$ctxt, df$orph35class, df$commonPMD))
  order2mean <- tapply(df$betas, df$order, mean)
  order2mean
}, mc.cores=28)
allctxt.mean <- do.call(cbind, allcat.mean)
colnames(allctxt.mean) <- gsub('.rda','',allsamples)
save(allctxt.mean, file='/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/cgcontext.hg19/allWGBS.ctxt.rda')
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/cgcontext.hg19/allWGBS.ctxt.rda')
samples <- read.csv2('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/samplesheet.txt', header=T, sep='\t', row.names='samplename', stringsAsFactors=F)
pdf('~/gallery/2017_02_28_allsample_ctxt.pdf', width=80, height=10)
WHeatmap(allctxt.mean[, rownames(samples)], cmp=CMPar(dmin=0.3,dmax=0.9,stop.points=c('#005AA0','white','#CC0000')), yticklabels=T, xticklabels=T, xticklabels.n=ncol(allctxt.mean)) + WCustomize(mar.bottom=0.5)
dev.off()

## just mean heatmap, all 10 context
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/cgcontext.hg19/ctxtextra.rda')
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/archived_Huy/2015_10_12_Huy_segments/commonPMDanno.rda')
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/annotations/hg19.cpg.cgi.rda')
allsamples <- list.files('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/rda/', '.rda')
ctxt1 <- mergeOv(mergeOv(ctxt, common.pmd), cgi)
allcat.mean <- mclapply(allsamples, function(fn) {
  cat(fn,'\n')
  sample <- gsub('.rda','',fn)
  load(paste0('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/rda/', fn))
  cancer <- gr
  cancer <- mergeOv(cancer, ctxt1)
  cancer <- cancer[cancer$cgi == 'NonCGI' & cancer$commonPMD!='Neither']
  df <- as.data.frame(mcols(cancer))
  df$orph35class <- factor(df$orph35class, rev(levels(df$orph35class)))
  df$order <- droplevels(interaction(df$ctxtextra, df$orph35class, df$commonPMD))
  order2mean <- tapply(df$betas, df$order, mean)
  order2mean
}, mc.cores=28)
allctxt.mean <- do.call(cbind, allcat.mean)
colnames(allctxt.mean) <- gsub('.rda','',allsamples)
save(allctxt.mean, file='/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/cgcontext.hg19/allWGBS.ctxtextra.rda')
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/cgcontext.hg19/allWGBS.ctxtextra.rda')
samples <- read.csv2('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/samplesheet.txt', header=T, sep='\t', row.names='samplename', stringsAsFactors=F)
pdf('~/gallery/2017_03_05_allsample_ctxtextra.pdf', width=80, height=20)
WHeatmap(allctxt.mean[, rownames(samples)], cmp=CMPar(dmin=0.3,dmax=0.9,stop.points=c('#005AA0','white','#CC0000')), yticklabels=T, xticklabels=T, xticklabels.n=ncol(allctxt.mean)) + WCustomize(mar.bottom=0.5)
dev.off()

for (domain in c('comPMD','comHMD')) {
  for (ctxt in c('[0,1)','[1,2)','[2,3)','[3,Inf)')) {
    ctxtlv <- c('ACGT','ACGA','TCGA','ACGG','CCGA','ACGC','GCGA','CCGG','CCGC','GCGC')
    df <- allctxt.mean[grepl(paste0(ctxt,'.',domain), rownames(allctxt.mean), fixed=TRUE),]
    df <- melt(t(df))
    df$ctxt <- factor(substr(df$Var2,1,4), levels=ctxtlv)
    png(sprintf('~/gallery/2017_03_05_ctxt_%s_%s.png', domain, ctxt), width=1000, height=900, res=200)
    print(ggplot(df) + geom_line(aes(ctxt, value, group=Var1), size=0.1) + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) + geom_point(aes(ctxt, value), size=0.5) + xlab('Tetranucleotide Context') + ylab('Beta Value Average') + ggtitle(paste0(domain,'.',ctxt)))
    dev.off()
  }
}

## collating
## sort samplesheet by tile and then tile order
cd ~/gallery/2017_02_27_WGBS_violin_all
awk 'BEGIN{key="";a=""}NR>1{if (key!=$7 && key!="") {print "montage "a" -gravity South -tile 1x -geometry 20000x200+1 png:- | convert -trim - collated/"key".png"; a="";} key=$7; a=a" -label "key","$5","$1" "$1".png";}END{print "montage "a" -gravity South -tile 1x -geometry 20000x200+1 png:- | convert -trim - collated/"key".png";}' /secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/samplesheet.txt | while read f; do eval $f; done

###################
## PMD depth
###################
setwd('~/projects/laird-secondary/2016_12_26_TCGA_WGBS')
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/archived_Huy/2015_10_12_Huy_segments/commonPMDanno.rda')
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/cgcontext.hg19/ctxt.rda')
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/annotations/hg19.cpg.cgi.rda')
allfns <- list.files('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/rda/', '.rda')
depths <- mclapply(allfns, function(fn) {
  gc()
  load(paste0('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/rda/', fn))
  gr <- mergeOv(mergeOv(mergeOv(gr, ctxt), common.pmd), cgi)
  gr.pmd <- gr[gr$cgi == 'NonCGI' & gr$commonPMD=='comPMD' & gr$ctxt=='WCGW' & gr$orph35class=='[0,1)']
  gr.hmd <- gr[gr$cgi == 'NonCGI' & gr$commonPMD=='comHMD' & gr$ctxt=='WCGW' & gr$orph35class=='[0,1)']
  c(mean(gr.pmd$betas), mean(gr.hmd$betas))
}, mc.cores=10)
depths <- t(simplify2array(depths))
rownames(depths) <- gsub('.rda','',allfns)
colnames(depths) <- c('comPMD','comHMD')
write.table(depths, file='comPMDdepth.tsv', sep='\t')
save(depths, file='comPMDdepth.rda')

########################
## PMD depth vs purity
########################

setwd('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/')
depths <- read.table('comPMDdepth.tsv', sep='\t', header=T)
purity <- read.table('TCGA_WGBS_purity.tsv', sep='\t', header=T, row.names='sample')
pdf('~/gallery/2017_01_04_purity_vs_PMD_depth.pdf', width=4,height=3)
qplot(purity, depths, color=cancertype, data=data.frame(purity=purity$purity, depths=depths[rownames(purity),'x'], cancertype=sapply(strsplit(rownames(purity),'_'), function(x) x[2][[1]]))) + xlab('Purity') + ylab('common PMD depth')
dev.off()

##################
## GC content
##################
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/rda/TCGA_LUAD_7156.rda')
cancer <- gr
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/archived_Huy/2015_10_12_Huy_segments/commonPMDanno.rda')
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/cgcontext.hg19/ctxt.rda')

cancer <- mergeOv(mergeOv(cancer, ctxt), common.pmd)
cancer$k <- interaction(cancer$orph35class, cancer$gc35class)
cp <- cancer[cancer$commonPMD=='comPMD']
ch <- cancer[cancer$commonPMD=='comHMD']
cps <- as.data.frame(cp[cp$ctxt=='SCGS'])
cpw <- as.data.frame(cp[cp$ctxt=='SCGW'])
chs <- as.data.frame(ch[ch$ctxt=='SCGS'])
chw <- as.data.frame(ch[ch$ctxt=='SCGW'])
ggplot(data=cpw) + geom_violin(aes(orph35class, betas, fill=gc35class, linetype=NA), bw=0.05)
## ongoing

##################
## CGI
##################

a <- read.table('~/projects/laird-secondary/2016_12_26_TCGA_WGBS/annotations/hg19.cpg.cgi.bed',sep='\t',header=F,stringsAsFactors=F)
chrms <- c(paste0('chr',1:22),'chrX','chrY','chrM')
a <- a[a$V1 %in% chrms,]
cgi <- GRanges(seqnames=a$V1, ranges=IRanges(a$V2+1,a$V3), seqinfo=Seqinfo(chrms))
cgi$cgi <- as.factor(a$V4)
save(cgi, file='~/projects/laird-secondary/2016_12_26_TCGA_WGBS/annotations/hg19.cpg.cgi.rda')

##################
## centromere
##################

setwd('~/projects/laird-secondary/2016_12_26_TCGA_WGBS')
chrms <- c(paste0('chr',1:22),'chrX','chrY','chrM')
a <- read.table('~/references/hg19/annotation/hg19.centromere_3Mflank.bed',sep='\t',header=F,stringsAsFactors=F)
a <- a[a$V1 %in% chrms,]
centromere <- sort(GRanges(seqnames=a$V1, ranges=IRanges(a$V2+1,a$V3), seqinfo=Seqinfo(chrms)))
save(centromere, file='bin/centromere3Mflank.rda')

##################
## bins
##################

library(GenomicRanges)
setwd('~/projects/laird-secondary/2016_12_26_TCGA_WGBS')
load('bin/centromere3Mflank.rda')
chrms <- c(paste0('chr',1:22),'chrX','chrY','chrM')
a <- read.table('~/references/hg19/annotation/hg19.cpg.bed',sep='\t',header=F,stringsAsFactors=F)
a <- a[a$V1 %in% chrms,]
gr <- GRanges(seqnames=a$V1, ranges=IRanges(a$V2+1,a$V3), seqinfo=Seqinfo(chrms))
save(a, file='bin/allcpg.rda')

## multiscale bin-sizes
pws <- seq(2,7,by=0.1)
binsizes <- floor(10^pws)
lapply(seq_along(pws), function(i) {
  cat(pws[i],'\n')
  binsize <- binsizes[i]
  mcols(gr) <- NULL
  bin <- paste0(seqnames(gr), ':', sapply(floor(start(gr)/binsize), function(x) sprintf('%d-%d', x*binsize,(x+1)*binsize)))
  gr$bin <- factor(bin, levels=rle(bin)$values)
  gr$centromere <- lojOv(gr, centromere)
  bin <- gr
  save(bin, file=sprintf('bin/scale.%1.1f.hg19.cpg.bin.rda', pws[i]))
})

## manual bin-sizes
binsize <- 100000
mcols(gr) <- NULL
bin <- paste0(seqnames(gr), ':', sapply(floor(start(gr)/binsize), function(x) sprintf('%d-%d', x*binsize,(x+1)*binsize)))
gr$bin <- factor(bin, levels=rle(bin)$values)
gr$centromere <- lojOv(gr, centromere)
bin <- gr
save(bin, file='bin/hg19.cpg.100kb.bin.rda')

binsize <- 30000
mcols(gr) <- NULL
bin <- paste0(seqnames(gr), ':', sapply(floor(start(gr)/binsize), function(x) sprintf('%d-%d', x*binsize,(x+1)*binsize)))
gr$bin <- factor(bin, levels=rle(bin)$values)
gr$centromere <- lojOv(gr, centromere)
bin <- gr
save(bin, file='bin/hg19.cpg.30kb.bin.rda')

binsize <- 1000000
mcols(gr) <- NULL
bin <- paste0(seqnames(gr), ':', sapply(floor(start(gr)/binsize), function(x) sprintf('%d-%d', x*binsize,(x+1)*binsize)))
gr$bin <- factor(bin, levels=rle(bin)$values)
gr$centromere <- lojOv(gr, centromere)
bin <- gr
save(bin, file='bin/hg19.cpg.1mb.bin.rda')

binsize <- 10000
mcols(gr) <- NULL
bin <- paste0(seqnames(gr), ':', sapply(floor(start(gr)/binsize), function(x) sprintf('%d-%d', x*binsize,(x+1)*binsize)))
gr$bin <- factor(bin, levels=rle(bin)$values)
gr$centromere <- lojOv(gr, centromere)
bin <- gr
save(bin, file='bin/hg19.cpg.10kb.bin.rda')

binsize <- 1000
mcols(gr) <- NULL
bin <- paste0(seqnames(gr), ':', sapply(floor(start(gr)/binsize), function(x) sprintf('%d-%d', x*binsize,(x+1)*binsize)))
gr$bin <- factor(bin, levels=rle(bin)$values)
gr$centromere <- lojOv(gr, centromere)
bin <- gr
save(bin, file='bin/hg19.cpg.1kb.bin.rda')

binsize <- 10000000
mcols(gr) <- NULL
bin <- paste0(seqnames(gr), ':', sapply(floor(start(gr)/binsize), function(x) sprintf('%d-%d', x*binsize,(x+1)*binsize)))
gr$bin <- factor(bin, levels=rle(bin)$values)
gr$centromere <- lojOv(gr, centromere)
bin <- gr
save(bin, file='bin/hg19.cpg.10mb.bin.rda')

###########################
## 100kb bin mean context
###########################

load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/rda/TCGA_COAD_3158.rda')
cancer <- gr
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/archived_Huy/2015_10_12_Huy_segments/commonPMDanno.rda')
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/cgcontext.hg19/ctxt.rda')
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/bin/hg19.cpg.100kb.bin.rda')
bin100kb <- gr

cancer <- mergeOv(mergeOv(mergeOv(cancer, ctxt), common.pmd), bin100kb)
prone <- cancer[cancer$ctxt == 'WCGW' & cancer$orph35class == '[0,1)']
resis <- cancer[cancer$ctxt == 'SCGS' & cancer$orph35class != '[0,1)']

##########################
## all prone probe mean
##########################

## load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/cgcontext.hg19/ctxt.rda')
## ctxt <- gr
## setwd('~/projects/laird-secondary/2016_12_26_TCGA_WGBS')
## basedir <- '/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/rda'
## allfns <- list.files(basedir,'.rda')

## for (binlen in c('100kb','1mb','1kb','10kb','10mb')) {
##   load(sprintf('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/bin/hg19.cpg.%s.bin.rda', binlen))
##   bin <- gr
##   betas <- do.call(cbind, mclapply(allfns, function(fn) {
##     load(paste0(basedir,'/',fn))
##     gr <- mergeOv(mergeOv(gr, ctxt), bin)
##     gr <- gr[gr$ctxt=='WCGW' & gr$orph35class=='[0,1)']
##     bin.mean <- tapply(gr$betas, gr$bin, mean)
##     bin.cnt <- tapply(gr$betas, gr$bin, length)
##     thres <- min(5, quantile(bin.cnt, probs=c(0.01), na.rm=T))
##     message(cat(binlen,':\t',thres,'\n'))
##     bin.mean[bin.cnt<thres] <- NA
##     bin.mean
##   }, mc.cores=40))
##   colnames(betas) <- gsub('.rda','',allfns)
##   save(betas, file=sprintf('allprone/bin_mean_%s.rda', binlen))
## }

load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/cgcontext.hg19/ctxt.rda')
load('~/projects/laird-secondary/2016_12_26_TCGA_WGBS/annotations/hg19.cpg.cgi.rda')
setwd('~/projects/laird-secondary/2016_12_26_TCGA_WGBS')
basedir <- '/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/rda'
allfns <- list.files(basedir,'.rda')
thres <- 1                              # minimum number of cpg per bin
for (binlen in '100kb') { #c('100kb','1mb','1kb','10kb','10mb'))
  gc()
  load(sprintf('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/bin/hg19.cpg.%s.bin.rda', binlen))
  betas <- do.call(cbind, mclapply(allfns, function(fn) {
    load(paste0(basedir,'/',fn))
    gr <- mergeOv(mergeOv(mergeOv(gr, ctxt), bin), cgi)
    gr <- gr[gr$ctxt=='WCGW' & gr$orph35class=='[0,1)' & gr$cgi=='NonCGI']
    bin.mean <- tapply(gr$betas, gr$bin, mean)
    bin.cnt <- tapply(gr$betas, gr$bin, length)
    bin.mean[bin.cnt<thres] <- NA
    bin.mean
  }, mc.cores=20))
  colnames(betas) <- gsub('.rda','',allfns)
  save(betas, file=sprintf('allprone/bin_mean_%s.rda', binlen))
}

## all scale bin average
setwd('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS')
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/cgcontext.hg19/ctxt.rda')
load('~/projects/laird-secondary/2016_12_26_TCGA_WGBS/annotations/hg19.cpg.cgi.rda')
basedir <- '/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/rda'
pws <- seq(2,7,by=0.1)
allfns <- list.files(basedir,'.rda')
ctxt.cgi <- mergeOv(ctxt, cgi)
thres <- 1
for (pw in pws) {
  gc()
  cat(pw,'\n')
  load(sprintf('bin/scale.%1.1f.hg19.cpg.bin.rda', pw))
  bin.ctxt.cgi <- mergeOv(ctxt.cgi, bin)
  betas <- do.call(cbind, mclapply(allfns, function(fn) {
    cat(fn,'\n')
    load(paste0(basedir,'/',fn))
    gc()
    gr <- mergeOv(gr, bin.ctxt.cgi)
    gr <- gr[gr$ctxt=='WCGW' & gr$orph35class=='[0,1)' & gr$cgi=='NonCGI']
    gr <- gr[!gr$centromere,]
    bin.mean <- tapply(gr$betas, gr$bin, mean)
    bin.cnt <- tapply(gr$betas, gr$bin, length)
    bin.mean[bin.cnt<thres] <- NA
    gr <- NULL
    bin.mean
  }, mc.cores=20))
  colnames(betas) <- gsub('.rda','',allfns)
  save(betas, file=sprintf('allprone/multiscale/bin_mean_scale_%1.1f.rda', pw))
  betas <- NULL
}

# for pw in $(seq 3 0.1 7); do pbsgen one -ppn 28 -memG 240 -hour 24 -dest pbs/binmean_scale_$pw.pbs "Rscript ~/wzprojects/2017_02_07_binmean_scale.R $pw" -submit; done

######################################################################################
## mutli-scale (obsolete, since this is nonoverlapping bin, we need overlapping bin)
######################################################################################

load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/cgcontext.hg19/ctxt.rda')
load('~/projects/laird-secondary/2016_12_26_TCGA_WGBS/annotations/hg19.cpg.cgi.rda')
chrms <- c(paste0('chr',1:22),'chrX','chrY','chrM')
setwd('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/')
load('bin/allcpg.rda')
prone <- mergeOv(ctxt, cgi)
prone <- prone[prone$ctxt=='WCGW' & prone$orph35class=='[0,1)' & prone$cgi=='NonCGI']
binsizes <- c(1000, 3000, seq(10000,90000, by=10000), seq(100000,900000, by=100000), seq(1000000,9000000, by=1000000), 10000000, 30000000)
allpronebin <- mclapply(binsizes, function(binsize) {
  cat(binsize, '\n')
  gr <- GRanges(seqnames=a$V1, ranges=IRanges(a$V2+1,a$V3), seqinfo=Seqinfo(chrms))
  bin <- paste0(seqnames(gr), ':', sapply(floor(start(gr)/binsize), function(x) sprintf('%d-%d', x*binsize,(x+1)*binsize)))
  gr$bin <- factor(bin, levels=rle(bin)$values)
  bin <- gr
  pronebin <- mergeOv(bin, prone)
  pronebin
}, mc.cores=40)
options(scipen=999)
names(allpronebin) <- binsizes
save(allpronebin, file='bin/allpronebin2.rda')

load('bin/allpronebin2.rda')
allfns <- list.files('rda/','.rda')
for (i in seq_along(allpronebin)) {
  gc()
  cat(names(allpronebin)[i], '\n')
  pronebin <- allpronebin[[i]]
  betas <- NULL
  betas <- do.call(cbind, mclapply(allfns, function(fn) {
    gc()
    load(paste0('rda/',fn))
    gr <- mergeOv(gr, pronebin)
    bin.mean <- tapply(gr$betas, gr$bin, mean)
    bin.mean
  }, mc.cores=20))
  colnames(betas) <- gsub('.rda','',allfns)
  save(betas, file=sprintf('allprone/multiscale_%s.rda', names(allpronebin)[i]))
  NULL
}

## skip 1k, which is too small
binsizes <- c(3000, seq(10000,90000, by=10000), seq(100000,900000, by=100000), seq(1000000,9000000, by=1000000), 10000000, 30000000)
allbetas <- lapply(binsizes, function(binsize) {
  load(sprintf('allprone/multiscale_%d.rda', binsize))
  betas
})
save(allbetas, file='allprone/mutliscale_allbetas.rda')

smallestbin <- rownames(allbetas[[1]])
size2smallbin <- mclapply(binsizes, function(binsize) {
  smallbinname <- sapply(smallestbin, function(bb) {
    binsplit <- strsplit(bb, '[:-]')[[1]]
    bid <- floor(as.integer(binsplit[2]) / binsize) * binsize
    sprintf('%s:%d-%d', binsplit[1], bid, bid + binsize)
  })
}, mc.cores=20)
save(size2smallbin, file='bin/size2smallbin.rda')

## all mean and all sd in multi-scale
betas <- allbetas[[1]]
alltcga.samples <- colnames(betas)[grep('TCGA_', colnames(betas))]
tumorsamples <- alltcga.samples[!grepl('_N', alltcga.samples)]
normalsamples <- alltcga.samples[grepl('_N', alltcga.samples)]

allmean <- lapply(allbetas, function(betas) apply(betas[,tumorsamples], 1, mean, na.rm=T))
allmean0 <- simplify2array(lapply(seq_along(allmean), function(i) allmean[[i]][size2smallbin[[i]]]))
options(scipen=999)
colnames(allmean0) <- as.character(binsizes)
allmean <- allmean0[!is.na(allmean0[,1]),]
save(allmean, file='allprone/mutliscale_allmean.rda')

allsd <- lapply(allbetas, function(betas) apply(betas[,tumorsamples], 1, sd, na.rm=T))
allsd0 <- simplify2array(lapply(seq_along(allsd), function(i) allsd[[i]][size2smallbin[[i]]]))
options(scipen=999)
colnames(allsd0) <- as.character(binsizes)
allsd <- allsd0[!is.na(allsd0[,1]),]
save(allsd, file='allprone/mutliscale_allsd.rda')

## plot multi-scale
setwd('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS')
load('allprone/mutliscale_allsd.rda')
chr16short <- rownames(allsd)[grep('chr16',rownames(allsd))]
chr16short <- chr16short[sapply(strsplit(chr16short,"[:-]"), function(x) as.integer(x[2])) < 35000000]
WHeatmap(t(allsd[chr16short,]), cmp=CMPar(stop.points=c('blue','yellow'), dmin=0, dmax=1), yticklabels=T)

load('allprone/mutliscale_allmean.rda')
chr16short <- rownames(allmean)[grep('chr16',rownames(allmean))]
chr16short <- chr16short[sapply(strsplit(chr16short,"[:-]"), function(x) as.integer(x[2])) < 35000000]
WHeatmap(t(allmean[chr16short,]), cmp=CMPar(stop.points=c('blue','yellow'), dmin=0, dmax=1), yticklabels=T)

################
## multiscale
################

## bed to rda
setwd('~/projects/laird-secondary/2016_12_26_TCGA_WGBS')
pws <- seq(4,7,0.1)
trash <- mclapply(pws, function(pw) {
  allfns <- list.files('windows/bed/', sprintf('scale_%1.1f_.*.bed', pw))
  betas <- do.call(cbind, lapply(allfns, function(fn) {
    a <- read.table(sprintf('windows/bed/%s', fn), stringsAsFactors=F)
    b <- setNames(as.numeric(a$V4), paste0(a$V1,':',a$V2))
    b
  }))
  colnames(betas) <- sub('.bed','',sub(sprintf('scale_%1.1f_', pw),'',allfns))
  save(betas, file=sprintf('windows/rda/betas_scale_%1.1f.rda', pw))
}, mc.cores=32)

## multi-scale sd and mean for both tumor and normal
setwd('~/projects/laird-secondary/2016_12_26_TCGA_WGBS')
load('allprone/bin_mean_100kb.rda')
alltcga.samples <- colnames(betas)[grep('TCGA_', colnames(betas))]
tumorsamples <- alltcga.samples[!grepl('_N', alltcga.samples)]
normalsamples <- alltcga.samples[grepl('_N', alltcga.samples)]

setwd('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS')
options(scipen=999)
pws <- seq(4,7,0.1)
mean.binmeanT <- do.call(cbind, mclapply(pws, function(pw) {
  load(sprintf('windows/rda/betas_scale_%1.1f.rda', pw))
  apply(betas[,tumorsamples], 1, mean, na.rm=T)
}, mc.cores=20))
colnames(mean.binmeanT) <- sprintf('%1.1f',pws)

sd.binmeanT <- do.call(cbind, mclapply(pws, function(pw) {
  load(sprintf('windows/rda/betas_scale_%1.1f.rda', pw))
  apply(betas[,tumorsamples], 1, sd, na.rm=T)
}, mc.cores=20))
colnames(sd.binmeanT) <- sprintf('%1.1f',pws)

mean.binmeanN <- do.call(cbind, mclapply(pws, function(pw) {
  load(sprintf('windows/rda/betas_scale_%1.1f.rda', pw))
  apply(betas[,normalsamples], 1, mean, na.rm=T)
}, mc.cores=20))
colnames(mean.binmeanN) <- sprintf('%1.1f',pws)

sd.binmeanN <- do.call(cbind, mclapply(pws, function(pw) {
  load(sprintf('windows/rda/betas_scale_%1.1f.rda', pw))
  apply(betas[,normalsamples], 1, sd, na.rm=T)
}, mc.cores=20))
colnames(sd.binmeanN) <- sprintf('%1.1f',pws)

save(mean.binmeanT, sd.binmeanT, mean.binmeanN, sd.binmeanN, file='/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/windows/binmean.rda')

## plot
setwd('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS')
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/windows/binmean.rda')
chr16 <- rownames(sd.binmeanT)[grep('chr16', rownames(sd.binmeanT))]
chr16short <- chr16[sapply(strsplit(chr16,"[:-]"), function(x) as.integer(x[2])) < 35000000]

dt <- t(sd.binmeanT[chr16short,ncol(sd.binmeanT):1])
rownames(dt) <- sprintf('%1.0f',10^as.numeric(rownames(dt)))
pdf('~/gallery/2017_01_10_multiscale_legend.pdf', width=5, height=5)
WHeatmap(dt[,1:3], cmp=CMPar(stop.points=c('blue','yellow'), dmin=0, dmax=0.2), yticklabels=T, name='a', yticklabels.n=nrow(dt)) + WLegendV('a', RightOf()) + WCustomize(mar.left=0.1)
dev.off()

pdf('~/gallery/2017_01_10_multiscale_legend2.pdf', width=5, height=5)
WHeatmap(dt[,1:3], cmp=CMPar(stop.points=c('#91bfdb','#ffffbf','#fc8d59'), dmin=0, dmax=0.2), yticklabels=T, name='a', yticklabels.n=nrow(dt)) + WLegendV('a', RightOf()) + WCustomize(mar.left=0.1)
dev.off()

dt <- t(sd.binmeanT[chr16short,ncol(sd.binmeanT):1])
rownames(dt) <- sprintf('%1.0f',10^as.numeric(rownames(dt)))
png('~/gallery/2017_01_10_multiscale_sd_tumor.png', width=800, height=200, res=1000)
WHeatmap(dt, cmp=CMPar(stop.points=c('#91bfdb','#ffffbf','#fc8d59'), dmin=0.0, dmax=0.2), name='a') + WCustomize(mar.left=0, mar.top=0, mar.right=0, mar.bottom=0)
dev.off()

dt <- t(sd.binmeanN[chr16short,ncol(sd.binmeanN):1])
rownames(dt) <- sprintf('%1.0f',10^as.numeric(rownames(dt)))
png('~/gallery/2017_01_10_multiscale_sd_normal.png', width=800, height=200, res=1000)
WHeatmap(dt, cmp=CMPar(stop.points=c('#91bfdb','#ffffbf','#fc8d59'), dmin=0, dmax=0.1), name='a') + WCustomize(mar.left=0, mar.top=0, mar.right=0, mar.bottom=0)
dev.off()

dt <- t(mean.binmeanT[chr16short,ncol(mean.binmeanT):1])
rownames(dt) <- sprintf('%1.0f',10^as.numeric(rownames(dt)))
png('~/gallery/2017_01_10_multiscale_mean_tumor.png', width=800, height=200, res=1000)
WHeatmap(dt, cmp=CMPar(stop.points=c('blue','yellow'), dmin=0, dmax=1), name='a') + WCustomize(mar.left=0, mar.top=0, mar.right=0, mar.bottom=0)
dev.off()

dt <- t(mean.binmeanN[chr16short,ncol(mean.binmeanN):1])
rownames(dt) <- sprintf('%1.0f',10^as.numeric(rownames(dt)))
png('~/gallery/2017_01_10_multiscale_mean_normal.png', width=800, height=200, res=1000)
WHeatmap(dt, cmp=CMPar(stop.points=c('blue','yellow'), dmin=0, dmax=1), name='a') + WCustomize(mar.left=0, mar.top=0, mar.right=0, mar.bottom=0)
dev.off()

## more multi-scale sd and mean
samples <- read.csv2('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/samplesheet.txt', header=T, sep='\t', row.names='samplename', stringsAsFactors=F)
samples <- list(
  'cellline'=samples[samples$tumor=='CL',],
  'fetal'=samples[samples$isFetal == 'Fetal' & samples$ok == 'OK',], # no germcells and PGC
  'normal'=samples[samples$plotorder >= 200 & samples$plotorder <= 300,],
  'blood'=samples[samples$plotorder >= 300 & samples$plotorder <= 400,],
  'tumor'=samples[samples$plotorder >= 500,],
  'bloodtumor'=samples[samples$plotorder >= 500 & !grepl('TCGA_',rownames(samples)),],
  'TCGAtumor'=samples[samples$plotorder >= 500 & samples$source %in% c('methbase','TCGA'),])

setwd('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS')
options(scipen=999)
pws <- seq(4,7,0.1)
mean.binmeanExtra <- list()
for (cname in names(samples)) {
  cat(cname, dim(samples[[cname]]), '\n')
  mean.binmeanExtra[[cname]] <- do.call(cbind, mclapply(pws, function(pw) {
    load(sprintf('windows/rda/betas_scale_%1.1f.rda', pw))
    apply(betas[,rownames(samples[[cname]])], 1, mean, na.rm=T)
  }, mc.cores=20))
  colnames(mean.binmeanExtra[[cname]]) <- sprintf('%1.1f',pws)
}
sd.binmeanExtra <- list()
for (cname in names(samples)) {
  cat(cname, dim(samples[[cname]]), '\n')
  sd.binmeanExtra[[cname]] <- do.call(cbind, mclapply(pws, function(pw) {
    load(sprintf('windows/rda/betas_scale_%1.1f.rda', pw))
    apply(betas[,rownames(samples[[cname]])], 1, sd, na.rm=T)
  }, mc.cores=20))
  colnames(sd.binmeanExtra[[cname]]) <- sprintf('%1.1f',pws)
}
save(mean.binmeanExtra, sd.binmeanExtra, file='/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/windows/binmeanExtra.rda')

load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/windows/binmeanExtra.rda')
betas1 <- mean.binmeanExtra[[1]]
chr16short <- rownames(betas1)[grep('chr16',rownames(betas1))]
chr16short <- chr16short[sapply(strsplit(chr16short,"[:-]"), function(x) as.integer(x[2])) < 35000000]
for (cname in names(mean.binmeanExtra)) {
  dt <- t(mean.binmeanExtra[[cname]][chr16short,ncol(mean.binmeanExtra[[cname]]):1])
  rownames(dt) <- sprintf('%1.0f',10^as.numeric(rownames(dt)))
  png(sprintf('~/gallery/2017_01_10_multiscale_mean_%s.png', cname), width=800, height=200, res=1000)
  print(WHeatmap(dt, cmp=CMPar(stop.points=c('blue','yellow'), dmin=0, dmax=1), name='a') + WCustomize(mar.left=0, mar.top=0, mar.right=0, mar.bottom=0))
  dev.off()
}
sd.limits <- list(
  'cellline'=c(0,0.3),
  'fetal'=c(0,0.09),
  'normal'=c(0,0.15),
  'blood'=c(0,0.11),
  'tumor'=c(0,0.21),
  'bloodtumor'=c(0,0.21))
for (cname in names(sd.binmeanExtra)) {
  dt <- t(sd.binmeanExtra[[cname]][chr16short,ncol(sd.binmeanExtra[[cname]]):1])
  rownames(dt) <- sprintf('%1.0f',10^as.numeric(rownames(dt)))
  png(sprintf('~/gallery/2017_01_10_multiscale_sd_%s.png', cname), width=800, height=200, res=1000)
  print(WHeatmap(dt, cmp=CMPar(stop.points=c('#91bfdb','#ffffbf','#fc8d59'), dmin=sd.limits[[cname]][1], dmax=sd.limits[[cname]][2]), name='a') + WCustomize(mar.left=0, mar.top=0, mar.right=0, mar.bottom=0))
  dev.off()
}

## 100 kb bin SD
for (cname in names(sd.binmeanExtra)) {
  pdf(sprintf('~/gallery/2017_01_10_100kb_sd_%s.pdf', cname), width=7, height=1)
  par(mar=c(1,1,1,1))
  barplot(sd.binmeanExtra[[cname]][chr16short,'5.0'], ylim=sd.limits[[cname]])
  dev.off()
}

## tumor vs others, composite heatmap
library(RColorBrewer)
library(MASS)
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
sd.limits <- list(
  'cellline'=c(0,0.3),
  'fetal'=c(0,0.12),
  'normal'=c(0,0.18),
  'blood'=c(0,0.15),
  'tumor'=c(0,0.25),
  'bloodtumor'=c(0,0.25),
  'TCGAtumor'=c(0,0.25))
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/windows/binmeanExtra.rda')
options(stringsAsFactors=T)             # IMPORTANT FOR stat_density2d TO WORK!!!
# 100kb
for (cname in names(sd.binmeanExtra)) {
  df <- data.frame(x=sd.binmeanExtra[[cname]][,'5.0'], y=sd.binmeanExtra[['TCGAtumor']][,'5.0'])
  df <- df[sapply(strsplit(rownames(df), ':'), function(x) as.integer(x[2])) %% 100000 == 0,]
  df <- df[complete.cases(df),]
  cat('\n\n', cname,'\n')
  print(cor.test(df$x, df$y))
  ## pdf(sprintf('~/gallery/2017_03_02_100kb_%s_mean_sd_2d_1.pdf', cname), width=5.5, height=4)
  ## print(ggplot(df[df$x<0.25&df$y<0.25,], aes(x,y)) + stat_density2d(aes(fill=..level..), n=50, geom='polygon') + scale_fill_gradientn(colors=jet.colors(50)) + ylim(0,0.3) + xlim(sd.limits[[cname]][1], sd.limits[[cname]][2]) + theme(panel.background=element_rect(fill=jet.colors(50)[1])) + xlab(sprintf('Standard Deviation (%s)', cname)) + ylab('Standard Deviation (TCGA Tumor)') + theme(text=element_text(size=18, face='bold')))
  ## dev.off()

  ## pdf(sprintf('~/gallery/2017_03_02_100kb_%s_mean_sd_2d_2.pdf', cname), width=4, height=1.5)
  ## print(ggplot(data.frame(x=df$x)) + geom_histogram(aes(x), bins=50, fill='#CC0000') + xlab('Standard Deviation') + xlim(sd.limits[[cname]][1], sd.limits[[cname]][2]) + theme(text=element_text(size=18, face='bold')))
  ## dev.off()
}

# 1mb
for (cname in names(sd.binmeanExtra)) {
  df <- data.frame(x=sd.binmeanExtra[[cname]][,'6.0'], y=sd.binmeanExtra[['TCGAtumor']][,'6.0'])
  df <- df[sapply(strsplit(rownames(df), ':'), function(x) as.integer(x[2])) %% 1000000 == 0,]
  df <- df[complete.cases(df),]
  pdf(sprintf('~/gallery/2017_03_02_1mb_%s_mean_sd_2d_1.pdf', cname), width=5.5, height=4)
  print(ggplot(df[df$x<0.25&df$y<0.25,], aes(x,y)) + stat_density2d(aes(fill=..level..), n=50, geom='polygon') + scale_fill_gradientn(colors=jet.colors(50)) + ylim(0,0.3) + xlim(sd.limits[[cname]][1], sd.limits[[cname]][2]) + theme(panel.background=element_rect(fill=jet.colors(50)[1])) + xlab(sprintf('Standard Deviation (%s)', cname)) + ylab('Standard Deviation (TCGA Tumor)') + theme(text=element_text(size=18, face='bold')))
  dev.off()

  pdf(sprintf('~/gallery/2017_03_02_1mb_%s_mean_sd_2d_2.pdf', cname), width=4, height=1.5)
  print(ggplot(data.frame(x=df$x)) + geom_histogram(aes(x), bins=50, fill='#CC0000') + xlab('Standard Deviation') + xlim(sd.limits[[cname]][1], sd.limits[[cname]][2]) + theme(text=element_text(size=18, face='bold')))
  dev.off()
}

## interaction count, normal: 0.083, tumor: 0.125
cname <- 'normal'
df <- data.frame(x=sd.binmeanExtra[[cname]][,'6.0'], y=sd.binmeanExtra[['TCGAtumor']][,'6.0'])
df <- df[sapply(strsplit(rownames(df), ':'), function(x) as.integer(x[2])) %% 1000000 == 0,]
df <- df[complete.cases(df),]
a <- interaction(df$x > 0.083, df$y > 0.125)
table(a)
table(a) / length(a)
aa <- table(a)
aa[4] / (aa[4] + aa[2] + aa[3])

## FALSE.FALSE  TRUE.FALSE  FALSE.TRUE   TRUE.TRUE
##       8657        1931        1354       15048

## FALSE.FALSE  TRUE.FALSE  FALSE.TRUE   TRUE.TRUE
## 0.32074843  0.07154502  0.05016673  0.55753983

## aa[4] / (aa[4] + aa[2] + aa[3])
## 0.8208149

#############################
## sample-based multi-scale
#############################

pws <- seq(4,7,0.1)
setwd('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS')
allbetas <- lapply(pws, function(pw) {
  cat(pw, '\n')
  load(sprintf('windows/rda/betas_scale_%1.1f.rda', pw))
  betas
})

samples <- colnames(allbetas[[1]])
trash <- mclapply(samples, function(sample) {
  cat(sample, '\n')
  b1 <- do.call(cbind, lapply(allbetas, function(betas) {
    betas[,sample]
  }))
  colnames(b1) <- sprintf('%1.1f', pws)
  betas <- b1
  save(betas, file=sprintf('windows/samplebased/%s.rda', sample))
}, mc.cores=30)

## plot brain

brains <- c('Ziller_515_fetal_brain', 'methbase_Brain_Human_FetalCerebCortex', 'TCGA_GBM_1788', 'TCGA_GBM_1460', 'TCGA_GBM_0128', 'TCGA_GBM_1454', 'TCGA_GBM_3477', 'TCGA_GBM_1401', 'methbase_Brain_Human_DorsPrefrontNeuron53Yr', 'methbase_Brain_Human_DorsPrefrontNeuronMale55Yr', 'Ziller_1020_Neuroepithelial_236461', 'Ziller_1019_Neuroectoderm_236460', 'Ziller_GSM916050_hippocampus_middle_292', 'Ziller_GSM1204462_1209_mid_frontal_cortex_281166', 'Ziller_1021_Radial_Glia_VZ_236462', 'Ziller_GSM1204461_1208_mid_frontal_cortex_281165', 'Ziller_GSM1204460_1205_mid_frontal_cortex_281162', 'Ziller_GSM1204459_1204_mid_frontal_cortex_281161', 'REMC_E071_Brain_Hippocampus_Middle', 'methbase_Zeng_2012_Human_PreFrontCortex', 'methbase_Brain_Human_DorsPrefrontNonNeuron53Yr', 'Ziller_246_hippocampus_middle', 'methbase_Brain_Human_DorsPrefrontNonNeuronMale55Yr', 'methbase_Brain_Human_MidFrontGyr16Yr', 'REMC_E070_Brain_Germinal_Matrix', 'methbase_Brain_Human_MidFrontGyr12Yr', 'Ziller_483_substantia_nigra', 'methbase_Brain_Human_MidFrontGyr5Yr', 'methbase_Brain_Human_MidFrontGyr35Day', 'methbase_Brain_Human_MidFrontGyr2Yr', 'methbase_Brain_Human_FrontCortexFemale64Yr')

setwd('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS')

for (brain in brains) {
  load(sprintf('windows/samplebased/%s.rda', brain))
  chr16 <- rownames(betas)[grep('chr16', rownames(betas))]
  chr16short <- chr16[sapply(strsplit(chr16,"[:-]"), function(x) as.integer(x[2])) < 35000000]

  dt <- t(betas[chr16short,])
  rownames(dt) <- sprintf('%1.0f',10^as.numeric(rownames(dt)))
  png(sprintf('~/gallery/2017_01_10_brain_multiscale_%s.png', brain), width=800, height=200, res=1000)
  print(WHeatmap(dt, cmp=CMPar(stop.points=c('blue','yellow'), dmin=0, dmax=1), name='a') + WCustomize(mar.left=0, mar.top=0, mar.right=0, mar.bottom=0) + WLabel(brain, BottomRightOf('a', just=c('right','bottom')), color='red', fontsize=1))
  dev.off()
}

## just chromosome 16 short arm
setwd('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS')
allfns <- list.files('windows/samplebased/','.rda')
for (sfn in allfns) {
  sname <- sub('.rda','',sfn)
  load(sprintf('windows/samplebased/%s', sfn))
  chr16 <- rownames(betas)[grep('chr16', rownames(betas))]
  target <- chr16[sapply(strsplit(chr16,"[:-]"), function(x) as.integer(x[2])) < 35000000]
  cat(sname, '\n')
  dt <- t(betas[target, ncol(betas):1])
  rownames(dt) <- sprintf('%1.0f',10^as.numeric(rownames(dt)))
  dir.create(sprintf('~/gallery/2017_01_13_allsample_multiscale/chr16short'), showWarnings = FALSE)
  png(sprintf('~/gallery/2017_01_13_allsample_multiscale/chr16short/%s.png', sname), width=1000, height=200, res=1000)
  print(WHeatmap(dt, cmp=CMPar(stop.points=c('blue','yellow'), dmin=0, dmax=1), name='a') + WCustomize(mar.left=0, mar.top=0, mar.right=0, mar.bottom=0))
  ##+ WLabel(paste(strsplit(gsub("([[:alnum:]_-]{40})", "\\1 ", sname), " ")[[1]], collapse='\n'), BottomRightOf('a', just=c('right','bottom')), color='red', fontsize=2))
  dev.off()
}


samples <- c('GSE63818_EmbryonicHeart_lowcov','GSE63818_EmbryonicBrain_lowcov','GSE49828_EmbryonicLiver_lowcov','GSE63818_PGC_10W_F','GSE63818_PGC_10W_M','GSE63818_PGC_11_to_13W_M','GSE63818_PGC_11W_F','GSE63818_PGC_17W_F','GSE63818_PGC_19W_M','GSE63818_PGC_7W_M','GSE63818_FetalGonadalSoma_7to9W','GSE63393_Fetal_PGC_F_D113','GSE63393_Fetal_PGC_F_D57to67','GSE63393_Fetal_PGC_M_D59to137_lowcov','Okae_2014_oocytes_lowcov','GSE49828_ICM_lowcov','Okae_2014_blastocyst')
allfns <- paste0(samples,'.rda')

## all sample multiscale plot
setwd('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS')
allfns <- list.files('windows/samplebased/','.rda')
chrms <- c(paste0('chr',1:2))
chr2len <- read.table('~/references/hg19/hg19.fa.fai', sep='\t', row.names=1)[,1,drop=FALSE]
for (sfn in allfns) {
  sname <- sub('.rda','',sfn)
  load(sprintf('windows/samplebased/%s', sfn))
  cat(sname, '\n')
  gc()
  for (chrm in chrms) {
    target <- rownames(betas)[grep(paste0(chrm,":"), rownames(betas))]
    ## chr16short <- chr16[sapply(strsplit(chr16,"[:-]"), function(x) as.integer(x[2])) < 35000000]
    dt <- t(betas[target, ncol(betas):1])
    rownames(dt) <- sprintf('%1.0f',10^as.numeric(rownames(dt)))
    dir.create(sprintf('~/gallery/2017_01_13_allsample_multiscale/%s', chrm), showWarnings = FALSE)
    png(sprintf('~/gallery/2017_01_13_allsample_multiscale/%s/%s.png', chrm, sname), width=2000*chr2len[chrm,]/1e8, height=400, res=2000*chr2len[chrm,]/1e8)
    print(WHeatmap(dt, cmp=CMPar(stop.points=c('blue','yellow'), dmin=0, dmax=1), name='a') + WCustomize(mar.left=0, mar.top=0, mar.right=0, mar.bottom=0))
    ##+ WLabel(paste(strsplit(gsub("([[:alnum:]_-]{40})", "\\1 ", sname), " ")[[1]], collapse='\n'), BottomRightOf('a', just=c('right','bottom')), color='red', fontsize=2))
    dev.off()
  }
}

## imagemagick, montage all
## montage BP_*.png -tile 1x -geometry 800x100+1+1 tile.png
mkdir -p merged/chr16short
awk -v chrm="chr16short" 'BEGIN{key="";a=""}NR>1{if (key!=$6 && key!="") {print "montage "a" -transpose -gravity South -tile x1 -geometry 100x800+1 png:- | convert -trim - merged/"chrm"/"key".png"; a="";} key=$6; a=a" -label "$4" "chrm"/"$1".png";}END{print "montage "a" -transpose -gravity South -tile x1 -geometry 100x800+1 png:- | convert -trim - merged/"chrm"/"key".png"}' /secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/samplesheet.txt | while read f; do eval $f; done

## resize
cd ~/Dropbox/presentation/gallery/2017-01/raw/2017_01_13_allsample_multiscale/resize
parallel -j 8 'mkdir -p {}; for f in ../{}/*.png; do convert $f -resize 3000x90\! {}/$(basename $f); done' ::: chr{1..22} "chr16short"

## montage
## for chrm in chr{1..22} "chr16short"; do
doawk() { awk -v chrm=$1 'BEGIN{key="";a=""}NR>1{if (key!=$6 && key!="") {print "montage "a" -gravity South -tile 1x -geometry 3000x30+1 png:- | convert -trim - merged/"chrm"/"key".png"; a="";} key=$6; a=a" -label "key","$4","$1" "chrm"/"$1".png";}END{print "montage "a" -gravity South -tile 1x -geometry 3000x30+1 png:- | convert -trim - merged/"chrm"/"key".png";}' $2; }
export -f doawk
parallel -j 8 'mkdir -p merged/{}; doawk {} /secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/samplesheet.txt | while read f; do eval $f; done;' ::: chr{1..22} chr16short

## montage all deep PMD
doawk() { awk -v chrm=$1 'BEGIN{key="";delete all[0];a="";}(NR>1&&$2<0.5){if(key!=$6 && key!="") {all[length(all)+1]=a; a="";} key=$6; a=a" -label "key","$4","$1" "chrm"/"$1".png";}END{all[length(all)+1]=a; printf("montage "); for(i=1;i<=length(all);++i) printf("%s ", all[i]); printf(" -gravity South -tile 1x -geometry 3000x90+1 png:- | convert -trim - merged/%s/deepPMD.png\n", chrm)}' $2; }
export -f doawk
parallel -j 8 'doawk {} /secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/samplesheet.txt | while read f; do eval "$f"; done;' ::: chr16short chr{1..22}

## example output command
## montage -label Imflammatory chr16short/BP_cord_blood_S0018A52_inflammatory_macrophage_DiseaseFree.png -label Imflammatory chr16short/BP_cord_blood_S007SK51_inflammatory_macrophage_DiseaseFree.png -label Imflammatory chr16short/BP_venous_blood_S0022I53_inflammatory_macrophage_DiseaseFree.png -label Imflammatory chr16short/BP_venous_blood_S001S753_inflammatory_macrophage_DiseaseFree.png -label Imflammatory chr16short/BP_venous_blood_S00H6O51_inflammatory_macrophage_DiseaseFree.png -label Imflammatory chr16short/BP_venous_blood_S001MJ51_inflammatory_macrophage_DiseaseFree.png -label AlternativelyActivated chr16short/BP_cord_blood_S00C1H51_alternatively_activated_macrophage_DiseaseFree.png -label AlternativelyActivated chr16short/BP_venous_blood_S00BS451_alternatively_activated_macrophage_DiseaseFree.png -label AlternativelyActivated chr16short/BP_venous_blood_S00FTN51_alternatively_activated_macrophage_DiseaseFree.png -label AlternativelyActivated chr16short/BP_cord_blood_S00E8W51_alternatively_activated_macrophage_DiseaseFree.png -label AlternativelyActivated chr16short/BP_venous_blood_S0062252_alternatively_activated_macrophage_DiseaseFree.png -label AlternativelyActivated chr16short/BP_venous_blood_S006VI53_alternatively_activated_macrophage_DiseaseFree.png -label macrophage chr16short/BP_venous_blood_S00V3BN1_macrophage_DiseaseFree.png -label macrophage chr16short/BP_venous_blood_S00V49N1_macrophage_DiseaseFree.png -label macrophage chr16short/BP_venous_blood_C005VG51_macrophage_DiseaseFree.png -label macrophage chr16short/BP_cord_blood_S00BHQ51_macrophage_DiseaseFree.png -label macrophage chr16short/BP_venous_blood_S0022I51_macrophage_DiseaseFree.png -label macrophage chr16short/BP_venous_blood_S001S751_macrophage_DiseaseFree.png -label macrophage chr16short/BP_cord_blood_S00DVR51_macrophage_DiseaseFree.png -label macrophage chr16short/BP_venous_blood_S0039051_macrophage_DiseaseFree.png -transpose -gravity South -tile x1 -geometry 100x800+1 png:- | convert -trim - chr16short/merged/macrophage.png

## montage -label Inflammatory chr16short/BP_cord_blood_S0018A52_inflammatory_macrophage_DiseaseFree.png -label Inflammatory chr16short/BP_cord_blood_S007SK51_inflammatory_macrophage_DiseaseFree.png -label Inflammatory chr16short/BP_venous_blood_S0022I53_inflammatory_macrophage_DiseaseFree.png -label Inflammatory chr16short/BP_venous_blood_S001S753_inflammatory_macrophage_DiseaseFree.png -label Inflammatory chr16short/BP_venous_blood_S00H6O51_inflammatory_macrophage_DiseaseFree.png -label Inflammatory chr16short/BP_venous_blood_S001MJ51_inflammatory_macrophage_DiseaseFree.png -label AlternativelyActivated chr16short/BP_cord_blood_S00C1H51_alternatively_activated_macrophage_DiseaseFree.png -label AlternativelyActivated chr16short/BP_venous_blood_S00BS451_alternatively_activated_macrophage_DiseaseFree.png -label AlternativelyActivated chr16short/BP_venous_blood_S00FTN51_alternatively_activated_macrophage_DiseaseFree.png -label AlternativelyActivated chr16short/BP_cord_blood_S00E8W51_alternatively_activated_macrophage_DiseaseFree.png -label AlternativelyActivated chr16short/BP_venous_blood_S0062252_alternatively_activated_macrophage_DiseaseFree.png -label AlternativelyActivated chr16short/BP_venous_blood_S006VI53_alternatively_activated_macrophage_DiseaseFree.png -label macrophage chr16short/BP_venous_blood_S00V3BN1_macrophage_DiseaseFree.png -label macrophage chr16short/BP_venous_blood_S00V49N1_macrophage_DiseaseFree.png -label macrophage chr16short/BP_venous_blood_C005VG51_macrophage_DiseaseFree.png -label macrophage chr16short/BP_cord_blood_S00BHQ51_macrophage_DiseaseFree.png -label macrophage chr16short/BP_venous_blood_S0022I51_macrophage_DiseaseFree.png -label macrophage chr16short/BP_venous_blood_S001S751_macrophage_DiseaseFree.png -label macrophage chr16short/BP_cord_blood_S00DVR51_macrophage_DiseaseFree.png -label macrophage chr16short/BP_venous_blood_S0039051_macrophage_DiseaseFree.png -gravity South -geometry 1000x60+1 -tile 1x png:- | convert -trim - 1.png

##########################
## SD of 100kb mean
##########################
setwd('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS')
load('allprone/bin_mean_100kb.rda')
alltcga.samples <- colnames(betas)[grep('TCGA_', colnames(betas))]
tumorsamples <- alltcga.samples[!grepl('_N', alltcga.samples)]
normalsamples <- alltcga.samples[grepl('_N', alltcga.samples)]
sd.binmeanT <- na.omit(apply(betas[,tumorsamples], 1, sd, na.rm=T))
sd.binmeanN <- na.omit(apply(betas[,normalsamples], 1, sd, na.rm=T))
pdf('~/gallery/2017_01_01_tumor_100kb_mean_sd.pdf', width=6, height=3)
ggplot(wzbind(sd.binmeanN, sd.binmeanT, names=c('Normal','Tumor'))) + geom_histogram(aes(x, fill=cat), position='dodge', bins=50) + scale_fill_manual(values=c('#CC0000','#005AA0')) + xlim(0,0.25) + xlab('Standard Deviation') + ylab('# 100kb bins') + theme(text=element_text(size=18, face='bold'))
dev.off()

# percentage of PMDs
# plot(ecdf(sd.binmeanT))
# abline(v=c(0.1,0.12,0.15,0.17))
# > 1-ecdf(sd.binmeanT)(c(0.1,0.12,0.15,0.17))
# [1] 0.7374250 0.6333803 0.4862212 0.3719164

library(RColorBrewer)
library(MASS)

ovNames <- intersect(names(sd.binmeanN), names(sd.binmeanT))
df <- data.frame(x=sd.binmeanT[ovNames], y=sd.binmeanN[ovNames])
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

## Gaussian mixture
library("mixtools")
xT <- sd.binmeanT[ovNames]
xN <- sd.binmeanN[ovNames]
aa <- (!is.na(xT)) & (!is.na(xN))
xT <- xT[aa]
xN <- xN[aa]
mixT <- normalmixEM(xT, k=2)
mixN <- normalmixEM(xN, k=2)
Tiscomp1 <- apply(mixT$posterior, 1, function(x) x[1] > x[2])
## Niscomp1 <- apply(mixN$posterior, 1, function(x) x[1] > x[2])
## Tiscomp1 <- xT < 0.125
Niscomp1 <- xN < 0.056
a <- table(interaction(Tiscomp1, Niscomp1))
a
## FALSE.FALSE  TRUE.FALSE  FALSE.TRUE   TRUE.TRUE
##       15092        2596        1934        7361
## a / sum(a)
## FALSE.FALSE  TRUE.FALSE  FALSE.TRUE   TRUE.TRUE
##  0.55931512  0.09620872  0.07167476  0.27280139
## make sure comp1 is low in methylation
plot(density(xN))
lines(density(xN[Niscomp1]), col='blue')
plot(density(xT))
lines(density(xT[Tiscomp1]), col='blue')
(max(xN[Niscomp1]) + min(xN[!Niscomp1])) / 2
## [1] 0.04685495
(max(xT[Tiscomp1]) + min(xT[!Tiscomp1])) / 2
## [1] 0.1203376

options(stringsAsFactors=T)              # IMPORTANT FOR stat_density2d TO WORK!!!
pdf('~/gallery/2017_01_01_100kb_mean_sd_2d_1.pdf', width=5.5, height=4)
ggplot(df[df$x<0.25&df$y<0.25,], aes(x,y)) + stat_density2d(aes(fill=..level..), n=50, geom='polygon') + scale_fill_gradientn(colors=jet.colors(50)) + xlim(0,0.25) + ylim(0,0.25) + theme(panel.background=element_rect(fill=jet.colors(50)[1])) + xlab('Standard Deviation (Tumor)') + ylab('Standard Deviation (Normal)') + theme(text=element_text(size=18, face='bold'))
dev.off()

pdf('~/gallery/2017_01_01_100kb_mean_sd_2d_2.pdf', width=4, height=1.5)
ggplot(data.frame(x=df$x)) + geom_histogram(aes(x), bins=50, fill='#CC0000') + xlab('Standard Deviation') + xlim(0,0.25) + theme(text=element_text(size=18, face='bold'))
dev.off()

pdf('~/gallery/2017_01_01_100kb_mean_sd_2d_3.pdf', width=4, height=1.5)
ggplot(data.frame(x=df$y)) + geom_histogram(aes(x), bins=50, fill='#005AA0') + xlab('Standard Deviation') + xlim(0,0.25) + theme(text=element_text(size=18, face='bold'))
dev.off()

mean.binmeanT <- na.omit(apply(betas[,tumorsamples], 1, mean, na.rm=T))
pdf('~/gallery/2017_01_10_mean_vs_sd.pdf', width=3.5, height=3.5)
par(font.lab=2, font.axis=2)
par(mar=c(5,5,1,1))
smoothScatter(mean.binmeanT, sd.binmeanT, xlab='Sample Mean', ylab='Sample SD', main='', colramp=colorRampPalette(c('white','black')), nrpoints=0, cex.lab=1.3, cex.axis=1.3)
dev.off()

#################################
## sample meta data
#################################

setwd('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS')
samples <- read.csv2('samplesheet.txt', header=T, sep='\t', row.names='samplename', stringsAsFactors=F)
## save(samples, file='samples.rda')


## TCGA PMD depth ordering match
setwd('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS')
samples <- read.csv2('samplesheet.txt', header=T, sep='\t', row.names='samplename', stringsAsFactors=F)
a <- samples[grep('TCGA',rownames(samples)),]
aa <- data.frame(depth=as.numeric(a$PMDdepth), sample=rownames(a))
aa$cancertype <- NA
aa$cancertype[grep('_N',aa$sample)] <- paste0('normal_',sapply(strsplit(aa$sample[grep('_N',aa$sample)],'\\_'), function(x) x[2]))
aa$cancertype[!grepl('_N',aa$sample)] <- sapply(strsplit(aa$sample[!grepl('_N',aa$sample)],'\\_'), function(x) x[2])
ggplot(aa) + geom_jitter(aes(cancertype, depth)) + theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
## THE RESULT DOESN'T SUPPORT CONCLUSION

###############################
## rescaling
###############################

setwd('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS')
load('allprone/bin_mean_100kb.rda')

## scaling plot
chr16short <- rownames(betas)[grep('chr16',rownames(betas))]
chr16short <- chr16short[sapply(strsplit(chr16short,"[:-]"), function(x) as.integer(x[2])) < 35000000]
samples <- read.csv2('samplesheet.txt', header=T, sep='\t', row.names='samplename', stringsAsFactors=F)
rsamples <- samples[samples$tumor != 'CL' & samples$isFetal == 'Adult' & samples$source != 'BLUEPRINT' & samples$ok=='OK' & !(rownames(samples) %in% c('methbase_Heyn_2012_Human_BCell-ICF')),]

alltcga.samples <- colnames(betas)[grep('TCGA_', colnames(betas))]
tumorsamples <- alltcga.samples[!grepl('_N', alltcga.samples)]
normalsamples <- alltcga.samples[grepl('_N', alltcga.samples)]
sd.binmeanT <- na.omit(apply(betas[,tumorsamples], 1, sd, na.rm=T))
mean.binmeanT <- na.omit(apply(betas[,tumorsamples], 1, mean, na.rm=T))
sd.binmeanN <- na.omit(apply(betas[,normalsamples], 1, sd, na.rm=T))

pdf('~/gallery/2017_02_20_binmean_sample_mean_not_bimodal.pdf',width=6,height=6)
par(mfrow=c(2,2), mar=c(5,4,3,1), oma=c(2,2,1,1), lwd=1.5)
plot(density(mean.binmeanN, bw=0.0005), main='Cross-sample Mean (normal)')
plot(density(mean.binmeanT, bw=0.0005), main='Cross-sample Mean (tumor)')
plot(density(na.omit(sd.binmeanN)), main='Cross-sample SD (normal)')
plot(density(na.omit(sd.binmeanT)), main='Cross-sample SD (tumor)')
dev.off()

## redefine common PMD using TCGA tumor sd > 0.18
comPMD <- names(na.omit(sd.binmeanT[sd.binmeanT >= 0.18]))
quintiles <- apply(betas[comPMD,], 2, function(x) quantile(x, probs=c(0.2,0.8), na.rm=T))
## quintiles[1,] <- pmin(quintiles[2,]-0.12, quintiles[1,])

## output comPMD depth by bin mean
## depth <- apply(betas[comPMD,],2,mean,na.rm=T)
## write.table(comPMDdepth, file='comPMDdepth_bybinmean.tsv', sep='\t')

betas.scaled <- do.call(cbind, lapply(colnames(betas),function(sn) {
  b <- betas[,sn]; qt <- quintiles[,sn];
  b[!is.na(b) & b<qt[1]] <- qt[1];
  b[!is.na(b) & b>qt[2]] <- qt[2];
  bb <- !is.na(b) & b>=qt[1] & b<=qt[2]
  b[bb] <- (b[bb]-qt[1]) / (qt[2]-qt[1]);
  b
}))
colnames(betas.scaled) <- colnames(betas)

pdf('~/gallery/2017_01_04_rescaling.pdf', width=7, height=4)
WHeatmap(t(betas[chr16short, rownames(rsamples)]), cmp=CMPar(stop.points=c('blue','yellow'), dmin=0, dmax=1), name='h1') + WColorBarV(rsamples$source=='TCGA', RightOf()) + WColorBarV(rsamples$subtype2, cmp=CMPar(label2color=setNames(as.character(rsamples$color), rsamples$subtype2)), RightOf(), name='a') + WColorBarV(rsamples$tumor=='N', RightOf()) + WLegendV('a', BottomRightOf(h.pad=0.03), height=0.05) + WHeatmap(t(betas.scaled[chr16short, rownames(rsamples)]), Beneath('h1'), cmp=CMPar(stop.points=c('blue','yellow'), dmin=0, dmax=1))
dev.off()

pdf('~/gallery/2017_01_04_tumor_sd_binmean.pdf', width=7, height=1)
par(mar=c(1,1,1,1))
barplot(sd.binmeanT[chr16short], col='black')
dev.off()

pdf('~/gallery/2017_01_04_normal_sd_binmean.pdf', width=7, height=1)
par(mar=c(1,1,1,1))
barplot(sd.binmeanN[chr16short], col='black', ylim=c(0,0.1))
dev.off()

pdf('~/gallery/2017_01_04_segment.pdf', width=7,height=1)
par(mar=c(1,1,1,1))
WColorBarH(as.character(c(0, cumsum(seg0[-1] != seg0[-length(seg0)])) %% 2), cmp=CMPar(label2color=setNames(c('#CC0000','#005AA0'),c('0','1'))))
dev.off()

plot(density(na.omit(betas[,'TCGA_COAD_A00R'])))
plot(density(na.omit(betas[,'TCGA_COAD_N3158'])))
plot(density(na.omit(betas[,'methbase_Berman_2012_Human_ColonNormal'])))
plot(density(na.omit(betas[,'TCGA_COAD_3158'])))

## raw data
for (sname in c('TCGA_COAD_A00R', 'TCGA_COAD_N3158', 'methbase_Berman_2012_Human_ColonNormal', 'TCGA_COAD_3158')) {
  pdf(sprintf('~/gallery/2017_01_04_density_%s.pdf', sname), width=4, height=3)
  par(lwd=3)
  plot(density(na.omit(betas[comPMD,sname]), bw=0.03), col='#005AA0', main=sname, cex.axis=1.5, cex.lab=1.5)
  lines(density(na.omit(betas[,sname]),bw=0.03), col='#DAA753', cex.axis=1.5, cex.lab=1.5)
  dev.off()
}

## scaled version
for (sname in c('TCGA_COAD_A00R', 'TCGA_COAD_N3158', 'methbase_Berman_2012_Human_ColonNormal', 'TCGA_COAD_3158')) {
  pdf(sprintf('~/gallery/2017_01_04_scaled_density_%s.pdf', sname), width=4, height=3)
  par(lwd=3)
  plot(density(na.omit(betas.scaled[,sname]), bw=0.03), col='#DAA753', main=sname, cex.axis=1.5, cex.lab=1.5)
  lines(density(na.omit(betas.scaled[comPMD,sname]),bw=0.03), col='#005AA0', cex.axis=1.5, cex.lab=1.5)
  dev.off()
}

###############################
## rescaling V2, middle
###############################

setwd('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS')
load('allprone/bin_mean_100kb.rda')

## scaling plot
chr16short <- rownames(betas)[grep('chr16',rownames(betas))]
chr16short <- chr16short[sapply(strsplit(chr16short,"[:-]"), function(x) as.integer(x[2])) < 35000000]

betas.scaled <- apply(betas, 2, function(b1) { ifelse(b1 >= min(median(b1, na.rm=T)-0.03, quantile(b1, 0.5, na.rm=T)), 1, 0); })
samples <- read.csv2('samplesheet.txt', header=T, sep='\t', row.names='samplename', stringsAsFactors=F)
rsamples <- samples[samples$tumor != 'CL' & samples$isFetal == 'Adult' & samples$source != 'BLUEPRINT' & samples$ok=='OK' & !(rownames(samples) %in% c('methbase_Heyn_2012_Human_BCell-ICF')),]

pdf('~/gallery/2017_01_09_rescaling.pdf', width=7, height=4)
WHeatmap(t(betas[chr16short, rownames(rsamples)]), cmp=CMPar(stop.points=c('blue','yellow'), dmin=0, dmax=1), name='h1') + WColorBarV(rsamples$source=='TCGA', RightOf()) + WColorBarV(rsamples$subtype2, cmp=CMPar(label2color=setNames(as.character(rsamples$color), rsamples$subtype2)), RightOf(), name='a') + WColorBarV(rsamples$tumor=='N', RightOf()) + WLegendV('a', BottomRightOf(h.pad=0.03), height=0.05) + WHeatmap(t(betas.scaled[chr16short, rownames(rsamples)]), Beneath('h1'), cmp=CMPar(stop.points=c('blue','yellow'), dmin=0, dmax=1))
dev.off()

############################
## rescaling V3 just TCGA
############################

###############################
## rescaling
###############################

setwd('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS')
load('allprone/bin_mean_100kb.rda')

## scaling plot
chr16short <- rownames(betas)[grep('chr16',rownames(betas))]
chr16short <- chr16short[sapply(strsplit(chr16short,"[:-]"), function(x) as.integer(x[2])) < 35000000]
samples <- read.csv2('samplesheet.txt', header=T, sep='\t', row.names='samplename', stringsAsFactors=F)
rsamples <- samples[samples$tumor != 'CL' & samples$isFetal == 'Adult' & samples$source != 'BLUEPRINT' & samples$ok=='OK' & !(rownames(samples) %in% c('methbase_Heyn_2012_Human_BCell-ICF')),]

alltcga.samples <- colnames(betas)[grep('TCGA_', colnames(betas))]
tumorsamples <- alltcga.samples[!grepl('_N', alltcga.samples)]
normalsamples <- alltcga.samples[grepl('_N', alltcga.samples)]
sd.binmeanT <- na.omit(apply(betas[,tumorsamples], 1, sd, na.rm=T))
mean.binmeanT <- na.omit(apply(betas[,tumorsamples], 1, mean, na.rm=T))
sd.binmeanN <- na.omit(apply(betas[,normalsamples], 1, sd, na.rm=T))

## redefine common PMD using TCGA tumor sd > 0.18
comPMD <- names(na.omit(sd.binmeanT[sd.binmeanT >= 0.18]))
quintiles <- apply(betas[comPMD,], 2, function(x) quantile(x, probs=c(0.2,0.8), na.rm=T))
## quintiles[1,] <- pmin(quintiles[2,]-0.12, quintiles[1,])

## output comPMD depth by bin mean
## depth <- apply(betas[comPMD,],2,mean,na.rm=T)
## write.table(comPMDdepth, file='comPMDdepth_bybinmean.tsv', sep='\t')

betas.scaled <- do.call(cbind, lapply(colnames(betas),function(sn) {
  b <- betas[,sn]; qt <- quintiles[,sn];
  b[!is.na(b) & b<qt[1]] <- qt[1];
  b[!is.na(b) & b>qt[2]] <- qt[2];
  bb <- !is.na(b) & b>=qt[1] & b<=qt[2]
  b[bb] <- (b[bb]-qt[1]) / (qt[2]-qt[1]);
  b
}))
colnames(betas.scaled) <- colnames(betas)

normalsamples.manual <- c(
  "TCGA_BRCA_NA0CE", "TCGA_UCEC_NA1CI", "TCGA_LUSC_N2722", "TCGA_LUAD_N6148",
  "TCGA_STAD_N6452", "TCGA_COAD_N3158", 'methbase_Berman_2012_Human_ColonNormal',
  "TCGA_READ_N2689", "TCGA_BLCA_NA20V")
tumorsamples.manual <- c(
  "TCGA_GBM_1788",  "TCGA_GBM_1460",  "TCGA_GBM_0128",  "TCGA_GBM_3477", "TCGA_GBM_1401", "TCGA_GBM_1454",
  "TCGA_BRCA_A0YG", "TCGA_BRCA_A0CE", "TCGA_BRCA_A07I", "TCGA_BRCA_A15H", "TCGA_BRCA_A04X", "TCGA_UCEC_A0G2", "TCGA_UCEC_A1CI", "TCGA_UCEC_A0K6", "TCGA_UCEC_A05J", "TCGA_UCEC_A1CK", "TCGA_LUSC_2695", "TCGA_LUSC_1078", "TCGA_LUSC_2722", "TCGA_LUSC_2600", "TCGA_LUAD_6148", "TCGA_LUAD_6215", "TCGA_LUAD_7156", "TCGA_LUAD_4630", "TCGA_LUAD_6840", "TCGA_COAD_3158", "methbase_Berman_2012_Human_ColonCancer", "TCGA_COAD_A00R", "TCGA_READ_3593", "TCGA_READ_2689", "TCGA_STAD_6452", "TCGA_STAD_6519", "TCGA_STAD_6177", "TCGA_STAD_5730", "TCGA_BLCA_A2LA", "TCGA_BLCA_A13J", "TCGA_BLCA_A20V", "TCGA_BLCA_A1AG", "TCGA_BLCA_A2HQ", "TCGA_BLCA_A1AA")
  
rsamplestcga <- rsamples[c(normalsamples.manual, tumorsamples.manual),]
pdf('~/gallery/2017_02_17_rescaling_TCGA.pdf', width=7, height=4)
WHeatmap(t(betas[chr16short, rownames(rsamplestcga)]), cmp=CMPar(stop.points=c('blue','yellow'), dmin=0, dmax=1), name='h1') + WColorBarV(rsamplestcga$source=='TCGA', RightOf()) + WColorBarV(rsamplestcga$subtype2, cmp=CMPar(label2color=setNames(as.character(rsamplestcga$color), rsamplestcga$subtype2)), RightOf(), name='a') + WColorBarV(rsamplestcga$tumor=='N', RightOf()) + WLegendV('a', BottomRightOf(h.pad=0.03), height=0.05) + WHeatmap(t(betas.scaled[chr16short, rownames(rsamplestcga)]), Beneath('h1'), cmp=CMPar(stop.points=c('blue','yellow'), dmin=0, dmax=1))
dev.off()

pdf('~/gallery/2017_02_17_tumor_sd_binmean_justTCGA.pdf', width=7, height=1)
par(mar=c(1,1,1,1))
barplot(sd.binmeanT[chr16short], col='black')
dev.off()

pdf('~/gallery/2017_02_17_normal_sd_binmean_justTCGA.pdf', width=7, height=1)
par(mar=c(1,1,1,1))
barplot(sd.binmeanN[chr16short], col='black', ylim=c(0,0.1))
dev.off()

################################
## all sample heatmap chr16
################################
setwd('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS')
load('allprone/bin_mean_100kb.rda')
chr16short <- rownames(betas)[grep('chr16',rownames(betas))]
chr16short <- chr16short[sapply(strsplit(chr16short,"[:-]"), function(x) as.integer(x[2])) < 35000000]

## all samples
pdf('~/gallery/2017_01_03_allsample_100kbbin_prone.pdf', colormodel='cmyk', width=8, height=6)
WHeatmap(t(betas[chr16short, rownames(samples)]), cmp=CMPar(stop.points=c('blue','yellow'), dmin=0, dmax=1)) + WColorBarV(samples$source, RightOf()) + WColorBarV(rownames(samples), cmp=CMPar(label2color=setNames(as.character(samples$color), rownames(samples))), RightOf())
dev.off()

## cell line
samples <- read.csv2('samplesheet.txt', header=T, sep='\t', row.names='samplename', stringsAsFactors=F)
clsamples <- samples[samples$tumor=='CL',]
dim(clsamples)
pdf('~/gallery/2017_01_03_Cellline_100kbbin_prone.pdf', colormodel='cmyk', width=8, height=2)
WHeatmap(t(betas[chr16short, rownames(clsamples)]), cmp=CMPar(stop.points=c('blue','yellow'), dmin=0, dmax=1)) + WColorBarV(clsamples$subtype, cmp=CMPar(label2color=setNames(clsamples$color2, clsamples$subtype)), RightOf(), name='a') + WLegendV('a',BottomRightOf(h.pad=0.03), height=0.1)
dev.off()
pdf('~/gallery/2017_01_03_Cellline_100kbbin_prone_scaled.pdf', colormodel='cmyk', width=8, height=2)
WHeatmap(t(betas.scaled[chr16short, rownames(clsamples)]), cmp=CMPar(stop.points=c('blue','yellow'), dmin=0, dmax=1)) + WColorBarV(clsamples$subtype, cmp=CMPar(label2color=setNames(clsamples$color2, clsamples$subtype)), RightOf(), name='a') + WLegendV('a',BottomRightOf(h.pad=0.03), height=0.1)
dev.off()

## fetal
samples <- read.csv2('samplesheet.txt', header=T, sep='\t', row.names='samplename', stringsAsFactors=F)
fsamples <- samples[(samples$isFetal %in% c('Other','Fetal')) & samples$ok == 'OK',]
dim(fsamples)
pdf('~/gallery/2017_01_03_normal_fetal_100kbbin_prone.pdf', colormodel='cmyk', width=8, height=2)
WHeatmap(t(betas[chr16short, rownames(fsamples)]), cmp=CMPar(stop.points=c('blue','yellow'), dmin=0, dmax=1)) + WColorBarV(fsamples$subtype, cmp=CMPar(label2color=setNames(fsamples$color2, fsamples$subtype)), RightOf(), name='a') + WLegendV('a',BottomRightOf(h.pad=0.03), height=0.1)
dev.off()
pdf('~/gallery/2017_01_03_normal_fetal_100kbbin_prone_scaled.pdf', colormodel='cmyk', width=8, height=2)
WHeatmap(t(betas.scaled[chr16short, rownames(fsamples)]), cmp=CMPar(stop.points=c('blue','yellow'), dmin=0, dmax=1)) + WColorBarV(fsamples$tissue, cmp=CMPar(label2color=setNames(fsamples$color, fsamples$tissue)), RightOf(), name='a') + WLegendV('a',BottomRightOf(h.pad=0.03), height=0.1)
dev.off()

## postnatal normal non-blood
samples <- read.csv2('samplesheet.txt', header=T, sep='\t', row.names='samplename', stringsAsFactors=F)
nsamples <- samples[samples$plotorder >= 200 & samples$plotorder <= 300,]
dim(nsamples)
pdf('~/gallery/2017_01_03_normal_postnatal_100kbbin_prone.pdf', colormodel='cmyk', width=8, height=3)
WHeatmap(t(betas[chr16short, rownames(nsamples)]), cmp=CMPar(stop.points=c('blue','yellow'), dmin=0, dmax=1)) + WColorBarV(nsamples$tissue, cmp=CMPar(label2color=setNames(nsamples$color, nsamples$tissue)), RightOf(), name='a') + WLegendV('a',BottomRightOf(h.pad=0.03), height=0.1)
dev.off()
pdf('~/gallery/2017_01_03_normal_postnatal_100kbbin_prone_scaled.pdf', colormodel='cmyk', width=8, height=3)
WHeatmap(t(betas.scaled[chr16short, rownames(nsamples)]), cmp=CMPar(stop.points=c('blue','yellow'), dmin=0, dmax=1)) + WColorBarV(nsamples$tissue, cmp=CMPar(label2color=setNames(nsamples$color, nsamples$tissue)), RightOf(), name='a') + WLegendV('a',BottomRightOf(h.pad=0.03), height=0.1)
dev.off()

## postnatal normal blood
samples <- read.csv2('samplesheet.txt', header=T, sep='\t', row.names='samplename', stringsAsFactors=F)
bsamples <- samples[samples$plotorder >= 300 & samples$plotorder <= 400,]
dim(bsamples)
pdf('~/gallery/2017_01_03_blood_100kbbin_prone.pdf', colormodel='cmyk', width=8, height=3)
WHeatmap(t(betas[chr16short, rownames(bsamples)]), cmp=CMPar(stop.points=c('blue','yellow'), dmin=0, dmax=1)) + WColorBarV(bsamples$subtype, cmp=CMPar(label2color=setNames(bsamples$color, bsamples$subtype)), RightOf(), name='a') + WLegendV('a',BottomRightOf(h.pad=0.03), height=0.05)
dev.off()
pdf('~/gallery/2017_01_03_blood_100kbbin_prone_scaled.pdf', colormodel='cmyk', width=8, height=3)
WHeatmap(t(betas.scaled[chr16short, rownames(bsamples)]), cmp=CMPar(stop.points=c('blue','yellow'), dmin=0, dmax=1)) + WColorBarV(bsamples$subtype, cmp=CMPar(label2color=setNames(bsamples$color, bsamples$subtype)), RightOf(), name='a') + WLegendV('a',BottomRightOf(h.pad=0.03), height=0.05)
dev.off()

## tumor
samples <- read.csv2('samplesheet.txt', header=T, sep='\t', row.names='samplename', stringsAsFactors=F)
tsamples <- samples[samples$plotorder >= 500,]
dim(tsamples)
pdf('~/gallery/2017_01_03_tumor_100kbbin_prone.pdf', colormodel='cmyk', width=8, height=3)
WHeatmap(t(betas[chr16short, rownames(tsamples)]), cmp=CMPar(stop.points=c('blue','yellow'), dmin=0, dmax=1)) + WColorBarV(tsamples$tumor, cmp=CMPar(label2color=setNames(tsamples$color, tsamples$tumor)), RightOf(), name='a') + WLegendV('a',BottomRightOf(h.pad=0.03), height=0.05)
dev.off()
pdf('~/gallery/2017_01_03_tumor_100kbbin_prone_scaled.pdf', colormodel='cmyk', width=8, height=3)
WHeatmap(t(betas.scaled[chr16short, rownames(tsamples)]), cmp=CMPar(stop.points=c('blue','yellow'), dmin=0, dmax=1)) + WColorBarV(tsamples$tumor, cmp=CMPar(label2color=setNames(tsamples$color, tsamples$tumor)), RightOf(), name='a') + WLegendV('a',BottomRightOf(h.pad=0.03), height=0.05)
dev.off()

###############################################
## cluster within category (not very useful)
## I ended up just sorting by PMD depth
###############################################
betas16 <- t(betas[chr16short,])
samples$consecat <- c(0,cumsum(samples$color[-1] != samples$color[-length(samples$color)]))
a <- split(rownames(samples), samples$consecat)
aa <- lapply(a, function(sns) if (length(sns) > 1) rownames(row.cluster(betas16[sns,])$mat) else sns)
newsampleorder <- do.call(c, aa)
samples <- samples[newsampleorder,]

#############################################
## HM450
#############################################
library(GenomicRanges)
chrms <- c(paste0('chr',1:22),'chrX','chrY','chrM')
a <- read.table('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450probes.bed',sep='\t',header=F,stringsAsFactors=F)
a <- a[a$V1 %in% chrms,]
gr <- GRanges(seqnames=a$V1, ranges=IRanges(a$V2+1,a$V3), seqinfo=Seqinfo(chrms))
names(gr) <- a$V4
mcols(gr) <- data.frame(gene=as.character(a$V5), orphan35=a$V6, gc35=a$V7, ctxt=a$V8, ctxtclass=a$V10, pmd=a$V11, cgi=a$V12)
probes <- gr
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/cgcontext.hg19/ctxt_w_sdbinmeanT.rda')
probes <- mergeOv(probes, ctxt)

## Huy's common PMD set
com.pmd.probes <- names(probes[probes$pmd=='comPMD' & probes$orphan35==0 & probes$cgi=='nonCGI' & probes$ctxtclass=='WCGW'])
com.hmd.probes <- names(probes[probes$pmd=='comHMD' & probes$orphan35==0 & probes$cgi=='nonCGI' & probes$ctxtclass=='WCGW'])
length(com.pmd.probes)
## [1] 1026
length(com.hmd.probes)
## [1] 9768

## SD-based common PMD set
com.pmd.probes0 <- unique(names(probes[probes$sd.binmeanT >=0.15 & probes$orphan35==0 & probes$cgi=='nonCGI' & probes$ctxtclass=='WCGW']))
com.hmd.probes0 <- unique(names(probes[probes$sd.binmeanT <=0.10 & probes$orphan35==0 & probes$cgi=='nonCGI' & probes$ctxtclass=='WCGW']))
## remove TCGA normal unmethylated probes
load('/primary/projects/laird/projects/2016_04_05_TCGA_pancan_renormalization/450k/normal.betas.rda')
pmd.unmeth.cnt <- apply(normal.betas[com.pmd.probes0, ], 1, function(x) sum(x <= 0.2, na.rm=T))
hmd.unmeth.cnt <- apply(normal.betas[com.hmd.probes0, ], 1, function(x) sum(x <= 0.2, na.rm=T))
com.pmd.probes <- com.pmd.probes0[pmd.unmeth.cnt < 20]
com.hmd.probes <- com.hmd.probes0[hmd.unmeth.cnt < 20]
length(com.pmd.probes)
## [1] 6214
length(com.hmd.probes)
## [1] 9040
save(com.pmd.probes, com.hmd.probes, probes, file='/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450probes.rda')
plot(density(na.omit(normal.betas[com.hmd.probes,])))
lines(density(na.omit(normal.betas[com.pmd.probes,])), col='blue')

## overlapping Ben's list gives too few probes
## Ben's stringent list
a <- read.table('hm450/groupC-Mincvg3-covminprctilelow5-covmaxprctilehigh90-minCov0.02-highcov-hg19.bed', sep='\t', header=F, stringsAsFactors=F)
a <- a[a$V1 %in% chrms,]
gr <- GRanges(seqnames=a$V1, ranges=IRanges(a$V2+1,a$V3), seqinfo=Seqinfo(chrms))
BenListStringent <- gr

## Ben's relaxed list
a <- read.table('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/Mincvg3-covminprctilelow10-covmaxprctilehigh75-minCov0.02-highcov-hg19.bed', sep='\t', header=F, stringsAsFactors=F)
a <- a[a$V1 %in% chrms,]
gr <- GRanges(seqnames=a$V1, ranges=IRanges(a$V2+1,a$V3+1), seqinfo=Seqinfo(chrms))
BenListRelaxed <- gr
probesBLrelax <- mergeOv(probes, BenListRelaxed)
com.pmd.probesB <- names(probesBLrelax[probesBLrelax$orphan35==0 & probesBLrelax$cgi=='nonCGI' & probesBLrelax$ctxtclass=='WCGW'])
com.hmd.probesB <- names(probes[probes$pmd=='comHMD' & probes$orphan35==0 & probes$cgi=='nonCGI' & probes$ctxtclass=='WCGW'])
length(com.pmd.probesB)
## [1] 418
length(com.hmd.probesB)
## [1] 9768
save(com.pmd.probesB, com.hmd.probesB, probes, probesBLrelax, file='/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450probes_Ben.rda')

###########################
## single cpg K36 analysis
###########################

readsignal <- function(fn) {
  a <- read.table(fn)
  setNames(a$V4,paste0(a$V1,":",a$V2,"-",a$V3))
}

## H1
df <- data.frame(
  repliseq = readsignal('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/REMC_H3K36me3/H1_repliseq_singlecpg.bed'),
  dname = readsignal('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/REMC_H3K36me3/dname_singlecpg/E003_H1_Cell_Line.bed'),
  k36 = readsignal('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/REMC_H3K36me3/K36_singlecpg/E003.H3K36me3.bed'),
  genebody = readsignal('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/REMC_H3K36me3/genebody_singlecpg.bed'))
df$seg <- paste0(df$k36, '.', df$genebody)
dfc <- data.frame(
  dname = tapply(df$dname, df$seg, mean),
  repliseq = tapply(df$repliseq, df$seg, mean),
  k36 = tapply(df$k36, df$seg, function(x) x[1]),
  genebody = tapply(df$genebody, df$seg, function(x) x[1]),
  cnt=table(df$seg))
save(df, dfc, file='/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/REMC_H3K36me3/H1.rda')

df <- df[complete.cases(df),]
pdf('~/gallery/2017_04_05_smoothScatterK36_H1.pdf', width=8, height=8)
par(mfrow=c(2,2))
a <- df[df$k36>0 & df$genebody>0,]
smoothScatter(-a$repliseq, a$dname, nrpoints=0, colramp=colorRampPalette(c('white','black')), xlab='RepliChip', ylab='Solo-WCGW CpG\nMethylation', main=sprintf('K36+/Gene+ (%d cpgs)', dim(a)[1]), xlim=c(-2,1), ylim=c(0,1))

a <- df[df$k36>0 & df$genebody<0,]
smoothScatter(-a$repliseq, a$dname, nrpoints=0, colramp=colorRampPalette(c('white','black')), xlab='RepliChip', ylab='Solo-WCGW CpG\nMethylation', main=sprintf('K36+/Gene- (%d cpgs)', dim(a)[1]), xlim=c(-2,1), ylim=c(0,1))

a <- df[df$k36<0 & df$genebody>1,]
smoothScatter(-a$repliseq, a$dname, nrpoints=0, colramp=colorRampPalette(c('white','black')), xlab='RepliChip', ylab='Solo-WCGW CpG\nMethylation', main=sprintf('K36-/Gene+ (%d cpgs)', dim(a)[1]), xlim=c(-2,1), ylim=c(0,1))

a <- df[df$k36<0 & df$genebody<0,]
smoothScatter(-a$repliseq, a$dname, nrpoints=0, colramp=colorRampPalette(c('white','black')), xlab='RepliChip', ylab='Solo-WCGW CpG\nMethylation', main=sprintf('K36-/Gene- (%d cpgs)', dim(a)[1]), xlim=c(-2,1), ylim=c(0,1))
dev.off()

load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/REMC_H3K36me3/H1.rda')
pdf('~/gallery/2017_04_05_smoothScatterK36_seg_H1_repliseq.pdf', width=8, height=8)
par(mfrow=c(2,2))
a <- dfc[dfc$k36>0 & dfc$genebody>0,]
plot(-a$repliseq, a$dname, cex=0.1, pch=16, xlab='RepliChip', ylab='Solo-WCGW CpG\nMethylation', main=sprintf('K36+/Gene+ (%d segments)', dim(a)[1]), xlim=c(-2,1), ylim=c(0,1))

a <- dfc[dfc$k36>0 & dfc$genebody<0,]
plot(-a$repliseq, a$dname, cex=0.1, pch=16, xlab='RepliChip', ylab='Solo-WCGW CpG\nMethylation', main=sprintf('K36+/Gene- (%d segments)', dim(a)[1]), xlim=c(-2,1), ylim=c(0,1))

a <- dfc[dfc$k36<0 & dfc$genebody>0,]
plot(-a$repliseq, a$dname, cex=0.1, pch=16, xlab='RepliChip', ylab='Solo-WCGW CpG\nMethylation', main=sprintf('K36-/Gene+ (%d segments)', dim(a)[1]), xlim=c(-2,1), ylim=c(0,1))

a <- dfc[dfc$k36<0 & dfc$genebody<0,]
plot(-a$repliseq, a$dname, cex=0.1, pch=16, xlab='RepliChip', ylab='Solo-WCGW CpG\nMethylation', main=sprintf('K36-/Gene- (%d segments)', dim(a)[1]), xlim=c(-2,1), ylim=c(0,1))
dev.off()

## IMR90
bb <- read.table('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/REMC_H3K36me3/IMR90_K36_genebody.bed')
df <- data.frame(
  repliseq = as.numeric(readsignal('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/REMC_H3K36me3/IMR90_repliseq_singlecpg.bed')),
  dname = as.numeric(readsignal('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/REMC_H3K36me3/dname_singlecpg/E017_IMR90_Cell_Line.bed')),
  genebody = bb$V4,
  k36 = bb$V5,
  seglength = bb$V6)
## k36 = readsignal('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/REMC_H3K36me3/K36_singlecpg/E017.H3K36me3.bed'),
## genebody = readsignal('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/REMC_H3K36me3/genebody_singlecpg.bed'))
df <- df[complete.cases(df),]
df$seg <- paste0(df$k36, '.', df$genebody)
dfc <- data.frame(
  dname = tapply(df$dname, df$seg, mean),
  repliseq = tapply(df$repliseq, df$seg, mean),
  k36 = tapply(df$k36, df$seg, function(x) x[1]),
  genebody = tapply(df$genebody, df$seg, function(x) x[1]),
  len = tapply(df$seglength, df$seg, function(x) x[1]),
  cnt=table(df$seg))
save(df, dfc, file='/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/REMC_H3K36me3/IMR90.rda')

pdf('~/gallery/2017_04_05_smoothScatterK36_IMR90.pdf', width=8, height=8)
par(mfrow=c(2,2))
a <- df[df$k36>0 & df$genebody>0,]
smoothScatter(90-a$repliseq, a$dname, nrpoints=0, colramp=colorRampPalette(c('white','black')), xlab='RepliSeq', ylab='Solo-WCGW CpG\nMethylation', main=sprintf('K36+/Gene+ (%d cpgs)', dim(a)[1]), ylim=c(0,1))

a <- df[df$k36>0 & df$genebody<0,]
smoothScatter(90-a$repliseq, a$dname, nrpoints=0, colramp=colorRampPalette(c('white','black')), xlab='RepliSeq', ylab='Solo-WCGW CpG\nMethylation', main=sprintf('K36+/Gene- (%d cpgs)', dim(a)[1]), ylim=c(0,1))

a <- df[df$k36<0 & df$genebody>0,]
smoothScatter(90-a$repliseq, a$dname, nrpoints=0, colramp=colorRampPalette(c('white','black')), xlab='RepliSeq', ylab='Solo-WCGW CpG\nMethylation', main=sprintf('K36-/Gene+ (%d cpgs)', dim(a)[1]), ylim=c(0,1))

a <- df[df$k36<0 & df$genebody<0,]
smoothScatter(90-a$repliseq, a$dname, nrpoints=0, colramp=colorRampPalette(c('white','black')), xlab='RepliSeq', ylab='Solo-WCGW CpG\nMethylation', main=sprintf('K36-/Gene- (%d cpgs)', dim(a)[1]), ylim=c(0,1))
dev.off()

pdf('~/gallery/2017_04_05_smoothScatterK36_seg_IMR90_repliseq.pdf', width=8, height=8)
par(mfrow=c(2,2))
a <- dfc[dfc$k36>0 & dfc$genebody>0,]
plot(90-a$repliseq, a$dname, cex=0.1, pch=16, xlab='RepliSeq', ylab='Solo-WCGW CpG\nMethylation', main=sprintf('K36+/Gene+ (%d segments)', dim(a)[1]), ylim=c(0,1))

a <- dfc[dfc$k36>0 & dfc$genebody<0,]
plot(90-a$repliseq, a$dname, cex=0.1, pch=16, xlab='RepliSeq', ylab='Solo-WCGW CpG\nMethylation', main=sprintf('K36+/Gene- (%d segments)', dim(a)[1]), ylim=c(0,1))

a <- dfc[dfc$k36<0 & dfc$genebody>0,]
plot(90-a$repliseq, a$dname, cex=0.1, pch=16, xlab='RepliSeq', ylab='Solo-WCGW CpG\nMethylation', main=sprintf('K36-/Gene+ (%d segments)', dim(a)[1]), ylim=c(0,1))

a <- dfc[dfc$k36<0 & dfc$genebody<0,]
plot(90-a$repliseq, a$dname, cex=0.1, pch=16, xlab='RepliSeq', ylab='Solo-WCGW CpG\nMethylation', main=sprintf('K36-/Gene- (%d segments)', dim(a)[1]), ylim=c(0,1))
dev.off()

pdf('~/gallery/2017_04_05_smoothScatterK36_seg_IMR90_repliseq_lengtheffect_early.pdf', width=8, height=8)
par(mfrow=c(2,2))
a <- dfc[dfc$k36>0 & dfc$genebody>0 & dfc$repliseq > 60,]
plot(log2(1+a$len), a$dname, cex=0.3, pch=16, xlab='log2 segment length', ylab='Solo-WCGW CpG\nMethylation', main=sprintf('K36+/Gene+ (%d segments)', dim(a)[1]), xlim=c(5,20), ylim=c(0,1))

a <- dfc[dfc$k36>0 & dfc$genebody<0 & dfc$repliseq > 60,]
plot(log2(1+a$len), a$dname, cex=0.3, pch=16, xlab='log2 segment length', ylab='Solo-WCGW CpG\nMethylation', main=sprintf('K36+/Gene- (%d segments)', dim(a)[1]), xlim=c(5,20), ylim=c(0,1))

a <- dfc[dfc$k36<0 & dfc$genebody>0 & dfc$repliseq > 60,]
plot(log2(1+a$len), a$dname, cex=0.3, pch=16, xlab='log2 segment length', ylab='Solo-WCGW CpG\nMethylation', main=sprintf('K36-/Gene+ (%d segments)', dim(a)[1]), xlim=c(5,20), ylim=c(0,1))

a <- dfc[dfc$k36<0 & dfc$genebody<0 & dfc$repliseq > 60,]
plot(log2(1+a$len), a$dname, cex=0.3, pch=16, xlab='log2 segment length', ylab='Solo-WCGW CpG\nMethylation', main=sprintf('K36-/Gene- (%d segments)', dim(a)[1]), xlim=c(5,20), ylim=c(0,1))
dev.off()


pdf('~/gallery/2017_04_05_smoothScatterK36_seg_IMR90_repliseq_lengtheffect_late.pdf', width=8, height=8)
par(mfrow=c(2,2))
a <- dfc[dfc$k36>0 & dfc$genebody>0 & dfc$repliseq < 30,]
plot(log2(1+a$len), a$dname, cex=0.3, pch=16, xlab='log2 segment length', ylab='Solo-WCGW CpG\nMethylation', main=sprintf('K36+/Gene+ (%d segments)', dim(a)[1]), xlim=c(5,20), ylim=c(0,1))

a <- dfc[dfc$k36>0 & dfc$genebody<0 & dfc$repliseq < 30,]
plot(log2(1+a$len), a$dname, cex=0.3, pch=16, xlab='log2 segment length', ylab='Solo-WCGW CpG\nMethylation', main=sprintf('K36+/Gene- (%d segments)', dim(a)[1]), xlim=c(5,20), ylim=c(0,1))

a <- dfc[dfc$k36<0 & dfc$genebody>0 & dfc$repliseq < 30,]
plot(log2(1+a$len), a$dname, cex=0.3, pch=16, xlab='log2 segment length', ylab='Solo-WCGW CpG\nMethylation', main=sprintf('K36-/Gene+ (%d segments)', dim(a)[1]), xlim=c(5,20), ylim=c(0,1))

a <- dfc[dfc$k36<0 & dfc$genebody<0 & dfc$repliseq < 30,]
plot(log2(1+a$len), a$dname, cex=0.3, pch=16, xlab='log2 segment length', ylab='Solo-WCGW CpG\nMethylation', main=sprintf('K36-/Gene- (%d segments)', dim(a)[1]), xlim=c(5,20), ylim=c(0,1))
dev.off()








dfc$cnt.bin = cut(dfc$cnt.Freq, breaks=c(quantile(dfc$cnt.Freq, probs = seq(0, 1, by = 0.20)))[-1], include.lowest=T)
for (lv in levels(dfc$cnt.bin)) {
  dfc1 <- dfc[dfc$cnt.bin == lv,]
  pdf(sprintf('~/gallery/2017_04_05_smoothScatterK36_seg_IMR90_repliseq_%s.pdf', lv), width=8, height=8)
  par(mfrow=c(2,2))
  a <- dfc1[dfc1$k36>0 & dfc1$genebody>0,]
  plot(90-a$repliseq, a$dname, cex=0.1, pch=16, xlab='RepliSeq', ylab='Solo-WCGW CpG\nMethylation', main=sprintf('K36+/Gene+ (%d segments)', dim(a)[1]), ylim=c(0,1))

  a <- dfc1[dfc1$k36>0 & dfc1$genebody<0,]
  plot(90-a$repliseq, a$dname, cex=0.1, pch=16, xlab='RepliSeq', ylab='Solo-WCGW CpG\nMethylation', main=sprintf('K36+/Gene- (%d segments)', dim(a)[1]), ylim=c(0,1))

  a <- dfc1[dfc1$k36<0 & dfc1$genebody>0,]
  plot(90-a$repliseq, a$dname, cex=0.1, pch=16, xlab='RepliSeq', ylab='Solo-WCGW CpG\nMethylation', main=sprintf('K36-/Gene+ (%d segments)', dim(a)[1]), ylim=c(0,1))

  a <- dfc1[dfc1$k36<0 & dfc1$genebody<0,]
  plot(90-a$repliseq, a$dname, cex=0.1, pch=16, xlab='RepliSeq', ylab='Solo-WCGW CpG\nMethylation', main=sprintf('K36-/Gene- (%d segments)', dim(a)[1]), ylim=c(0,1))
  dev.off()
}

## IMR90 lamin
df <- data.frame(
  lamin = as.numeric(readsignal('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/REMC_H3K36me3/IMR90_lamin_singlecpg.bed')),
  dname = as.numeric(readsignal('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/REMC_H3K36me3/dname_singlecpg/E017_IMR90_Cell_Line.bed')),
  k36 = readsignal('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/REMC_H3K36me3/K36_singlecpg/E017.H3K36me3.bed'),
  genebody = readsignal('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/REMC_H3K36me3/genebody_singlecpg.bed'))
df <- df[complete.cases(df),]
df$seg <- paste0(df$k36, '.', df$genebody)
dfc <- data.frame(
  dname = tapply(df$dname, df$seg, mean),
  lamin = tapply(df$lamin, df$seg, mean),
  k36 = tapply(df$k36, df$seg, function(x) x[1]),
  genebody = tapply(df$genebody, df$seg, function(x) x[1]))
save(df, dfc, file='/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/REMC_H3K36me3/IMR90_lamin.rda')

pdf('~/gallery/2017_04_05_smoothScatterK36_seg_IMR90_lamin.pdf', width=8, height=8)
par(mfrow=c(2,2))
a <- dfc[dfc$k36>0 & dfc$genebody>0,]
plot(a$lamin, a$dname, cex=0.1, pch=16, xlab='LaminB1', ylab='Solo-WCGW CpG\nMethylation', main=sprintf('K36+/Gene+ (%d segments)', dim(a)[1]), ylim=c(0,1), xlim=c(-3,3))

a <- dfc[dfc$k36>0 & dfc$genebody<0,]
plot(a$lamin, a$dname, cex=0.1, pch=16, xlab='LaminB1', ylab='Solo-WCGW CpG\nMethylation', main=sprintf('K36+/Gene- (%d segments)', dim(a)[1]), ylim=c(0,1), xlim=c(-3,3))

a <- dfc[dfc$k36<0 & dfc$genebody>0,]
plot(a$lamin, a$dname, cex=0.1, pch=16, xlab='LaminB1', ylab='Solo-WCGW CpG\nMethylation', main=sprintf('K36-/Gene+ (%d segments)', dim(a)[1]), ylim=c(0,1), xlim=c(-3,3))

a <- dfc[dfc$k36<0 & dfc$genebody<0,]
plot(a$lamin, a$dname, cex=0.1, pch=16, xlab='LaminB1', ylab='Solo-WCGW CpG\nMethylation', main=sprintf('K36-/Gene- (%d segments)', dim(a)[1]), ylim=c(0,1), xlim=c(-3,3))
dev.off()

###########################
## IMR90 (old obsolete)
###########################

#### 30kb-bin ####

setwd('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/final_nocentromere')
mnames <- c('DNase','H3K27ac','H3K27me3','H3K36me3','H3K4me1','H3K4me3','H3K9me3','laminB1','pronecpg','repliseq','expression','genedensity')
chrms <- c(paste0('chr',1:22),'chrX','chrY','chrM')
allgrs <- lapply(mnames, function(mname) {
  cat(mname, '\n')
  a <- read.table(sprintf('%s_30kbbin.bed', mname), header=F, sep='\t', stringsAsFactors=F)
  a <- a[a$V1 %in% chrms,]
  gr <- GRanges(seqnames=a$V1, ranges=IRanges(a$V2+1,a$V3), seqinfo=Seqinfo(chrms))
  gr$val <- a$V4
  colnames(mcols(gr)) <- c(mname)
  gr
})
names(allgrs) <- mnames
gr <- Reduce(mergeOv, allgrs)
gr0 <- gr

## clean data a bit, rescale outliers
gr$laminB1[gr$laminB1 < -3] <- -3
gr$laminB1[gr$laminB1 > 3] <- 3

gr$repliseq[gr$repliseq < 0] <- 0

gr$DNase[gr$DNase < 20000] <- 20000
gr$DNase[gr$DNase > 400000] <- 400000
gr$DNase <- log2(1+gr$DNase)

gr$H3K9me3 <- log2(1+gr$H3K9me3)
gr$H3K9me3[gr$H3K9me3 < 6.5] <- 6.5
gr$H3K9me3[gr$H3K9me3 > 11.5] <- 11.5

gr$H3K27me3 <- log2(1+gr$H3K27me3)
gr$H3K27me3[gr$H3K27me3 < 5] <- 5
gr$H3K27me3[gr$H3K27me3 > 12] <- 12

gr$H3K36me3 <- log2(1+gr$H3K36me3)
gr$H3K36me3[gr$H3K36me3 < 4.3] <- 4.3
gr$H3K36me3[gr$H3K36me3 > 12] <- 12

gr$H3K4me3 <- log2(1+gr$H3K4me3)
gr$H3K4me3[gr$H3K4me3 < 2] <- 2
gr$H3K4me3[gr$H3K4me3 > 15] <- 15

gr$H3K27ac <- log2(1+gr$H3K27ac)
gr$H3K27ac[gr$H3K27ac < 2] <- 2
gr$H3K27ac[gr$H3K27ac > 13] <- 13

gr$H3K4me1 <- log2(1+gr$H3K4me1)
gr$H3K4me1[gr$H3K4me1 < 3] <- 3
gr$H3K4me1[gr$H3K4me1 > 12] <- 12

gr$genedensity <- log2(1+gr$genedensity)
gr$genedensity[gr$genedensity < 5] <- 5

vv <- unique(gr$expression)
v2rank <- setNames(rank(vv), vv)
gr$expression <- v2rank[as.character(gr$expression)]

save(gr0, gr, file='/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/gr.rda')

panel.cor <- function(x, y, digits=2, prefix="", cex.cor) {
  usr <- par("usr"); on.exit(par(usr)) 
  par(usr = c(0, 1, 0, 1)) 
  rraw <- cor(x, y)
  r <- abs(rraw)
  txt <- format(c(rraw, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="") 
  if(missing(cex.cor)) cex <- 0.8/strwidth(txt) 
  
  test <- cor.test(x,y) 
  ## borrowed from printCoefmat
  Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                   symbols = c("***", "**", "*", ".", " ")) 
  
  text(0.5, 0.5, txt, cex = max(0.8, cex * r))
  ## text(.8, .8, Signif, cex=cex, col=2) 
}

panel.dens = function (x, ...) {
  par(new = TRUE)
  plot(density(x), col = "red", lwd = 2, main="", axes=FALSE, font.main=2)
}

df <- as.data.frame(mcols(gr))[,c('genedensity',"expression","DNase","H3K4me1","H3K27ac","H3K4me3","H3K27me3","H3K36me3","H3K9me3","repliseq","laminB1","pronecpg")]
pdf('~/gallery/2017_01_07_hypomethylation_unsupervised_clustering.pdf', width=10, height=10)
pairs(df, panel = function(...) smoothScatter(..., colramp=colorRampPalette(c('white','black')), nrpoints=0, add = TRUE), upper.panel=panel.cor, diag.panel=panel.dens)
dev.off()

png('~/gallery/2017_01_07_hypomethylation_unsupervised_clustering.png', width=2000, height=2000)
pairs(df, panel = function(...) smoothScatter(..., colramp=colorRampPalette(c('white','black')), nrpoints=0, add = TRUE), upper.panel=panel.cor, diag.panel=panel.dens)
dev.off()

## no log-transform
gr2 <- gr0
gr2$H3K9me3[gr2$H3K9me3 > 2000] <- 2000
gr2$H3K27me3[gr$H3K27me3 > 2000] <- 2000
gr2$laminB1[gr2$laminB1 > 3] <- 3
gr2$laminB1[gr2$laminB1 < -3] <- -3

df <- as.data.frame(mcols(gr2))[,c('genedensity',"expression","DNase","H3K4me1","H3K27ac","H3K4me3","H3K27me3","H3K36me3","H3K9me3","repliseq","laminB1","pronecpg")]
png('~/gallery/2017_01_30_hypomethylation_unsupervised_clustering_nolog.png', width=2000, height=2000)
pairs(df, panel = function(...) smoothScatter(..., colramp=colorRampPalette(c('white','black')), nrpoints=0, add = TRUE), upper.panel=panel.cor, diag.panel=panel.dens)
dev.off()

df <- as.data.frame(mcols(gr0))[,c('genedensity',"expression","DNase","H3K4me1","H3K27ac","H3K4me3","H3K27me3","H3K36me3","H3K9me3","repliseq","laminB1","pronecpg")]
png('~/gallery/2017_01_30_hypomethylation_unsupervised_clustering_raw.png', width=2000, height=2000)
pairs(df, panel = function(...) smoothScatter(..., colramp=colorRampPalette(c('white','black')), nrpoints=0, add = TRUE), upper.panel=panel.cor, diag.panel=panel.dens)
dev.off()

###############
## IMR90 HiC
###############
setwd('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/HiC/hIMR90/nij')
bins <- read.table('/primary/vari/genomicdata/genomes/hg18/annotation/hg18.windows40kb.binname.bed')$V1
chrms <- c(paste0('chr',1:22))
windowmap <- read.table('/primary/vari/genomicdata/genomes/hg18/annotation/hg18.windows40kb.to.hg19.windows30kb.bed')
map19to18 <- setNames(windowmap$V2, windowmap$V1)
map18to19 <- setNames(windowmap$V1, windowmap$V2)

for (chrm in chrms) {
  cat(chrm,'\n')
  ## chrm <- 'chr16'
  chrmbins <- bins[grep(paste0(chrm,':'), bins)]
  nij <- as.matrix(read.table(paste0('nij.',chrm)))
  colnames(nij) <- chrmbins
  rownames(nij) <- chrmbins
  
  chrmlevels <- map18to19[grep(paste0(chrm,':'),map18to19)]
  chrmlevels <- unname(chrmlevels[unname(!c('a',chrmlevels[-length(chrmlevels)]) != chrmlevels)])
  chrmbins19 <- factor(map18to19[chrmbins], levels=chrmlevels)
  nij1 <- apply(nij, 1, function(x) tapply(x, chrmbins19, mean, na.rm=T))
  nij2 <- apply(nij1, 1, function(x) tapply(x, chrmbins19, mean, na.rm=T))
  nij <- nij2
  save(nij, file=paste0('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/HiC/hIMR90/nij/',chrm,'.rda'))
}

load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/HiC/hIMR90/nij/chr16.rda')
chr16short <- colnames(nij)
chr16short <- chr16short[sapply(strsplit(chr16short,"[:-]"), function(x) as.integer(x[2])) < 35000000]

png('~/gallery/2017_01_29_IMR90_HiC_chr16.png', width=600, height=600)
WHeatmap(log2(1+nij[chr16short, chr16short]), cmp=CMPar(dmin=0, dmax=4)) + WCustomize(mar.left=0, mar.top=0, mar.right=0, mar.bottom=0)
dev.off()

pdf('~/gallery/2017_01_29_IMR90_HiC_chr16_legend.pdf', width=2, height=2)
WColorBarH(seq(0,4,0.1)) + WCustomize(mar.left=0, mar.top=0, mar.right=0, mar.bottom=0)
dev.off()

###########################################
## IMR90 multiscale, nonoverlapping bins
###########################################

library(plyr)
convertrds <- function(dname, oname) {
  setwd(dname)
  pws <- sprintf('%1.1f', seq(2,7,0.1))
  a <- lapply(pws, function(pw) {
    a <- read.table(paste0('scale_',pw,'.bed'))
    score <- setNames(as.numeric(a$V6), paste0(a$V1,':',a$V4,'-',a$V5))
    score
  })
  measure <- Reduce(function(x,y) t(rbind.fill.matrix(t(x), t(y))), a)
  colnames(measure) <- pws
  saveRDS(measure, file=paste0('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/multiscale/',oname,'.rds'))
}

## pronecpg
convertrds('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/pronecpg/final/', 'pronecpg')
convertrds('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/pronecpg/cnt/', 'pronecpgcnt')

## genedensity
convertrds('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/genebody/final/', 'genebodyoverlap')
convertrds('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/exons/final/', 'exonoverlap')

## rnaseq
convertrds('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/RNAseq/final/', 'rnaseq')

## CTCF
convertrds('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/CTCF/final/', 'ctcf')

## Lamin B1
convertrds('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/LaminB1/final/', 'laminb1')
ob <- readRDS('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/multiscale/laminb1.rds')
ob[ob < -3] <- NA
ob[ob > 3] <- NA
saveRDS(ob, file='/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/multiscale/laminb1.rds')

## DNase-seq
convertrds('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/DNase/final/', 'dnase')

## RepliSeq
convertrds('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/repliseq/final/wavesmooth', 'repliseq')
for (phase in c('G1b','S1','S2','S3','S4','G2')) {
  cat(phase,'\n')
  convertrds(paste0('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/repliseq/final/',phase), paste0('repliseq',phase))
}

## histone marks
for (mark in c('H3K9me3','H3K27ac','H3K27me3','H3K4me3','H3K4me1','H3K36me3')) {
  cat(mark,'\n')
  convertrds(paste0('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/', mark, '/final/'), mark)
}

## H3K36me3ov
convertrds(paste0('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/H3K36me3/finalcnt/'), 'H3K36me3ov')

#######################################
## plot IMR90 correlation heatmap
#######################################
features <- c('ctcf', 'H3K27ac', 'H3K36me3', 'H3K4me3', 'laminb1', 'rnaseq', 'dnase', 'H3K27me3', 'H3K4me1', 'H3K9me3', 'pronecpg', 'pronecpgcnt', 'repliseq', 'repliseqG1b', 'repliseqG2', 'repliseqS1', 'repliseqS2', 'repliseqS3', 'repliseqS4', 'genebodyoverlap', 'exonoverlap', 'H3K36me3ov')
toofficial <- setNames(features, features)
toofficial['ctcf'] <- 'CTCF'
toofficial['dnase'] <- 'DNase-Seq'
toofficial['rnaseq'] <- 'RNA-Seq'
toofficial['laminb1'] <- 'Lamin B1'
toofficial['pronecpg'] <- 'Solo-WCGW CpG\nMethylation'
toofficial <- sub('repliseq','RepliSeq ',toofficial)

signals <- mclapply(features, function(oname) {
  readRDS(paste0('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/multiscale/',oname,'.rds'))
}, mc.cores=20)
names(signals) <- features
save(signals, file='/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/multiscale/signals.rda')

####### Spearman's correlation coefficients #######
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/multiscale/signals.rda')
getspearmans <- function(o1, o2) {
  pws <- sprintf('%1.1f', seq(2,7,0.1))
  multicors <- sapply(1:51, function(i) {
    sapply(1:51, function(j) {
      a <- o1[,i]
      b <- o2[,j]
      cor(a, b, use='na.or.complete', method='spearman')
    })
  })
  colnames(multicors) <- pws
  rownames(multicors) <- pws
  multicors
}

corfeatures <- c('H3K27ac', 'H3K36me3', 'H3K4me3', 'laminb1', 'rnaseq', 'dnase', 'H3K27me3', 'H3K4me1', 'H3K9me3', 'pronecpg', 'repliseq', 'genebodyoverlap')
allcors <- mclapply(signals[corfeatures], function(sig1) {
  a <- mclapply(signals[corfeatures], function(sig2) {
    getspearmans(sig1, sig2)
  }, mc.cores=6)
  names(a) <- corfeatures
  a
}, mc.cores=6)
names(allcors) <- corfeatures

## repliseq data can have 6 missing correlation coefficients
save(allcors, file='/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/allspearman.rda')

##################################
### correlation heatmap by marks
##################################

load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/allspearman.rda')
pws <- seq(2,7,0.1)

maxrowind <- sapply(allcors, function(allcors1) {
  sapply(allcors1, function(cors) {
    ind <- which(abs(cors) == max(abs(cors)), arr.ind=T)
    if (length(ind) < 10)
      ind[1]
    else
      1
  })
})
maxcolind <- sapply(allcors, function(allcors1) {
  sapply(allcors1, function(cors) {
    ind <- which(abs(cors) == max(abs(cors)), arr.ind=T)
    if (length(ind) < 10)
      ind[2]
    else
      1
  })
})
maxcor <- sapply(allcors, function(allcors1) {
  sapply(allcors1, function(cors) {
    max(abs(cors))
  })
})

nsmp <- dim(maxcolind)[1]
maxcor <- both.cluster(maxcor)$mat
maxcolind <- maxcolind[rownames(maxcor), colnames(maxcor)]
maxrowind <- maxrowind[rownames(maxcor), colnames(maxcor)]

a <- WHeatmap(maxcor, name='heat', xticklabels=TRUE, yticklabels=TRUE, cmp=CMPar(brewer.name='Reds', dmin=0, dmax=1))
for (i in 1:nsmp) {
  for (j in 1:nsmp) {
    if (i < j) {
      a <- a + WLabel(bquote(10^.(pws[maxrowind[i,j]])*", "*10^.(pws[maxcolind[i,j]])),
                      WPosition((i-0.5)/nsmp,1-(j-0.5)/nsmp, 'heat', just=c('center','center')))
    } else if (i > j) {
      a <- a + WLabel(sprintf('%1.2f', maxcor[i,j]),
                      WPosition((i-0.5)/nsmp,1-(j-0.5)/nsmp, 'heat', just=c('center','center')))
    }
  }
}
pdf('~/gallery/2017_02_02_allscale_best.pdf', width=15, height=10)
print(a)
dev.off()

##################################
### correlation heatmap by scale
##################################

## subfeatures <- c('pronecpg', 'laminb1', 'repliseqG2',
##                  'H3K36me3', 'H3K4me3', 'repliseqG1b', 
##                  'repliseqS1', 'H3K4me1', 'H3K27ac', 'dnase',
##                   'H3K27me3', 'repliseqS2', 'rnaseq',
##                  'repliseqS3', 'repliseqS4', 'H3K9me3')

## for (pw in c('2.0','3.0','4.0','5.0','6.0','7.0')) {
##   cors1 <- sapply(allcors, function(allcors1) {
##     sapply(allcors1, function(cors) {
##       abs(cors[pw,pw])
##     })
##   })
##   pdf(sprintf('~/gallery/2017_02_04_allcors_scale_%s.pdf', pw), width=6, height=6)
##   print(WHeatmap(cors1[subfeatures, subfeatures], xticklabels=T, yticklabels=T, cmp=CMPar(dmin=0,dmax=1)))
##   dev.off()
## }

## with clustering
## for (pw in c('2.0','3.0','4.0','5.0','6.0','7.0')) {
##   cors1 <- sapply(allcors, function(allcors1) {
##     sapply(allcors1, function(cors) {
##       abs(cors[pw,pw])
##     })
##   })
##   clus <- column.cluster(cors1)
##   colnames(clus$mat)
##   feats <- colnames(clus$mat)
##   pdf(sprintf('~/gallery/2017_04_03_allcors_scale_100kb_dendro_%s.pdf', pw), width=6, height=6)
##   print(WHeatmap(cors1[feats,feats], xticklabels=T, yticklabels=T, cmp=CMPar(dmin=0,dmax=1)) + WDendrogram(clus$column.clust, TopOf(), facing='bottom'))
##   dev.off()
## }

## single cpg
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/singlecpg/allcors.rda')
cors1 <- abs(allcors)
clus <- column.cluster(cors1)
feats <- colnames(clus$mat)
pdf('~/gallery/2017_04_09_allcors_scale_dendro_singlecpg.pdf', width=6, height=6)
WHeatmap(cors1[feats, feats], xticklabels=T, yticklabels=T, cmp=CMPar(dmin=0,dmax=1)) + WDendrogram(clus$column.clust, TopOf(), facing='bottom')
dev.off()
subfeatures <- c('solowcgw','repliseq','H3K9me3','laminb1','H3K36me3','genebody', 'H3K27me3','rnaseq','H3K4me1', 'H3K27ac', 'dnaseq', 'H3K4me3')
pdf('~/gallery/2017_04_09_allcors_scale_all_singlecpg.pdf', width=6, height=6)
print(WHeatmap(cors1[subfeatures, subfeatures], xticklabels=T, yticklabels=T, cmp=CMPar(stop.points=c('white','red'), dmin=0,dmax=1)))
dev.off()

load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/allspearman.rda')

## 2.0
pw <- '2.0'
cors1 <- sapply(allcors, function(allcors1) {
  sapply(allcors1, function(cors) {
    abs(cors[pw,pw])
  })
})
clus <- column.cluster(cors1)
colnames(clus$mat)
feats <- colnames(clus$mat)
pdf(sprintf('~/gallery/2017_04_09_allcors_scale_dendro_%s.pdf', pw), width=6, height=6)
print(WHeatmap(cors1[feats,feats], xticklabels=T, yticklabels=T, cmp=CMPar(dmin=0,dmax=1)) + WDendrogram(clus$column.clust, TopOf(), facing='bottom'))
dev.off()
subfeatures <- c('pronecpg','repliseq','H3K36me3','genebodyoverlap','rnaseq','laminb1','H3K9me3','H3K27me3','H3K4me1', 'H3K27ac', 'H3K4me3', 'dnase')
pdf(sprintf('~/gallery/2017_04_09_allcors_scale_all_%s.pdf', pw), width=6, height=6)
print(WHeatmap(cors1[subfeatures, subfeatures], xticklabels=T, yticklabels=T, cmp=CMPar(stop.points=c('white','red'), dmin=0,dmax=1)))
dev.off()

## 3.0
pw <- '3.0'
cors1 <- sapply(allcors, function(allcors1) {
  sapply(allcors1, function(cors) {
    abs(cors[pw,pw])
  })
})
clus <- column.cluster(cors1)
colnames(clus$mat)
feats <- colnames(clus$mat)
pdf(sprintf('~/gallery/2017_04_09_allcors_scale_dendro_%s.pdf', pw), width=6, height=6)
print(WHeatmap(cors1[feats,feats], xticklabels=T, yticklabels=T, cmp=CMPar(dmin=0,dmax=1)) + WDendrogram(clus$column.clust, TopOf(), facing='bottom'))
dev.off()
subfeatures <- c('pronecpg','repliseq','laminb1','H3K9me3','H3K36me3','rnaseq','genebodyoverlap','H3K27me3','H3K4me1', 'H3K27ac', 'dnase', 'H3K4me3')
pdf(sprintf('~/gallery/2017_04_09_allcors_scale_all_%s.pdf', pw), width=6, height=6)
print(WHeatmap(cors1[subfeatures, subfeatures], xticklabels=T, yticklabels=T, cmp=CMPar(stop.points=c('white','red'), dmin=0,dmax=1)))
dev.off()

## 5.0
pw <- '5.0'
cors1 <- sapply(allcors, function(allcors1) {
  sapply(allcors1, function(cors) {
    abs(cors[pw,pw])
  })
})
clus <- column.cluster(cors1)
colnames(clus$mat)
feats <- colnames(clus$mat)
pdf(sprintf('~/gallery/2017_04_09_allcors_scale_dendro_%s.pdf', pw), width=6, height=6)
print(WHeatmap(cors1[feats,feats], xticklabels=T, yticklabels=T, cmp=CMPar(dmin=0,dmax=1)) + WDendrogram(clus$column.clust, TopOf(), facing='bottom'))
dev.off()
subfeatures <- c('pronecpg','laminb1','repliseq','dnase','H3K27ac','H3K4me1','H3K36me3','rnaseq','H3K4me3','H3K9me3','genebodyoverlap','H3K27me3')
pdf(sprintf('~/gallery/2017_04_09_allcors_scale_all_%s.pdf', pw), width=6, height=6)
print(WHeatmap(cors1[subfeatures, subfeatures], xticklabels=T, yticklabels=T, cmp=CMPar(stop.points=c('white','red'), dmin=0,dmax=1)))
dev.off()

############# singlecpg plus multiscale ###########

load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/singlecpg/allcors.rda')
features1 <- c("repliseq", "H3K36me3", "H3K9me3", "genebody", "laminb1", "H3K4me1", "dnaseq", "H3K27me3", "H3K4me3", "H3K27ac", "rnaseq")
singlecpg <- allcors['solowcgw', features1]
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/allspearman.rda')
features2 <- c('repliseq', "H3K36me3", "H3K9me3", "genebodyoverlap", "laminb1", "H3K4me1", "dnase", "H3K27me3", "H3K4me3", "H3K27ac", "rnaseq")
pws <- sprintf('%1.1f', seq(2,7,0.1))
solowcgwcors <- do.call(cbind, lapply(pws, function(pw) sapply(features2, function(x) allcors$pronecpg[[x]][pw,pw])))
colnames(solowcgwcors) <- pws
solowcgwcors <- cbind(single=singlecpg, solowcgwcors)
pdf('~/gallery/2017_04_11_solowcgwcors_scale_all.pdf', width=6, height=3)
print(WHeatmap(abs(solowcgwcors), xticklabels=T, yticklabels=T, cmp=CMPar(stop.points=c('white','red'), dmin=0,dmax=1)))
dev.off()

load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/multiscale/signals.rda')
chr16short <- rownames(signals$laminb1)[grep('chr16',rownames(signals$laminb1))]
chr16short <- chr16short[sapply(strsplit(chr16short,"[:-]"), function(x) as.integer(x[2])) < 35000000]
nml <- function(x) {
  (na.omit(x) - min(x, na.rm=T)) / (max(x, na.rm=T) - min(x, na.rm=T))
}

pdf('~/gallery/2017_04_11_explore_correlation.pdf', width=6, height=8)
par(mfrow=c(5,1), mar=c(1,1,3,1))
pw <- '2.0'
chr16short1 <- chr16short[!is.na(signals$pronecpg[chr16short, pw]) & !is.na(signals$repliseq[chr16short, pw]) & !is.na(signals$laminb1[chr16short, pw])]
plot(nml(signals$pronecpg[chr16short1, pw]), type='l', ylab='signal', main='100bp')
lines(nml(signals$repliseq[chr16short1, pw]), type='l', col='blue')
lines(1-nml(signals$laminb1[chr16short1, pw]), type='l', col='red')
## lines(nml(signals$repliseqG2[chr16short, pw]), type='l', col='green')
pw <- '3.0'
chr16short1 <- chr16short[!is.na(signals$pronecpg[chr16short, pw]) & !is.na(signals$repliseq[chr16short, pw]) & !is.na(signals$laminb1[chr16short, pw])]
plot(nml(signals$pronecpg[chr16short1, pw]), type='l', ylab='signal', main='1kbp')
lines(nml(signals$repliseq[chr16short1, pw]), type='l', col='blue')
lines(1-nml(signals$laminb1[chr16short1, pw]), type='l', col='red')
## lines(nml(signals$repliseqG2[chr16short, pw]), type='l', col='green')
pw <- '4.0'
chr16short1 <- chr16short[!is.na(signals$pronecpg[chr16short, pw]) & !is.na(signals$repliseq[chr16short, pw]) & !is.na(signals$laminb1[chr16short, pw])]
plot(nml(signals$pronecpg[chr16short1, pw]), type='l', ylab='signal', main='10kbp')
lines(nml(signals$repliseq[chr16short1, pw]), type='l', col='blue')
lines(1-nml(signals$laminb1[chr16short1, pw]), type='l', col='red')
## lines(nml(signals$repliseqG2[chr16short, pw]), type='l', col='green')
pw <- '5.0'
chr16short1 <- chr16short[!is.na(signals$pronecpg[chr16short, pw]) & !is.na(signals$repliseq[chr16short, pw]) & !is.na(signals$laminb1[chr16short, pw])]
plot(nml(signals$pronecpg[chr16short1, pw]), type='l', ylab='signal', main='100kbp')
lines(nml(signals$repliseq[chr16short1, pw]), type='l', col='blue')
lines(1-nml(signals$laminb1[chr16short1, pw]), type='l', col='red')
## lines(nml(signals$repliseqG2[chr16short, pw]), type='l', col='green')
pw <- '6.0'
chr16short1 <- chr16short[!is.na(signals$pronecpg[chr16short, pw]) & !is.na(signals$repliseq[chr16short, pw]) & !is.na(signals$laminb1[chr16short, pw])]
plot(nml(signals$pronecpg[chr16short1, pw]), type='l', ylab='signal', main='1mbp')
lines(nml(signals$repliseq[chr16short1, pw]), type='l', col='blue')
lines(1-nml(signals$laminb1[chr16short1, pw]), type='l', col='red')
## lines(nml(signals$repliseqG2[chr16short, pw]), type='l', col='green')
dev.off()

## ## without clustering
## ## subfeatures <- c('pronecpg','laminb1','repliseq','dnase','H3K27ac','H3K4me1','H3K27me3','rnaseq','H3K36me3','H3K4me3','H3K9me3')
## subfeatures <- c('pronecpg','repliseq','H3K36me3','laminb1','dnase','H3K27ac','H3K4me1','H3K27me3','rnaseq','H3K4me3','H3K9me3')
## for (pw in c('2.0','3.0','4.0','5.0','6.0','7.0')) {
##   cors1 <- sapply(allcors, function(allcors1) {
##     sapply(allcors1, function(cors) {
##       cors[pw,pw]
##     })
##   })
##   pdf(sprintf('~/gallery/2017_04_03_allcors_scale_%s.pdf', pw), width=6, height=6)
##   print(WHeatmap(cors1[subfeatures, subfeatures], xticklabels=T, yticklabels=T, cmp=CMPar(stop.points=c('blue','white','red'), dmin=-1,dmax=1)))
##   dev.off()
## }

## load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/allspearman.rda')
## pw <- '7.0'
## cors1 <- sapply(allcors, function(allcors1) {
##   sapply(allcors1, function(cors) {
##     abs(cors[pw,pw])
##   })
## })
## clus1 <- both.cluster(cors1[subfeatures, subfeatures])
## pdf('~/gallery/2017_02_02_1MB_dendrogram.pdf')
## colnames(clus1$mat) <- toofficial[colnames(clus1$mat)]
## WHeatmap(clus1$mat, xticklabels=T) + WDendrogram(clus1$column.clust, TopOf(height=0.5))
## dev.off()

#########################################
######## 1Mb scatter plots of signals
#########################################

## signals1[signals1[,'repliseqG2']>1000000,'repliseqG2'] <- NA
## signals1[signals1[,'dnase']<5000000,'dnase'] <- NA
## signals1[signals1[,'rnaseq']>18000000,'rnaseq'] <- NA
## signals1[signals1[,'repliseqS4']>5e5,'repliseqS4'] <- NA

## ## png
## for (ft1 in 'pronecpg') {
##   for (ft2 in subfeatures) {
##     png(sprintf('~/gallery/2017_02_04_1mb_allcors_png/%s_vs_%s.png', ft1, ft2), width=400, height=600)
##     par(mar=c(0,0,0,0), lwd=2, cex.lab=1.3, cex.axis=1.3, font.lab=2, font.axis=2)
##     cc <- cor(signals1[,ft2], signals1[,ft1], use='complete', method='spearman')
##     smoothScatter(signals1[,ft2], signals1[,ft1], colramp=colorRampPalette(c('white','black')), nrpoints=0, xlab=ft2, ylab=ft1, main='', xaxt='n', yaxt='n')
##     dev.off()
##   }
## }

## pdf
## pw <- '7.0'
## 
## signals1 <- do.call(cbind, lapply(signals, function(s1) s1[,pw]))
## for (ft1 in 'pronecpg') {
##   for (ft2 in subfeatures) {
##     pdf(sprintf('~/gallery/2017_02_04_10mb_allcors/%s_vs_%s.pdf', ft1, ft2), width=3.7, height=5)
##     par(mar=c(5,5,4,2), lwd=3, cex.lab=2, cex.axis=1.7, cex.main=2, font.main=2, font.lab=2, font.axis=2)
##     cc <- cor(signals1[,ft2], signals1[,ft1], use='complete', method='spearman')
##     ## smoothScatter(signals1[,ft2], signals1[,ft1], colramp=colorRampPalette(c('white','black')), nrpoints=0, xlab=ft2, ylab=ft1, main=sprintf('Spearman: %1.3f', cc))
##     plot(signals1[,ft2], signals1[,ft1], pch=20, xlab=toofficial[ft2], ylab=toofficial[ft1], main=sprintf("Spearman's\nr = %1.3f", cc))
##     dev.off()
##   }
## }

load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/multiscale/signals.rda')
subfeatures <- c('H3K27ac', 'H3K36me3', 'H3K4me3', 'laminb1', 'rnaseq', 'dnase', 'H3K27me3', 'H3K4me1', 'H3K9me3', 'pronecpg', 'repliseq', 'genebodyoverlap')
toofficial <- setNames(subfeatures, subfeatures)
toofficial['ctcf'] <- 'CTCF'
toofficial['dnase'] <- 'DNase-Seq'
toofficial['rnaseq'] <- 'RNA-Seq'
toofficial['laminb1'] <- 'Lamin B1'
toofficial['pronecpg'] <- 'Solo-WCGW CpG\nMethylation'
toofficial <- sub('repliseq','RepliSeq ',toofficial)

pw <- '5.0'
signals1 <- do.call(cbind, lapply(signals, function(s1) s1[,pw]))
signals1[,'laminb1'][signals1[,'laminb1'] > 2] <- NA
signals1[,'H3K4me3'][signals1[,'H3K4me3'] > 15000] <- NA
signals1[,'dnase'][signals1[,'dnase'] > 3e6] <- NA
signals1[,'dnase'][signals1[,'dnase'] < 0.5e6] <- NA
signals1[,'H3K4me1'][signals1[,'H3K4me1'] > 6e4] <- NA
signals1[,'H3K27ac'][signals1[,'H3K27ac'] > 6e4] <- NA
signals1[,'rnaseq'] <- log2(signals1[,'rnaseq'])
signals1[,'rnaseq'][signals1[,'rnaseq'] > 21] <- NA
for (ft1 in 'pronecpg') {
  for (ft2 in subfeatures) {
    pdf(sprintf('~/gallery/2017_04_06_100kb_allcors/%s_vs_%s.pdf', ft1, ft2), width=3.9, height=4)
    par(mar=c(5,7,4,2), lwd=3, cex.lab=2, cex.axis=1.7, cex.main=2, font.main=2, font.lab=2, font.axis=2)
    cc <- cor(signals1[,ft2], signals1[,ft1], use='complete', method='spearman')
    ## smoothScatter(signals1[,ft2], signals1[,ft1], colramp=colorRampPalette(c('white','black')), nrpoints=0, xlab=ft2, ylab=ft1, main=sprintf('Spearman: %1.3f', cc))
    coef <- sprintf("= %1.3f", cc)
    plot(signals1[,ft2], signals1[,ft1], pch=16, xlab=toofficial[ft2], ylab=toofficial[ft1], main=bquote(rho~.(coef)), cex=0.1)
    dev.off()
  }
}

pw <- '4.5'
signals1 <- do.call(cbind, lapply(signals, function(s1) s1[,pw]))
signals1[,'laminb1'][signals1[,'laminb1'] > 2.5] <- NA
signals1[,'rnaseq'][signals1[,'rnaseq'] > 2e5] <- NA
signals1[,'dnase'][signals1[,'dnase'] < 1e5] <- NA
signals1[,'dnase'][signals1[,'dnase'] >1.2e6] <- NA
signals1[,'H3K4me3'][signals1[,'H3K4me3'] >2e4] <- NA
signals1[,'H3K27ac'][signals1[,'H3K27ac'] >2e4] <- NA
for (ft1 in 'pronecpg') {
  for (ft2 in subfeatures) {
    pdf(sprintf('~/gallery/2017_04_06_30kbb_allcors/%s_vs_%s.pdf', ft1, ft2), , width=3.9, height=4)
    par(mar=c(5,7,4,2), lwd=3, cex.lab=2, cex.axis=1.7, cex.main=2, font.main=2, font.lab=2, font.axis=2)
    cc <- cor(signals1[,ft2], signals1[,ft1], use='complete', method='spearman')
    coef <- sprintf("= %1.3f", cc)
    plot(signals1[,ft2], signals1[,ft1], pch=16, xlab=toofficial[ft2], ylab=toofficial[ft1], main=bquote(rho~.(coef)), cex=0.1)
    dev.off()
  }
}

#################################
## single CpG IMR90 all features
#################################

beds <- list.files('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/singlecpg','.bed')
features <- sub('.bed','',beds)
signals <- lapply(beds, function(bed) {
  a <- read.table(paste0('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/singlecpg/', bed))
  as.numeric(a$V4)
})
names(signals) <- features
signals$laminb1[!is.na(signals$laminb1) & signals$laminb1 < -30] <- -30
signals$laminb1[!is.na(signals$laminb1) & signals$laminb1 > 30] <- 30
save(signals, file='/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/singlecpg/signals.rda')

load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/singlecpg/signals.rda')
allcors <- simplify2array(mclapply(1:length(signals), function(i) {
  cat(i,'\n');
  simplify2array(mclapply(1:length(signals), function(j) {
    cor(signals[[i]], signals[[j]], use='complete')
  }, mc.cores=6))
}, mc.cores=6))
colnames(allcors) <- names(signals)
rownames(allcors) <- names(signals)
save(allcors, file='/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/singlecpg/allcors.rda')

###############################
## Meta-gene analysis (IMR90)
###############################

df <- read.table('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/gene_structure/IMR90.genic.solowcgw_withK36repliseq.bed', col.names=c('chrm','beg','end','genename','val','reg','beta','k36','repliseq'))
df$binrepliseq <- cut(df$repliseq, breaks=c(-Inf,40,60,75,Inf), include.lowest=T)
## quantile(df$repliseq, 0:4/4, na.rm=T)
df$x <- NA
df$x[df$reg == -1] <- as.integer(cut(df$val[df$reg == -1], seq(-10000, 0, length=30), include.lowest=T))
df$x[df$reg == 0] <- 30 + as.integer(cut(df$val[df$reg == 0], 0:30/30, include.lowest=T))
df$x[df$reg == 1] <- 60 + as.integer(cut(df$val[df$reg == 1], seq(0, 10000, length=30), include.lowest=T))

## > quantile(df$k36, probs=0:5/5, na.rm=T)
##          0%         20%         40%         60%         80%        100%
## 0.000102191 0.049545500 0.227815000 0.461296000 0.730113000 1.000000000
## > a <- cut(df$k36, breaks=quantile(df$k36, probs=0:5/5, na.rm=T), include.lowest=T)
## > levels(a)
## [1] "[0.000102,0.0495]" "(0.0495,0.228]"    "(0.228,0.461]"
## [4] "(0.461,0.73]"      "(0.73,1]"

pdf('~/gallery/2017_04_07_metageneplot_k36neg.pdf', width=5, height=3)
df0 <- df[df$k36==0,]
df0b <- melt(tapply(df0$beta, list(df0$x, df0$binrepliseq), mean))
print(ggplot(df0b) + geom_line(aes(as.integer(Var1), value, color=Var2)) + ylim(0,1) + ggtitle(sprintf('%d cpgs, %d genes', dim(df0)[1], length(unique(df0$genename)))) + xlab('Position') + ylab('Solo-WCGW CpG\nMethylation'))
dev.off()

## > table(df0$binrepliseq)
## [-Inf,40]   (40,60]   (60,75] (75, Inf]
##    458070    381251    341502    236631
## > sum(table(df0$binrepliseq))
## [1] 1417454

pdf('~/gallery/2017_04_07_metageneplot_k36pos.pdf', width=5, height=3)
df0 <- df[df$k36==1,]
df0b <- melt(tapply(df0$beta, list(df0$x, df0$binrepliseq), mean))
print(ggplot(df0b) + geom_line(aes(as.integer(Var1), value, color=Var2)) + ylim(0,1) + ggtitle(sprintf('%d cpgs, %d genes', dim(df0)[1], length(unique(df0$genename)))) + xlab('Position') + ylab('Solo-WCGW CpG\nMethylation'))
dev.off()

## > table(df0$binrepliseq)
## [-Inf,40]   (40,60]   (60,75] (75, Inf]
##     33428     84243    143393    177168
## > sum(table(df0$binrepliseq))
## [1] 438232

## intergenic part
###################
df <- read.table('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/gene_structure/IMR90.intergenic.solowcgw_withK36repliseq.bed', col.names=c('chrm','beg','end','solowcgw','k36','repliseq'))
df$binrepliseq <- cut(df$repliseq, breaks=c(-Inf,40,60,75,Inf), include.lowest=T)

## quantile(df$repliseq, 0:4/4, na.rm=T)
pdf('~/gallery/2017_04_07_intergenic_k36neg.pdf', width=3.5, height=3)
df0 <- df[df$k36==0,]
ggplot(df0) + geom_violin(aes(binrepliseq, solowcgw), bw=0.03, draw_quantiles=c(0.5)) + xlab('RepliSeq') + ggtitle(sprintf('%d cpgs', dim(df0)[1])) + stat_summary(aes(binrepliseq, solowcgw), fun.y=mean, geom="point", shape=23, size=2)
dev.off()

pdf('~/gallery/2017_04_07_intergenic_k36pos.pdf', width=3.5, height=3)
df0 <- df[df$k36==1,]
ggplot(df0) + geom_violin(aes(binrepliseq, solowcgw), bw=0.07, draw_quantiles=c(0.5)) + xlab('RepliSeq')+ ggtitle(sprintf('%d cpgs', dim(df0)[1])) + stat_summary(aes(binrepliseq, solowcgw), fun.y=mean, geom="point", shape=23, size=2)
dev.off()

###############################
## Meta-gene analysis (H1)
###############################

df <- read.table('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/REMC_H3K36me3/H1.genic.solowcgw_withK36repliseq.bed', col.names=c('chrm','beg','end','genename','val','reg','beta','k36','repliseq'))

quantile(na.omit(df$repliseq), 0:4/4)
##         0%        25%        50%        75%       100%
## -1.3777900 -0.4983145  0.3943540  1.1546600  1.9422000

df$binrepliseq <- cut(df$repliseq, breaks=c(-Inf, -0.4983145,  0.3943540,  1.1546600, Inf), include.lowest=T)
df$x <- NA
df$x[df$reg == -1] <- as.integer(cut(df$val[df$reg == -1], seq(-10000, 0, length=30), include.lowest=T))
df$x[df$reg == 0] <- 30 + as.integer(cut(df$val[df$reg == 0], 0:30/30, include.lowest=T))
df$x[df$reg == 1] <- 60 + as.integer(cut(df$val[df$reg == 1], seq(0, 10000, length=30), include.lowest=T))

pdf('~/gallery/2017_04_10_H1_metageneplot_k36neg.pdf', width=5, height=3)
df0 <- df[df$k36==0,]
df0b <- melt(tapply(df0$beta, list(df0$x, df0$binrepliseq), mean))
print(ggplot(df0b) + geom_line(aes(as.integer(Var1), value, color=Var2)) + ylim(0.5,1) + ggtitle(sprintf('%d cpgs, %d genes', dim(df0)[1], length(unique(df0$genename)))) + xlab('Position') + ylab('Solo-WCGW CpG\nMethylation'))
dev.off()

table(df0$binrepliseq)
sum(table(df0$binrepliseq))
## > table(df0$binrepliseq)
##  [-Inf,-0.498] (-0.498,0.394]   (0.394,1.15]    (1.15, Inf]
##         434928         394777         348725         307385
## > sum(table(df0$binrepliseq))
## [1] 1485815

pdf('~/gallery/2017_04_10_H1_metageneplot_k36pos.pdf', width=5, height=3)
df0 <- df[df$k36==1,]
df0b <- melt(tapply(df0$beta, list(df0$x, df0$binrepliseq), mean))
print(ggplot(df0b) + geom_line(aes(as.integer(Var1), value, color=Var2)) + ylim(0.5,1) + ggtitle(sprintf('%d cpgs, %d genes', dim(df0)[1], length(unique(df0$genename)))) + xlab('Position') + ylab('Solo-WCGW CpG\nMethylation'))
dev.off()

table(df0$binrepliseq)
sum(table(df0$binrepliseq))
## > table(df0$binrepliseq)
##  [-Inf,-0.498] (-0.498,0.394]   (0.394,1.15]    (1.15, Inf]
##          23561          63713         109771         151096
## > sum(table(df0$binrepliseq))
## [1] 348141

## intergenic part
###################
df <- read.table('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/REMC_H3K36me3/H1.intergenic.solowcgw_withK36repliseq.bed', col.names=c('chrm','beg','end','solowcgw','k36','repliseq'))
df$binrepliseq <- cut(df$repliseq, breaks=c(-Inf, -0.4983145,  0.3943540,  1.1546600, Inf), include.lowest=T)

pdf('~/gallery/2017_04_10_H1_intergenic_k36neg.pdf', width=3.5, height=3)
df0 <- df[df$k36==0,]
ggplot(df0) + geom_violin(aes(binrepliseq, solowcgw), bw=0.03, draw_quantiles=c(0.5)) + xlab('RepliSeq') + ggtitle(sprintf('%d cpgs', dim(df0)[1])) + stat_summary(aes(binrepliseq, solowcgw), fun.y=mean, geom="point", shape=23, size=2)
dev.off()

table(df0$binrepliseq)
sum(table(df0$binrepliseq))
## > table(df0$binrepliseq)
##  [-Inf,-0.498] (-0.498,0.394]   (0.394,1.15]    (1.15, Inf]
##         722527         280391         178180         129762
## > sum(table(df0$binrepliseq))
## [1] 1310860

pdf('~/gallery/2017_04_10_H1_intergenic_k36pos.pdf', width=3.5, height=3)
df0 <- df[df$k36==1,]
ggplot(df0) + geom_violin(aes(binrepliseq, solowcgw), bw=0.07, draw_quantiles=c(0.5)) + xlab('RepliSeq') + ggtitle(sprintf('%d cpgs', dim(df0)[1])) + stat_summary(aes(binrepliseq, solowcgw), fun.y=mean, geom="point", shape=23, size=2)
dev.off()

table(df0$binrepliseq)
sum(table(df0$binrepliseq))
## > table(df0$binrepliseq)
##  [-Inf,-0.498] (-0.498,0.394]   (0.394,1.15]    (1.15, Inf]
##           8474           4827           4962           5260
## > sum(table(df0$binrepliseq))
## [1] 23523

#############################################
## K36me3 positive/negative
#############################################
## load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/gr.rda')
## pdf('~/gallery/2017_01_30_H3K36me3_dist.pdf', width=3.5, height=3)
## par(mar=c(5,5,1,1))
## plot(density(gr0$H3K36me3), xlab='#H3K36me3 tags / 30kb', ylab='Density', main='')
## abline(v=300, lty='dashed')
## dev.off()

## ## 30kb, stringent
## set.seed(1)
## k36neg <- sample(gr0[gr0$H3K36me3<100], 1000)
## k36pos <- sample(gr0[gr0$H3K36me3>=500], 1000)
## pdf('~/gallery/2017_02_05_lamin_vs_meth_H3K36me3_stringent.pdf', width=6, height=6.5)
## par(mfrow=c(2,2), lwd=3, cex.lab=2, cex.axis=1.7, cex.main=2, font.main=2, font.lab=2, font.axis=2, mar=c(5,5,5,1))
## pcex <- 0.5
## ppch <- 16
## plot(k36neg$laminB1, k36neg$pronecpg, xlim=c(-3,3), ylim=c(0.1,1), pch=ppch, cex=pcex, main='K36 negative', xlab='LaminB1', ylab='HP-CpG methylation')
## plot(k36pos$laminB1, k36pos$pronecpg, xlim=c(-3,3), ylim=c(0.1,1), pch=ppch, cex=pcex, main='K36 positive', xlab='LaminB1', ylab='HP-CpG methylation')
## plot(k36neg$repliseq, k36neg$pronecpg, xlim=c(0,90), ylim=c(0.1,1), pch=ppch, cex=pcex, main='K36 negative', xlab='RepliSeq', ylab='HP-CpG methylation')
## plot(k36pos$repliseq, k36pos$pronecpg, xlim=c(0,90), ylim=c(0.1,1), pch=ppch, cex=pcex, main='K36 positive', xlab='RepliSeq', ylab='HP-CpG methylation')
## dev.off()

## ## 30kb, relaxed
## set.seed(1)
## k36neg <- sample(gr0[gr0$H3K36me3<500], 1000)
## k36pos <- sample(gr0[gr0$H3K36me3>=500], 1000)
## pdf('~/gallery/2017_02_05_lamin_vs_meth_H3K36me3.pdf', width=7.5, height=4.2)
## par(mfrow=c(1,2), lwd=3, cex.lab=2, cex.axis=1.7, cex.main=2, font.main=2, font.lab=2, font.axis=2, mar=c(5,8,5,0.2), oma=c(0.3,0.2,0.3,1))
## pcex <- 0.5
## ppch <- 16
## plot(k36neg$laminB1, k36neg$pronecpg, xlim=c(-3,3), ylim=c(0.1,1), pch=ppch, cex=pcex, main='K36 negative', xlab='LaminB1', ylab='Solo-WCGW CpG\nMethylation')
## plot(k36pos$laminB1, k36pos$pronecpg, xlim=c(-3,3), ylim=c(0.1,1), pch=ppch, cex=pcex, main='K36 positive', xlab='LaminB1', ylab='Solo-WCGW CpG\nMethylation')
## ## plot(k36neg$repliseqG2, k36neg$pronecpg, xlim=c(0,90), ylim=c(0.1,1), pch=ppch, cex=pcex, main='K36 negative', xlab='RepliSeq', ylab='HP-CpG methylation')
## ## plot(k36pos$repliseqG2, k36pos$pronecpg, xlim=c(0,90), ylim=c(0.1,1), pch=ppch, cex=pcex, main='K36 positive', xlab='RepliSeq', ylab='HP-CpG methylation')
## dev.off()

## ## pw == 4.4
## load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/multiscale/signals.rda')
## set.seed(1)
## gr0 <- do.call(cbind, lapply(signals, function(x) x[,'4.4']))
## k36neg <- gr0[gr0[,'H3K36me3']<500,]
## k36pos <- gr0[gr0[,'H3K36me3']>=500,]
## k36neg <- k36neg[!is.na(rownames(k36neg)),]
## k36pos <- k36pos[!is.na(rownames(k36pos)),]
## k36neg <- k36neg[sample(1:nrow(k36neg), 2000),]
## k36pos <- k36pos[sample(1:nrow(k36pos), 2000),]
## pdf('~/gallery/2017_02_18_lamin_vs_meth_H3K36me3_4.4.pdf', width=5, height=6.5)
## par(mfrow=c(2,2), lwd=3, cex.lab=2, cex.axis=1.7, cex.main=2, font.main=2, font.lab=2, font.axis=2, mar=c(5,5,5,1))
## pcex <- 0.5
## ppch <- 16
## plot(k36neg[,'laminb1'], k36neg[,'pronecpg'], xlim=c(-3,3), ylim=c(0.1,1), pch=ppch, cex=pcex, main='K36 negative', xlab='LaminB1', ylab='HP-CpG methylation')
## plot(k36pos[,'laminb1'], k36pos[,'pronecpg'], xlim=c(-3,3), ylim=c(0.1,1), pch=ppch, cex=pcex, main='K36 positive', xlab='LaminB1', ylab='HP-CpG methylation')
## plot(k36neg[,'repliseqG2'], k36neg[,'pronecpg'], xlim=c(0,17000), ylim=c(0.1,1), pch=ppch, cex=pcex, main='K36 negative', xlab='RepliSeqG2', ylab='HP-CpG methylation')
## plot(k36pos[,'repliseqG2'], k36pos[,'pronecpg'], xlim=c(0,17000), ylim=c(0.1,1), pch=ppch, cex=pcex, main='K36 positive', xlab='RepliSeqG2', ylab='HP-CpG methylation')
## dev.off()

## ## 1MB
## pdf('~/gallery/2017_02_18_H3K36me3_1mb.pdf', width=4, height=2.6)
## plot(density(na.omit(gr0[,'H3K36me3'])))
## dev.off()

## load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/multiscale/signals.rda')
## set.seed(1)
## gr0 <- do.call(cbind, lapply(signals, function(x) x[,'6.0']))
## k36neg <- gr0[gr0[,'H3K36me3']<20000,]
## k36pos <- gr0[gr0[,'H3K36me3']>=20000,]
## pdf('~/gallery/2017_02_18_lamin_vs_meth_H3K36me3_1mb_relaxed_20k.pdf', width=6.5, height=6.9)
## par(mfrow=c(2,2), lwd=3, cex.lab=2, cex.axis=1.7, cex.main=2, font.main=2, font.lab=2, font.axis=2, mar=c(5,8,4,1))
## pcex <- 0.5
## ppch <- 16
## plot(k36neg[,'laminb1'], k36neg[,'pronecpg'], xlim=c(-3,3), ylim=c(0.1,1), pch=ppch, cex=pcex, main='K36 negative', xlab='LaminB1', ylab='Solo-WCGW CpG\nMethylation')
## plot(k36pos[,'laminb1'], k36pos[,'pronecpg'], xlim=c(-3,3), ylim=c(0.1,1), pch=ppch, cex=pcex, main='K36 positive', xlab='LaminB1', ylab='Solo-WCGW CpG\nMethylation')
## plot(k36neg[,'repliseqG2'], k36neg[,'pronecpg'], xlim=c(0,600000), ylim=c(0.1,1), pch=ppch, cex=pcex, main='K36 negative', xlab='RepliSeqG2', ylab='Solo-WCGW CpG\nMethylation')
## plot(k36pos[,'repliseqG2'], k36pos[,'pronecpg'], xlim=c(0,600000), ylim=c(0.1,1), pch=ppch, cex=pcex, main='K36 positive', xlab='RepliSeqG2', ylab='Solo-WCGW CpG\nMethylation')
## dev.off()

## ## 1MB for repliseq and lamin, 30kb for solo-wcgw
## load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/multiscale/signals.rda')
## a0 <- data.frame(h3k36me3=signals$H3K36me3[,'4.4'], pronecpg=signals$pronecpg[,'4.4'], laminb1=signals$laminb1[,'6.0'], repliseq=signals$repliseqG2[,'6.0'], pronecpgcnt=signals$pronecpgcnt[,'4.4'])
## a0 <- a0[complete.cases(a0) & a0$pronecpgcnt >= 10,]
## pdf('~/gallery/2017_03_31_lamin_vs_meth_H3K36me3_mixed2.pdf', width=10, height=10)
## par(mfrow=c(4,4), lwd=3, cex.lab=2, cex.axis=1.7, cex.main=2, font.main=1.5, font.lab=2, font.axis=2, mar=c(5,8,4,1))
## pcex <- 0.5
## ppch <- 16
## for (thres in c(50, 100, 200, 500)) {
##   k36neg <- a0[a0$h3k36me3 < thres,]
##   k36pos <- a0[a0$h3k36me3 >= thres,]
##   plot(k36neg[,'laminb1'], k36neg[,'pronecpg'], xlim=c(-2,2), ylim=c(0.1,1), pch=ppch, cex=pcex, main=sprintf('K36 (<%d)', thres), xlab='LaminB1', ylab='Solo-WCGW CpG\nMethylation')
##   plot(k36pos[,'laminb1'], k36pos[,'pronecpg'], xlim=c(-2,2), ylim=c(0.1,1), pch=ppch, cex=pcex, main=sprintf('K36 (>=%d)', thres), xlab='LaminB1', ylab='Solo-WCGW CpG\nMethylation')
##   plot(k36neg[,'repliseq'], k36neg[,'pronecpg'], xlim=c(0,90), ylim=c(0.1,1), pch=ppch, cex=pcex, main=sprintf('K36 (<%d)', thres), xlab='RepliSeq', ylab='Solo-WCGW CpG\nMethylation')
##   plot(k36pos[,'repliseq'], k36pos[,'pronecpg'], xlim=c(0,90), ylim=c(0.1,1), pch=ppch, cex=pcex, main=sprintf('K36 (>=%d)', thres), xlab='RepliSeq', ylab='Solo-WCGW CpG\nMethylation')
## }
## dev.off()


## a0 <- data.frame(h3k36me3=signals$H3K36me3[,'4.4'], pronecpg=signals$pronecpg[,'4.4'], laminb1=signals$laminb1[,'6.0'], repliseq=signals$repliseqG2[,'6.0'], pronecpgcnt=signals$pronecpgcnt[,'4.4'])
## a0 <- a0[complete.cases(a0) & a0$pronecpgcnt >= 10,]
## pdf('~/gallery/2017_03_31_lamin_vs_meth_H3K36me3_mixed2.pdf', width=10, height=10)
## par(mfrow=c(4,4), lwd=3, cex.lab=2, cex.axis=1.7, cex.main=2, font.main=1.5, font.lab=2, font.axis=2, mar=c(5,8,4,1))
## pcex <- 0.5
## ppch <- 16
## for (thres in c(50, 100, 200, 500)) {
##   k36neg <- a0[a0$h3k36me3 < thres,]
##   k36pos <- a0[a0$h3k36me3 >= thres,]
##   plot(k36neg[,'laminb1'], k36neg[,'pronecpg'], xlim=c(-2,2), ylim=c(0.1,1), pch=ppch, cex=pcex, main=sprintf('K36 (<%d)', thres), xlab='LaminB1', ylab='Solo-WCGW CpG\nMethylation')
##   plot(k36pos[,'laminb1'], k36pos[,'pronecpg'], xlim=c(-2,2), ylim=c(0.1,1), pch=ppch, cex=pcex, main=sprintf('K36 (>=%d)', thres), xlab='LaminB1', ylab='Solo-WCGW CpG\nMethylation')
##   plot(k36neg[,'repliseq'], k36neg[,'pronecpg'], xlim=c(0,600000), ylim=c(0.1,1), pch=ppch, cex=pcex, main=sprintf('K36 (<%d)', thres), xlab='RepliSeq', ylab='Solo-WCGW CpG\nMethylation')
##   plot(k36pos[,'repliseq'], k36pos[,'pronecpg'], xlim=c(0,600000), ylim=c(0.1,1), pch=ppch, cex=pcex, main=sprintf('K36 (>=%d)', thres), xlab='RepliSeq', ylab='Solo-WCGW CpG\nMethylation')
## }
## dev.off()

## a <- signals$pronecpg[,c('4.4','6.0')]
## a <- a[complete.cases(a),]
## plot(a[,1],a[,2], xlab='30kb solo-wcgw', ylab='1mb solo-wcgw')

## a <- data.frame(pronecpg=signals$pronecpg[,'4.4'], repliseq=signals$repliseqG2[,'6.0'])
## a <- a[complete.cases(a),]
## plot(a[,2],a[,1], xlim=c(0,600000), xlab='1mb repliseqG2', ylab='30kb solo-wcgw')

## a <- data.frame(pronecpg=signals$pronecpg[,'6.0'], repliseq=signals$repliseqG2[,'6.0'])
## a <- a[complete.cases(a),]
## plot(a[,2],a[,1], xlim=c(0,600000), xlab='1mb repliseqG2', ylab='1mb solo-wcgw')

## a <- data.frame(pronecpg=signals$pronecpg[,'4.4'], repliseq=signals$repliseqG2[,'5.0'], genebody=signals$genebodyoverlap[,'4.4'], h3k36=signals$H3K36me3[,'4.4'], laminb1=signals$laminb1[,'5.0'])
## a <- a[complete.cases(a),]
## ggplot(a, aes(repliseq, pronecpg, color=genebody)) + geom_point() + xlim(0,600000)
## plot(a[,2],a[,1], xlim=c(0,600000), xlab='1mb repliseqG2', ylab='30kb solo-wcgw')

## ## 100kb repliseq G2 vs 30kb pronecpg
## a <- data.frame(pronecpg=signals$pronecpg[,'4.4'], repliseqg2=signals$repliseqG2[,'5.0'], repliseq=signals$repliseq[,'5.0'], genebody=signals$genebodyoverlap[,'4.4'], h3k36=signals$H3K36me3[,'4.4'], laminb1=signals$laminb1[,'5.0'], exon=signals$exonoverlap[,'4.4'], cnt=signals$pronecpgcnt[,'4.4'])
## a <- a[complete.cases(a) & a$cnt >10 ,]
## jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
## pdf('~/gallery/2017_03_31_lamin_vs_cpg.pdf', width=7, height=5)
## ggplot(a, aes(laminb1, pronecpg, color=log2(1+pmax(3,h3k36)))) + geom_point(size=0.2) + scale_color_gradientn(colours = jet.colors(7))
## dev.off()

## pdf('~/gallery/2017_03_31_repliseq_vs_cpg.pdf', width=7, height=5)
## ggplot(a, aes(repliseq, pronecpg, color=log2(1+pmax(3,h3k36)))) + geom_point(size=0.2) + scale_color_gradientn(colours = jet.colors(7))
## dev.off()

## pdf('~/gallery/2017_03_31_repliseqg2_vs_cpg.pdf', width=7, height=5)
## ggplot(a, aes(repliseqg2, pronecpg, color=log2(1+pmax(3,h3k36)))) + geom_point(size=0.2) + scale_color_gradientn(colours = jet.colors(7)) + xlim(0,60000)
## dev.off()

## pdf('~/gallery/2017_03_31_repliseqg2_vs_repliseq.pdf', width=7, height=5)
## ggplot(a, aes(repliseqg2, repliseq, color=log2(1+pmax(3,h3k36)))) + geom_point(size=0.2) + scale_color_gradientn(colours = jet.colors(7)) + xlim(0,60000)
## dev.off()

## ## 1mb repliseq G2 vs 30kb pronecpg
## a <- data.frame(pronecpg=signals$pronecpg[,'4.4'], repliseqg2=signals$repliseqG2[,'6.0'], repliseq=signals$repliseq[,'6.0'], genebody=signals$genebodyoverlap[,'4.4'], h3k36=signals$H3K36me3[,'4.4'], laminb1=signals$laminb1[,'6.0'], exon=signals$exonoverlap[,'4.4'])
## a <- a[complete.cases(a),]
## jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
## pdf('~/gallery/2017_03_31_1mb_lamin_vs_cpg.pdf', width=7, height=5)
## ggplot(a, aes(laminb1, pronecpg, color=log2(1+pmax(3,h3k36)))) + geom_point(size=2) + scale_color_gradientn(colours = jet.colors(7))
## dev.off()

## pdf('~/gallery/2017_03_31_1mb_repliseq_vs_cpg.pdf', width=7, height=5)
## ggplot(a, aes(repliseq, pronecpg, color=log2(1+pmax(3,h3k36)))) + geom_point(size=2) + scale_color_gradientn(colours = jet.colors(7))
## dev.off()

## pdf('~/gallery/2017_03_31_1mb_repliseqg2_vs_cpg.pdf', width=7, height=5)
## ggplot(a, aes(repliseqg2, pronecpg, color=log2(1+pmax(3,h3k36)))) + geom_point(size=2) + scale_color_gradientn(colours = jet.colors(7)) + xlim(0,60000)
## dev.off()

## pdf('~/gallery/2017_03_31_1mb_repliseqg2_vs_repliseq.pdf', width=7, height=5)
## ggplot(a, aes(repliseqg2, repliseq, color=log2(1+pmax(3,h3k36)))) + geom_point(size=2) + scale_color_gradientn(colours = jet.colors(7)) + xlim(0,60000)
## dev.off()

## 30kb bin two color version
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/multiscale/signals.rda')
a0 <- data.frame(pronecpg=signals$pronecpg[,'4.5'],
                 repliseq=signals$repliseq[,'4.5'],
                 h3k36=signals$H3K36me3[,'4.5'],
                 laminb1=signals$laminb1[,'4.5'],
                 k36ov = signals$H3K36me3ov[,'4.5'],
                 geneov = signals$genebodyoverlap[,'4.5'],
                 cnt=signals$pronecpgcnt[,'4.5'])
a <- a0[complete.cases(a0) & a0$cnt >= 10,]
a$k36 <- ifelse(a$geneov > 20000 & a$k36ov > 20000, 'red', 'blue')
a$k36 <- ifelse(a$k36 == 'blue' & !(a$geneov == 0 & a$k36ov == 0), NA, a$k36)
a$k36[a$geneov > 25000 & a$k36ov == 0] <- 'green'
table(a$k36)
## FALSE  TRUE
## 14048  4380
a1 <- a[!is.na(a),]
pdf('~/gallery/2017_04_03_laminb1_vs_solowcgw_imr90.pdf', width=4, height=3.7)
par(mar=c(5,7,1,1), lwd=3, cex.lab=2, cex.axis=2, font.lab=2, font.axis=2)
plot(a1$laminb1, a1$pronecpg, col=a1$k36, pch=16, cex=0.18, xlab='Lamin B1', ylab='Solo-WCGW\nMethylation')
dev.off()

pdf('~/gallery/2017_04_03_laminb1_vs_solowcgw_imr90_sbs.pdf', width=8, height=3.7)
par(mar=c(5,7,1,1), lwd=3, cex.lab=2, cex.axis=2, font.lab=2, font.axis=2, mfrow=c(1,3))
a2 <- a1[a1$k36=='red',]
plot(a2$laminb1, a2$pronecpg, col=a2$k36, pch=16, cex=0.18, xlab='Lamin B1', ylab='Solo-WCGW\nMethylation')
a2 <- a1[a1$k36=='blue',]
plot(a2$laminb1, a2$pronecpg, col=a2$k36, pch=16, cex=0.18, xlab='Lamin B1', ylab='Solo-WCGW\nMethylation')
a2 <- a1[a1$k36=='green',]
plot(a2$laminb1, a2$pronecpg, col=a2$k36, pch=16, cex=0.18, xlab='Lamin B1', ylab='Solo-WCGW\nMethylation')
dev.off()

pdf('~/gallery/2017_04_03_repliseq_vs_solowcgw_imr90.pdf', width=4, height=3.7)
par(mar=c(5,7,1,1), lwd=3, cex.lab=2, cex.axis=2, font.lab=2, font.axis=2)
plot(100-a1$repliseq, a1$pronecpg, col=a1$k36, pch=16, cex=0.18, xlab='RepliSeq', ylab='Solo-WCGW\nMethylation')
dev.off()

pdf('~/gallery/2017_04_03_repliseq_vs_solowcgw_imr90_sbs.pdf', width=8, height=3.7)
par(mar=c(5,7,1,1), lwd=3, cex.lab=2, cex.axis=2, font.lab=2, font.axis=2, mfrow=c(1,3))
a2 <- a1[a1$k36=='red',]
plot(100-a2$repliseq, a2$pronecpg, col=a2$k36, pch=16, cex=0.18, xlab='Lamin B1', ylab='Solo-WCGW\nMethylation')
a2 <- a1[a1$k36=='blue',]
plot(100-a2$repliseq, a2$pronecpg, col=a2$k36, pch=16, cex=0.18, xlab='Lamin B1', ylab='Solo-WCGW\nMethylation')
a2 <- a1[a1$k36=='green',]
plot(100-a2$repliseq, a2$pronecpg, col=a2$k36, pch=16, cex=0.18, xlab='Lamin B1', ylab='Solo-WCGW\nMethylation')
dev.off()

## ggplot(a) + geom_point(aes(laminb1, pronecpg, color=k36), size=0.5)

## ggplot(a, aes(laminb1, pronecpg, color=log2(1+genebody))) + geom_point(size=0.2) + scale_color_gradient2(low='blue',mid='white',high='red',midpoint=9)
## ggplot(a, aes(laminb1, pronecpg, color=log2(1+exon))) + geom_point(size=0.2) + scale_color_gradient2(low='blue',mid='white',high='red',midpoint=9)
## hist(a$genebody)
## ggplot(a[a$genebody>24000,], aes(laminb1, pronecpg)) + geom_point(size=0.2)
## ggplot(a[a$genebody>500,], aes(laminb1, pronecpg)) + geom_point(size=0.2)
## ggplot(a, aes(genebody, h3k36)) + geom_point(size=0.5)

## pdf('~/gallery/2017_01_30_lamin_vs_meth_H3K36me3_smooth.pdf', width=6.5, height=6.5)
## par(mfrow=c(2,2))
## par(mar=c(5,5,4,1))
## par(cex.lab=1.4, cex.axis=1.2, cex.main=1.5, font.lab=2, font.axis=2, lwd=1.5)
## pcex <- 0.2
## smoothScatter(gr0[gr0$H3K36me3<300]$laminB1, gr0[gr0$H3K36me3<300]$pronecpg, xlim=c(-3,3), pch=20, cex=pcex, main='K36 negative', xlab='LaminB1', ylab='HP-CpG methylation', colramp=colorRampPalette(c('white','black')), nrpoints=0)
## smoothScatter(gr0[gr0$H3K36me3>=300]$laminB1, gr0[gr0$H3K36me3>=300]$pronecpg, xlim=c(-3,3), pch=20, cex=pcex, main='K36 positive', xlab='LaminB1', ylab='HP-CpG methylation', colramp=colorRampPalette(c('white','black')), nrpoints=0)
## smoothScatter(gr0[gr0$H3K36me3<300]$repliseq, gr0[gr0$H3K36me3<300]$pronecpg, pch=20, cex=pcex, main='K36 negative', ylab='HP-CpG methylation', xlab='RepliSeq', colramp=colorRampPalette(c('white','black')), nrpoints=0)
## smoothScatter(gr0[gr0$H3K36me3>=300]$repliseq, gr0[gr0$H3K36me3>=300]$pronecpg, pch=20, cex=pcex, main='K36 positive', ylab='HP-CpG methylation', xlab='RepliSeq', colramp=colorRampPalette(c('white','black')), nrpoints=0)
## dev.off()

###################################################
## other REMC samples doing K36 analysis as IMR90
###################################################

Enumbers <- c('E003', 'E005', 'E007', 'E011', 'E013', 'E017', 'E022', 'E050', 'E054', 'E065', 'E070', 'E079', 'E085', 'E095', 'E097', 'E100', 'E105', 'E109', 'E113','E004', 'E006', 'E008', 'E012', 'E016', 'E021', 'E024', 'E053', 'E058', 'E066', 'E071', 'E084', 'E094', 'E096', 'E098', 'E104', 'E106', 'E112')

library(plyr)
convertrds <- function(dname, oname) {
  setwd(dname)
  pws <- sprintf('%1.1f', seq(2,7,0.1))
  a <- lapply(pws, function(pw) {
    a <- read.table(paste0('scale_',pw,'.bed'))
    score <- setNames(as.numeric(a$V6), paste0(a$V1,':',a$V4,'-',a$V5))
    score
  })
  measure <- Reduce(function(x,y) t(rbind.fill.matrix(t(x), t(y))), a)
  colnames(measure) <- pws
  saveRDS(measure, file=paste0('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/REMC_H3K36me3/multiscale/',oname,'.rds'))
}

## K36
for (Enumber in Enumbers) {
  convertrds(sprintf('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/REMC_H3K36me3/K36_finalcnt/%s/', Enumber), sprintf('%s.K36', Enumber))
}

## pronecpg
trash <- mclapply(Enumbers, function(Enumber) {
  convertrds(sprintf('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/REMC_H3K36me3/dname_final/%s/', Enumber), sprintf('%s.dname', Enumber))
}, mc.cores=20)

## replichip for H1 and H9
convertrds('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/REMC_H3K36me3/ENCODE_replichip/H1_final/ENCFF000KUF/', 'replichip.H1.ENCFF000KUF')
convertrds('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/REMC_H3K36me3/ENCODE_replichip/H1_final/ENCFF000KUG/', 'replichip.H1.ENCFF000KUG')
convertrds('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/REMC_H3K36me3/ENCODE_replichip/H1_final/ENCFF000KUH/', 'replichip.H1.ENCFF000KUH')
convertrds('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/REMC_H3K36me3/ENCODE_replichip/H9_final/ENCFF000KUN/', 'replichip.H9')

## sdbinmean100kb
convertrds('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/REMC_H3K36me3/binmeansd/binmeansd100kb_final/', 'sdbinmean100kb')

## all signals
ffiles <- list.files('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/REMC_H3K36me3/multiscale/','.rds')
signals <- mclapply(ffiles, function(oname) { readRDS(paste0('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/REMC_H3K36me3/multiscale/',oname)); }, mc.cores=20)
names(signals) <- sub('.rds','',ffiles)
signals$genebodyov <- readRDS('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/multiscale/genebodyoverlap.rds')
signals$pronecpgcnt <- readRDS('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/multiscale/pronecpgcnt.rds')
save(signals, file='/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/REMC_H3K36me3/binmeansd/binmeansd100kb_final/signals.rda')

## plot H1
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/REMC_H3K36me3/binmeansd/binmeansd100kb_final/signals.rda')
a0 <- data.frame(repliseq=signals$replichip.H1.ENCFF000KUF[,'4.5'],
                 dname=signals$E003.dname[,'4.5'],
                 k36ov=signals$E003.K36[,'4.5'],
                 geneov=signals$genebodyov[,'4.5'],
                 cnt=signals$pronecpgcnt[,'4.5'])
a <- a0[complete.cases(a0) & a0$cnt >= 10,]
a$k36 <- ifelse(a$geneov > 20000 & a$k36ov > 20000, TRUE, FALSE)
a$k36 <- ifelse(!a$k36 & !(a$geneov == 0 & a$k36ov == 0), NA, a$k36)
a1 <- a[!is.na(a),]
png('~/gallery/2017_04_03_H1_repliseq_vs_methylation.png', width=900, height=800)
par(mar=c(10,13,5,1), lwd=5, cex.lab=3, cex.axis=3, font.lab=2, font.axis=2, cex.main=2, mgp=c(7,3,0))
plot(-a1$repliseq, a1$dname, col=ifelse(a1$k36,'red','blue'), pch=16, cex=0.5, xlab='RepliChIP', ylab='Solo-WCGW\nMethylation', xaxt='n', yaxt='n', ylim=c(0,1))
axis(1, at=c(-2,-1,0,1,2), labels=c(-2,-1,0,1,2), lwd=3)
axis(2, at=1:5/5, labels=1:5/5, lwd=3, las=2)
legend('bottomleft', legend=paste0(c('Active Gene Bodies (', 'Inactive Intergenic ('), rev(table(a$k36)), ')'), col=c('red','blue'), pch=16, cex=3)
dev.off()

pdf('~/gallery/2017_04_03_H1_repliseq_vs_methylation.pdf', width=4, height=3.7)
par(mar=c(5,7,1,1), lwd=3, cex.lab=2, cex.axis=2, font.lab=2, font.axis=2)
## par(mar=c(10,13,5,1), lwd=5, cex.lab=3, cex.axis=3, font.lab=2, font.axis=2, cex.main=2, mgp=c(7,3,0))
plot(-a1$repliseq, a1$dname, col=ifelse(a1$k36,'red','blue'), pch=16, cex=0.18, xlab='RepliChIP', ylab='Solo-WCGW\nMethylation', xaxt='n', yaxt='n', ylim=c(0,1))
axis(1, at=c(-2,-1,0,1,2), labels=c(-2,-1,0,1,2), lwd=3)
axis(2, at=1:5/5, labels=1:5/5, lwd=3, las=2)
legend('bottomleft', legend=paste0(c('Active Gene Bodies (', 'Inactive Intergenic ('), rev(table(a$k36)), ')'), col=c('red','blue'), pch=16, cex=0.5)
dev.off()

## plot H9
a0 <- data.frame(repliseq=signals$replichip.H9[,'4.5'],
                 dname=signals$E008.dname[,'4.5'],
                 k36ov=signals$E008.K36[,'4.5'],
                 geneov=signals$genebodyov[,'4.5'],
                 cnt=signals$pronecpgcnt[,'4.5'])
a <- a0[complete.cases(a0) & a0$cnt >= 10,]
a$k36 <- ifelse(a$geneov > 20000 & a$k36ov > 20000, TRUE, FALSE)
a$k36 <- ifelse(!a$k36 & !(a$geneov == 0 & a$k36ov == 0), NA, a$k36)
a1 <- a[!is.na(a),]
png('~/gallery/2017_04_03_H9_repliseq_vs_methylation.png', width=900, height=800)
par(mar=c(10,13,5,1), lwd=5, cex.lab=3, cex.axis=3, font.lab=2, font.axis=2, cex.main=2, mgp=c(7,3,0))
plot(-a1$repliseq, a1$dname, col=ifelse(a1$k36,'red','blue'), pch=16, cex=0.5, xlab='RepliChIP', ylab='Solo-WCGW\nMethylation', xaxt='n', yaxt='n', ylim=c(0,1))
axis(1, at=c(-2,-1,0,1,2), labels=c(-2,-1,0,1,2), lwd=3)
axis(2, at=1:5/5, labels=1:5/5, lwd=3, las=2)
legend('bottomleft', legend=paste0(c('Active Gene Bodies (', 'Inactive Intergenic ('), rev(table(a$k36)), ')'), col=c('red','blue'), pch=16, cex=3)
dev.off()

## plot other REMC samples
Enumbers <- c('E003', 'E005', 'E007', 'E011', 'E013', 'E017', 'E022', 'E050', 'E054', 'E065', 'E070', 'E079', 'E085', 'E095', 'E097', 'E100', 'E105', 'E109', 'E113','E004', 'E006', 'E008', 'E012', 'E016', 'E021', 'E024', 'E053', 'E058', 'E066', 'E071', 'E084', 'E094', 'E096', 'E098', 'E104', 'E106', 'E112')
## mapping Enumber to name
df <- read.table('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/REMC_H3K36me3/samples')
Enumber2name <- setNames(df$V2, df$V1)
for (Enumber in Enumbers) {
  a0 <- data.frame(sdbinmeanT=signals$sdbinmean100kb[,'4.5'],
                   dname=signals[[paste0(Enumber,'.dname')]][,'4.5'],
                   k36ov=signals[[paste0(Enumber,'.K36')]][,'4.5'],
                   geneov=signals$genebodyov[,'4.5'],
                   cnt=signals$pronecpgcnt[,'4.5'])
  a <- a0[complete.cases(a0) & a0$cnt >= 10,]
  a$k36 <- ifelse(a$geneov > 20000 & a$k36ov > 20000, TRUE, FALSE)
  a$k36 <- ifelse(!a$k36 & !(a$geneov == 0 & a$k36ov == 0), NA, a$k36)
  table(a$k36)
  a1 <- a[!is.na(a),]
  png(sprintf('~/gallery/2017_04_03_REMC_%s_binmeanSD_vs_methylation.png', Enumber), width=800, height=800)
  par(mar=c(10,11,5,1), lwd=5, cex.lab=3, cex.axis=3, font.lab=2, font.axis=2, cex.main=2, mgp=c(5,1,0))
  plot(a1$sdbinmeanT, a1$dname, col=ifelse(a1$k36,'red','blue'), pch=16, cex=0.8, xlab='SD bin average\n(Tumor)', ylab='Solo-WCGW\nMethylation', xaxt='n', yaxt='n', ylim=c(0,1), main=sprintf('%s\n%s', Enumber, Enumber2name[Enumber]))
  axis(1, at=c(-2,-1,0,1,2), labels=c(-2,-1,0,1,2), lwd=3)
  axis(2, at=1:5/5, labels=1:5/5, lwd=3, las=2)
  legend('bottomleft', legend=paste0(c('Active Gene Bodies (', 'Inactive Intergenic ('), rev(table(a$k36)), ')'), col=c('red','blue'), pch=16, cex=3)
  dev.off()
}

#############################
## match cell line with REMC
#############################

library(plyr)
convertrds <- function(dname) {
  setwd(dname)
  pws <- sprintf('%1.1f', seq(2,7,0.1))
  a <- lapply(pws, function(pw) {
    a <- read.table(paste0('scale_',pw,'.bed'))
    score <- setNames(as.numeric(a$V6), paste0(a$V1,':',a$V4,'-',a$V5))
    score
  })
  measure <- Reduce(function(x,y) t(rbind.fill.matrix(t(x), t(y))), a)
  colnames(measure) <- pws
  measure
}

gsms <- c('GSM923439', 'GSM923440', 'GSM923441', 'GSM923442', 'GSM923443', 'GSM923444', 'GSM923445', 'GSM923446', 'GSM923447', 'GSM923448', 'GSM923449', 'GSM923450', 'GSM923451', 'GSM923452', 'GSM923453')
gsmsnames <- c('GSM923439_GM12812', 'GSM923440_GM12801', 'GSM923441_SKNSH', 'GSM923442_MCF7', 'GSM923443_GM06990', 'GSM923444_BJ', 'GSM923445_NHEK', 'GSM923446_HepG2', 'GSM923447_IMR90', 'GSM923448_K562', 'GSM923449_HeLaS3', 'GSM923450_GM12813', 'GSM923451_GM12878', 'GSM923452_HUVEC', 'GSM923453_BG02ES')
signals1 <- lapply(gsms, function(gsm) convertrds(paste0('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/REMC_H3K36me3/ENCODE_repliseq/final/', gsm)))
names(signals1) <- gsmsnames

gsmsnames <- c(
  ## 'GSM923439_GM12812'='',
  ## 'GSM923440_GM12801'=,
  ## 'GSM923450_GM12813'=,
  ## 'GSM923444_BJ'=,
  ## 'GSM923447_IMR90'=,
  ## 'GSM923448_K562'='',
  ## 'GSM923453_BG02ES'=''
  'GSM923443_GM06990'='BP_venous_blood_S001JP51_CD38negative_naive_B_cell_DiseaseFree',
  'GSM923451_GM12878'='BP_venous_blood_csMBC_NC11_41_class_switched_memory_B_cell_DiseaseFree',
  'GSM923441_SKNSH'='REMC_E071_Brain_Hippocampus_Middle',
  'GSM923442_MCF7'='REMC_UCSF-UBC.Breast_Myoepithelial_Cells.RM066',
  'GSM923445_NHEK'='REMC_E058_Penis_Foreskin_Keratinocyte_Primary_Cells_skin03',
  'GSM923446_HepG2'='REMC_E066_Adult_Liver',
  'GSM923449_HeLaS3'='REMC_E097_Ovary',
  'GSM923452_HUVEC'='BP_cord_blood_S00BJM51_endothelial_cell_of_umbilical_vein_proliferating_DiseaseFree')
signals2 <- lapply(gsmsnames, function(gsm) convertrds(paste0('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/REMC_H3K36me3/ENCODE_repliseq/final_dname/', gsm)))
names(signals2) <- gsmsnames
signals3 <- c(signals1, signals2)
signals3$GM12878.K36 <- convertrds('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/REMC_H3K36me3/ENCODE_repliseq/final_GM12878_K36/')
signals3 <- signals
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/REMC_H3K36me3/binmeansd/binmeansd100kb_final/signals.rda')
signals3$genebodyov <- signals$genebodyov
signals3$pronecpgcnt <- signals$pronecpgcnt
signals <- signals3
save(signals, file='/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/REMC_H3K36me3/ENCODE_repliseq/multiscale/signals.rda')


a0 <- data.frame(repliseq=signals$GSM923451_GM12878[,'4.5'],
                 solowcgw=signals$BP_venous_blood_S001JP51_CD38negative_naive_B_cell_DiseaseFree[,'4.5'],
                 k36ov=signals$GM12878.K36[,'4.5'],
                 geneov=signals$genebodyov[,'4.5'],
                 cnt=signals$pronecpgcnt[,'4.5'])
a <- a0[complete.cases(a0) & a0$cnt >= 10,]
a$k36 <- ifelse(a$geneov > 20000 & a$k36ov > 20000, TRUE, FALSE)
a$k36 <- ifelse(!a$k36 & !(a$geneov == 0 & a$k36ov == 0), NA, a$k36)
table(a$k36)
a1 <- a[!is.na(a),]
png('~/gallery/2017_04_03_naiveB_GM12878.png', width=800, height=800)
par(mar=c(10,11,5,1), lwd=5, cex.lab=3, cex.axis=3, font.lab=2, font.axis=2, cex.main=2, mgp=c(5,1,0))
plot(-a1$repliseq, a1$solowcgw, col=ifelse(a1$k36,'red','blue'), pch=16, cex=0.3, xlab='RepliSeq', ylab='Solo-WCGW\nMethylation', yaxt='n', ylim=c(0,1))
## axis(1, at=c(-2,-1,0,1,2), labels=c(-2,-1,0,1,2), lwd=3)
axis(2, at=0:5/5, labels=0:5/5, lwd=3, las=2)
legend('bottomleft', legend=paste0(c('Active Gene Bodies (', 'Inactive Intergenic ('), rev(table(a$k36)), ')'), col=c('red','blue'), pch=16, cex=3)
dev.off()

a0 <- data.frame(repliseq=signals$GSM923451_GM12878[,'4.5'],
                 solowcgw=signals$BP_venous_blood_csMBC_NC11_41_class_switched_memory_B_cell_DiseaseFree[,'4.5'],
                 k36ov=signals$GM12878.K36[,'4.5'],
                 geneov=signals$genebodyov[,'4.5'],
                 cnt=signals$pronecpgcnt[,'4.5'])
a <- a0[complete.cases(a0) & a0$cnt >= 10,]
a$k36 <- ifelse(a$geneov > 20000 & a$k36ov > 20000, TRUE, FALSE)
a$k36 <- ifelse(!a$k36 & !(a$geneov == 0 & a$k36ov == 0), NA, a$k36)
table(a$k36)
a1 <- a[!is.na(a),]
png('~/gallery/2017_04_03_memoryB_GM12878.png', width=800, height=800)
par(mar=c(10,11,5,1), lwd=5, cex.lab=3, cex.axis=3, font.lab=2, font.axis=2, cex.main=2, mgp=c(5,1,0))
plot(-a1$repliseq, a1$solowcgw, col=ifelse(a1$k36,'red','blue'), pch=16, cex=0.3, xlab='RepliSeq', ylab='Solo-WCGW\nMethylation', yaxt='n', ylim=c(0,1))
## axis(1, at=c(-2,-1,0,1,2), labels=c(-2,-1,0,1,2), lwd=3)
axis(2, at=0:5/5, labels=0:5/5, lwd=3, las=2)
legend('bottomleft', legend=paste0(c('Active Gene Bodies (', 'Inactive Intergenic ('), rev(table(a$k36)), ')'), col=c('red','blue'), pch=16, cex=3)
dev.off()

#############################
## chr16 short arm IMR90
#############################
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/gr.rda')
gr1 <- gr0[seqnames(gr0) == 'chr16' & start(gr0) < 35000000]
pdf('~/gallery/2017_01_28_IMR90_allmarks.pdf', width=10, height=10)
par(mfrow=c(12,1))
par(mar=c(0,5,0,0), oma=c(3,3,3,3))
plot(gr1$pronecpg, type='l', ylab='Prone CpG')
plot(gr1$laminB1, type='l', ylab='Lamin B1')
plot(gr1$repliseq, type='l', ylab='Replication')
plot(gr1$H3K27me3, type='l', ylab='H3K27me3')
plot(gr1$H3K36me3, type='l', ylab='H3K36me3')
plot(gr1$DNase, type='l', ylab='DNase')
plot(gr1$H3K9me3, type='l', ylab='H3K9me3')
plot(gr1$H3K4me3, type='l', ylab='H3K4me3')
plot(gr1$H3K27ac, type='l', ylab='H3K27ac')
plot(gr1$H3K4me1, type='l', ylab='H3K4me1')
plot(gr1$expression, type='l', ylab='expression')
plot(gr1$genedensity, type='l', ylab='GeneDensity')
dev.off()

###########################################
## IMR90 multiscale, overlapping bins
###########################################

####### load data and merge
library(plyr)
convertrds <- function(dname, oname) {
  setwd(dname)
  pws <- sprintf('%1.1f', seq(2,7,0.1))
  a <- lapply(pws, function(pw) {
    a <- read.table(paste0('scale_',pw,'.bed'))
    score <- setNames(as.numeric(a$V6), paste0(a$V1,':',a$V4,'-',a$V5))
    score
  })
  measure <- Reduce(function(x,y) t(rbind.fill.matrix(t(x), t(y))), a)
  colnames(measure) <- pws
  saveRDS(measure, file=paste0('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/multiscaleOv/',oname,'.rds'))
}

## pronecpg
convertrds('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/pronecpg/finalOv/', 'pronecpg')

## rnaseq
convertrds('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/RNAseq/finalOv/', 'rnaseq')

## CTCF
convertrds('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/CTCF/finalOv/', 'ctcf')

## Lamin B1
convertrds('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/LaminB1/finalOv/', 'laminb1')
ob <- readRDS('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/multiscaleOv/laminb1.rds')
ob[ob < -3] <- NA
ob[ob > 3] <- NA
saveRDS(ob, file='/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/multiscaleOv/laminb1.rds')

## DNase-seq
convertrds('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/DNase/finalOv/', 'dnase')

## RepliSeq
convertrds('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/repliseq/finalOv/wavelet', 'repliseq')
## for (phase in c('G1b','S1','S2','S3','S4','G2')) {
##   cat(phase,'\n')
##   convertrds(paste0('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/repliseq/finalOv/',phase), paste0('repliseq',phase))
## }

## histone marks
for (mark in c('H3K9me3','H3K27ac','H3K27me3','H3K4me3','H3K4me1','H3K36me3')) {
  cat(mark,'\n')
  convertrds(paste0('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/', mark, '/finalOv/'), mark)
}

## merge
## features <- c('ctcf', 'H3K27ac', 'H3K36me3', 'H3K4me3', 'laminb1', 'rnaseq', 'dnase', 'H3K27me3', 'H3K4me1', 'H3K9me3', 'pronecpg', 'repliseqG1b', 'repliseqG2', 'repliseqS1', 'repliseqS2', 'repliseqS3', 'repliseqS4')
## subfeatures <- c('pronecpg', 'laminb1', 'repliseqG2',
##                  'H3K36me3', 'H3K4me3', 'repliseqG1b', 
##                  'repliseqS1', 'H3K4me1', 'H3K27ac', 'dnase',
##                   'H3K27me3', 'repliseqS2', 'rnaseq',
##                  'repliseqS3', 'repliseqS4', 'H3K9me3')
features <-  c('dnase', 'H3K27ac', 'H3K27me3', 'H3K36me3', 'H3K4me1', 'H3K4me3', 'H3K9me3', 'laminb1', 'pronecpg', 'repliseq', 'rnaseq')
toofficial <- setNames(features, features)
toofficial['ctcf'] <- 'CTCF'
toofficial['dnase'] <- 'DNase-Seq'
toofficial['rnaseq'] <- 'RNA-Seq'
toofficial['laminb1'] <- 'Lamin B1'
toofficial['pronecpg'] <- 'HP-CpG methylation'
toofficial <- sub('repliseq','RepliSeq ',toofficial)

signals <- lapply(features, function(oname) {
  readRDS(paste0('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/multiscaleOv/',oname,'.rds'))
})
names(signals) <- features
save(signals, features, toofficial, file='/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/multiscaleOv/signals.rda')

###### get spearman's correlation

load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/multiscaleOv/signals.rda')

getspearmans <- function(o1, o2) {
  pws <- sprintf('%1.1f', seq(2,7,0.1))
  multicors <- sapply(1:51, function(i) {
    sapply(1:51, function(j) {
      a <- o1[,i]
      b <- o2[,j]
      cor(a, b, use='na.or.complete', method='spearman')
    })
  })
  colnames(multicors) <- pws
  rownames(multicors) <- pws
  multicors
}

allcors <- mclapply(signals, function(sig1) {
  a <- mclapply(signals, function(sig2) {
    getspearmans(sig1, sig2)
  }, mc.cores=6)
  names(a) <- features
  a
}, mc.cores=6)
names(allcors) <- features
save(allcors, features, toofficial, file='/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/multiscaleOv/allspearman.rda')

###################################
##### plot multiscale correlation2
###################################

## ##### plot multiscale correlation

## library(RColorBrewer)
## plotcors <- function(cors, fn, xlab, ylab='HP-CpG methylation') {
##   pdf(fn, width=4.4, height=3)
##   ## png(fn, width=1000, height=800)
##   par(mar=c(5,5,2,1), lwd=1)
##   filled.contour(x=1:51, y=1:51, z=t(cors), zlim=c(-1,1),
##                  color.palette=colorRampPalette(rev(brewer.pal(11, 'RdBu'))),
##                  plot.axes={
##                    contour(x=1:51, y=1:51, z=t(cors), nlevels = 7, drawlabels = FALSE, axes = FALSE, frame.plot = FFALSE, add = TRUE);
##                    axis(1, which(seq(2,7,0.1) %in% c(2:7)), sapply(2:7, function(i) parse(text=paste0('bold(10^',i,')'))));
##                    axis(2, which(seq(2,7,0.1) %in% c(2:7)), sapply(2:7, function(i) parse(text=paste0('bold(10^',i,')'))));
##                  }, cex.axis=2, xlab=xlab, ylab=ylab, cex.lab=1.3, font.lab=2, font.axis=2, lwd=2, nlevels=20)
##   dev.off()
## }

## for (ft1 in features) {
##   for (ft2 in 'pronecpg') {
##     plotcors(allcors[[ft1]][[ft2]], sprintf('~/gallery/2017_02_02_allcors_multiscale/%s_vs_%s.pdf', ft2, ft1), toofficial[ft1], toofficial[ft2])
##   }
## }

load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/multiscaleOv/allspearman.rda')
library(RColorBrewer)
plotcors <- function(cors, fn, xlab, ylab='Solo-WCGW CpG Methylation') {
  pdf(fn, width=4.4, height=3)
  ## png(fn, width=1000, height=800)
  par(mar=c(5,5,2,1), lwd=1)
  filled.contour(x=1:51, y=1:51, z=t(cors), zlim=c(-1,1),
                 color.palette=colorRampPalette(rev(brewer.pal(11, 'RdBu'))),
                 plot.axes={
                   contour(x=1:51, y=1:51, z=t(cors), nlevels = 7, drawlabels = FALSE, axes = FALSE, frame.plot = FFALSE, add = TRUE);
                   axis(1, which(seq(2,7,0.1) %in% c(2:7)), sapply(2:7, function(i) parse(text=paste0('bold(10^',i,')'))));
                   axis(2, which(seq(2,7,0.1) %in% c(2:7)), sapply(2:7, function(i) parse(text=paste0('bold(10^',i,')'))));
                 }, cex.axis=2, xlab=xlab, ylab=ylab, cex.lab=1.3, font.lab=2, font.axis=2, lwd=2, nlevels=20)
  dev.off()
}

features <- c('pronecpg','laminb1','repliseq','dnase','H3K27ac','H3K4me1','H3K27me3','rnaseq','H3K36me3','H3K4me3','H3K9me3')
toofficial <- setNames(features, features)
toofficial['ctcf'] <- 'CTCF'
toofficial['dnase'] <- 'DNase-Seq'
toofficial['rnaseq'] <- 'RNA-Seq'
toofficial['laminb1'] <- 'Lamin B1'
toofficial['pronecpg'] <- 'Solo-WCGW CpG\nMethylation'
toofficial <- sub('repliseq','RepliSeq ',toofficial)

for (ft1 in features) {
  for (ft2 in 'pronecpg') {
    plotcors(allcors[[ft1]][[ft2]], sprintf('~/gallery/2017_04_04_allcors_multiscale/%s_vs_%s.pdf', ft2, ft1), toofficial[ft1], toofficial[ft2])
  }
}

########################################################
### multiscale heatmap plot, from overlapping bin data
########################################################
## laminB1 
pws <- seq(4,7,0.1)
betas <- do.call(cbind, lapply(pws, function(pw) {
  a <- read.table(sprintf('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/LaminB1/multiscale_ov/scale_%1.1f.bed', pw), sep='\t')
  setNames(as.numeric(a$V6), paste0(a$V1,':',a$V4))
}))
colnames(betas) <- pws
save(betas, file='/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/LaminB1/multiscale_ov/multiscale.rda')

chr16short <- rownames(betas)[grep('chr16',rownames(betas))]
chr16short <- chr16short[sapply(strsplit(chr16short,"[:-]"), function(x) as.integer(x[2])) < 35000000]
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/LaminB1/multiscale_ov/multiscale.rda')
## qts <- quantile(betas, c(0.01,0.99), na.rm=T)
## betas[is.na(betas)] <- -2.5
betas[betas < -2.5] <- -2.5
betas[betas > 2.5] <- 2.5
png('~/gallery/2017_02_18_IMR90_laminB1_multiscale.png', width=2000, height=400)
WHeatmap(t(betas[chr16short,]), cmp=CMPar(stop.points=c('#fff7bc','#d95f0e'))) + WCustomize(mar.left=0, mar.top=0, mar.right=0, mar.bottom=0)
dev.off()
pdf('~/gallery/2017_02_18_IMR90_laminB1_multiscale_legend.pdf', width=2, height=2)
WColorBarH(seq(-2.5,2.5,0.2), cmp=CMPar(stop.points=c('#fff7bc','#d95f0e')))
dev.off()
pdf('~/gallery/2017_02_18_IMR90_betas_legend.pdf', width=2, height=2)
WColorBarH(seq(0,1,0.01), cmp=CMPar(stop.points=c('blue','yellow')))
dev.off()

## LOAD DATA
marks <- c('H3K36me3', 'H3K27ac','H3K4me1','H3K27me3','repliseq','H3K36me3','H3K4me3','H3K9me3','genebody', 'DNase', 'RNAseq', 'exon', 'intron')
## marks <- c('intergenic')
for (mark in marks) {
  cat(mark,'\n')
  pws <- seq(4,7,0.1)
  betas <- do.call(cbind, mclapply(pws, function(pw) {
    a <- read.table(sprintf('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/%s/multiscale_ov/scale_%1.1f.bed', mark, pw), sep='\t')
    setNames(as.numeric(a$V6), paste0(a$V1,':',a$V4))
  }, mc.cores=3))
  colnames(betas) <- pws
  save(betas, file=sprintf('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/%s/multiscale_ov/multiscale.rda', mark))
}

## histone marks
marks <- c('H3K36me3', 'H3K27ac','H3K4me1','H3K27me3','H3K36me3','H3K4me3','H3K9me3')
for (mark in marks) {
  load(sprintf('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/%s/multiscale_ov/multiscale.rda', mark))
  a <- do.call(cbind, lapply(colnames(betas), function(x) betas[,x]))
  colnames(a) <- colnames(betas)
  betas <- a
  chr16short <- rownames(betas)[grep('chr16',rownames(betas))]
  chr16short <- chr16short[sapply(strsplit(chr16short,"[:-]"), function(x) as.integer(x[2])) < 35000000]
  png(sprintf('~/gallery/2017_04_04_IMR90_%s_multiscale.png', mark), width=2000, height=400)
  print(WHeatmap(t(betas[chr16short,]), cmp=CMPar(stop.points=c('#FFFFFF','#7D8A2E'))) + WCustomize(mar.left=0, mar.top=0, mar.right=0, mar.bottom=0))
  dev.off()
}

## repliseq
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/repliseq/multiscale_ov/multiscale.rda')
a <- do.call(cbind, lapply(colnames(betas), function(x) betas[,x]))
colnames(a) <- colnames(betas)
betas <- a
chr16short <- rownames(betas)[grep('chr16',rownames(betas))]
chr16short <- chr16short[sapply(strsplit(chr16short,"[:-]"), function(x) as.integer(x[2])) < 35000000]
png(sprintf('~/gallery/2017_04_04_IMR90_repliseq_multiscale.png', mark), width=2000, height=400)
print(WHeatmap(t(betas[chr16short,]), cmp=CMPar(stop.points=c('#FFFFFF','#332532'))) + WCustomize(mar.left=0, mar.top=0, mar.right=0, mar.bottom=0))
dev.off()

pdf('~/gallery/2017_04_04_IMR90_histone_multiscale_legend.pdf', width=2, height=2)
WColorBarH(seq(0,1,0.1), cmp=CMPar(stop.points=c('#FFFFFF','#7D8A2E')))
dev.off()

pdf('~/gallery/2017_04_04_IMR90_repliseq_multiscale_legend.pdf', width=2, height=2)
WColorBarH(seq(0,90,5), cmp=CMPar(stop.points=c('#FFFFFF','#332532')))
dev.off()

pdf('~/gallery/2017_04_04_IMR90_rnaseq_multiscale_legend.pdf', width=2, height=2)
WColorBarH(seq(0,90,5), cmp=CMPar(stop.points=c('#FFFFFF','#CC0000')))
dev.off()

## rnaseq, log2(rpkm)
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/RNAseq/multiscale_ov/multiscale.rda')
a <- do.call(cbind, lapply(colnames(betas), function(x) betas[,x]))
colnames(a) <- colnames(betas)
betas <- a
chr16short <- rownames(betas)[grep('chr16',rownames(betas))]
chr16short <- chr16short[sapply(strsplit(chr16short,"[:-]"), function(x) as.integer(x[2])) < 35000000]
png(sprintf('~/gallery/2017_04_04_IMR90_rnaseq_multiscale.png', mark), width=2000, height=400)
print(WHeatmap(t(log2(1+betas[chr16short,])), cmp=CMPar(stop.points=c('#FFFFFF','#CC0000','red'))) + WCustomize(mar.left=0, mar.top=0, mar.right=0, mar.bottom=0))
dev.off()

## dnaseq
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/DNase/multiscale_ov/multiscale.rda')
a <- do.call(cbind, lapply(colnames(betas), function(x) betas[,x]))
colnames(a) <- colnames(betas)
betas <- a
chr16short <- rownames(betas)[grep('chr16',rownames(betas))]
chr16short <- chr16short[sapply(strsplit(chr16short,"[:-]"), function(x) as.integer(x[2])) < 35000000]
png(sprintf('~/gallery/2017_04_04_IMR90_dnaseq_multiscale.png', mark), width=2000, height=400)
print(WHeatmap(t(betas[chr16short,]), cmp=CMPar(stop.points=c('#FFFFFF','black','black'))) + WCustomize(mar.left=0, mar.top=0, mar.right=0, mar.bottom=0))
dev.off()

## genebody
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/genebody/multiscale_ov/multiscale.rda')
a <- do.call(cbind, lapply(colnames(betas), function(x) betas[,x]))
colnames(a) <- colnames(betas)
betas <- a
chr16short <- rownames(betas)[grep('chr16',rownames(betas))]
chr16short <- chr16short[sapply(strsplit(chr16short,"[:-]"), function(x) as.integer(x[2])) < 35000000]
png(sprintf('~/gallery/2017_04_04_IMR90_genebody_multiscale.png', mark), width=2000, height=400)
print(WHeatmap(t(betas[chr16short,]), cmp=CMPar(stop.points=c('#FFFFFF','#7D8A2E'))) + WCustomize(mar.left=0, mar.top=0, mar.right=0, mar.bottom=0))
dev.off()

## exon
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/exon/multiscale_ov/multiscale.rda')
a <- do.call(cbind, lapply(colnames(betas), function(x) betas[,x]))
colnames(a) <- colnames(betas)
betas <- a
chr16short <- rownames(betas)[grep('chr16',rownames(betas))]
chr16short <- chr16short[sapply(strsplit(chr16short,"[:-]"), function(x) as.integer(x[2])) < 35000000]
png(sprintf('~/gallery/2017_04_04_IMR90_exon_multiscale.png', mark), width=2000, height=400)
print(WHeatmap(t(betas[chr16short,]), cmp=CMPar(stop.points=c('#FFFFFF','black','black'))) + WCustomize(mar.left=0, mar.top=0, mar.right=0, mar.bottom=0))
dev.off()

## intron
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/intron/multiscale_ov/multiscale.rda')
a <- do.call(cbind, lapply(colnames(betas), function(x) betas[,x]))
colnames(a) <- colnames(betas)
betas <- a
chr16short <- rownames(betas)[grep('chr16',rownames(betas))]
chr16short <- chr16short[sapply(strsplit(chr16short,"[:-]"), function(x) as.integer(x[2])) < 35000000]
png(sprintf('~/gallery/2017_04_04_IMR90_intron_multiscale.png', mark), width=2000, height=400)
print(WHeatmap(t(betas[chr16short,]), cmp=CMPar(stop.points=c('#FFFFFF','black'))) + WCustomize(mar.left=0, mar.top=0, mar.right=0, mar.bottom=0))
dev.off()

## intergenic
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/intergenic/multiscale_ov/multiscale.rda')
a <- do.call(cbind, lapply(colnames(betas), function(x) betas[,x]))
colnames(a) <- colnames(betas)
betas <- a
chr16short <- rownames(betas)[grep('chr16',rownames(betas))]
chr16short <- chr16short[sapply(strsplit(chr16short,"[:-]"), function(x) as.integer(x[2])) < 35000000]
png(sprintf('~/gallery/2017_04_04_IMR90_intergenic_multiscale.png', mark), width=2000, height=400)
print(WHeatmap(t(betas[chr16short,]), cmp=CMPar(stop.points=c('#FFFFFF','black'))) + WCustomize(mar.left=0, mar.top=0, mar.right=0, mar.bottom=0))
dev.off()

##################
## domains
##################

setwd('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/')

chrms <- c(paste0('chr',1:22),'chrX','chrY','chrM')

a <- read.table('domain/window100kb_HiCBcompartment.bed', stringsAsFactors=F)
a <- a[a$V1 %in% chrms,]
hicb <- setNames(a$V4, paste0(a$V1,":",a$V2,"-",a$V3))

a <- read.table('domain/window100kb_laminB.bed', stringsAsFactors=F)
a <- a[a$V1 %in% chrms,]
laminb <- setNames(a$V4, paste0(a$V1,":",a$V2,"-",a$V3))

a <- read.table('domain/window100kb_repliseq_imr90.bed', stringsAsFactors=F)
a <- a[a$V1 %in% chrms,]
repliseq <- setNames(a$V4, paste0(a$V1,":",a$V2,"-",a$V3))

a <- read.table('domain/window100kb_h3k36.bed', stringsAsFactors=F)
a <- a[a$V1 %in% chrms,]
h3k36 <- setNames(a$V4, paste0(a$V1,":",a$V2,"-",a$V3))

a <- read.table('domain/window100kb_h3k9.bed', stringsAsFactors=F)
a <- a[a$V1 %in% chrms,]
h3k9 <- setNames(a$V4, paste0(a$V1,":",a$V2,"-",a$V3))

save(hicb, laminb, repliseq, h3k9, h3k36, file='domain/window100kb.rda')

load('allprone/bin_mean_100kb.rda')
chr16short <- rownames(betas)[grep('chr16',rownames(betas))]
chr16short <- chr16short[sapply(strsplit(chr16short,"[:-]"), function(x) as.integer(x[2])) < 35000000]

WColorBarH(hicb[chr16short], cmp=CMPar(label2color=c('A'='blue','B'='yellow'))) + WColorBarH(laminb[chr16short], Beneath(), cmp=CMPar(stop.points=c('blue','yellow'))) + WColorBarH(repliseq, Beneath(), cmp=CMPar(stop.points=c('blue','yellow'))) + WColorBarH(log2(1+h3k9),Beneath(), cmp=CMPar(stop.points=c('blue','yellow'))) + WColorBarH(log2(1+h3k36), Beneath(), cmp=CMPar(stop.points=c('blue','yellow')))

###########################
## IMR90 K27 vs prone cpg
###########################

setwd('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90')

pws <- sprintf('%1.1f', seq(4,7,0.1))
betas <- lapply(pws, function(pw) {
  a <- read.table(paste0('pronecpg/scale_', pw, '.bed'), stringsAsFactor=F)
  setNames(as.numeric(a$V4), paste0(a$V1, ':', a$V2, '-', a$V3))
})

dbetas <- lapply(betas, function(b) {
  dd <- abs(b[-c(1,2)] - b[-c(length(b),length(b)-1)])
  c(0,dd,0)
})

K27me3 <- lapply(pws, function(pw) {
  a <- read.table(paste0('H3K27me3/scale_', pw, '_merged.bed'), stringsAsFactor=F)
  setNames(as.numeric(a$V4), paste0(a$V1, ':', a$V2, '-', a$V3))
})

pws <- sprintf('%1.1f', seq(4,7,0.1))
names(betas) <- pws
names(dbetas) <- pws
names(K27me3) <- pws

save(betas, dbetas, K27me3, file='dname_vs_k27me3.rda')
save(betas, dbetas, K27me3, file='/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/dname_vs_k27me3.rda')

## plot
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/dname_vs_k27me3.rda')

a <- dbetas[[4.1]]
b <- log2(1+K27me3[[4.1]])
smoothScatter(a,b)

pdf('~/gallery/2017_01_10_dname_vs_k27me3_scatter.pdf')
smoothScatter(betas[[4.1]], log2(1+K27me3[[4.1]]), xlab='Betas', ylab='K27me3')
dev.off()
pdf('~/gallery/2017_01_10_dname_vs_k27me3_scatter_delta.pdf')
smoothScatter(dbetas[[4.1]], log2(1+K27me3[[4.1]]), xlab='Delta Betas', ylab='K27me3', xlim=c(0,0.4))
dev.off()

for (scale in pws) {
  betas1 <- betas[[scale]]
  dbetas1 <- dbetas[[scale]]
  k27me3 <- K27me3[[scale]]
  chr16 <- names(betas1)[grep('chr16', names(betas1))]
  pdf(sprintf('~/gallery/2017_01_10_dname_vs_k27me3_scale_%s.pdf', scale), width=18, height=4)
  par(mar=c(5,5,1,1))
  plot(betas1[chr16], type='l', col='blue', ylim=c(0,3.5), xlab='Genomic Location', ylab='Signal')
  lines(dbetas1[chr16]+1, type='l', col='red')
  lines(k27me3[chr16]/max(k27me3[chr16])+2, type='l', col='green')
  legend('top', c('K27me3', 'DNAme', 'delta_DNAme'), lty=1, col=c('green','blue','red'))
  dev.off()
}


#######################################
## all pair-wise sample correlation
#######################################

setwd('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/')
load('allprone/bin_mean_100kb.rda')
allcors <- simplify2array(lapply(1:ncol(betas), function(i) {
  cat(i,'\n');
  simplify2array(mclapply(1:ncol(betas), function(j) {
    cor(betas[,i], betas[,j], use='complete')
  }))
}))
colnames(allcors) <- colnames(betas)
rownames(allcors) <- colnames(betas)
save(allcors, file='tissuespecific/allcors.rda')

## cormat <- allcors[
setwd('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/')
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/tissuespecific/allcors.rda')
samples <- read.csv2('samplesheet.txt', header=T, sep='\t', row.names='samplename', stringsAsFactors=F)
## tumorsamples <- rownames(samples[!(samples$tumor %in% c('N','CL')) & !(samples$tissue %in% c('Blood', 'Other')) & samples$isFetal == 'Adult' & samples$source != 'BLUEPRINT',])
## normalsamples <- rownames(samples[samples$tumor %in% c('N') & !(samples$tissue %in% c('Blood', 'Other')) & samples$isFetal == 'Adult' & samples$source != 'BLUEPRINT',])
isamples <- rownames(samples[!(samples$tumor %in% c('CL')) & !(samples$tissue %in% c('Blood', 'Other')) & samples$isFetal == 'Adult' & samples$source != 'BLUEPRINT',])
cormat <- allcors[isamples, isamples]
cormat <- both.cluster(cormat)$mat
WHeatmap(cormat, name='a', xticklabels=T, yticklabels=T) + WColorBarV(rownames(cormat), RightOf('a'), cmp=CMPar(label2color=setNames(samples[rownames(cormat), 'color3'], rownames(cormat)))) + WColorBarH(colnames(cormat), TopOf('a'), cmp=CMPar(label2color=setNames(samples[colnames(cormat),'color3'], colnames(cormat))))

## all pairwise correlation of scaled betas
setwd('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/')
load('allprone/bin_mean_100kb.rda')
alltcga.samples <- colnames(betas)[grep('TCGA_', colnames(betas))]
tumorsamples <- alltcga.samples[!grepl('_N', alltcga.samples)]
sd.binmeanT <- na.omit(apply(betas[,tumorsamples], 1, sd, na.rm=T))
comPMD <- names(na.omit(sd.binmeanT[sd.binmeanT >= 0.18]))
quintiles <- apply(betas[comPMD,], 2, function(x) quantile(x, probs=c(0.2,0.8), na.rm=T))
quintiles[1,] <- pmin(quintiles[2,]-0.12, quintiles[1,])
betas.scaled <- do.call(cbind, lapply(colnames(betas),function(sn) {
  b <- betas[,sn]; qt <- quintiles[,sn];
  b[!is.na(b) & b<qt[1]] <- qt[1];
  b[!is.na(b) & b>qt[2]] <- qt[2];
  bb <- !is.na(b) & b>=qt[1] & b<=qt[2]
  b[bb] <- (b[bb]-qt[1]) / (qt[2]-qt[1]);
  b
}))
colnames(betas.scaled) <- colnames(betas)
allcors <- simplify2array(lapply(1:ncol(betas.scaled), function(i) {
  cat(i,'\n');
  simplify2array(mclapply(1:ncol(betas.scaled), function(j) {
    cor(betas.scaled[,i], betas.scaled[,j], use='complete')
  }))
}))
colnames(allcors) <- colnames(betas.scaled)
rownames(allcors) <- colnames(betas.scaled)
save(allcors, file='tissuespecific/allcors.scaled.rda')

setwd('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/')
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/tissuespecific/allcors.rda')
samples <- read.csv2('samplesheet.txt', header=T, sep='\t', row.names='samplename', stringsAsFactors=F)
## tumorsamples <- rownames(samples[!(samples$tumor %in% c('N','CL')) & !(samples$tissue %in% c('Blood', 'Other')) & samples$isFetal == 'Adult' & samples$source != 'BLUEPRINT',])
## normalsamples <- rownames(samples[samples$tumor %in% c('N') & !(samples$tissue %in% c('Blood', 'Other')) & samples$isFetal == 'Adult' & samples$source != 'BLUEPRINT',])
samples <- read.csv2('samplesheet.txt', header=T, sep='\t', row.names='samplename', stringsAsFactors=F)
isamples <- rownames(samples[!(samples$tumor %in% c('CL')) & !(samples$tissue %in% c('Blood', 'Other')) & samples$isFetal == 'Adult' & samples$source != 'BLUEPRINT' & !(samples$source %in% 'methbase' & samples$tissue=='Brain'),])
cormat <- allcors[isamples, isamples]
## cormat <- both.cluster(cormat)$mat
WHeatmap(cormat, name='a', xticklabels=T, yticklabels=T, cmp=CMPar(stop.points=c('DarkBlue','yellow'), dmin=0.1)) + WColorBarV(rownames(cormat), RightOf('a'), cmp=CMPar(label2color=setNames(samples[rownames(cormat), 'color3'], rownames(cormat)))) + WColorBarH(colnames(cormat), TopOf('a'), cmp=CMPar(label2color=setNames(samples[colnames(cormat),'color3'], colnames(cormat))))

###########################################################
## all pairwise correlation of samples, at different scales
###########################################################

## all raw correlation
pws <- seq(3,5,by=0.1)
for (pw in pws) {
  cat('Now processing', pw,'\n')
  load(sprintf('allprone/multiscale/bin_mean_scale_%1.1f.rda', pw))
  allcors <- simplify2array(lapply(1:ncol(betas), function(i) {
    cat(i,' ');
    simplify2array(mclapply(1:ncol(betas), function(j) {
      ## cat(j,'\n')
      cor(betas[,i], betas[,j], use='complete', method='spearman')
    }))
  }))
  colnames(allcors) <- colnames(betas)
  rownames(allcors) <- colnames(betas)
  save(allcors, file=sprintf('tissuespecific/multiscale/allcors.scale.%1.1f.rda', pw))
  cat('\n')
}

## tumor specific variance
setwd('~/projects/laird-secondary/2016_12_26_TCGA_WGBS')
pws <- seq(3,7,by=0.1)
mclapply(pws, function(pw) {
  cat('Now processing', pw,'\n')
  load(sprintf('allprone/multiscale/bin_mean_scale_%1.1f.rda', pw))
  alltcga.samples <- colnames(betas)[grep('TCGA_', colnames(betas))]
  tumorsamples <- alltcga.samples[!grepl('_N', alltcga.samples)]
  ## normalsamples <- alltcga.samples[grepl('_N', alltcga.samples)]
  sd.binmeanT <- na.omit(apply(betas[,tumorsamples], 1, sd, na.rm=T))
  save(sd.binmeanT, file=sprintf('tissuespecific/multiscale/alltumorsd.scale.%1.1f.rda', pw))
  cat('\n')
}, mc.cores=10)

#####################################################
## all pairwise correlation of bins, at 100kb scale
#####################################################
setwd('~/projects/laird-secondary/2016_12_26_TCGA_WGBS')
load('allprone/multiscale/bin_mean_scale_5.0.rda')

smoothScatter(betas[,'TCGA_UCEC_A1CI'], betas[,'TCGA_READ_2689'])
smoothScatter(betas[,'TCGA_READ_3593'], betas[,'TCGA_READ_2689'])

panel.cor <- function(x, y, digits=2, prefix="", cex.cor) {
  usr <- par("usr"); on.exit(par(usr)) 
  par(usr = c(0, 1, 0, 1)) 
  rraw <- cor(x, y)
  r <- abs(rraw)
  txt <- format(c(rraw, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="") 
  if(missing(cex.cor)) cex <- 0.8/strwidth(txt) 
  test <- cor(x,y, use='complete', method='spearman')
  text(0.5, 0.5, txt, cex = max(0.8, cex * r))
}

panel.dens = function (x, ...) {
  par(new = TRUE)
  plot(density(na.omit(x)), col = "red", lwd = 2, main="", axes=FALSE, font.main=2)
}

df <- betas[, c(
  'TCGA_LUAD_4630','TCGA_LUAD_6148','TCGA_LUAD_6215','TCGA_LUAD_N6148',
  'TCGA_READ_N2689','TCGA_COAD_3158','TCGA_READ_2689','TCGA_READ_3593')]
## df <- apply(df, 2, function(x) rank(x))
png('~/gallery/2017_02_11_COAD_LUAD_clustering.png', width=500, height=500)
pairs(df, panel = function(...)
  smoothScatter(..., colramp=colorRampPalette(c('white','black')), nrpoints=0, add = TRUE), upper.panel=panel.cor, diag.panel=panel.dens)
dev.off()

between.normal <- abs(rank(betas[,'TCGA_LUAD_N6148']) - rank(betas[,'TCGA_READ_N2689']))
between.tumor <- abs(rank(betas[,'TCGA_READ_2689']) - rank(betas[,'TCGA_READ_N2689']))

head(sort(between.normal - between.tumor), n=20) # most telomere/centromere, tumor is too low!
tail(sort(between.normal - between.tumor), n=20) # most in middle of each arm, REAL!

betas <- t(betas[complete.cases(betas),])
pc <- prcomp(betas)


load('allprone/bin_mean_100kb.rda')
chr16 <- rownames(betas)[grep('chr16',rownames(betas))]

load('allprone/multiscale/bin_mean_scale_5.0.rda')
par(mfrow=c(8,1), mar=c(1,5,3,1))
plot(betas[chr16,'TCGA_READ_2689'], col='blue', type='l', main='READ_2689 vs READ_N2689 raw', ylab='')
lines(betas[chr16,'TCGA_READ_N2689'], col='red')
plot(rank(betas[chr16,'TCGA_READ_2689']), col='blue', type='l', main='READ_2689 vs READ_N2689', ylab='')
lines(rank(betas[chr16,'TCGA_READ_N2689']), col='red')
plot(rank(betas[chr16,'TCGA_READ_2689'])-rank(betas[chr16,'TCGA_READ_N2689']), type='l', main='READ_2689 - READ_N2689', ylab='')
abline(h=0, lty='dashed')
plot(rank(betas[chr16,'TCGA_LUAD_N6148']), col='darkgreen', type='l', main='LUAD_N6148 VS READ_N2689', ylab='')
lines(rank(betas[chr16,'TCGA_READ_N2689']), col='red')
plot(rank(betas[chr16,'TCGA_LUAD_N6148'])-rank(betas[chr16,'TCGA_READ_N2689']), type='l', main='LUAD_N6148 - READ_N2689', ylab='')
abline(h=0, lty='dashed')
plot(betas[chr16,'TCGA_LUAD_6148'], col='blue', type='l', main='LUAD_6148 vs LUAD_N6148 raw', ylab='')
lines(betas[chr16,'TCGA_LUAD_N6148'], col='red')
plot(rank(betas[chr16,'TCGA_LUAD_6148']), col='purple', type='l', main='LUAD_6148 VS LUAD_N6148', ylab='')
lines(rank(betas[chr16,'TCGA_LUAD_N6148']), col='darkgreen')
plot(rank(betas[chr16,'TCGA_LUAD_6148'])-rank(betas[chr16,'TCGA_LUAD_N6148']), type='l', main='LUAD_6148 - READ_N6148', ylab='')
abline(h=0, lty='dashed')

par(mfrow=c(7,1), mar=c(1,5,3,1))
plot(betas[chr16,'TCGA_READ_2689'], col='blue', type='l', main='READ_2689 vs READ_N2689 raw', ylab='')
lines(betas[chr16,'TCGA_READ_N2689'], col='red')
plot(rank(betas[chr16,'TCGA_READ_2689']), col='blue', type='l', main='READ_2689 vs READ_N2689', ylab='')
lines(rank(betas[chr16,'TCGA_READ_N2689']), col='red')
plot(rank(betas[chr16,'TCGA_READ_2689'])-rank(betas[chr16,'TCGA_READ_N2689']), type='l', main='READ_2689 - READ_N2689', ylab='')
abline(h=0, lty='dashed')
plot(betas[chr16,'TCGA_COAD_3158'], col='blue', type='l', main='COAD_3158 vs COAD_N3158 raw', ylab='')
lines(betas[chr16,'TCGA_COAD_N3158'], col='red')
plot(rank(betas[chr16,'TCGA_COAD_3158']), col='blue', type='l', main='COAD_3158 vs COAD_N3158', ylab='')
lines(rank(betas[chr16,'TCGA_COAD_N3158']), col='red')
plot(rank(betas[chr16,'TCGA_COAD_3158'])-rank(betas[chr16,'TCGA_COAD_N3158']), type='l', main='COAD_3158 - COAD_N3158', ylab='')
abline(h=0, lty='dashed')
plot(betas[chr16,'TCGA_READ_2689'], col='blue', type='l', main='READ_2689 vs READ_3593 raw', ylab='')
lines(betas[chr16,'TCGA_READ_3593'], col='darkgreen')

###############################
## PCA analysis
###############################
load('allprone/bin_mean_100kb.rda')
betas <- t(betas[complete.cases(betas),])
pc <- prcomp(betas)

depths <- read.table('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/comPMDdepth.tsv', sep='\t', header=T)
for (i in 1:5) {
  pdf(sprintf('~/gallery/2017_02_06_PCA_commonPMD_PC%d.pdf', i), width=4, height=3.7)
  par(mar=c(5,5,1,1))
  plot(pc$x[rownames(depths),sprintf('PC%d',i)], depths$x, col=samples[rownames(depths),]$color, pch=16, xlab=sprintf('PC%d',i), ylab='common PMD methylation')
  dev.off()
}

## remove PC1
pcx <- pc$x
pcx[,'PC1'] <- 0
a <- pcx %*% t(pc$rotation)
betas.noPC1 <- apply(a, 1, function(x) x+pc$center)
save(betas.noPC1, file='/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/tissuespecific/betas.noPC1.rda')

## compute correlation without PC1
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/tissuespecific/betas.noPC1.rda')
allcors <- simplify2array(lapply(1:ncol(betas.noPC1), function(i) {
  cat(i,'\n');
  simplify2array(mclapply(1:ncol(betas.noPC1), function(j) {
    cor(betas.noPC1[,i], betas.noPC1[,j], use='complete', method='spearman')
  }, mc.cores=8))
}))
colnames(allcors) <- colnames(betas.noPC1)
rownames(allcors) <- colnames(betas.noPC1)
save(allcors, file='/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/tissuespecific/allcors.noPC1.rda')

## plot correlation
setwd('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/')
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/tissuespecific/allcors.noPC1.rda')
samples <- read.csv2('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/samplesheet.txt', header=T, sep='\t', row.names='samplename', stringsAsFactors=F)
isamples <- rownames(samples[!(samples$tumor %in% c('CL')) & !(samples$tissue %in% c('Blood', 'Other')) & samples$isFetal == 'Adult' & samples$source != 'BLUEPRINT' & !(samples$source %in% 'methbase' & samples$tissue=='Brain'),])
cormat <- allcors[isamples, isamples]
cormat <- both.cluster(cormat)$mat
WHeatmap(cormat, name='a', xticklabels=T, yticklabels=T, cmp=CMPar(stop.points=c('DarkBlue','yellow'), dmin=0.5)) + WColorBarV(rownames(cormat), RightOf('a'), cmp=CMPar(label2color=setNames(samples[rownames(cormat), 'color3'], rownames(cormat)))) + WColorBarH(colnames(cormat), TopOf('a'), cmp=CMPar(label2color=setNames(samples[colnames(cormat),'color3'], colnames(cormat))))

isamples <- rownames(samples[samples$source == 'BLUEPRINT',])
cormat <- allcors[isamples, isamples]
cormat <- both.cluster(cormat)$mat
pdf('~/gallery/2017_02_07_allcorrelation_heatmap_noPC1.pdf', width=10, height=10)
WHeatmap(cormat, name='a', xticklabels=T, yticklabels=T, xticklabels.n=50, yticklabels.n=50, cmp=CMPar(stop.points=c('DarkBlue','yellow'), dmin=0.5)) + WColorBarV(rownames(cormat), RightOf('a'), cmp=CMPar(label2color=setNames(samples[rownames(cormat), 'color3'], rownames(cormat)))) + WColorBarH(colnames(cormat), TopOf('a'), cmp=CMPar(label2color=setNames(samples[colnames(cormat),'color3'], colnames(cormat))))
dev.off()

################################
## BLUEPRINT blood vs age
################################

setwd('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/')
samples <- read.table('samplesheetBP.txt', sep='\t', header=T, quote="", comment.char='')
samples <- samples[samples$BP_DONOR_HEALTH_STATUS == 'Healthy',]
samples1 <- samples[samples$BP_CELL_TYPE %in% c('alternatively activated macrophage', 'mature neutrophil', 'inflammatory macrophage'), ]
pdf('~/gallery/2017_01_17_BLUEPRINT_age_vs_pmd.pdf', width=8, height=5)
ggplot(samples1, aes(AGE_MIDDLE, PMDdepth, color=BP_CELL_TYPE)) + geom_point() + geom_smooth(method='lm', se=FALSE) + xlab('Age') + ylab('PMD methylation') + scale_color_discrete(guide=guide_legend(title='Cell Type'))
dev.off()

################################
## SATB2
################################
setwd('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS')
load('allprone/bin_mean_100kb.rda')
bincoords <- as.data.frame(do.call(rbind, strsplit(rownames(betas),"[:-]")))
colnames(bincoords) <- c('chrm','beg','end')
bincoords$ID <- rownames(betas)
bincoords$beg <- as.integer(bincoords$beg)
bincoords$end <- as.integer(bincoords$end)
satb2 <- subset(bincoords, chrm=="chr2" & beg > 197420463 & end < 201920243)$ID

samples <- read.csv2('samplesheet.txt', header=T, sep='\t', row.names='samplename', stringsAsFactors=F)
tsamples <- subset(samples, plotorder >= 500 & PMDdepth < 0.5)

tsamplesn <- c(
  'TCGA_COAD_N3158', 'TCGA_READ_N2689', 'REMC_E106_Sigmoid_Colon', 'methbase_Berman_2012_Human_ColonNormal', 'Ziller_157_REMC_19_colonic_mucosa',
  'TCGA_UCEC_NA1CI', 'TCGA_BRCA_NA0CE', 'TCGA_BLCA_NA20V', 'TCGA_LUAD_N6148', 'TCGA_LUSC_N2722',
  'TCGA_LUAD_6215',
  "TCGA_LUAD_7156",
  "TCGA_LUAD_4630",
  "TCGA_LUAD_6840",
  "TCGA_LUSC_1078",
  "TCGA_LUSC_2600",
  "TCGA_LUSC_2722",
  "TCGA_BRCA_A0YG",
  "TCGA_BRCA_A04X",
  "TCGA_BRCA_A07I",
  "TCGA_BRCA_A15H",
  'BP_venous_blood_S00B1LU1_NA_ChronicLymphocyticLeukemia',
  'BP_venous_blood_S00GPRA1_NA_ChronicLymphocyticLeukemia',
  'BP_venous_blood_S01FG4A1_NA_ChronicLymphocyticLeukemia',
  'BP_venous_blood_S00B2JA1_NA_ChronicLymphocyticLeukemia',
  'BP_venous_blood_S00AYXU1_NA_ChronicLymphocyticLeukemia',
  'TCGA_BLCA_A2LA',
  'TCGA_BLCA_A20V',
  'TCGA_BLCA_A13J',
  'TCGA_BLCA_A2HQ',
  'TCGA_BLCA_A1AG',
  'TCGA_BLCA_A1AA',
  "TCGA_COAD_A00R",
  "TCGA_COAD_3158",
  "Ziller_GSM1204465_1120_Colon_Primary_Tumor",
  "methbase_Berman_2012_Human_ColonCancer",
  "TCGA_READ_3593",
  "TCGA_READ_2689"
  )
pdf('~/gallery/2017_02_10_SATB2_100kb.pdf', width=6, height=3.5)
WHeatmap(t(betas[satb2,tsamplesn]), yticklabels=T, yticklabels.n=nrow(tsamples), cmp=CMPar(stop.points=c('blue','yellow'), dmin=0, dmax=1))
dev.off()


#######################
## Fetal and Embryo
#######################
setwd('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS')

### not much use below
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/archived_Huy/2015_10_12_Huy_segments/commonPMDanno.rda')
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/cgcontext.hg19/ctxt.rda')
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/annotations/hg19.cpg.cgi.rda')
domain.betas <- mclapply(samples, function(sample) {
  gc()
  load(paste0('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/rda/', sample, '.rda'))
  gr <- mergeOv(mergeOv(mergeOv(gr, ctxt), common.pmd), cgi)
  gr.pmd <- gr[gr$cgi == 'NonCGI' & gr$commonPMD=='comPMD' & gr$ctxt=='WCGW' & gr$orph35class=='[0,1)']
  gr.hmd <- gr[gr$cgi == 'NonCGI' & gr$commonPMD=='comHMD' & gr$ctxt=='WCGW' & gr$orph35class=='[0,1)']
  list(pmd=gr.pmd$betas, hmd=gr.hmd$betas)
}, mc.cores=28)
names(domain.betas) <- samples
domain.betas <- do.call(rbind, lapply(names(domain.betas), function(x) rbind(data.frame(betas=domain.betas[[x]]$pmd, domain='pmd',sample=x), data.frame(betas=domain.betas[[x]]$hmd, domain='hmd',sample=x))))
save(domain.betas, file='/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/fetalbetas.rda')

load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/fetalbetas.rda')

### early development plots

samples <- c('GSE49828_ICM_lowcov', 'methbase_Lister_2009_Human_H1ESC', 'methbase_Lister_2011_Human_H9ESC',
             'GSE63818_EmbryonicHeart_lowcov', 'GSE63818_EmbryonicBrain_lowcov', 'GSE49828_EmbryonicLiver_lowcov',
             'REMC_BI.Fetal_Muscle_Leg.UW_H24996', 'REMC_BI.Fetal_Thymus.UW_H24943', 'REMC_E084_Fetal_Intestine_Large', 'Ziller_515_fetal_brain', 'Ziller_1119_fetal_heart', 'Ziller_1238_fThymus_304082', 'Ziller_1243_fMuscle_Leg_304083', 'Ziller_1244_fAdrenal_304084',
             'REMC_E095_Left_Ventricle', 'REMC_E104_Right_Atrium', 'REMC_E105_Right_Ventricle', 'REMC_E100_Psoas_Muscle',
             'Ziller_liver_162_REMC3_null', 'REMC_E066_Adult_Liver', 'Ziller_GSM916049_liver_161_REMC2',
             'REMC_E112_Thymus',
             'methbase_Brain_Human_FrontCortexFemale64Yr',
             'methbase_Brain_Human_MidFrontGyr16Yr',
             'methbase_Brain_Human_MidFrontGyr12Yr',
             'methbase_Brain_Human_MidFrontGyr5Yr',
             'methbase_Brain_Human_MidFrontGyr35Day',
             'methbase_Brain_Human_MidFrontGyr2Yr',
             'REMC_E071_Brain_Hippocampus_Middle',
             'REMC_E070_Brain_Germinal_Matrix',
             'TCGA_GBM_1401', 'TCGA_BRCA_A07I',
             'methbase_Berman_2012_Human_ColonNormal',
             'methbase_Berman_2012_Human_ColonCancer',
             'REMC_UCSD.Adrenal_Gland.STL003'
             )
sample <- 'methbase_Lister_2009_Human_H1ESC'
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/windows/binmean.rda')
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/allprone/multiscale/bin_mean_scale_5.0.rda')
rownames(betas) <- sapply(strsplit(rownames(betas),'\\-'), function(x) x[1])
chr16 <- rownames(betas)[grep('chr16',rownames(betas))]
chr16 <- chr16[chr16 != 'chr16:0']

pdf('~/gallery/2017_02_19_early_dev_scatter.pdf', width=11, height=19)
par(mfrow=c(9,4), mar=c(1,5,5,1), oma=c(4,4,2,2), lwd=2.5, font.lab=1, font.axis=2, cex.lab=0.5, cex.axis=2, cex.main=2, font.main=2)
for (sample in samples) {
  a <- (!is.na(betas[chr16, sample]) & !is.na(mean.binmeanT[chr16, '5.0']))
  rho <- cor(mean.binmeanT[chr16[a], '5.0'], betas[chr16[a], sample], method='spearman')
  plot(rank(mean.binmeanT[chr16[a], '5.0']), rank(betas[chr16[a], sample]), ylab=sample, main=sprintf('rho: %1.3f', rho), pch=16)
}
dev.off()
for (sample in samples) {
  a <- (!is.na(betas[chr16, sample]) & !is.na(mean.binmeanT[chr16, '5.0']))
  rho <- cor(mean.binmeanT[chr16[a], '5.0'], betas[chr16[a], sample], method='spearman')
  png(sprintf('~/gallery/2017_02_19_earlydev/%s.png', sample), width=2000, height=2000, res=500)
  par(mar=c(0,0,0,0))
  plot(rank(mean.binmeanT[chr16[a], '5.0']), rank(betas[chr16[a], sample]), pch=16, cex=1.6)
  dev.off()
}

## violin plot
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/archived_Huy/2015_10_12_Huy_segments/commonPMDanno.rda')
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/cgcontext.hg19/ctxt.rda')
cpgs <- ctxt
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/cgcontext.hg19/ctxt_w_sdbinmeanT.rda')
cpgs <- mergeOv(cpgs, ctxt)
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/annotations/hg19.cpg.cgi.rda')
allsamples <- list.files('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/rda/', '.rda')
ctxt1 <- mergeOv(mergeOv(cpgs, common.pmd), cgi)
ctxt <- ctxt1[ctxt1$ctxt == 'WCGW' & ctxt1$orph35class == '[0,1)' & ctxt1$cgi == 'NonCGI']
pmd <- ctxt[ctxt$sd.binmeanT >= 0.15]
hmd <- ctxt[ctxt$sd.binmeanT <= 0.10]
sampledata <- do.call(rbind, lapply(samples, function(sample) {
  load(paste0('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/rda/', sample, '.rda'))
  cancer <- gr
  cancerhmd <- mergeOv(cancer, hmd)
  cancerpmd <- mergeOv(cancer, pmd)
  df <- rbind(data.frame(betas=cancerpmd$betas, pmd=TRUE),
              data.frame(betas=cancerhmd$betas, pmd=FALSE))
  df$sample <- sample
  df
}))
order2m <- tapply(sampledata$betas, interaction(sampledata$pmd, sampledata$sample), mean, na.rm=T)
sampledata$mbetas <- order2m[interaction(sampledata$pmd, sampledata$sample)]
save(sampledata, file='/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/fetalsampledata.rda')
sampledata$mbetas <- pmax(sampledata$mbetas, 0.6)
sampledata$mbetas <- pmin(sampledata$mbetas, 0.8)
pdf('~/gallery/2017_02_28_fetal_violin.pdf', width=15, height=5)
ggplot(sampledata, aes(interaction(pmd, sample), betas)) + geom_violin(aes(fill=mbetas), scale='area', draw_quantiles=0.5, bw=0.1, alpha=0.8, size=0.3) + theme(axis.text.x=element_text(size=10, angle=90, vjust=0.5, hjust=1)) + ylab('Beta Value') + scale_fill_gradient2(high='#CC0000',low='#005AA0', mid='white', midpoint=0.7, limits=c(0.6,0.8))
dev.off()

## smoothScatter of all chromosomes
## jet.colors <- colorRampPalette(c("#00007F", "#00007F", "#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00"))
jet.colors <-  colorRampPalette(c("#00007F", "#00007F", "#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
ovnames <- intersect(rownames(betas), rownames(mean.binmeanT))
pdf('~/gallery/2017_02_19_early_dev_scatter_smooth_wholegenome.pdf', width=11, height=19)
par(mfrow=c(9,4), mar=c(1,5,5,1), oma=c(4,4,2,2), lwd=2.5, font.lab=1, font.axis=2, cex.lab=0.5, cex.axis=2, cex.main=2, font.main=2)
for (sample in samples) {
  a <- (!is.na(betas[ovnames, sample]) & !is.na(mean.binmeanT[ovnames, '5.0']))
  rho <- cor(mean.binmeanT[ovnames[a], '5.0'], betas[ovnames[a], sample], method='spearman')
  smoothScatter(rank(mean.binmeanT[ovnames[a], '5.0']), rank(betas[ovnames[a], sample]), ylab=sample, main=sprintf('rho: %1.3f', rho), pch=16, nrpoints=0, colramp=jet.colors)
}
dev.off()
for (sample in samples) {
  a <- (!is.na(betas[ovnames, sample]) & !is.na(mean.binmeanT[ovnames, '5.0']))
  rho <- cor(mean.binmeanT[ovnames[a], '5.0'], betas[ovnames[a], sample], method='spearman')
  png(sprintf('~/gallery/2017_02_19_earlydev_wholegenome/%s.png', sample), width=2000, height=2000, res=500)
  par(mar=c(0,0,0,0))
  plot(rank(mean.binmeanT[ovnames[a], '5.0']), rank(betas[ovnames[a], sample]), pch=16, cex=0.2)
  dev.off()
}
for (sample in samples) {
  a <- (!is.na(betas[ovnames, sample]) & !is.na(mean.binmeanT[ovnames, '5.0']))
  rho <- cor(mean.binmeanT[ovnames[a], '5.0'], betas[ovnames[a], sample], method='spearman')
  png(sprintf('~/gallery/2017_02_19_earlydev_smooth/%s.png', sample), width=2000, height=2000, res=500)
  par(mar=c(0,0,0,0))
  smoothScatter(rank(mean.binmeanT[ovnames[a], '5.0']), rank(betas[ovnames[a], sample]), ylab=sample, main=sprintf('rho: %1.3f', rho), pch=16, nrpoints=0, colramp=jet.colors)
  dev.off()
}
