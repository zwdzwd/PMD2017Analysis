########################################
## WGBS-matched hm450 sequence context
########################################
WGBSsamples <- c("TCGA_GBM_1788",  "TCGA_GBM_1460",  "TCGA_GBM_0128",  "TCGA_GBM_3477", "TCGA_GBM_1401", "TCGA_GBM_1454", "TCGA_BRCA_A0YG", "TCGA_BRCA_A0CE", "TCGA_BRCA_A07I", "TCGA_BRCA_A15H", "TCGA_BRCA_A04X", "TCGA_UCEC_A0G2", "TCGA_UCEC_A1CI", "TCGA_UCEC_A0K6", "TCGA_UCEC_A05J", "TCGA_UCEC_A1CK", "TCGA_LUSC_2695", "TCGA_LUSC_1078", "TCGA_LUSC_2722", "TCGA_LUSC_2600", "TCGA_LUAD_6148", "TCGA_LUAD_6215", "TCGA_LUAD_7156", "TCGA_LUAD_4630", "TCGA_LUAD_6840", "TCGA_COAD_3158", "TCGA_COAD_A00R", "TCGA_READ_3593", "TCGA_READ_2689", "TCGA_STAD_6452", "TCGA_STAD_6519", "TCGA_STAD_6177", "TCGA_STAD_5730", "TCGA_BLCA_A2LA", "TCGA_BLCA_A13J", "TCGA_BLCA_A20V", "TCGA_BLCA_A1AG", "TCGA_BLCA_A2HQ", "TCGA_BLCA_A1AA")
cancers <- sapply(strsplit(WGBSsamples, '_'), function(x) x[2])
patients <- sapply(strsplit(WGBSsamples, '_'), function(x) x[3])
load('/primary/projects/laird/projects/2016_04_05_TCGA_pancan_renormalization/450k/merged_mapping.rda')
load('/primary/projects/laird/projects/2016_04_05_TCGA_pancan_renormalization/450k/normal.betas.rda')
hm450.samples <- merged.mapping[substr(merged.mapping$barcode,9,12) %in% patients,]
hm450.samples$patients <- substr(hm450.samples$barcode,9,12)
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450probes.rda')
probes$orphan35class <- cut(probes$orphan35, breaks=c(0,1,2,3,Inf), right=F)
probes$orphan35class <- factor(probes$orphan35class, rev(levels(probes$orphan35class)))
probes <- probes[probes$pmd != 'Neither']
probes <- probes[probes$cgi == 'nonCGI']
pn <- names(probes)
unmeth.cnt <- apply(normal.betas[pn, ], 1, function(x) sum(x <= 0.2, na.rm=T))
pn <- pn[unmeth.cnt < 3]
probes <- probes[pn]
probeclass <- interaction(probes$ctxtclass, probes$orphan35class, probes$pmd)
probename <- names(probes)

for (i in 1:nrow(hm450.samples)) {
  fn <- hm450.samples[i,'idat']
  cat(fn, '\n')
  load(paste0('/primary/projects/laird/projects/2016_04_05_TCGA_pancan_renormalization/450k/betas/', fn, '.rda'))
  df <- data.frame(beta=betas[probename], order=probeclass)
  order2mean <- tapply(df$beta, df$order, mean, na.rm=T)
  df$mbetas <- order2mean[df$order]
  df$mbetas[df$mbetas > 0.9] <- 0.9
  df$mbetas[df$mbetas < 0.3] <- 0.3
  png(sprintf('~/gallery/2017_02_27_hm450_violin/%s_%s_%s.png', hm450.samples[i,'cancer'], hm450.samples[i,'patients'], hm450.samples[i,'barcode']), width=1200, height=200)
  print(ggplot(aes(order, beta), data=df) + geom_violin(aes(fill=mbetas), draw_quantiles=c(0.5), bw=0.05, size=0.8) + theme(axis.text.x = element_text(size=0, angle=90, vjust=0.5, hjust=1), plot.title = element_text(size=40)) + scale_fill_gradient2(high='#CC0000',low='#005AA0', mid='white', midpoint=0.5, limits=c(0.3,0.9)) + xlab('') + ggtitle(paste0(hm450.samples[i,'cancer'], '_', hm450.samples[i,'patients'], ': ', hm450.samples[i,'barcode'])))
  dev.off()
}

###############
## blood
###############

load('~/projects/hs-tcga/2015_03_18_tumor_purity/dataset/blood.rda')
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450probes.rda')
load('/Users/wandingzhou/projects/hs-tcga/data/2015_04_10_sorted_cell_population/sample_age.rda')

aa <- age[age$class=='myeloid',]
a <- do.call(rbind, lapply(aa$sname, function(sn) rbind(data.frame(betas = betas[com.pmd.probes, sn], sn = as.character(sn), pmd='PMD', stringsAsFactors=F), data.frame(betas = betas[com.hmd.probes, sn], sn=as.character(sn), pmd='HMD', stringsAsFactors=F))))
a$pmd <- factor(a$pmd, levels=c('HMD','PMD'))
a$sn <- factor(a$sn, levels=as.character(aa$sname))

pdf('~/gallery/2017_01_04_blood_age_myeloid.pdf', width=11, height=5)
ggplot(a) + annotate("rect", fill = "grey", alpha = 0.3, xmin = c(7,19)-0.5, xmax = c(12,24)+0.5, ymin = -Inf, ymax = Inf) + geom_violin(aes(sn, betas, fill=pmd, color='white'), scale='width', draw_quantiles=c(0.5)) + scale_color_manual(values=c('white')) + scale_fill_manual(values=c('#CC0000','#005AA0')) + scale_x_discrete(labels=paste0(aa$sample,'/',aa$age)) + theme(text=element_text(size=18, face='bold'), axis.text.x = element_text(size=15, angle=90, hjust=1, vjust=0.5), axis.text.y = element_text(size=17)) + ylab('Beta Value') + ggtitle('myeloid')
dev.off()

aa <- age[age$class=='lymphoid',]
a <- do.call(rbind, lapply(aa$sname, function(sn) rbind(data.frame(betas = betas[com.pmd.probes, sn], sn = as.character(sn), pmd='PMD', stringsAsFactors=F), data.frame(betas = betas[com.hmd.probes, sn], sn=as.character(sn), pmd='HMD', stringsAsFactors=F))))
a$pmd <- factor(a$pmd, levels=c('HMD','PMD'))
a$sn <- factor(a$sn, levels=as.character(aa$sname))

pdf('~/gallery/2017_01_04_blood_age_lymphoid.pdf', width=11, height=5)
ggplot(a) + annotate("rect", fill = "grey", alpha = 0.3, xmin = c(7,19)-0.5, xmax = c(12,24)+0.5, ymin = -Inf, ymax = Inf) + geom_violin(aes(sn, betas, fill=pmd, color='white'), scale='width', draw_quantiles=c(0.5)) + scale_color_manual(values=c('white')) + scale_fill_manual(values=c('#CC0000','#005AA0')) + scale_x_discrete(labels=paste0(aa$sample,'/',aa$age)) + theme(text=element_text(size=18, face='bold'), axis.text.x = element_text(size=15, angle=90, hjust=1, vjust=0.5), axis.text.y = element_text(size=17)) + ylab('Beta Value') + ggtitle('lymphoid')
dev.off()

## merged scatter plot
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/Reinus_Blood/betas.rda')
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450probes.rda')
samples$pmd <- apply(betas[com.pmd.probes,], 2, mean, na.rm=T)[rownames(samples)]
samples$hmd <- apply(betas[com.hmd.probes,], 2, mean, na.rm=T)[rownames(samples)]
pdf('~/gallery/2017_01_13_blood_age_vs_pmd.pdf', width=5.5, height=3)
ggplot(samples, aes(age, pmd, color=sample)) + geom_jitter() + xlab('Age') + ylab('PMD depth') + scale_color_discrete(guide=guide_legend(title="Cell Types")) + theme(text=element_text(face='bold', size=18), axis.text = element_text(colour='black', size=16), panel.border = element_rect(linetype = "solid", colour = "black", size=1.2))
dev.off()

## axis.line.x = element_line(colour = "black", size=1.2), axis.line.y = element_line(colour = "black", size=1.2), 

## by cell types
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/Reinus_Blood/betas.rda')
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450probes.rda')
samples$pmd <- apply(betas[com.pmd.probes,], 2, mean, na.rm=T)[rownames(samples)]
samples$hmd <- apply(betas[com.hmd.probes,], 2, mean, na.rm=T)[rownames(samples)]
celltype2samples <- split(samples, samples$sample)
for (sn in names(celltype2samples)) {
  pdf(sprintf('~/gallery/2017_01_14_blood_age_vs_pmd_%s.pdf', sn), width=3.5, height=3)
  print(ggplot(celltype2samples[[sn]], aes(age, pmd)) + geom_smooth(se=FALSE, linetype='solid', size=2, method='lm') + geom_point(size=1) + xlab('Age') + ylab('PMD methylation') + ylim(0,1) + ggtitle(sn) + theme(text=element_text(face='bold', size=18), axis.text = element_text(colour='black', size=16)))
  dev.off()
}
for (sn in names(celltype2samples)) {
  cat(sn,'\n'); cat(dim(celltype2samples[[sn]]),'\n');
  print(cor.test(celltype2samples[[sn]]$age, celltype2samples[[sn]]$pmd));
  cat(lm(pmd~age, data=celltype2samples[[sn]])$coefficients['age'], '\n');
}

###########################################
## Renauer paper, sorted blood from PBMC
###########################################
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/GSE61195/betas.rda')
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450probes.rda')
samples$pmd <- apply(betas[com.pmd.probes,], 2, mean, na.rm=T)[rownames(samples)]
samples$hmd <- apply(betas[com.hmd.probes,], 2, mean, na.rm=T)[rownames(samples)]
celltype2samples <- split(samples, samples$sourceName)
for (sn in names(celltype2samples)) {
  pdf(sprintf('~/gallery/2017_01_14_blood_age_vs_pmd_Tcell_%s.pdf', sn), width=3.5, height=3)
  print(ggplot(celltype2samples[[sn]], aes(age, pmd)) + geom_smooth(se=FALSE, linetype='solid', size=2, method='lm') + geom_point(size=1) + xlab('Age') + ylab('PMD methylation') + ylim(0,1) + ggtitle(sn) + theme(text=element_text(face='bold', size=18), axis.text = element_text(colour='black', size=16)))
  dev.off()
}

## Tserel paper CD4 CD8 blood
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/GSE59065/betas.rda')
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450probes.rda')
samples$pmd <- apply(betas[com.pmd.probes,], 2, mean, na.rm=T)[rownames(samples)]
samples$hmd <- apply(betas[com.hmd.probes,], 2, mean, na.rm=T)[rownames(samples)]
celltype2samples <- split(samples, samples[['cell/tissue type']])
for (sn in names(celltype2samples)) {
  pdf(sprintf('~/gallery/2017_01_14_blood_age_vs_pmd_Tserel_%s.pdf', sn), width=3.5, height=3)
  print(ggplot(celltype2samples[[sn]], aes(age, pmd)) + geom_smooth(, linetype='solid', size=2, method='lm') + geom_point(size=1) + xlab('Age') + ylab('PMD methylation') + ylim(0,1) + ggtitle(sn) + theme(text=element_text(face='bold', size=18), axis.text = element_text(colour='black', size=16)))
  dev.off()
}
for (sn in names(celltype2samples)) {
  cat(sn,'\n'); cat(dim(celltype2samples[[sn]]),'\n');
  print(cor.test(celltype2samples[[sn]]$age, celltype2samples[[sn]]$pmd));
  cat(lm(pmd~age, data=celltype2samples[[sn]])$coefficients['age'], '\n');
}


## Reynolds CD14 paper
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/GSE56046/samplespmd.rda')
pdf('~/gallery/2017_01_14_blood_age_vs_pmd_Reynold.pdf', width=3.5, height=3)
print(ggplot(samples, aes(age, pmd)) + geom_jitter(size=1) + geom_smooth(method='lm') + ylim(0,1) + ggtitle('CD14+') + xlab('Age') + ylab('PMD methylation') + theme(text=element_text(face='bold', size=18), axis.text = element_text(colour='black', size=16)))
dev.off()

lm(pmd~age, data=samples)$coefficients['age']

## Limbach CD4, CD8 Grave's diseease paper
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/GSE71955/betas.rda')
samples$pmd <- apply(betas[intersect(com.pmd.probes, rownames(betas)),], 2, mean, na.rm=T)[rownames(samples)]
samples$hmd <- apply(betas[intersect(com.hmd.probes, rownames(betas)),], 2, mean, na.rm=T)[rownames(samples)]
celltype2samples <- split(samples, samples[['celltype']])
for (sn in names(celltype2samples)) {
  pdf(sprintf('~/gallery/2017_01_14_blood_age_vs_pmd_Limbach_%s.pdf', sn), width=3.5, height=3)
  print(ggplot(celltype2samples[[sn]], aes(age, pmd)) + geom_smooth(, linetype='solid', size=2, method='lm') + geom_point(size=1) + xlab('Age') + ylab('PMD methylation') + ylim(0,1) + ggtitle(sn) + theme(text=element_text(face='bold', size=18), axis.text = element_text(colour='black', size=16)))
  dev.off()
}
for (sn in names(celltype2samples)) {
  cat(sn,'\n'); cat(dim(celltype2samples[[sn]]),'\n');
  print(cor.test(celltype2samples[[sn]]$age, celltype2samples[[sn]]$pmd));
  cat(lm(pmd~age, data=celltype2samples[[sn]])$coefficients['age'], '\n');
}


################################################
## GSE30870, newborn vs very old (nonagenarian)
################################################
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450probes.rda')
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/GSE30870/betas.rda')
samples$pmd <- apply(betas[intersect(com.pmd.probes, rownames(betas)),], 2, mean, na.rm=T)[rownames(samples)]
samples$hmd <- apply(betas[intersect(com.hmd.probes, rownames(betas)),], 2, mean, na.rm=T)[rownames(samples)]
samples$age <- gsub('.years','',samples$age)
samples$age[samples$age == 'Newborn'] <- 0
samples$age <- as.numeric(samples$age)
pdf('~/gallery/2017_01_17_PBMNC_newborn_vs_old.pdf', width=4, height=4)
ggplot(samples, aes(age, pmd)) + geom_point() + geom_smooth(method='lm') + xlab('Age') + ylab('PMD methylation') + ggtitle('GSE30870')
dev.off()
pdf('~/gallery/2017_01_17_PBMNC_newborn_vs_old_boxplot.pdf', width=3.5, height=3.7)
ggplot(samples, aes(sourceName, pmd)) + geom_boxplot() + xlab('') + ylab('PMD methylation') + ggtitle('PBMC (GSE30870)') + theme(axis.text.x=element_text(angle=20, hjust=1))
dev.off()

## co-plot with liver
samples0 <- data.frame(stage=samples$sourceName, pmd=samples$pmd)
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/GSE61278/betas.rda')
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450probes.rda')
samples$pmd <- apply(betas[intersect(com.pmd.probes, rownames(betas)),], 2, mean, na.rm=T)[rownames(samples)]
samples$hmd <- apply(betas[intersect(com.hmd.probes, rownames(betas)),], 2, mean, na.rm=T)[rownames(samples)]
samples0 <- rbind(samples0, data.frame(
  stage=samples[samples$group!='Accident','agegroup'],
  pmd=samples[samples$group!='Accident','pmd']))
samples0$stage <- factor(samples0$stage, levels=c('Newborns','Nonagenarians','Fetal','Adult'))
pdf('~/gallery/2017_01_17_PBMNC_newborn_vs_old_boxplot_wLiver.pdf', width=3.5, height=3.7)
ggplot(samples0, aes(stage, pmd)) + geom_boxplot() + xlab('') + ylab('PMD methylation') + theme(axis.text.x=element_text(angle=20, hjust=1))
dev.off()

## with Ben's probe set
samples$pmd <- apply(betas[intersect(com.pmd.probesB, rownames(betas)),], 2, mean, na.rm=T)[rownames(samples)]
samples$hmd <- apply(betas[intersect(com.hmd.probesB, rownames(betas)),], 2, mean, na.rm=T)[rownames(samples)]
samples$age <- gsub('.years','',samples$age)
samples$age[samples$age == 'Newborn'] <- 0
samples$age <- as.numeric(samples$age)
pdf('~/gallery/2017_01_17_PBMNC_newborn_vs_old_boxplotB.pdf', width=3.5, height=3.7)
ggplot(samples, aes(sourceName, pmd)) + geom_boxplot() + xlab('') + ylab('PMD methylation') + ggtitle('PBMC (GSE30870)') + theme(axis.text.x=element_text(angle=20, hjust=1))
dev.off()

##########################
## colon and brain
##########################

## brain and age
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/Sorted_Brain/betas.rda')
samples$pmd <- apply(betas[intersect(com.pmd.probes, rownames(betas)),], 2, mean, na.rm=T)[rownames(samples)]
samples$hmd <- apply(betas[intersect(com.hmd.probes, rownames(betas)),], 2, mean, na.rm=T)[rownames(samples)]
celltype2samples <- split(samples, samples[['celltype']])
for (sn in names(celltype2samples)) {
  pdf(sprintf('~/gallery/2017_01_14_brain_age_vs_pmd_%s.pdf', sn), width=3.5, height=3)
  print(ggplot(celltype2samples[[sn]], aes(age, pmd)) + geom_smooth(, linetype='solid', size=2, method='lm') + geom_point(size=1) + xlab('Age') + ylab('PMD methylation') + ylim(0,1) + ggtitle(sn) + theme(text=element_text(face='bold', size=18), axis.text = element_text(colour='black', size=16)))
  dev.off()
}
for (sn in names(celltype2samples)) {
  cat(sn,'\n'); cat(dim(celltype2samples[[sn]]),'\n');
  print(cor.test(celltype2samples[[sn]]$age, celltype2samples[[sn]]$pmd));
  cat(lm(pmd~age, data=celltype2samples[[6]])$coefficients['age'], '\n');
}

## colon mucosa
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/GSE32146/betas.rda')
samples$pmd <- apply(betas[intersect(com.pmd.probes, rownames(betas)),], 2, mean, na.rm=T)[rownames(samples)]
samples$hmd <- apply(betas[intersect(com.hmd.probes, rownames(betas)),], 2, mean, na.rm=T)[rownames(samples)]
pdf('~/gallery/2017_01_13_colonic_mucosa_age_vs_pmd.pdf', width=3.5, height=3)
print(ggplot(samples, aes(age, pmd)) + geom_smooth(, linetype='solid', size=2, method='lm') + geom_point(size=1) + xlab('Age') + ylab('PMD methylation') + ylim(0,1) + ggtitle('Colonic Mucosa') + theme(text=element_text(face='bold', size=18), axis.text = element_text(colour='black', size=16)))
dev.off()

############################
## old woman
############################

load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/GSE64491/betas.rda')
samples$pmd <- apply(betas[com.pmd.probes,], 2, mean, na.rm=T)[rownames(samples)]
samples$hmd <- apply(betas[com.hmd.probes,], 2, mean, na.rm=T)[rownames(samples)]
pdf('~/gallery/2017_01_17_old_woman.pdf', width=8, height=5)
ggplot(samples, aes(tissue, pmd)) + geom_point() + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) + xlab('Tissue') + ylab('PMD methylation')
dev.off()

pdf('~/gallery/2017_01_17_old_woman_hmdvspmd.pdf', width=9, height=5)
ggplot(samples, aes(hmd, pmd, color=tissue)) + geom_point() + geom_abline(intercept=0, slope=1) + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) + xlab('HMD methylation') + ylab('PMD methylation')
dev.off()

## liver and pericardium are the two biggest hits

## check liver
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/GSE61278/betas.rda')
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450probes.rda')
samples$pmd <- apply(betas[intersect(com.pmd.probes, rownames(betas)),], 2, mean, na.rm=T)[rownames(samples)]
samples$hmd <- apply(betas[intersect(com.hmd.probes, rownames(betas)),], 2, mean, na.rm=T)[rownames(samples)]
pdf('~/gallery/2017_01_17_liver_age_group.pdf', width=5, height=2.85)
ggplot(samples, aes(hmd, pmd, color=agegroup)) + geom_point() + xlab('HMD Methylation') + ylab('PMD Methylation')
dev.off()

## liver vs age
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/GSE61258/betas.rda')
samples$age <- as.numeric(samples$age)
samples$pmd <- apply(betas[com.pmd.probes,], 2, mean, na.rm=T)[rownames(samples)]
samples$hmd <- apply(betas[com.hmd.probes,], 2, mean, na.rm=T)[rownames(samples)]
pdf('~/gallery/2017_01_17_liver_vs_age.pdf', width=6, height=4)
ggplot(samples, aes(hmd, pmd, color=diseasestatus)) + geom_point() + geom_abline(slope=1, intercept=0) + xlab('HMD methylation') + ylab('PMD methylation')
dev.off()

## bone vs age
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/GSE64490/betas.rda')
samples$pmd <- apply(betas[com.pmd.probes,], 2, mean, na.rm=T)[rownames(samples)]
samples$hmd <- apply(betas[com.hmd.probes,], 2, mean, na.rm=T)[rownames(samples)]
## no correlation with age, at least from this plot

## progeria
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/GSE42865/betas.rda')
samples$pmd <- apply(betas[com.pmd.probes,], 2, mean, na.rm=T)[rownames(samples)]
samples$hmd <- apply(betas[com.hmd.probes,], 2, mean, na.rm=T)[rownames(samples)]
pdf('~/gallery/2017_01_17_progeria_pmd_vs_hmd.pdf', width=8, height=4)
ggplot(samples, aes(hmd, pmd, color=diseasestate)) + geom_point() + geom_abline(slope=1, intercept=0) + xlab('HMD methylation') + ylab('PMD methylation')
dev.off()

##################
## merged blood
##################

load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450probes.rda')
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/Reinus_Blood/betas.rda')
samples$pmd <- apply(betas[com.pmd.probes,], 2, mean, na.rm=T)[rownames(samples)]
samples$hmd <- apply(betas[com.hmd.probes,], 2, mean, na.rm=T)[rownames(samples)]
samples0 <- data.frame(celltype=samples$sample, class=samples$class, pmd=samples$pmd, hmd=samples$hmd, age=samples$age, study='GSE35069')

## load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/GSE61195/betas.rda')
## samples$pmd <- apply(betas[com.pmd.probes,], 2, mean, na.rm=T)[rownames(samples)]
## samples$hmd <- apply(betas[com.hmd.probes,], 2, mean, na.rm=T)[rownames(samples)]
## samples <- samples[samples$celltype %in% c('TCRαβ+.CD4+.T.cells', 'TCRαβ+.CD8+.T.cells'),]
## samples$celltype[samples$celltype=='TCRαβ+.CD4+.T.cells'] <- 'CD4.Tcell'
## samples$celltype[samples$celltype=='TCRαβ+.CD8+.T.cells'] <- 'CD8.Tcell'
## samples0 <- rbind(samples0, data.frame(celltype=samples$celltype, class='lymphoid', pmd=samples$pmd, hmd=samples$hmd, age=samples$age, study='GSE61195'))

## print correlation coefficient
celltype2samples <- split(samples, samples$sample)
for (sn in names(celltype2samples)) {
  cat(sn,'\n'); cat(dim(celltype2samples[[sn]]),'\n');
  print(cor.test(celltype2samples[[sn]]$age, celltype2samples[[sn]]$pmd));
  cat(lm(pmd~age, data=celltype2samples[[sn]])$coefficients['age'], '\n');
}

load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/GSE59065/betas.rda')
samples$pmd <- apply(betas[com.pmd.probes,], 2, mean, na.rm=T)[rownames(samples)]
samples$hmd <- apply(betas[com.hmd.probes,], 2, mean, na.rm=T)[rownames(samples)]
colnames(samples)[7] <- 'celltype'
samples <- samples[samples$celltype %in% c('CD4','CD8'),]
samples$celltype[samples$celltype == 'CD4'] <- 'CD4.Tcell'
samples$celltype[samples$celltype == 'CD8'] <- 'CD8.Tcell'
samples0 <- rbind(samples0, data.frame(celltype=samples$celltype, class='lymphoid', pmd=samples$pmd, hmd=samples$hmd, age=samples$age, study='GSE56065'))

## print correlation coefficient
celltype2samples <- split(samples, samples$celltype)
for (sn in names(celltype2samples)) {
  cat(sn,'\n'); cat(dim(celltype2samples[[sn]]),'\n');
  print(cor.test(celltype2samples[[sn]]$age, celltype2samples[[sn]]$pmd));
  cat(lm(pmd~age, data=celltype2samples[[sn]])$coefficients['age'], '\n');
}

load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/GSE71955/betas.rda')
samples$pmd <- apply(betas[intersect(com.pmd.probes, rownames(betas)),], 2, mean, na.rm=T)[rownames(samples)]
samples$hmd <- apply(betas[intersect(com.hmd.probes, rownames(betas)),], 2, mean, na.rm=T)[rownames(samples)]
samples <- samples[samples$diagnosis == 'Healthy',]
samples$celltype[samples$celltype == 'CD4.T.cells'] <- 'CD4.Tcell'
samples$celltype[samples$celltype == 'CD8.T.cells'] <- 'CD8.Tcell'
samples0 <- rbind(samples0, data.frame(celltype=samples$celltype, class='lymphoid', pmd=samples$pmd, hmd=samples$hmd, age=samples$age, study='GSE71955'))

## print correlation coefficient
celltype2samples <- split(samples, samples$celltype)
for (sn in names(celltype2samples)) {
  cat(sn,'\n'); cat(dim(celltype2samples[[sn]]),'\n');
  print(cor.test(celltype2samples[[sn]]$age, celltype2samples[[sn]]$pmd));
  cat(lm(pmd~age, data=celltype2samples[[sn]])$coefficients['age'], '\n');
}

## on hpc
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/GSE56046/betas.rda')
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450probes.rda')
samples$pmd <- apply(betas[com.pmd.probes,], 2, mean, na.rm=T)[rownames(samples)]
samples$hmd <- apply(betas[com.hmd.probes,], 2, mean, na.rm=T)[rownames(samples)]
save(samples, file='/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/GSE56046/samplespmd.rda')

load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/GSE56046/samplespmd.rda')
samples0 <- rbind(samples0, data.frame(celltype='CD14.Monocytes', class='myeloid', pmd=samples$pmd, hmd=samples$hmd, age=samples$age, study='GSE56046'))

## print correlation coefficient
cat(dim(samples),'\n');
print(cor.test(samples$age, samples$pmd));
cat(lm(pmd~age, data=samples)$coefficients['age'], '\n');

samples0$ss <- samples0$study
samples0$study[samples0$study=='GSE35069'] <- 'S1'
## samples0$study[samples0$study=='GSE61195'] <- 'S2'
samples0$study[samples0$study=='GSE56065'] <- 'S2'
samples0$study[samples0$study=='GSE71955'] <- 'S3'
samples0$study[samples0$study=='GSE56046'] <- 'S4'

samples.l <- subset(samples0, class=='lymphoid')
samples.l$source <- interaction(samples.l$study, samples.l$celltype)

samples.m <- subset(samples0, class=='myeloid')
samples.m$source <- interaction(samples.m$study, samples.m$celltype)
samples.m <- rbind(samples.m[samples.m$study=='S4',], samples.m[samples.m$study!='S4',])


pdf('~/gallery/2017_01_16_blood_age_vs_pmd_lymphoid.pdf', width=6, height=4.5)
ggplot(samples.l, aes(age, pmd, color=source)) + geom_point() + geom_smooth(method='lm', se=FALSE, size=0.5) + ylim(0.5,0.85) + xlim(20,85) + scale_color_manual(values=c('S1.CD8.Tcell'='#FF530D', 'S2.CD8.Tcell'='#FF0000', 'S3.CD8.Tcell'='#E80C7A', 'S1.CD4.Tcell'='#95AB63', 'S2.CD4.Tcell'='#33A02C', 'S3.CD4.Tcell'='#468966', 'S1.CD56.NKcell'='#3498DB', 'S1.CD19.Bcell'='#DB9E36'), guide=guide_legend(title='Source')) + scale_shape_discrete(guide=guide_legend(title='Cell Type')) + xlab('Age') + ylab('PMD methylation')
dev.off()

pdf('~/gallery/2017_01_16_blood_age_vs_pmd_myeloid.pdf', width=6, height=4.5)
ggplot(samples.m, aes(age, pmd, color=source)) + geom_point(size=0.5) + geom_smooth(method='lm', se=FALSE, size=0.5) + ylim(0.5,0.85) + xlim(20,85) + scale_color_discrete(guide=guide_legend(title='Source')) + scale_shape_discrete(guide=guide_legend(title='Cell Type')) + xlab('Age') + ylab('PMD methylation')
dev.off()

##############################################
## merged blood Ben's high variance probes
##############################################

## Ben's high variance probes
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450probesB.rda')
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/Reinus_Blood/betas.rda')
samples$pmd <- apply(betas[com.pmd.probesB,], 2, mean, na.rm=T)[rownames(samples)]
samples$hmd <- apply(betas[com.hmd.probesB,], 2, mean, na.rm=T)[rownames(samples)]
samples0 <- data.frame(celltype=samples$sample, class=samples$class, pmd=samples$pmd, hmd=samples$hmd, age=samples$age, study='GSE35069')

## print correlation coefficient
celltype2samples <- split(samples, samples$sample)
for (sn in names(celltype2samples)) {
  cat(sn,'\n'); cat(dim(celltype2samples[[sn]]),'\n');
  print(cor.test(celltype2samples[[sn]]$age, celltype2samples[[sn]]$pmd));
  cat(lm(pmd~age, data=celltype2samples[[sn]])$coefficients['age'], '\n');
}

load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/GSE59065/betas.rda')
samples$pmd <- apply(betas[com.pmd.probesB,], 2, mean, na.rm=T)[rownames(samples)]
samples$hmd <- apply(betas[com.hmd.probesB,], 2, mean, na.rm=T)[rownames(samples)]
colnames(samples)[7] <- 'celltype'
samples <- samples[samples$celltype %in% c('CD4','CD8'),]
samples$celltype[samples$celltype == 'CD4'] <- 'CD4.Tcell'
samples$celltype[samples$celltype == 'CD8'] <- 'CD8.Tcell'
samples0 <- rbind(samples0, data.frame(celltype=samples$celltype, class='lymphoid', pmd=samples$pmd, hmd=samples$hmd, age=samples$age, study='GSE59065'))

## print correlation coefficient
celltype2samples <- split(samples, samples$celltype)
for (sn in names(celltype2samples)) {
  cat(sn,'\n'); cat(dim(celltype2samples[[sn]]),'\n');
  print(cor.test(celltype2samples[[sn]]$age, celltype2samples[[sn]]$pmd));
  cat(lm(pmd~age, data=celltype2samples[[sn]])$coefficients['age'], '\n');
}

load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/GSE71955/betas.rda')
samples$pmd <- apply(betas[intersect(com.pmd.probesB, rownames(betas)),], 2, mean, na.rm=T)[rownames(samples)]
samples$hmd <- apply(betas[intersect(com.hmd.probesB, rownames(betas)),], 2, mean, na.rm=T)[rownames(samples)]
samples <- samples[samples$diagnosis == 'Healthy',]
samples$celltype[samples$celltype == 'CD4.T.cells'] <- 'CD4.Tcell'
samples$celltype[samples$celltype == 'CD8.T.cells'] <- 'CD8.Tcell'
samples0 <- rbind(samples0, data.frame(celltype=samples$celltype, class='lymphoid', pmd=samples$pmd, hmd=samples$hmd, age=samples$age, study='GSE71955'))

## print correlation coefficient
celltype2samples <- split(samples, samples$celltype)
for (sn in names(celltype2samples)) {
  cat(sn,'\n'); cat(dim(celltype2samples[[sn]]),'\n');
  print(cor.test(celltype2samples[[sn]]$age, celltype2samples[[sn]]$pmd));
  cat(lm(pmd~age, data=celltype2samples[[sn]])$coefficients['age'], '\n');
}

load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/GSE56046/samplespmdB.rda')
samples <- data.frame(celltype='CD14.Monocytes', class='myeloid', pmd=samples$pmd, hmd=samples$hmd, age=samples$age, study='GSE56046')
samples0 <- rbind(samples0, samples)

## print correlation coefficient
cat(sn,'\n'); cat(dim(samples),'\n');
print(cor.test(samples$age, samples$pmd));
cat(lm(pmd~age, data=samples)$coefficients['age'], '\n');

samples0$ss <- samples0$study
samples0$study[samples0$study=='GSE35069'] <- 'S1'
## samples0$study[samples0$study=='GSE61195'] <- 'S2'
samples0$study[samples0$study=='GSE59065'] <- 'S2'
samples0$study[samples0$study=='GSE71955'] <- 'S3'
samples0$study[samples0$study=='GSE56046'] <- 'S4'

samples.l <- subset(samples0, class=='lymphoid')
samples.l$source <- interaction(samples.l$study, samples.l$celltype)

samples.m <- subset(samples0, class=='myeloid')
samples.m$source <- interaction(samples.m$study, samples.m$celltype)
samples.m <- rbind(samples.m[samples.m$study=='S4',], samples.m[samples.m$study!='S4',])

pdf('~/gallery/2017_01_16_blood_age_vs_pmd_lymphoidB.pdf', width=6, height=4.5)
ggplot(samples.l, aes(age, pmd, color=source)) + geom_point() + geom_smooth(method='lm', se=FALSE, size=0.5) + ylim(0.4,0.85) + xlim(20,85) + scale_color_manual(values=c('S1.CD8.Tcell'='#FF530D', 'S2.CD8.Tcell'='#FF0000', 'S3.CD8.Tcell'='#E80C7A', 'S1.CD4.Tcell'='#95AB63', 'S2.CD4.Tcell'='#33A02C', 'S3.CD4.Tcell'='#468966', 'S1.CD56.NKcell'='#3498DB', 'S1.CD19.Bcell'='#DB9E36'), guide=guide_legend(title='Source')) + scale_shape_discrete(guide=guide_legend(title='Cell Type')) + xlab('Age') + ylab('PMD methylation')
dev.off()

pdf('~/gallery/2017_01_16_blood_age_vs_pmd_myeloidB.pdf', width=6, height=4.5)
ggplot(samples.m, aes(age, pmd, color=source)) + geom_point(size=1) + geom_smooth(method='lm', se=FALSE, size=0.5) + ylim(0.4,0.85) + xlim(20,85) + scale_color_discrete(guide=guide_legend(title='Source')) + scale_shape_discrete(guide=guide_legend(title='Cell Type')) + xlab('Age') + ylab('PMD methylation')
dev.off()

########################
## skin and fibroblast
########################

load('~/projects/hs-tcga/2015_03_18_tumor_purity/dataset/fibroblast_wagner.rda')
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450probes.rda')
a <- do.call(rbind, lapply(colnames(betas), function(sn) rbind(
  data.frame(betas = betas[com.pmd.probes, sn], sn = as.character(sn), pmd='PMD', stringsAsFactors=F),
  data.frame(betas = betas[com.hmd.probes, sn], sn = as.character(sn), pmd='HMD', stringsAsFactors=F))))
a$pmd <- factor(a$pmd, levels=c('HMD','PMD'))
ap <- a[a$pmd=='PMD',]
depths <- tapply(ap$betas, ap$sn, mean, na.rm=T)
sampleorder <- names(sort(depths))
a$cn <- factor(a$sn, levels=sampleorder)

pdf('~/gallery/2017_01_05_fibroblast.pdf', width=20, height=5)
ggplot(a) + annotate("rect", fill = "grey", alpha = 0.3, xmin = 1:(length(sampleorder)/2)*2-0.5, xmax = 1:(length(sampleorder)/2)*2+0.5, ymin = -Inf, ymax = Inf) + geom_violin(aes(sn, betas, fill=pmd, color='white'), scale='area', draw_quantiles=c(0.5)) + scale_color_manual(values=c('white')) + scale_fill_manual(values=c('#CC0000','#005AA0')) + scale_x_discrete(labels=sapply(strsplit(sampleorder, '_'), function(x) paste0(x[3],'_',x[4]))) + theme(text=element_text(size=17, face='bold'), axis.text.x = element_text(size=15, angle=90, hjust=1, vjust=0.5), axis.text.y = element_text(size=15)) + xlab('sample') + ylab('Beta Value') + ggtitle('Fibroblast (Wagner et al., GSM1257669)')
dev.off()

adf <- as.data.frame(tapply(a$betas, list(a$sn, a$pmd), mean, na.rm=T))
pdf('~/gallery/2017_01_05_fibroblast_scatter.pdf', width=3.5, height=3)
ggplot(aes(HMD, PMD), data=adf) + geom_point(size=1) + geom_abline(aes(slope=1, intercept=0), linetype='dashed') + xlim(0.2,0.7) + ylim(0.2,0.7) + xlab('HMD methylation') + ylab('PMD methylation')
dev.off()

depths <- adf
depths$tissue <- 'fibroblast'
save(depths, file='/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/fibroblast.rda')

## peri-umbilical punch skin, GSE90124
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/GSE90124/betas.rda')
samples$pmd <- apply(betas[intersect(com.pmd.probes, rownames(betas)),], 2, mean, na.rm=T)[rownames(samples)]
samples$hmd <- apply(betas[intersect(com.hmd.probes, rownames(betas)),], 2, mean, na.rm=T)[rownames(samples)]
pdf('~/gallery/2017_01_17_skin_pmd_vs_age.pdf', width=5, height=4)
ggplot(samples, aes(age.at.biopsy, pmd, color=as.factor(smoking))) + geom_point() + geom_smooth(method='lm',se=FALSE) + scale_color_discrete(guide=guide_legend(title='Smoking Status'))
dev.off()

## skin
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/GSE51954/betas.rda')
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450probes.rda')
## samples$age[samples$age=='>90'] <- 90
## samples$age <- as.numeric(samples$age)
samples$pmd <- apply(betas[com.pmd.probes,], 2, mean, na.rm=T)[rownames(samples)]
samples$hmd <- apply(betas[com.hmd.probes,], 2, mean, na.rm=T)[rownames(samples)]
pdf('~/gallery/2017_01_17_skin_age_vs_pmd.pdf', width=7.0, height=4.5)
ggplot(samples, aes(age, pmd, color=sun)) + geom_point() + geom_smooth(method='lm') + facet_wrap(~celltype) + xlab('Age') + ylab('PMD methylation') + ggtitle('GSE51954')
dev.off()

cor.test(samples[samples$celltype == 'dermis' & samples$sun == 'protected','pmd'], samples[samples$celltype == 'dermis' & samples$sun == 'protected', 'age'])
cat(lm(pmd~age, data=samples[samples$celltype == 'dermis' & samples$sun == 'protected',])$coefficients['age'], '\n');

cor.test(samples[samples$celltype == 'dermis' & samples$sun == 'exposed','pmd'], samples[samples$celltype == 'dermis' & samples$sun == 'exposed', 'age'])
cat(lm(pmd~age, data=samples[samples$celltype == 'dermis' & samples$sun == 'exposed',])$coefficients['age'], '\n');

cor.test(samples[samples$celltype == 'epidermis' & samples$sun == 'protected','pmd'], samples[samples$celltype == 'epidermis' & samples$sun == 'protected', 'age'])
cat(lm(pmd~age, data=samples[samples$celltype == 'epidermis' & samples$sun == 'protected',])$coefficients['age'], '\n');

cor.test(samples[samples$celltype == 'epidermis' & samples$sun == 'exposed','pmd'], samples[samples$celltype == 'epidermis' & samples$sun == 'exposed', 'age'])
cat(lm(pmd~age, data=samples[samples$celltype == 'epidermis' & samples$sun == 'exposed',])$coefficients['age'], '\n');

## Ben's list
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/GSE51954/betas.rda')
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450probesB.rda')
## samples$age[samples$age=='>90'] <- 90
## samples$age <- as.numeric(samples$age)
samples$pmd <- apply(betas[com.pmd.probesB,], 2, mean, na.rm=T)[rownames(samples)]
samples$hmd <- apply(betas[com.hmd.probesB,], 2, mean, na.rm=T)[rownames(samples)]
pdf('~/gallery/2017_01_17_skin_age_vs_pmdB.pdf', width=8, height=5)
ggplot(samples, aes(age, pmd, color=sun)) + geom_point() + geom_smooth(method='lm') + facet_wrap(~celltype) + xlab('Age') + ylab('PMD methylation') + ggtitle('GSE51954')
dev.off()

cor.test(samples[samples$celltype == 'dermis' & samples$sun == 'protected','pmd'], samples[samples$celltype == 'dermis' & samples$sun == 'protected', 'age'])

cor.test(samples[samples$celltype == 'dermis' & samples$sun == 'exposed','pmd'], samples[samples$celltype == 'dermis' & samples$sun == 'exposed', 'age'])

cor.test(samples[samples$celltype == 'epidermis' & samples$sun == 'protected','pmd'], samples[samples$celltype == 'epidermis' & samples$sun == 'protected', 'age'])

cor.test(samples[samples$celltype == 'epidermis' & samples$sun == 'exposed','pmd'], samples[samples$celltype == 'epidermis' & samples$sun == 'exposed', 'age'])

###############################
## fetal development 4 tissue
###############################
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/GSE56515/betas.rda')
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450probes.rda')
samples$pmd <- apply(betas[intersect(com.pmd.probes, rownames(betas)),], 2, mean, na.rm=T)[rownames(samples)]
samples$hmd <- apply(betas[intersect(com.hmd.probes, rownames(betas)),], 2, mean, na.rm=T)[rownames(samples)]
pdf('~/gallery/2017_02_26_fetal_4tissue.pdf', width=4.5, height=3.5)
ggplot(samples, aes(time.point, pmd, color=tissue)) + geom_point() + stat_smooth(method='lm') + xlab('Time Past Gestation') + ylab('PMD Methylation')
dev.off()

celltype2samples <- split(samples, samples$tissue)
for (sn in names(celltype2samples)) {
  cat(sn,'\n'); cat(dim(celltype2samples[[sn]]),'\n');
  print(cor.test(celltype2samples[[sn]]$time.point, celltype2samples[[sn]]$pmd));
  cat(lm(pmd~time.point, data=celltype2samples[[sn]])$coefficients['time.point'], '\n');
}


############################
## MSI data
############################
setwd('/primary/projects/laird/projects/2016_04_05_TCGA_pancan_renormalization/450k/msi')
msi <- read.table('nm.4191-S3.txt', header=T, sep='\t')
msi <- msi[,c(1,7,8)]
colnames(msi) <- c('sample', 'msi', 'cancer')
save(msi, file='/primary/projects/laird/projects/2016_04_05_TCGA_pancan_renormalization/450k/msi/msi.rda')

############################
## all TCGA
############################

setwd('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS')
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450probes.rda')
allfns <- list.files('/primary/projects/laird/projects/2016_04_05_TCGA_pancan_renormalization/450k/betas')
depths <- mclapply(allfns, function(fn) {
  load(paste0('/primary/projects/laird/projects/2016_04_05_TCGA_pancan_renormalization/450k/betas/', fn))
  c(mean(betas[com.pmd.probes], na.rm=T), mean(betas[com.hmd.probes], na.rm=T))
}, mc.cores=20)
depths <- t(simplify2array(depths))
rownames(depths) <- gsub('.rda','',allfns)
colnames(depths) <- c('commonPMD','commonHMD')
save(depths, file='/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/allTCGA.rda')

## ben's high var list
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450probesB.rda')
allfns <- list.files('/primary/projects/laird/projects/2016_04_05_TCGA_pancan_renormalization/450k/betas')
depths <- mclapply(allfns, function(fn) {
  load(paste0('/primary/projects/laird/projects/2016_04_05_TCGA_pancan_renormalization/450k/betas/', fn))
  c(mean(betas[com.pmd.probesB], na.rm=T), mean(betas[com.hmd.probesB], na.rm=T))
}, mc.cores=20)
depths <- t(simplify2array(depths))
rownames(depths) <- gsub('.rda','',allfns)
colnames(depths) <- c('commonPMD','commonHMD')
save(depths, file='/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/allTCGA-B.rda')

merged.mapping <- read.table('/primary/projects/laird/projects/2016_04_05_TCGA_pancan_renormalization/450k/merged_mapping', sep='\t', stringsAsFactors = F, col.names = c('cancer','barcode','idat'))
merged.mapping.normal <- merged.mapping[substr(merged.mapping$barcode, 14,14) == '1',]
merged.mapping.tumor <- merged.mapping[substr(merged.mapping$barcode, 14,14) == '0',]
idat2barcode <- setNames(merged.mapping.normal$barcode, merged.mapping.normal$idat)
barcode2cancer <- setNames(merged.mapping.normal$cancer, merged.mapping.normal$barcode)
barcode2idat <- setNames(merged.mapping.normal$idat, merged.mapping.normal$barcode)
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/allTCGA.rda')

## normal
merged.mapping.normal <- cbind(merged.mapping.normal, depths[merged.mapping.normal$idat,])
merged.mapping.normal$tissue <- 'Others'
merged.mapping.normal$tissue[merged.mapping.normal$cancer == 'ESCA'] <- 'Esophagus (ESCA)'
merged.mapping.normal$tissue[merged.mapping.normal$cancer == 'LIHC'] <- 'Liver (LIHC)'
merged.mapping.normal$tissue[merged.mapping.normal$cancer == 'BRCA'] <- 'Breast (BRCA)'
merged.mapping.normal$tissue[merged.mapping.normal$cancer == 'PRAD'] <- 'Prostate (PRAD)'
merged.mapping.normal$tissue[merged.mapping.normal$cancer == 'CHOL'] <- 'Bile duct (CHOL)'
pdf('~/gallery/2017_01_05_TCGA_allnormal_scatter.pdf', width=5, height=2.5)
ggplot(aes(commonHMD, commonPMD), data=merged.mapping.normal) + geom_point(aes(color=tissue), size=0.8, shape=21) + geom_abline(intercept=0, slope=1, linetype='dashed') + scale_color_manual(values=setNames(c('#FF0000','#924965','#1F78B6','#FF7F00','#33A02C','#999999'), c('Esophagus (ESCA)','Liver (LIHC)', 'Breast (BRCA)', 'Prostate (PRAD)', 'Bile duct (CHOL)', 'Others'))) + xlab('common HMD') + ylab('common PMD') + theme(text=element_text(size=15, face='bold')) + xlim(0.35,0.9) + ylim(0.35,0.9)
dev.off()

## tumor
cancer2mapping <- split(merged.mapping, merged.mapping$cancer)
for (cancer in names(cancer2mapping)) {
  mm <- cancer2mapping[[cancer]]
  mm <- cbind(mm, depths[mm$idat,])
  mm <- mm[substr(mm$barcode, 14,14) %in% c('0','1'),]
  mm$tumor <- ifelse(substr(mm$barcode, 14,14)=='0', 'Tumor', 'Adjacent Normal')
  pdf(sprintf('~/gallery/2017_01_05_TCGA_tumor_scatter_%s.pdf', cancer), width=3, height=2.5)
  print(ggplot(aes(commonHMD, commonPMD), data=mm) + geom_point(aes(color=tumor), fill=NA, size=0.6, alpha=0.8, shape=21) + geom_abline(intercept=0, slope=1, linetype='dashed') + xlab('common HMD') + ylab('common PMD') + theme(text=element_text(size=15, face='bold')) + xlim(0,0.9) + ylim(0,0.9) + scale_color_manual(values=setNames(c('#CC0000','#005AA0'), c('Adjacent Normal', 'Tumor')), guide=FALSE) + ggtitle(cancer))
  dev.off()
}

## mark subtypes
pancan.subtypes <- read.table('/Users/wandingzhou/projects/hs-tcga/data/2015_03_23_Hui_annotation/Tumor subtypes for PanCanPathways.txt', header=T)
load('~/projects/hs-tcga/data/2015_03_23_Hui_annotation/subtypes.rda')

## TGCT
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/allTCGA.rda')
samples <- cbind(merged.mapping.tumor, depths[merged.mapping.tumor$idat,])
samples.TGCT <- samples[(samples$cancer == 'TGCT'),]
a <- read.csv('/Users/wandingzhou/projects/hs-tcga/data/2015_03_23_Hui_annotation/TGCT.Histology.KH.csv', header=T, check.names=F)
samples.TGCT$subtype <- a$wzsubtype[match(substr(samples.TGCT$barcode,1,15), a[[2]])]
png('~/gallery/2017_03_08_TGCT_subtypes.png', width=1400, height=600, res=200)
ggplot(samples.TGCT) + geom_point(aes(commonHMD, commonPMD, color=subtype), size=0.5) + geom_abline(slope=1, intercept=0, lty='dashed') + xlab('Common HMD') + ylab('Common PMD') + ggtitle('TGCT')
dev.off()

## BRCA
samples.BRCA <- samples[samples$cancer == 'BRCA', ]
samples.BRCA <- cbind(samples.BRCA, subtypes[match(samples.BRCA$barcode, subtypes$id_TCGA),])
png('~/gallery/2017_03_08_BRCA_subtypes.png', width=900, height=600, res=200)
ggplot(samples.BRCA) + geom_point(aes(commonHMD, commonPMD, color=subtype), size=0.5) + geom_abline(slope=1, intercept=0, lty='dashed') + xlab('Common HMD') + ylab('Common PMD') + ggtitle('BRCA')
dev.off()
png('~/gallery/2017_03_08_BRCA_pam50.png', width=900, height=600, res=200)
ggplot(samples.BRCA) + geom_point(aes(commonHMD, commonPMD, color=pam50), size=0.5) + geom_abline(slope=1, intercept=0, lty='dashed') + xlab('Common HMD') + ylab('Common PMD') + ggtitle('BRCA')
dev.off()

## COAD
samples.COAD <- samples[samples$cancer == 'COAD', ]
samples.COAD$subtype <- pancan.subtypes$SUBTYPE[match(substr(samples.COAD$barcode,1,15), pancan.subtypes$SAMPLE_BARCODE)]
png('~/gallery/2017_03_08_COAD_subtypes.png', width=900, height=600, res=200)
ggplot(samples.COAD) + geom_point(aes(commonHMD, commonPMD, color=subtype), size=0.5) + geom_abline(slope=1, intercept=0, lty='dashed') + xlab('Common HMD') + ylab('Common PMD') + xlim(0.2,0.8) + ylim(0.2,0.8) + ggtitle('COAD')
dev.off()

## READ
samples.READ <- samples[samples$cancer == 'READ', ]
samples.READ$subtype <- pancan.subtypes$SUBTYPE[match(substr(samples.READ$barcode,1,15), pancan.subtypes$SAMPLE_BARCODE)]
png('~/gallery/2017_03_08_READ_subtypes.png', width=900, height=600, res=200)
ggplot(samples.READ) + geom_point(aes(commonHMD, commonPMD, color=subtype), size=0.5) + geom_abline(slope=1, intercept=0, lty='dashed') + xlab('Common HMD') + ylab('Common PMD') + xlim(0.2,0.85) + ylim(0.2,0.85) + ggtitle('READ')
dev.off()

## GBM
samples.GBM <- samples[samples$cancer == 'GBM', ]
samples.GBM <- cbind(samples.GBM, subtypes[match(samples.GBM$barcode, subtypes$id_TCGA),])
png('~/gallery/2017_03_08_GBM_subtypes.png', width=900, height=600, res=200)
ggplot(samples.GBM) + geom_point(aes(commonHMD, commonPMD, color=cancerlabel), size=0.5) + geom_abline(slope=1, intercept=0, lty='dashed') + xlab('Common HMD') + ylab('Common PMD') + xlim(0.2,0.85) + ylim(0.2,0.85) + ggtitle('GBM')
dev.off()

## UCEC
samples.UCEC <- samples[samples$cancer == 'UCEC', ]
samples.UCEC$subtype <- pancan.subtypes$SUBTYPE[match(substr(samples.UCEC$barcode,1,15), pancan.subtypes$SAMPLE_BARCODE)]
png('~/gallery/2017_03_08_UCEC_subtypes.png', width=900, height=600, res=200)
ggplot(samples.UCEC) + geom_point(aes(commonHMD, commonPMD, color=subtype), size=0.5) + geom_abline(slope=1, intercept=0, lty='dashed') + xlab('Common HMD') + ylab('Common PMD') + xlim(0.2,0.85) + ylim(0.2,0.85) + ggtitle('UCEC')
dev.off()

## CESC
samples.CESC <- samples[samples$cancer == 'CESC', ]
samples.CESC$subtype <- pancan.subtypes$SUBTYPE[match(substr(samples.CESC$barcode,1,15), pancan.subtypes$SAMPLE_BARCODE)]
samples.CESC <- samples.CESC[samples.CESC$subtype!="<NA>",]
png('~/gallery/2017_03_08_CESC_subtypes.png', width=1100, height=600, res=200)
ggplot(samples.CESC) + geom_point(aes(commonHMD, commonPMD, color=subtype), size=0.5) + geom_abline(slope=1, intercept=0, lty='dashed') + xlab('Common HMD') + ylab('Common PMD') + xlim(0.2,0.85) + ylim(0.2,0.85) + ggtitle('CESC')
dev.off()

## LGG
samples.LGG <- samples[samples$cancer == 'LGG', ]
samples.LGG$subtype <- pancan.subtypes$SUBTYPE[match(substr(samples.LGG$barcode,1,15), pancan.subtypes$SAMPLE_BARCODE)]
samples.LGG <- samples.LGG[samples.LGG$subtype!="<NA>",]
png('~/gallery/2017_03_08_LGG_subtypes.png', width=1100, height=600, res=200)
ggplot(samples.LGG) + geom_point(aes(commonHMD, commonPMD, color=subtype), size=0.1) + geom_abline(slope=1, intercept=0, lty='dashed') + xlab('Common HMD') + ylab('Common PMD') + xlim(0.45,0.85) + ylim(0.45,0.85) + ggtitle('LGG')
dev.off()

## STAD
samples.STAD <- samples[samples$cancer == 'STAD', ]
samples.STAD$subtype <- pancan.subtypes$SUBTYPE[match(substr(samples.STAD$barcode,1,15), pancan.subtypes$SAMPLE_BARCODE)]
samples.STAD <- samples.STAD[samples.STAD$subtype!="<NA>",]
png('~/gallery/2017_03_08_STAD_subtypes.png', width=900, height=600, res=200)
ggplot(samples.STAD) + geom_point(aes(commonHMD, commonPMD, color=subtype), size=0.8) + geom_abline(slope=1, intercept=0, lty='dashed') + xlab('Common HMD') + ylab('Common PMD') + xlim(0.2,0.85) + ylim(0.2,0.85) + ggtitle('STAD')
dev.off()

## landscape plot
merged.mapping <- read.table('/primary/projects/laird/projects/2016_04_05_TCGA_pancan_renormalization/450k/merged_mapping', sep='\t', stringsAsFactors = F, col.names = c('cancer','barcode','idat'))
merged.mapping.tumor <- merged.mapping[substr(merged.mapping$barcode, 14,14) == '0',]
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/allTCGA.rda')
samples <- cbind(merged.mapping.tumor, depths[merged.mapping.tumor$idat,])
samples$delta <- samples$commonPMD - samples$commonHMD
samples1 <- do.call(rbind, lapply(split(samples, samples$cancer), function(x) {
  xx <- x[order(x$delta),]
  xx$order <- seq_along(x$delta) / length(x$delta)
  xx$meandelta <- mean(x$delta, na.rm=T)
  xx
}))
deltamean <- tapply(samples$delta, samples$cancer, mean)
orderedtumors <- names(sort(deltamean))
samples1$cancer <- factor(samples1$cancer, level=orderedtumors)
samples1$deltacap <- pmax(-0.3, samples1$delta)
samples1$deltacap <- pmin(0, samples1$deltacap)
hlinedata <- data.frame(mean=deltamean, cancer=factor(names(deltamean), level=orderedtumors))
pdf('~/gallery/2017_02_26_landscapeplot_pmd_depth.pdf', width=17, height=2.7)
ggplot(samples1) + geom_point(aes(order, delta, color=deltacap), size=1.5) + geom_hline(aes(yintercept=mean), hlinedata) + geom_hline(yintercept=0, lty='dashed') + facet_wrap(~cancer, ncol=33) + scale_color_gradient2(high='#CC0000',low='#005AA0',mid='#AAAAAA',midpoint=-0.1, limits=c(-0.0, -0.3)) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), strip.background = element_blank(), panel.spacing=unit(0.1,'lines'), strip.text.x = element_text(size = 17, face='bold', colour = "black", angle = 90, vjust=0)) + ylab('PMD - HMD')
dev.off()

## landscape plot - normal
merged.mapping <- read.table('/primary/projects/laird/projects/2016_04_05_TCGA_pancan_renormalization/450k/merged_mapping', sep='\t', stringsAsFactors = F, col.names = c('cancer','barcode','idat'))
merged.mapping.normal <- merged.mapping[substr(merged.mapping$barcode, 14,14) == '1',]
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/allTCGA.rda')
samples <- cbind(merged.mapping.normal, depths[merged.mapping.normal$idat,])
samples$delta <- samples$commonPMD - samples$commonHMD
samples1 <- do.call(rbind, lapply(split(samples, samples$cancer), function(x) {
  xx <- x[order(x$delta),]
  xx$order <- seq_along(x$delta) / length(x$delta)
  xx$meandelta <- mean(x$delta, na.rm=T)
  xx
}))
deltamean <- tapply(samples$delta, samples$cancer, mean)
orderedtumors <- names(sort(deltamean))
samples1$cancer <- factor(samples1$cancer, level=orderedtumors)
samples1$deltacap <- pmax(-0.3, samples1$delta)
samples1$deltacap <- pmin(0, samples1$deltacap)
hlinedata <- data.frame(mean=deltamean, cancer=factor(names(deltamean), level=orderedtumors))
pdf('~/gallery/2017_02_26_landscapeplot_pmd_depth_normal.pdf', width=13, height=2.7)
ggplot(samples1) + geom_point(aes(order, delta, color=deltacap), size=1.5) + geom_hline(aes(yintercept=mean), hlinedata) + geom_hline(yintercept=0, lty='dashed') + facet_wrap(~cancer, ncol=33) + scale_color_gradient2(high='#CC0000',low='#005AA0',mid='#AAAAAA',midpoint=-0.1, limits=c(-0.0, -0.3)) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), strip.background = element_blank(), panel.spacing=unit(0.1,'lines'), strip.text.x = element_text(size = 17, face='bold', colour = "black", angle = 90, vjust=0)) + ylab('PMD - HMD')
dev.off()

## tumor vs normal by cancer types
samples <- cbind(merged.mapping.tumor, depths[merged.mapping.tumor$idat,])
samples$delta <- samples$commonPMD - samples$commonHMD
samples1 <- do.call(rbind, lapply(split(samples, samples$cancer), function(x) {
  xx <- x[order(x$delta),]
  xx$order <- seq_along(x$delta) / length(x$delta)
  xx$meandelta <- mean(x$delta, na.rm=T)
  xx
}))
deltameanT <- tapply(samples$delta, samples$cancer, mean)
sdT <- tapply(samples$delta, samples$cancer, sd)
merged.mapping.normal <- merged.mapping[substr(merged.mapping$barcode, 14,14) == '1',]
samples <- cbind(merged.mapping.normal, depths[merged.mapping.normal$idat,])
samples$delta <- samples$commonPMD - samples$commonHMD
samples1 <- do.call(rbind, lapply(split(samples, samples$cancer), function(x) {
  xx <- x[order(x$delta),]
  xx$order <- seq_along(x$delta) / length(x$delta)
  xx$meandelta <- mean(x$delta, na.rm=T)
  xx
}))
deltameanN <- tapply(samples$delta, samples$cancer, mean)
sdN <- tapply(samples$delta, samples$cancer, sd)
df <- data.frame(tumor = deltameanT[names(deltameanN)], normal = deltameanN, tumosd = sdT[names(deltameanN)], normsd = sdN[names(deltameanN)], cancer=names(deltameanN))
pdf('~/gallery/2017_03_30_normal_tumor_PMD_scatter.pdf', width=7, height=5)
ggplot(df, aes(normal, tumor)) + geom_point(aes(color=cancer), size=2) + geom_errorbar(aes(ymin = tumor-tumosd, ymax = tumor+tumosd, color=cancer)) + geom_errorbarh(aes(xmin = normal-normsd, xmax = normal+normsd, color=cancer))
dev.off()
pdf('~/gallery/2017_03_30_normal_tumor_PMD_scatter2.pdf', width=7, height=5)
ggplot(df, aes(normal, tumor)) + geom_point(aes(color=cancer), size=2) + geom_errorbar(aes(ymin = tumor-tumosd, ymax = tumor+tumosd, color=cancer)) + geom_errorbarh(aes(xmin = normal-normsd, xmax = normal+normsd, color=cancer)) + xlim(-0.35,0.02) + ylim(-0.35,0.02)
dev.off()

## some specific cancer type
samples <- merged.mapping
samples <- cbind(samples, depths[samples$idat,])
## KIRP
load('~/projects/hs-tcga/data/2015_03_23_Hui_annotation/subtypes.rda')
subtypes$subtype[match(samples$barcode, subtypes$id_TCGA)]

## cell lines
mmc <- merged.mapping[substr(merged.mapping$barcode, 14,14) == '2',]
mmc <- cbind(mmc, depths[mmc$idat,])
mmc$batch <- substr(mmc$barcode,9,12)
pdf('~/gallery/2017_01_05_TCGA_cellline_scatter.pdf', width=4, height=3)
ggplot(aes(commonHMD, commonPMD), data=mmc) + geom_point(aes(color=batch), fill=NA, size=0.6, alpha=0.8, shape=21) + geom_abline(intercept=0, slope=1, linetype='dashed') + xlab('common HMD') + ylab('common PMD') + theme(text=element_text(size=15, face='bold')) + xlim(0.25,0.65) + ylim(0.25,0.65) + ggtitle('Cell Lines')
dev.off()

pdf('~/gallery/2017_01_05_TCGA_cellline_scatter_zoom.pdf', width=4, height=3)
ggplot(aes(commonHMD, commonPMD), data=mmc) + geom_point(aes(color=batch), fill=NA, size=0.6, alpha=0.8, shape=21) + geom_abline(intercept=0, slope=1, linetype='dashed') + xlab('common HMD') + ylab('common PMD') + theme(text=element_text(size=15, face='bold')) + ggtitle('Cell Lines')
dev.off()

## batch information
load('/Users/wandingzhou/projects/hs-tcga/data/2016_06_20_TCGA_xls/HM450_packaging/batch.rda')
mmc$batch2 <- batch[mmc$idat]
mmcn <- mmc[!is.na(mmc$batch2),]
pdf('~/gallery/2017_01_05_TCGA_cellline_scatter_zoom_batch.pdf', width=6, height=2)
ggplot(mmcn) + geom_point(aes(batch2, commonPMD, color=batch), size=0.8) + geom_line(aes(batch2, commonPMD, color=batch), alpha=0.7) + xlab('Batch') + ylab('common PMD')
dev.off()

## mutation rates
load('/primary/projects/laird/projects/2016_04_05_TCGA_pancan_renormalization/450k/merged_mapping.rda')
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/allTCGA.rda')
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/mutationrate/mutrate.rda')
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/purity/purity.rda')
samples <- merged.mapping.tumor
samples2 <- cbind(samples, mutrate[match(substr(samples$barcode,1,20), substr(mutrate$name,1,20)),])
samples3 <- cbind(samples2, depths[samples2$idat,])
samples4 <- cbind(samples3,purity[substr(samples3$barcode, 1, 15),])
## samples4$purity <- samples4$purity
samples4$correctedPMD <- pmax(0, (samples4$commonPMD + samples4$commonHMD * samples4$purity - samples4$commonHMD) / samples4$purity)
samples <- samples4
save(samples, file='/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/TCGAtumorsamples.rda')

## purity vs corrected and uncorrected pmd
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/TCGAtumorsamples.rda')
pdf('~/gallery/2017_01_15_purity_vs_pmd.pdf', width=4, height=4)
par(mar=c(5,5,1,1))
smoothScatter(samples$purity, samples$commonPMD, nrpoints=0, colramp=colorRampPalette(c('white','black')), xlab='Purity', ylab='PMD methylation', font.lab=2, font.axis=2, cex.axis=1.3, cex.lab=1.5)
dev.off()

pdf('~/gallery/2017_01_15_purity_vs_pmd_corr.pdf', width=4, height=4)
par(mar=c(5,5,1,1))
smoothScatter(samples$purity, samples$correctedPMD, nrpoints=0, colramp=colorRampPalette(c('white','black')), xlab='Purity', ylab='PMD methylation', font.lab=2, font.axis=2, cex.axis=1.3, cex.lab=1.5)
dev.off()

########################
## mutation rate
########################

## mutation rate vs purity corrected pmd
## load('/Users/wandingzhou/projects/hs-tcga/data/2015_03_23_Hui_annotation/sampleAnnotSubWB20130619.rda')
pdf('~/gallery/2017_01_15_mutationrate_vs_pmd_corr.pdf', width=4.5, height=4)
par(mar=c(5,5,1,3))
smoothScatter(samples$correctedPMD, log10(samples$rate), nrpoints=0, colramp=colorRampPalette(c('white','black')), xlab='Purity-corrected PMD methylation', ylab='Mutation Rate', ylim=c(-7,-4), yaxt='n', font.lab=2, font.axis=2, cex.axis=1.3, cex.lab=1.5)
axis(2, -7:-4, sapply(-7:-4, function(x) parse(text=paste0('10^',x))), font.axis=2, cex.axis=1.3)
dev.off()

pdf('~/gallery/2017_01_15_mutationrate_vs_pmd_corr_scatter.pdf', width=4.5, height=4)
par(mar=c(5,5,1,3))
plot(samples$correctedPMD, log10(samples$rate), col=ifelse(substr(samples$barcode,1,12) %in% rownames(msi[msi$msi=='MSI-H',]), 'red', 'black') , xlab='Purity-corrected PMD methylation', ylab='Mutation Rate', ylim=c(-7,-4), yaxt='n', font.lab=2, font.axis=2, cex.axis=1.3, cex.lab=1.5, pch=20, cex=0.3)
axis(2, -7:-4, sapply(-7:-4, function(x) parse(text=paste0('10^',x))), font.axis=2, cex.axis=1.3)
legend('bottomleft', c('Normal','MSI'), c('black','red'))
dev.off()

## contrast MSI cases
load('/primary/projects/laird/projects/2016_04_05_TCGA_pancan_renormalization/450k/msi/msi.rda')
pdf('~/gallery/2017_01_15_mutationrate_vs_pmd_corr_scatter_msi.pdf', width=4.8, height=4.5)
par(mar=c(5,5,3,4))
samples.msi <- samples[substr(samples$barcode, 1, 12) %in% rownames(msi),]
plot(samples.msi$correctedPMD, log10(samples.msi$rate), col=ifelse(substr(samples.msi$barcode,1,12) %in% rownames(msi[msi$msi=='MSI-H',]), 'red', 'black') , xlab='Purity-corrected PMD methylation', ylab='Mutation Rate', ylim=c(-7,-3.5), yaxt='n', font.lab=2, font.axis=2, cex.axis=1.3, cex.lab=1.5, pch=20, cex=0.3, main='Restricted to 18 tumor types with MSI annot.')
axis(2, -7:-4, sapply(-7:-4, function(x) parse(text=paste0('10^',x))), font.axis=2, cex.axis=1.3)
legend('bottomleft', c('Normal','MSI'), col=c('black','red'), pch=20)
dev.off()

## contrast POLE and APOBEC cases
load('/Users/wandingzhou/projects/hs-tcga/data/2015_06_17_TCGA_mutations/mutations.rda')
samples.mut <- samples[substr(samples$barcode, 1, 12) %in% mutations$patient,]
nonsilent.pole <- nonsilent[nonsilent$gene == 'POLE',]
pdf('~/gallery/2017_01_15_mutationrate_vs_pmd_corr_scatter_POLE.pdf', width=4.8, height=4.5)
par(mar=c(5,5,3,4))
plot(samples.mut$correctedPMD, log10(samples.mut$rate), col=ifelse(substr(samples.mut$barcode,1,12) %in% nonsilent.pole$patient, 'red', 'black') , xlab='Purity-corrected PMD methylation', ylab='Mutation Rate', ylim=c(-7,-3.5), yaxt='n', font.lab=2, font.axis=2, cex.axis=1.3, cex.lab=1.5, pch=20, cex=0.3, main='Restricted to samples with mutation annot.')
axis(2, -7:-4, sapply(-7:-4, function(x) parse(text=paste0('10^',x))), font.axis=2, cex.axis=1.3)
legend('bottomleft', c('Normal','POLE mutation'), col=c('black','red'), pch=20)
dev.off()

nonsilent.apobec <- nonsilent[grep('APOBEC',nonsilent$gene),]
pdf('~/gallery/2017_01_15_mutationrate_vs_pmd_corr_scatter_APOBEC.pdf', width=4.8, height=4.5)
par(mar=c(5,5,3,4))
plot(samples.mut$correctedPMD, log10(samples.mut$rate), col=ifelse(substr(samples.mut$barcode,1,12) %in% nonsilent.apobec$patient, 'red', 'black') , xlab='Purity-corrected PMD methylation', ylab='Mutation Rate', ylim=c(-7,-3.5), yaxt='n', font.lab=2, font.axis=2, cex.axis=1.3, cex.lab=1.5, pch=20, cex=0.3, main='Restricted to samples with mutation annot.')
axis(2, -7:-4, sapply(-7:-4, function(x) parse(text=paste0('10^',x))), font.axis=2, cex.axis=1.3)
legend('bottomleft', c('Normal','APOBEC mutation'), col=c('black','red'), pch=20)
dev.off()

## filter purity, masked POLE, APOBEC and MSI, plot mutation rate vs PMD (raw)
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/TCGAtumorsamples.rda')
load('/Users/wandingzhou/projects/hs-tcga/data/2015_06_17_TCGA_mutations/mutations.rda')
load('/primary/projects/laird/projects/2016_04_05_TCGA_pancan_renormalization/450k/msi/msi.rda')
nonsilent.pole <- nonsilent[nonsilent$gene == 'POLE',]
nonsilent.apobec <- nonsilent[grep('APOBEC',nonsilent$gene),]
samples$msi <- msi[substr(samples$barcode, 1, 12),'msi']
samples$pole <- substr(samples$barcode,1,12) %in% nonsilent.pole$patient
samples$apobec <- substr(samples$barcode,1,12) %in% nonsilent.apobec$patient

## new version
s1 <- samples[samples$purity >= 0.6,]
s1 <- s1[(is.na(s1$msi) | s1$msi == 'MSI-H') & !s1$pole & !s1$apobec,]
cor.test(s1$commonPMD, s1$rate, use='complete', method='spearman')
pdf('~/gallery/2017_01_18_mutationrate_vs_pmd_scatter_highpurity.pdf', width=4.8, height=4.2)
par(mar=c(5,5,3,4))
plot(s1$commonPMD, log10(s1$rate), xlab='PMD methylation', ylab='Mutation Rate', ylim=c(-7,-3.5), yaxt='n', font.lab=2, font.axis=2, cex.axis=1.3, cex.lab=1.5, pch=20, cex=0.3)
axis(2, -7:-4, sapply(-7:-4, function(x) parse(text=paste0('10^',x))), font.axis=2, cex.axis=1.3)
dev.off()

## newer vesion
s1 <- samples[samples$purity >= 0.7,]
cor.test(s1$commonPMD, s1$rate, use='complete', method='spearman')
pdf('~/gallery/2017_02_26_mutationrate_vs_pmd_scatter_highpurity.pdf', width=4.8, height=4.2)
par(mar=c(5,5,3,4))
ggplot(s1, aes(commonPMD, log10(rate))) + geom_point(size=0.5) + stat_smooth(method='lm') + xlab('PMD methylation') + ylab('Mutation Rate') + ylim(c(-7,-3.5))
dev.off()
s2 <- s1[(is.na(s1$msi) | s1$msi == 'MSI-H') & !s1$pole & !s1$apobec,]
cor.test(s2$commonPMD, s2$rate, use='complete', method='spearman')
pdf('~/gallery/2017_02_26_mutationrate_vs_pmd_scatter_highpurity_nomsi.pdf', width=4.8, height=4.2)
ggplot(s2, aes(commonPMD, log10(rate))) + geom_point(size=0.5) + stat_smooth(method='lm') + xlab('PMD methylation') + ylab('Mutation Rate') + ylim(c(-7,-3.5))
dev.off()

thresholds <- seq(0.1, 0.95, by=0.02)
allpurity.cors <- lapply(thresholds, function(thres) {
  s1 <- samples
  s1$msi <- msi[substr(s1$barcode, 1, 12),'msi']
  s1$pole <- substr(s1$barcode,1,12) %in% nonsilent.pole$patient
  s1$apobec <- substr(s1$barcode,1,12) %in% nonsilent.apobec$patient
  s1 <- s1[s1$purity >= thres,]
  s1 <- s1[(is.na(s1$msi) | s1$msi == 'MSI-H') & !s1$pole & !s1$apobec,]
  list(cor=cor.test(s1$commonPMD, s1$rate, use='complete', method='spearman'), n=nrow(s1))
})
pdf('~/gallery/2017_02_13_pmd_mutationrate_correlation_by_purity.pdf', width=3, height=4)
par(mfrow=c(3,1), mar=c(2,5,1,1), lwd=2, oma=c(2,2,1,1), cex.lab=1.5, cex.axis=1.5)
plot(thresholds, sapply(allpurity.cors, function(x) x$cor$estimate), ylim=c(-1,0), xlab='Threshold', ylab="Spearman's Rho", type='o')
plot(thresholds, sapply(allpurity.cors, function(x) log10(x$cor$p.value)), xlab='Threshold', ylab="Log10 P-Value", type='o')
plot(thresholds, sapply(allpurity.cors, function(x) x$n), xlab='Threshold', ylab="Sample Size", type='o')
dev.off()

## by cancer type
cancer2samples <- split(samples, samples$cancer)
for (cancer in names(cancer2samples)) {
  s1 <- cancer2samples[[cancer]]
  pdf(sprintf('~/gallery/2017_01_15_mutationrate_vs_pmd_corr_%s.pdf', cancer), width=4, height=3.7)
  s1$rate <- pmax(s1$rate, 10^-7.5)
  print(ggplot(s1, aes(correctedPMD, rate)) + geom_point() + scale_y_log10(breaks=10^(-8:-4), labels=sapply(-8:-4, function(x) parse(text=paste0('10^',x))), limits=c(1e-8,1e-4)) + xlab('Purity-corrected\nPMD methylation') + ylab('Mutation Rate') + xlim(0,1) + ggtitle(cancer))
  dev.off()
}

###########################
## mutation rate vs pmd
###########################
## threshold by purity >=0.6, plot raw data, and exclude MSI, APOBEC and POLE
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/TCGAtumorsamples.rda')
load('/Users/wandingzhou/projects/hs-tcga/data/2015_06_17_TCGA_mutations/mutations.rda')
load('/primary/projects/laird/projects/2016_04_05_TCGA_pancan_renormalization/450k/msi/msi.rda')
nonsilent.pole <- nonsilent[nonsilent$gene == 'POLE',]
nonsilent.apobec <- nonsilent[grep('APOBEC',nonsilent$gene),]
samples$msi <- msi[substr(samples$barcode, 1, 12),'msi']
samples$pole <- substr(samples$barcode,1,12) %in% nonsilent.pole$patient
samples$apobec <- substr(samples$barcode,1,12) %in% nonsilent.apobec$patient
cancer2samples <- split(samples, samples$cancer)
for (cancer in names(cancer2samples)) {
  s1 <- cancer2samples[[cancer]]
  s1 <- s1[s1$purity >= 0.6,]
  s1 <- s1[(is.na(s1$msi) | s1$msi == 'MSI-H') & !s1$pole & !s1$apobec,]
  pdf(sprintf('~/gallery/2017_01_17_mutationrate_vs_pmd_%s.pdf', cancer), width=4, height=3.7)
  s1$rate <- pmax(s1$rate, 10^-7.5)
  print(ggplot(s1, aes(commonPMD, rate)) + geom_point() + scale_y_log10(breaks=10^(-8:-4), labels=sapply(-8:-4, function(x) parse(text=paste0('10^',x))), limits=c(1e-8,1e-4)) + xlab('PMD methylation') + ylab('Mutation Rate') + xlim(0,1) + ggtitle(cancer))
  dev.off()
}

## plot pmd vs purity
merged.mapping <- read.table('/primary/projects/laird/projects/2016_04_05_TCGA_pancan_renormalization/450k/merged_mapping', sep='\t', stringsAsFactors = F, col.names = c('cancer','barcode','idat'))
merged.mapping.normal <- merged.mapping[substr(merged.mapping$barcode, 14,14) == '1',]
merged.mapping.tumor <- merged.mapping[substr(merged.mapping$barcode, 14,14) == '0',]
merged.mapping.cellline <- merged.mapping[substr(merged.mapping$barcode, 14,14) == '2',]
save(merged.mapping, merged.mapping.normal, merged.mapping.tumor, merged.mapping.cellline, file='/primary/projects/laird/projects/2016_04_05_TCGA_pancan_renormalization/450k/merged_mapping.rda')

load('/primary/projects/laird/projects/2016_04_05_TCGA_pancan_renormalization/450k/merged_mapping.rda')
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/allTCGA.rda')
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/mutationrate/mutrate.rda')
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/purity/purity.rda')
samples <- merged.mapping
samples2 <- cbind(samples, rate=mutrate$rate[match(substr(samples$barcode,1,20), substr(mutrate$name,1,20))])
samples3 <- cbind(samples2, unname(depths[match(samples2$idat, rownames(depths)),]))
colnames(samples3)[c(5,6)] <- c('pmd','hmd')
samples4 <- cbind(samples3, purity=purity$purity[match(substr(samples3$barcode, 1, 15), rownames(purity))], ploidy=purity$ploidy[match(substr(samples3$barcode, 1, 15), rownames(purity))])
samples <- samples4
samples$disease <- ifelse(substr(samples$barcode, 14, 14) == '0', 'tumor', 'normal')
samples$disease[substr(samples$barcode, 14,14) == '2'] = 'cellline'
samples$purity[samples$disease=='normal'] <- 1 # purity of normals assumed to be 1.0
samples$pmdc <- pmax(0, (samples$pmd + samples$hmd * samples$purity - samples$hmd) / samples$purity)
load('/primary/projects/laird/projects/2016_04_05_TCGA_pancan_renormalization/450k/age.rda')
samples$age <- agetable$age_year[match(substr(samples$barcode,1,12), agetable$barcode0)] # add age
save(samples, file='/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/TCGAsamples.rda')

## hmd vs pmd
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/TCGAsamples.rda')
cancer2mapping <- split(samples, samples$cancer)
for (cancer in names(cancer2mapping)) {
  mm <- cancer2mapping[[cancer]]
  mm <- cbind(mm, depths[mm$idat,])
  mm <- mm[substr(mm$barcode, 14,14) %in% c('0','1'),]
  mm$tumor <- ifelse(substr(mm$barcode, 14,14)=='0', 'Tumor', 'Adjacent Normal')
  pdf(sprintf('~/gallery/2017_01_16_TCGA_tumor_scatter_%s_purity_corrected.pdf', cancer), width=3, height=2.5)
  print(ggplot(aes(hmd, pmdc), data=mm) + geom_point(aes(color=tumor), fill=NA, size=0.6, alpha=0.8, shape=21) + geom_abline(intercept=0, slope=1, linetype='dashed') + xlab('HMD methylation') + ylab('PMD methylation') + theme(text=element_text(size=15, face='bold')) + xlim(0,0.9) + ylim(0,0.9) + scale_color_manual(values=setNames(c('#CC0000','#005AA0'), c('Adjacent Normal', 'Tumor')), guide=FALSE) + ggtitle(cancer))
  dev.off()
}

## pmd vs age
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/TCGAsamples.rda')
cancer2mapping <- split(samples, samples$cancer)
for (cancer in names(cancer2mapping)) {
  mm <- cancer2mapping[[cancer]]
  mm <- cbind(mm, depths[mm$idat,])
  mm <- mm[substr(mm$barcode, 14,14) %in% c('0','1'),]
  mm$tumor <- ifelse(substr(mm$barcode, 14,14)=='0', 'Tumor', 'Adjacent Normal')
  pdf(sprintf('~/gallery/2017_01_16_TCGA_age_vs_pmd_purity_corrected_%s.pdf', cancer), width=3, height=2.5)
  print(ggplot(aes(age, pmdc), data=mm) + geom_point(aes(color=tumor), fill=NA, size=0.6, alpha=0.8, shape=21) + xlab('Age') + ylab('PMD methylation') + theme(text=element_text(size=15, face='bold')) + ylim(0,0.9) + scale_color_manual(values=setNames(c('#CC0000','#005AA0'), c('Adjacent Normal', 'Tumor')), guide=FALSE) + ggtitle(cancer))
  dev.off()
}

## TCGA normal age vs pmd
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/TCGAsamples.rda')
pdf('~/gallery/2017_03_01_TCGAnormal_age_vs_pmd.pdf', width=10, height=10)
ggplot(samples[samples$disease=='normal',]) + geom_point(aes(age, pmd-hmd)) + facet_wrap(~cancer)
dev.off()

bycancer <- split(samples, samples$cancer)
for(cancer in names(bycancer)) {
  df <- bycancer[[cancer]]
  df <- df[df$disease == 'normal',]
  if (sum(!is.na(df$age)) > 10) {
    cat('\n\n',cancer,'\n')
    print(cor.test(df$pmd, df$age))
  }
}

pdf('~/gallery/2017_03_01_TCGAtumor_age_vs_pmd.pdf', width=10, height=10)
ggplot(samples[samples$disease=='tumor',]) + geom_point(aes(age, pmd-hmd), size=0.5) + facet_wrap(~cancer)
dev.off()

###############################
## Hannum age related probes
###############################
negprobes <- c('cg05442902', 'cg22285878', 'cg09651136', 'cg20822990', 'cg06685111', 'cg20052760', 'cg02867102', 'cg16054275', 'cg00486113', 'cg22796704', 'cg02046143', 'cg04474832', 'cg08415592', 'cg10501210', 'cg13001142', 'cg19722847', 'cg04875128', 'cg06874016', 'cg19283806', 'cg14556683', 'cg03473532', 'cg08234504', 'cg01528542', 'cg00481951', 'cg22158769', 'cg25428494', 'cg16419235', 'cg07927379', 'cg09809672', 'cg23091758')
## 'ch.13.39564907R',
## 'ch.2.30415474F', 
posprobes <- c('cg23744638', 'cg02085953', 'cg22512670', 'cg22016779', 'cg24079702', 'cg07082267', 'cg07583137', 'cg07547549', 'cg07553761', 'cg25410668', 'cg25478614', 'cg22736354', 'cg22454769', 'cg23500537', 'cg00748589', 'cg23606718', 'cg21296230', 'cg03032497', 'cg18473521', 'cg06639320', 'cg08540945', 'cg06493994', 'cg04400972', 'cg02650266', 'cg03607117', 'cg14361627', 'cg16867657', 'cg04940570', 'cg04416734', 'cg19935065', 'cg06419846', 'cg07955995', 'cg11067179', 'cg21139312', 'cg20426994', 'cg14692377', 'cg22213242', 'cg08097417', 'cg03399905')
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450probes.rda')
a <- probes[negprobes]
a$orph35class <- cut(a$orphan35, breaks=c(0,1,2,3,Inf), right=F)
a <- probes[posprobes]
a$orph35class <- cut(a$orphan35, breaks=c(0,1,2,3,Inf), right=F)

##############################
## AML subtypes, PML-RARA
##############################
aml.matrix <- read.table('/Volumes/projects_secondary/laird/projects/2016_12_26_TCGA_WGBS/AML_matrix_0304_2hit.txt', sep='\t', header=T, check.names=F, row.names=1)
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/TCGAsamples.rda')
samples$delta <- samples$pmd-samples$hmd
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450probes.rda')
samples <- samples[order(samples$delta),]
sampleindex <- match(substr(samples$barcode,1,12), colnames(aml.matrix))
samples <- samples[!is.na(sampleindex),]
aml.matrix <- aml.matrix[,na.omit(sampleindex)]
idats <- sprintf('/primary/projects/laird/projects/2016_04_05_TCGA_pancan_renormalization/450k/betas/%s.rda', samples$idat)
amlbetas <- do.call(cbind, mclapply(idats, function(idat) {
  load(idat);
  betas[com.pmd.probes]
}, mc.cores=7))
colnames(amlbetas) <- samples$barcode
## remove NAs
amlbetas <- amlbetas[apply(amlbetas,1,function(x) sum(is.na(x)) < 10), ]
colnames(aml.matrix) <- colnames(amlbetas)
a <- read.csv('/Volumes/projects_secondary/laird/projects/2016_12_26_TCGA_WGBS/AML/aml200.csv')
WBC <- setNames(a$WBC, a$X)
aa <- both.cluster(amlbetas)$mat
png('~/gallery/2017_03_01_aml_pmd.png', width=1000, height=1000)
WHeatmap(aa) + WHeatmap(aml.matrix[,colnames(aa)], TopOf(height=1), yticklabels=T) + WColorBarH(log10(WBC[substr(colnames(aa),1,12)]), TopOf(height=0.1))
dev.off()

#############################
## LOLIPOP peripheral blood
#############################

setwd('/secondary/projects/jones/projects/2016_12_22_LOLIPOP')
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450probes.rda')
allfns <- list.files('betas')
depths <- mclapply(allfns, function(fn) {
  load(paste0('betas/', fn))
  c(mean(betas[com.pmd.probes], na.rm=T), mean(betas[com.hmd.probes], na.rm=T))
}, mc.cores=20)

depths <- t(simplify2array(depths))
rownames(depths) <- gsub('.rda','',allfns)
colnames(depths) <- c('commonPMD','commonHMD')
load('samplesheet.rda')
save(depths, samplesheet, file='/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/LOLIPOP_peripheral_blood_depth.rda')

## plot depth vs age
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/LOLIPOP_peripheral_blood_depth.rda')
mm <- cbind(samplesheet[match(rownames(depths), samplesheet$sampleID),],depths)
pdf('~/gallery/2017_01_06_LOLIPOP_age.pdf', width=6, height=4)
ggplot(mm, aes(age, commonPMD)) + geom_point(aes(color=gender), size=0.5) + geom_smooth() + xlab('Age') + ylab('common PMD depth') + ggtitle('LOLIPOP (peripheral blood)') + theme(text=element_text(size=15, face='bold'))
dev.off()

################################
## clean all TCGA mutations
################################
setwd('~/projects/hs-tcga/data/2015_06_17_TCGA_mutations')
mutations <- read.table('merged_maf', sep='\t', quote='', stringsAsFactors=FALSE)
colnames(mutations) <- c('cancer', 'gene', 'outcome', 'type', 'barcode_tumor', 'barcode_normal', 'genomicchange', 'transcript', 'proteinchange', 'source', 'chrm', 'beg','end', 'ref', 'alt')
mutations$patient <- substr(mutations$barcode_tumor,1,12)
nonsilent <- mutations[!(mutations$outcome %in% c('Silent',"3'UTR", "5'Flank", "5'UTR")),]
save(mutations, nonsilent, file='/Users/wandingzhou/projects/hs-tcga/data/2015_06_17_TCGA_mutations/mutations.rda')

################################
## passage number vs pmd depths
################################
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450probes.rda')
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/GSE31848/betas.rda')

a1 <- melt(betas[com.pmd.probes,])
a1$pmd <- 'pmd'
a2 <- melt(betas[com.hmd.probes,])
a2$pmd <- 'hmd'
a <- rbind(a1, a2)
a$passage <- samples[a$Var2, 'passagenumber']
a <- a[!is.na(a$passage) & !is.na(a$value),]
ggplot(data=a) + geom_violin(aes(x=as.factor(passage), y=value, fill=pmd))

load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/GSE37066/betas.rda')
samples$pmd <- apply(betas[intersect(com.pmd.probes, rownames(betas)),], 2, mean, na.rm=T)[rownames(samples)]
samples$hmd <- apply(betas[intersect(com.hmd.probes, rownames(betas)),], 2, mean, na.rm=T)[rownames(samples)]
pdf('~/gallery/2017_01_20_passage_number_MSC.pdf', width=5, height=3.5)
ggplot(samples, aes(passagenumber, pmd/hmd, color=irradiation)) + geom_point() + geom_smooth(method='lm') + ggtitle('MSC (GSE37066)') + ylab('PMD Me. / HMD Me.') + xlab('Passage Number')
dev.off()
pdf('~/gallery/2017_01_20_passage_number_MSC_pmd.pdf', width=5, height=3.5)
ggplot(samples, aes(passagenumber, pmd, color=irradiation)) + geom_point() + geom_smooth(method='lm') + ggtitle('MSC (GSE37066)') + ylab('PMD Methylation') + xlab('Passage Number')
dev.off()


load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/GSE40699/betas.rda')
a1 <- melt(betas[com.pmd.probes,])
a1$pmd <- 'pmd'
a2 <- melt(betas[com.hmd.probes,])
a2$pmd <- 'hmd'
a <- rbind(a1, a2)

makepmddf <- function(betas, samples) {
  a1 <- melt(betas[com.pmd.probes,])
  a1$pmd <- 'pmd'
  a2 <- melt(betas[com.hmd.probes,])
  a2$pmd <- 'hmd'
  a <- rbind(a1, a2)
  a3 <- do.call(cbind, lapply(colnames(samples), function(cn) {
    samples[a$Var2, cn]
  }))
  colnames(a3) <- colnames(samples)
  cbind(a, a3)
}

a$passage <- samples[a$Var2, 'passagenumber']
a$age <- samples[a$Var2, 'age']
a <- a[!is.na(a$passage) & !is.na(a$value),]

################################
## PMD vs age in fibroblast
################################

load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450probes.rda')
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/GSE77135/betas.rda')
samples$pmd <- apply(betas[com.pmd.probes,], 2, mean, na.rm=T)[rownames(samples)]
samples$hmd <- apply(betas[com.hmd.probes,], 2, mean, na.rm=T)[rownames(samples)]
pdf('~/gallery/2017_01_13_fibroblast_age.pdf', width=7, height=3)
ggplot(samples, aes(age, pmd, color=sourceName)) + geom_jitter() + xlab('Age') + ylab('PMD depth') + ggtitle('GSE77135')
dev.off()

pdf('~/gallery/2017_01_13_fibroblast_age_pmd_vs_hmd.pdf', width=4, height=3)
ggplot(samples, aes(hmd, pmd, color=age)) + geom_point() + geom_abline(slope=1,intercept=0,linetype='dashed') + xlim(0.45,0.7) + ylim(0.45,0.7) + ggtitle('GSE77135')
dev.off()

## ICF dataset
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/GSE68747/betas.rda')
samples$pmd <- apply(betas[com.pmd.probes,], 2, mean, na.rm=T)[rownames(samples)]
samples$hmd <- apply(betas[com.hmd.probes,], 2, mean, na.rm=T)[rownames(samples)]

pdf('~/gallery/2017_01_13_fibroblast_ICF.pdf', width=4, height=3)
ggplot(samples, aes(title, pmd, fill=diseasestate)) + geom_bar(stat='identity') + xlab('Samples') + ylab('PMD depth') + ggtitle('GSE68747') + theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
dev.off()

## Wagner fibroblast
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/GSE52025_Wagner_fibroblast/betas.rda')
samples$pmd <- apply(betas[com.pmd.probes,], 2, mean, na.rm=T)[rownames(samples)]
samples$hmd <- apply(betas[com.hmd.probes,], 2, mean, na.rm=T)[rownames(samples)]
pdf('~/gallery/2017_01_13_fibroblast_Wagner_pmd_vs_hmd.pdf', width=4,height=3)
ggplot(samples, aes(hmd, pmd, color=sex)) + geom_point() + geom_abline(slope=1, intercept=0, linetype='dashed') + xlim(0.3,0.75) + ylim(0.3,0.75) + ggtitle("GSE52025 (Wagner et al.)")
dev.off()

## breast fibroblast
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/GSE74877/betas.rda')
samples.melano <- samples[grep('epidermal.melanocyte', samples$title),]
samples.melano$pigmentation <- sapply(strsplit(samples.melano$title,'[\\.-]'), function(x) x[4])
df <- makepmddf(betas[,rownames(samples.melano)], samples.melano)
pdf('~/gallery/2017_01_13_melanocytes.pdf', width=8, height=3)
ggplot(df, aes(interaction(geo,factor(pigmentation, levels=c('light','medium','dark'))), value, fill=pmd)) + geom_violin(lwd=0.5, bw=0.05, draw_quantiles=0.5, alpha=0.7) + scale_fill_manual(values=c('#CC0000','#005AA0')) + scale_x_discrete(labels=c('Light','Light','Medium','Medium','Dark','Dark')) + xlab('Pigmentation') + ylab('DNA methylation') + ggtitle('GSE74877')
dev.off()

## cell line differentiation experiment
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/GSE85828/betas.rda')
samples$pmd <- apply(betas[com.pmd.probes,], 2, mean, na.rm=T)[rownames(samples)]
samples$hmd <- apply(betas[com.hmd.probes,], 2, mean, na.rm=T)[rownames(samples)]
pdf('~/gallery/2017_01_13_cellline_differentiation.pdf', width=5,height=3)
ggplot(samples, aes(hmd, pmd, color=differentiationstage)) + geom_point() + geom_abline(slope=1, intercept=0, linetype='dashed') + xlim(0.3,0.8) + ylim(0.3,0.8) + ggtitle("GSE85828 (PCBC)")
dev.off()

#############
## subtype
#############
subtypes <- read.table('/Users/wandingzhou/projects/hs-tcga/data/2015_03_23_Hui_annotation/subtypes.txt', header=T, sep='\t')
save(subtypes, file='/Users/wandingzhou/projects/hs-tcga/data/2015_03_23_Hui_annotation/subtypes.rda')

########################
## mutation analysis
########################

load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/mutations.rda')
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/allTCGA.rda')
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/purity/purity.rda')
load('/primary/projects/laird/projects/2016_04_05_TCGA_pancan_renormalization/450k/merged_mapping.rda')
## focus on high purity samples
tumor <- cbind(merged.mapping.tumor, purity=purity[match(substr(merged.mapping.tumor$barcode,1,15), purity$array),'purity'])
tumor <- cbind(tumor, depths[tumor$idat,])
genes <- unique(mutations$gene)
length(genes)

allgenedepth <- mclapply(genes, function(gene) {
  tumor[substr(tumor$barcode,1,12) %in% unique(mutations[mutations$gene == gene, 'patient']),]
}, mc.cores=28)
names(allgenedepth) <- genes
allgenen <- sapply(allgenedepth, function(x) dim(x)[1])
allgenedepth <- allgenedepth[allgenen > 0]
allgenedepth <- allgenedepth[names(allgenedepth) != 'Unknown']
names(allgenen) <- genes
allgenen <- sort(allgenen, decreasing=T)
save(allgenedepth, allgenen, tumor, file='/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/allmutationpmd.rda')

load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/allmutationpmd.rda')
allgenePMDmean <- sapply(allgenedepth, function(x) mean(x$commonPMD, na.rm=T))
a <- data.frame(PMD=allgenePMDmean, N=allgenen[names(allgenePMDmean)])
a <- a[order(a$PMD),]
allgenePMDmean <- a

genes <- c('NSD1', 'SETD2')
genes <- rownames(allgenePMDmean[allgenePMDmean$N>100,])[1:50]
genes <- rownames(allgenePMDmean[order(allgenePMDmean$N, decreasing=T),])[1:50]
genes <- c(
  'DNMT1','DNMT3A','DNMT3B','UHRF1','UHRF2',
  'NSD1', 'SETD2')

genes <- c('STAB2', 'EPHA5', 'NF1', 'TP53')

load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/allmutation.pmdtest.rda')
genes <- unique(rownames(which(genetests < 1e-5, arr.ind=T)))

genes <- 'IDH1'
for (gene in genes) {
  pdf(sprintf('~/gallery/2017_02_03_cancer_mutations/%s.pdf', gene), width=12, height=12)
  par(mfrow=c(6,6), mar=c(4,4,4,1), oma=c(3,3,3,3))
  for (cancer in unique(tumor$cancer)) {
    a <- tumor[tumor$cancer==cancer,]
    a$mutation <- a$barcode %in% allgenedepth[[gene]]$barcode
    a <- a[order(a$mutation),]
    plot(a$commonHMD, a$commonPMD, col=ifelse(a$mutation, '#CC0000', '#005AA0'), cex=0.5, xlim=c(0.0,0.85), ylim=c(0.0,0.85), main=cancer, xlab='HMD methylation', ylab='PMD methylation', font.main=1)
    abline(0,1, lty='dashed')
  }
  dev.off()

  ## pdf(sprintf('~/gallery/2017_02_03_cancer_mutations/%s_highpurity.pdf', gene), width=12, height=12)
  ## par(mfrow=c(6,6), mar=c(4,4,4,1), oma=c(3,3,3,3))
  ## for (cancer in unique(tumor$cancer)) {
  ##   a <- tumor[tumor$cancer==cancer & tumor$purity >= 0.7,]
  ##   a$mutation <- a$barcode %in% allgenedepth[[gene]]$barcode
  ##   a <- a[order(a$mutation),]
  ##   plot(a$commonHMD, a$commonPMD, col=ifelse(a$mutation, '#CC0000', '#005AA0'), cex=0.5, xlim=c(0,0.8), ylim=c(0,0.8), main=cancer, xlab='HMD methylation', ylab='PMD methylation', font.main=1)
  ##   abline(0,1, lty='dashed')
  ## }
  ## dev.off()
}

a <- tumor[tumor$cancer=='HNSC',]
a$mutation <- a$barcode %in% allgenedepth[['NSD1']]$barcode

############################
## mutation PMD regression
############################

## all gene mutation matrix
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/mutations.rda')
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/allTCGA.rda')
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/purity/purity.rda')
load('/primary/projects/laird/projects/2016_04_05_TCGA_pancan_renormalization/450k/merged_mapping.rda')
## focus on high purity samples
tumor <- cbind(merged.mapping.tumor, purity=purity[match(substr(merged.mapping.tumor$barcode,1,15), purity$array),'purity'])
tumor <- cbind(tumor, depths[tumor$idat,])
mutations <- mutations[mutations$outcome %in% c('Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del', 'In_Frame_Ins', 'Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation', 'Splice_Site', 'Translation_Start_Site'),]
mutations <- mutations[!grepl('^ENSG', mutations$gene),]
mutations <- mutations[!grepl('^LOC', mutations$gene),]
genes <- unique(mutations$gene)
genecnt <- table(mutations$gene)
genecnt <- genecnt[25:length(genecnt)]
sum(genecnt > 100)

a <- do.call(cbind, mclapply(names(genecnt), function(gene) substr(tumor$barcode,1,12) %in% mutations[mutations$gene == gene,'patient'], mc.cores=28))
colnames(a) <- names(genecnt)
tumor <- cbind(tumor, a)

load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/mutationrate/mutrate.rda')
mutation <- tumor[,7:ncol(tumor)]
basic <- tumor[,1:6]
basic$mutrate <- mutrate$rate_non[match(substr(basic$barcode,1,16), substr(mutrate$name,1,16))]
basic$delta <- basic$commonPMD - basic$commonHMD
tumor <- cbind(basic, mutation)
save(tumor, file='/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/allmutation.tumor.rda')

## t-test of mutant and wildtype in term of pmd depth
genes <- colnames(tumor)
genes <- genes[7:length(genes)]
bycancer <- split(tumor, tumor$cancer)
genetests <- simplify2array(mclapply(genes, function(gene) {
  sapply(bycancer, function(tt) {
    mut <- tt[tt[,gene],'commonPMD']
    wt <- tt[!tt[,gene],'commonPMD']
    if (length(mut) >= 5 && length(wt) >= 5) {
      rr <- t.test(mut, wt)
      if (rr$statistic > 0) {
        1
      } else {
        t.test(mut, wt)$p.value
      }
    } else {
      1
    }
  })
}, mc.cores=28))
colnames(genetests) <- genes
genetests <- t(genetests)
save(genetests, file='/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/allmutation.pmdtest.rda')

###########################
## all regression analysis
###########################

load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/allmutation.tumor.rda')
## tumor <- tumor[!is.na(tumor$mutrate) & tumor$mutrate < 1e-5,]
mutcount <- sapply(tumor[,9:ncol(tumor)], sum, na.rm=T)
tumor <- cbind(tumor[,1:8], tumor[,9:ncol(tumor)][,mutcount>=50])

## step1: by cancer test to select features within cancer types
bycancer <- split(tumor, tumor$cancer)
bycancer.pvals <- do.call(cbind, mclapply(bycancer, function(t1) {
  sapply(colnames(t1)[9:ncol(t1)], function(genename) {
    fitdata <- t1[,c('purity','delta',genename)]
    fit <- lm(delta~., data=fitdata)
    fitsumm <- summary(fit)
    coef <- fitsumm$coefficients
    if (dim(coef)[1] < 3) {
      1
    } else {
      coef[3,4]
    }
  })
}, mc.cores=28))

## merged test (no-interaction term) to select global features
fitdata <- tumor[,c('cancer','purity','delta',colnames(tumor)[9:ncol(tumor)])]
fit <- lm(delta~., data=fitdata)
fitsumm <- summary(fit)
coef <- fitsumm$coefficients
coef <- as.data.frame(coef[order(coef[,'Pr(>|t|)']),])
coef$padj <- p.adjust(coef[,'Pr(>|t|)'], method='BH')
coefglobal <- coef

save(bycancer.pvals, tumor, coefglobal, file='/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/allmutation.tumor.testbycancer.rda')

## step2: cancer-type mutation: classify mutations by the major cancer types it goes to
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/allmutation.tumor.testbycancer.rda')
genes <- sort(unique(rownames(which(bycancer.pvals<0.01, arr.ind=T))))
mutation0 <- tumor[,genes]
mutation0typed <- do.call(data.frame, lapply(mutation0, function(m1) {
  cnt <- tapply(m1, tumor$cancer, sum, na.rm=T)
  ct <- ifelse(tumor$cancer %in% names(cnt[cnt >= 30]), paste0('.',tumor$cancer), '.Other')
  interaction(ct, m1)
}))
fitdata <- cbind(tumor[,c('purity','cancer','delta')], mutation0typed)
fit <- lm(delta~., data=fitdata)
fitsumm <- summary(fit)
coef <- as.data.frame(fitsumm$coefficients)
coef <- coef[order(coef[,'Pr(>|t|)']),]
coef$padj <- p.adjust(coef[,'Pr(>|t|)'], method='BH')

save(fit, fitsumm, coef, file='/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/allmutation.tumor.testall3.rda')

## effect size
sgenes <- c('KRAS','IDH1','TP53','BRAF','EN2','FOXA1','NRAS','PRR22','MKI67','NSD1','LAMA5','PCDH20','NINL','ZNF536','KIF19','CIC','NUFIP2','RYR2','PCLO','MPP2','TNR','SHROOM2','HBD','PCDHGC5','RAPSN', 'SETD2', 'HOPX', 'KIF1A', 'SYNE1', 'ZNF208', 'SPOP', 'CASP9', 'KLHL1', 'ZNF286B', 'DYSF')
scancers <- c('ACC','BLCA','BRCA','CESC','ESCA','GBM','HNSC','KICH','KIRC','KIRP','LAML','LGG','LUAD','LUSC','PAAD','PRAD','SARC','SKCM','STAD','TGCT','THCA','THYM','UCEC','UCS')
mutcnts <- t(sapply(sgenes, function(gene) {
  sapply(scancers, function(cancer) {
    sum(tumor$cancer==cancer & tumor[[gene]])
  })
}))
pdf('~/gallery/2017_02_09_mutation_cnts.pdf', width=8, height=10)
WHeatmap(mutcnts, xticklabels=T, yticklabels=T, cmp=CMPar(stop.points=c('#A0A0A0','#000000')))
dev.off()

mutfracs <- t(sapply(sgenes, function(gene) {
  sapply(scancers, function(cancer) {
    sum(tumor$cancer==cancer & tumor[[gene]]) / sum(tumor$cancer==cancer)
  })
}))

pdf('~/gallery/2017_02_09_mutation_fracs.pdf', width=8, height=10)
WHeatmap(mutfracs, xticklabels=T, yticklabels=T, cmp=CMPar(stop.points=c('#A0A0A0','#000000')))
dev.off()

## genes <- apply(mutation, 2, sum)
## mutation0 <- mutation[,genes > 100]
## fitdata <- cbind(basic[,c('purity','cancer','commonPMD')], mutation0)
## fit <- lm(commonPMD~., data=fitdata)
## fitsumm <- summary(fit)
## coef <- as.data.frame(fitsumm$coefficients)
## coef <- coef[order(coef[,'Pr(>|t|)']),]
save(fit, fitsumm, coef, file='/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/allmutation.tumor.testall2.rda')

## (digression) multiple test correction
coef$padj <- p.adjust(coef[,4], method='bonferroni')
coefsig <- coef[coef[,4]<0.01,]
coefsig <- coefsig[rownames(coefsig) != '(Intercept)',]
coefsig <- coefsig[!(grepl('ENSG',rownames(coefsig))),]

## (digression) list plot
pdf('~/gallery/2017_02_08_pmd_regression.pdf', width=8, height=5)
par(mar=c(12,5,1,1), lwd=2)
plot(coefsig[,1], cex=log2(-log2(coefsig[,4])) / 2, xaxt='n', ylim=c(-0.3,0.3))
axis(1, at=1:nrow(coefsig), label=rownames(coefsig), las=2, lwd=2)
abline(h=0, lty='dashed')
dev.off()

## mutation-cancer plot
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/allmutation.tumor.testall2.rda')
coef$gene <- sapply(strsplit(rownames(coef),'\\.'), function(x) x[1])
coef$cancer <- sapply(strsplit(rownames(coef),'\\.'), function(x) x[2])
coef$mut <- sapply(strsplit(rownames(coef),'\\.'), function(x) as.logical(x[3]))
coef$cancer <- ifelse(grepl('^cancer', rownames(coef)), sub('cancer','',rownames(coef)), coef$cancer)
coef$gene[grepl('^cancer', rownames(coef)) & rownames(coef) != "purity"] <- 'cancertype'

## hits, genes and cancers
genes <- unique(coef[coef[,4] < 0.05,'gene'])
genes <- sort(genes[!(genes %in% c('(Intercept)','ENSG00000159247','C12orf62','purity','RP11'))], decreasing=T)
hits <- coef[coef$gene %in% genes,]
cancers <- unique(hits[hits$gene != 'cancertype',]$cancer)
cancers <- c(sort(cancers[cancers!='Other']), 'Other')
genes <- rev(hits$gene[!duplicated(hits$gene)])

## hits plot
sgenes <- rev(c('cancertype','KRAS','IDH1','TP53','BRAF','EN2','FOXA1','NRAS','PRR22','MKI67','NSD1','LAMA5','PCDH20','NINL','ZNF536','KIF19','CIC','NUFIP2','RYR2','PCLO','MPP2','TNR','SHROOM2','HBD','PCDHGC5','RAPSN', 'SETD2', 'HOPX', 'KIF1A', 'SYNE1', 'ZNF208', 'SPOP', 'CASP9', 'KLHL1', 'ZNF286B', 'DYSF'))
hits <- hits[hits$gene %in% sgenes,]
hits <- hits[hits$cancer %in% cancers,]
df <- data.frame(
  cancerindex=sapply(hits$cancer, function(x) which(x==cancers)))
df$geneindex <- sapply(hits$gene, function(x) which(x==sgenes))
df$estimate <- pmax(ifelse(is.na(hits$mut) | hits$mut, hits$Estimate, -hits$Estimate),-0.25)
df$pval <- pmax(1e-5, hits[,4])
pdf('~/gallery/2017_02_07_mutation_pmd_matrix.pdf', width=9, height=9)
## #C0C0C0
ggplot(data=df) + geom_point(aes(cancerindex, geneindex, size=-log2(pval), color=estimate)) + scale_color_gradient2(high='#CC0000',low='#005AA0',mid='#C4C4C4') + scale_x_continuous(breaks=seq_along(cancers)+1,labels=cancers) + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5), panel.grid.minor = element_line(colour="grey", size=0.5)) + scale_y_continuous(breaks=seq_along(sgenes),labels=sgenes)
dev.off()

pdf('~/gallery/2017_02_07_mutation_pmd_matrix_light.pdf', width=9, height=9)
## #C0C0C0
ggplot(data=df) + geom_point(aes(cancerindex, geneindex, size=-log2(pval), color=estimate)) + scale_color_gradient2(high='#FFE20A',low='#00BFFF',mid='#FFFFFF')+ scale_x_continuous(breaks=seq_along(cancers)+1,labels=cancers) + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5), panel.grid.minor = element_line(colour="grey", size=0.5)) + scale_y_continuous(breaks=seq_along(sgenes),labels=sgenes)
dev.off()


## regression of cancer types plot
## coef.cancer <- coef[is.na(coef$cancer),]
## coef.cancer <- coef.cancer[grepl('cancer', rownames(coef.cancer)),]
## rownames(coef.cancer) <- sub('cancer','',rownames(coef.cancer))
## coef.cancer <- coef.cancer[cancers,]
## pdf('~/gallery/2017_02_08_pmd_regression_cancertypes.pdf', width=8, height=3)
## par(mar=c(5,5,1,1), lwd=3, cex.axis=1.5)
## plot(coef.cancer[,1], cex=pmax(log2(-log2(coef.cancer[,4])) / 2, 0.2), xaxt='n', ylim=c(-0.4,0.4), ylab='', xlab='', col=ifelse(coef.cancer[,1]>0, '#CC0000', '#005AA0'), pch=19)
## axis(1, at=1:nrow(coef.cancer), label=rownames(coef.cancer), las=2, lwd=2)
## abline(h=0, lty='solid', lwd=1.5)
## dev.off()

##################################
## RNAseq regression correlation
##################################
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/allTCGA.rda')
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/purity/purity.rda')
load('/Users/wandingzhou/projects/hs-tcga/data/2015_04_30_TCGA_rnaseq/pancan.synapse.rnaseq.rda')
load('/primary/projects/laird/projects/2016_04_05_TCGA_pancan_renormalization/450k/merged_mapping.rda')

tumor <- cbind(merged.mapping.tumor, purity=purity[match(substr(merged.mapping.tumor$barcode,1,15), purity$array),'purity'])
gene.exp <- as.data.frame(t(gene.exp))
depths1 <- depths[match(merged.mapping.tumor$idat[match(substr(rownames(gene.exp),1,15), substr(merged.mapping.tumor$barcode,1,15))], rownames(depths)),]
rownames(depths1) <- rownames(gene.exp)
cancertype <- merged.mapping.tumor$cancer[match(substr(rownames(gene.exp),1,15), substr(merged.mapping.tumor$barcode,1,15))]
purity1 <- purity[match(substr(rownames(gene.exp),1,15), substr(rownames(purity),1,15)),'purity']
df <- cbind(data.frame(cancer=cancertype, purity=purity1), depths1, gene.exp)
df <- df[!is.na(df$commonPMD),]
df <- df[!is.na(df$purity),]
df <- df[,!grepl('\\?', colnames(df))]
save(df, file='/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/expression_correlation/allexp.rda')

exprange <- sapply(df[,5:ncol(df)], function(x) max(x, na.rm=T) - min(x, na.rm=T))
igenes <- names(exprange[exprange>=20])
fit <- lm(commonPMD~., data=df[c('cancer','purity','commonPMD',igenes)])

load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/expression_correlation/allexp.rda')
cor.all <- sapply(df[,5:ncol(df)], function(xx) {
  cor(xx, df$commonPMD, use='complete', method='spearman')
})

###### scatter plot ######
## negative correlation
## plot CDCA8 and CENP-A
pdf('~/gallery/2017_02_13_CDCA8_expression.pdf', width=10, height=10)
ggplot(df) + geom_point(aes(x=commonPMD, y=log2(get('CDCA8|55143'))), size=0.5) + facet_wrap(~cancer) + xlab('PMD methylation') + ylab('Log2 CDCA8 expression')
dev.off()
pdf('~/gallery/2017_02_13_CENPA_expression.pdf', width=10, height=10)
ggplot(df) + geom_point(aes(x=commonPMD, y=log2(get('CENPA|1058'))), size=0.5) + facet_wrap(~cancer) + xlab('PMD methylation') + ylab('Log2 CENPA expression')
dev.off()
## positive correlation
pdf('~/gallery/2017_02_13_CELF2_expression.pdf', width=10, height=10)
ggplot(df) + geom_point(aes(x=commonPMD, y=log2(get('CELF2|10659'))), size=0.5) + facet_wrap(~cancer) + xlab('PMD methylation') + ylab('Log2 CELF2 expression')
dev.off()
pdf('~/gallery/2017_02_13_CD37_expression.pdf', width=10, height=10)
ggplot(df) + geom_point(aes(x=commonPMD, y=log2(get('CD37|951'))), size=0.5) + facet_wrap(~cancer) + xlab('PMD methylation') + ylab('Log2 CD37 expression')
dev.off()
pdf('~/gallery/2017_02_13_CNRIP1_expression.pdf', width=10, height=10)
ggplot(df) + geom_point(aes(x=commonPMD, y=log2(get('CNRIP1|25927'))), size=0.5) + facet_wrap(~cancer) + xlab('PMD methylation') + ylab('Log2 CNRIP1 expression')
dev.off()

## high purity
pdf('~/gallery/2017_02_13_CDCA8_expression_highpurity.pdf', width=10, height=10)
ggplot(subset(df, purity>=0.7)) + geom_point(aes(x=commonPMD, y=log2(get('CDCA8|55143'))), size=0.5) + facet_wrap(~cancer) + xlab('PMD methylation') + ylab('Log2 CDCA8 expression')
dev.off()
pdf('~/gallery/2017_02_13_CENPA_expression_highpurity.pdf', width=10, height=10)
ggplot(subset(df, purity>=0.7)) + geom_point(aes(x=commonPMD, y=log2(get('CENPA|1058'))), size=0.5) + facet_wrap(~cancer) + xlab('PMD methylation') + ylab('Log2 CENPA expression')
dev.off()
## positive correlation
pdf('~/gallery/2017_02_13_CELF2_expression_highpurity.pdf', width=10, height=10)
ggplot(subset(df, purity>=0.7)) + geom_point(aes(x=commonPMD, y=log2(get('CELF2|10659'))), size=0.5) + facet_wrap(~cancer) + xlab('PMD methylation') + ylab('Log2 CELF2 expression')
dev.off()
pdf('~/gallery/2017_02_13_CD37_expression_highpurity.pdf', width=10, height=10)
ggplot(subset(df, purity>=0.7)) + geom_point(aes(x=commonPMD, y=log2(get('CD37|951'))), size=0.5) + facet_wrap(~cancer) + xlab('PMD methylation') + ylab('Log2 CD37 expression')
dev.off()
pdf('~/gallery/2017_02_13_CNRIP1_expression_highpurity.pdf', width=10, height=10)
ggplot(subset(df, purity>=0.7)) + geom_point(aes(x=commonPMD, y=log2(get('CNRIP1|25927'))), size=0.5) + facet_wrap(~cancer) + xlab('PMD methylation') + ylab('Log2 CNRIP1 expression')
dev.off()

load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/expression_correlation/allexp.rda')
df <- df[df$purity >= 0.6,]             #not quite consistent with the following threshold with 0.7
op <- options(warn=(-1))                      #suppress warning
cor.allcancertypes <- do.call(cbind, mclapply(split(df, df$cancer), function(df1) {
  sapply(df1[,5:ncol(df1)], function(xx) {
    cor(xx, df1$commonHMD-df1$commonPMD, use='na.or.complete', method='spearman')
  })
}, mc.cores=28))
cor.allcancertypes.pval <- do.call(cbind, mclapply(split(df, df$cancer), function(df1) {
  sapply(df1[,5:ncol(df1)], function(xx) {
    cor.test(xx, df1$commonHMD-df1$commonPMD, use='na.or.complete', method='spearman')$p.value
  })
}, mc.cores=28))
cor.allcancertypes.pval <- as.numeric(cor.allcancertypes.pval)
cor.allcancertypes.pval <- p.adjust(cor.allcancertypes.pval, method='BH') # BH-adjusted
dim(cor.allcancertypes.pval) <- dim(cor.allcancertypes)
dimnames(cor.allcancertypes.pval) <- dimnames(cor.allcancertypes)
options(op)
save(cor.allcancertypes, cor.allcancertypes.pval, file='/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/expression_correlation/allexpbycancer_highpurity.rda') #4749 tumors

###################
## DAVID analysis
###################
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/expression_correlation/allexpbycancer_highpurity.rda')
## cor.allcancertypes <- cor.allcancertypes[,colnames(cor.allcancertypes)!='OV']
## cor.allcancertypes.mean <- apply(cor.allcancertypes,1,function(x) mean(abs(x)))
top.pos <- sort(na.omit(apply(cor.allcancertypes,1,function(x) sum(x > 0.4))),decreasing=T)
top.pos <- data.frame(top.pos, gene=sapply(strsplit(names(top.pos),'\\|'), function(x) x[1]))
top.neg <- sort(na.omit(apply(cor.allcancertypes,1,function(x) sum(x < -0.4))),decreasing=T)
top.neg <- data.frame(top.neg, gene=sapply(strsplit(names(top.neg),'\\|'), function(x) x[1]))
write.table(top.pos, file='/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/expression_correlation/allexpbycancer_highpurity_top_pos.tsv', sep='\t', col.names=NA)
write.table(top.neg, file='/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/expression_correlation/allexpbycancer_highpurity_top_neg.tsv', sep='\t', col.names=NA)
## cd /secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/expression_correlation
## awk 'NR==FNR{a[$1]=$5;}FNR>1 && NR!=FNR{if(a[$3]) b=a[$3]; else b="nonperiodic"; print $1,$2,$3,b}' cellcycle.txt allexpbycancer_highpurity_top_neg.tsv >allexpbycancer_highpurity_top_neg_cellcycle.tsv
## awk 'NR==FNR{a[$1]=$5;}FNR>1 && NR!=FNR{if(a[$3]) b=a[$3]; else b="nonperiodic"; print $1,$2,$3,b}' cellcycle.txt allexpbycancer_highpurity_top_pos.tsv >allexpbycancer_highpurity_top_pos_cellcycle.tsv

df <- data.frame(
  BH.pval = c(9.58E-25, 6.72E-24, 3.38E-23, 2.61E-15, 9.23E-13, 5.01E-11),
  Fold.Enrichment = c(32.60708627, 16.13024673, 23.01901381, 47.09730313, 34.01349676, 39.22431866))
rownames(df) <- c('Mitosis', 'Cell cycle', 'Cell division', 'Sister Chromatid Cohesion', 'Centromere', 'Kinetochore')
df$logpval <- -log10(df$BH.pval)
pdf('~/gallery/2017_02_13_allexpbycancer_highpurity_top_pos_DAVID_pval.pdf', width=5, height=4.5)
barplot(df$logpval, names.arg=rownames(df))
dev.off()
pdf('~/gallery/2017_02_13_allexpbycancer_highpurity_top_pos_DAVID_FC.pdf', width=5, height=4.5)
barplot(df$Fold.Enrichment, names.arg=rownames(df), las=1)
dev.off()

df <- read.table('/Volumes/projects_secondary/laird/projects/2016_12_26_TCGA_WGBS/hm450/expression_correlation/allexpbycancer_highpurity_top_pos_DAVID_topFDR.txt', sep='\t', header=T)
df <- df[grep('GO:',df$Term),]
pdf('~/gallery/2017_02_20_allexpbycancer_highpurity_top_pos_DAVID_pval.pdf', width=5, height=4.5)
barplot(log10(df[1:20,]$FDR), names.arg=rownames(df[1:20,]$Term))
dev.off()
pdf('~/gallery/2017_02_20_allexpbycancer_highpurity_top_pos_DAVID_FC.pdf', width=5, height=4.5)
barplot(df[1:20,]$Fold.Enrichment, names.arg=df[1:20,]$Term, las=1, ylim=c(0,60))
dev.off()

#########################
## cell cycle analysis
#########################
## load correlation
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/expression_correlation/allexpbycancer_highpurity.rda')
cor.allcancertypes <- cor.allcancertypes[, colnames(cor.allcancertypes) != 'OV']

## load peak time
peaktime <- read.table('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/expression_correlation/allexpbycancer_highpurity_top_pos_cellcycle.tsv', header=F, col.names=c('id', 'pos0p4cnt', 'gene', 'peaktime'))
bypeak <- split(peaktime$id, peaktime$peaktime)
bypeak <- bypeak[!(names(bypeak) %in% c('non-periodic','nonperiodic'))]
allperiodic <- unique(do.call(c, bypeak))

## plot correlation by phase
aclus <- both.cluster(cor.allcancertypes[bypeak[['M']],])$mat
a <- WHeatmap(row.cluster(cor.allcancertypes[bypeak[['G1']],colnames(aclus)])$mat, cmp=CMPar(dmin=-0.5, dmax=0.5), xticklabels=T, xticklabel.side='top', yticklabels=T, yticklabels.n=length(bypeak[['G1']]))
for (phase in c('G1/S','S','G2','G2/M','M')) {
  a <- a + WHeatmap(row.cluster(cor.allcancertypes[bypeak[[phase]],colnames(aclus)])$mat, Beneath(pad=0.02), cmp=CMPar(dmin=-0.5, dmax=0.5), yticklabels=T, yticklabels.n=length(bypeak[[phase]]))
}
pdf('~/gallery/2017_02_14_expression_pmd_cor_bycancer.pdf', width=10, height=50)
print(a)
dev.off()

## plot correlation of all phases merged
mm <- both.cluster(cor.allcancertypes[allperiodic, ])
WHeatmap(mm$mat, cmp=CMPar(dmin=-0.5, dmax=0.5), xticklabels=T, xticklabel.side='top')

## output gene order (requested by Hui)
a <- do.call(rbind, lapply(c('G1','G1/S','S','G2','G2/M','M'), function(phase) {
  data.frame(gene=rownames(row.cluster(cor.allcancertypes[bypeak[[phase]],])$mat),
             mean=rowMeans(row.cluster(cor.allcancertypes[bypeak[[phase]],])$mat),
             phase=phase)
}))
rownames(a) <- sapply(strsplit(rownames(a), '\\|'), function(x) x[1])
write.table(a, file='~/gallery/2017_02_14_expression_by_pmd_geneorder.tsv', sep='\t', col.names=NA)

## replication timing
## awk 'NR>1' ~/gallery/2017_02_14_expression_by_pmd_geneorder.tsv | gsed 's/"//g' | cut -f1 | sort | awk 'NR==FNR{a[$1]}NR!=FNR && ($1 in a){print $7,$5,$6,$1;}' - ~/tools/transvari/transvar/transvar/transvar.download/hg19.ccds.txt.transvardb | sortbed | bedtools merge -i - -c 4 -o distinct | bedtools intersect -a - -b /secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/IMR90/repliseq/ENCFF000KUR.bed -sorted -wo | bedtools groupby -g 1-4 -c 8 -o mean >/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/expression_correlation/cellcycle_repliseq.tsv
a <- read.table('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/expression_correlation/cellcycle_repliseq.tsv', sep='\t')
genereplitime <- setNames(a$V5, a$V4)

## by phases
aclus <- both.cluster(cor.allcancertypes[bypeak[['M']],])$mat
a <- WHeatmap(row.cluster(cor.allcancertypes[bypeak[['G1']],colnames(aclus)])$mat, cmp=CMPar(dmin=-0.5, dmax=0.5), xticklabels=T, xticklabel.side='top', yticklabels=T, yticklabels.n=length(bypeak[['G1']]))
phase <- 'G1'
for (phase in c('G1/S','S','G2','G2/M','M')) {
  mm <- row.cluster(cor.allcancertypes[bypeak[[phase]],colnames(aclus)])$mat
  rownames(mm) <- sapply(strsplit(rownames(mm), '\\|'), function(x) x[1])
  a <- WHeatmap(mm, cmp=CMPar(dmin=-0.5, dmax=0.5), xticklabels=T, xticklabel.side='top', yticklabels=T, yticklabels.n=length(bypeak[[phase]])) + WColorBarV(genereplitime[rownames(mm)], RightOf())
  pdf(sprintf('~/gallery/2017_02_14_expression_pmd_cor_bycancer_phase_%s.pdf',sub('\\/','',phase)), width=10, height=15)
  print(a)
  dev.off()
}

#################################################################
### plot raw expression value heatmap (instead of correlation)
#################################################################
## load data
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/expression_correlation/allexp.rda')
df <- df[order(df$commonHMD - df$commonPMD), ] # order by HMD-PMD
## df <- df[complete.cases(df),]
## remove samples with more than 70% 0-expressed genes
zeroexpressioncnt <- apply(df, 1, function(x) sum(x[5:length(x)] < 1))
zeroexpressionperiodic <- apply(df[,allperiodic], 1, function(x) sum(x<1))

## by-phase correlation plot
df0 <- df[df$purity >= 0.7,]
cors <- sapply(df0[5:ncol(df0)], function(x) cor(x, df0$commonHMD-df0$commonPMD, use='complete', method='spearman'))

cors.df <- as.data.frame(cors)
cors.df$phase <- 'nonperiodic'
for (phase in c('G1','G1/S','S','G2','G2/M','M'))
  cors.df[bypeak[[phase]],'phase'] <- phase
cors.df$abscors <- abs(cors.df$cors)
cors.df$gene <- sapply(strsplit(as.character(rownames(cors.df)),'\\|'), function(x) x[1])
cors.df <- cors.df[!is.na(cors.df$abscors),]
cors.df1 <- cors.df[,c('gene','abscors')]
write.table(cors.df1, file='/Volumes/projects_secondary/laird/projects/2016_12_26_TCGA_WGBS/hm450/expression_correlation/2017_02_16_expression_pmd_cor.txt', quote=F, sep='\t', row.names=FALSE)

## fake GSEA
cors.df <- cors.df[order(cors.df$cors, decreasing=T),]
unitp <- 1/sum(cors.df[cors.df$phase!='nonperiodic','abscors'])
unitnp <- 1/sum(cors.df$phase=='nonperiodic')
cors.df$step <- 1
cors.df[cors.df$phase=='nonperiodic', 'step'] <- - unitnp
cors.df[cors.df$phase!='nonperiodic', 'step'] <- cors.df[cors.df$phase!='nonperiodic', 'abscors'] * unitp
cors.df$sumcors <- cumsum(cors.df$step)
pdf('~/gallery/2017_02_16_GSEA_expression_pmd_cor_periodic.pdf', width=4.5, height=4)
par(mar=c(4,4,1,1), mfrow=c(2,1), lwd=3, font.axis=2, cex.axis=1.8)
plot(cors.df$sumcors, type='l', ylab='ES(S)', xlab='Gene', lwd=3)
plot(cors.df$cors, type='l', ylab='Rho', col='red', lwd=3, ylim=c(-0.6, 0.6))
dev.off()

cors.df <- do.call(rbind, lapply(split(cors.df, cors.df$phase), function(x) {
  x <- x[order(x$cors),];
  x$order <- (1:nrow(x)) / nrow(x)
  x
}))
pdf('~/gallery/2017_02_16_expression_pmd_cor.pdf', width=5.5, height=3.5)
ggplot(cors.df) + geom_line(aes(order, cors, color=phase), size=1.2) + xlab('Ordered Genes') + ylab("Spearman's Rho") + guides(color=guide_legend(title='Phase'))
dev.off()

## version 1 rank-based
# heatmap
df0 <- cbind(df[,1:4], t(apply(df[,5:ncol(df)],1,rank))) # rank expression within sample
df0 <- df0[df0$purity >= 0.7,]
df0 <- df0[, c(1:4, 4+which(colMeans(df0[,5:ncol(df0)]) >= 10))] # threshold RPKM
a <- WColorBarH(df0$commonHMD-df0$commonPMD)
for (phase in c('G1','G1/S','S','G2','G2/M','M')) {
  df1 <- df0[,intersect(colnames(df0), bypeak[[phase]])]
  df11 <- apply(df1, 2, function(x) {
    (x - mean(x,na.rm=T)) / sd(x, na.rm=T)
  })
  ## colsd <- apply(df11, 2, function(x) sd(x, na.rm=T))
  ## browser()
  ## cat(colsd)
  a <- a + WHeatmap(row.cluster(t(df11))$mat, cmp=CMPar(stop.points=c('green', 'black', 'red'), dmin=-1.2, dmax=1.2), Beneath(pad=2), name=phase)
}
png('~/gallery/2017_02_14_expression_by_pmd_rank.png', width=3000, height=1000)
print(a)
dev.off()
pdf('~/gallery/2017_02_14_expression_by_pmd_rank_legend.pdf', width=3, height=3)
WColorBarH(1:10/10, cmp=CMPar(stop.points=c('green', 'black', 'red'), dmin=-1.2, dmax=1.2), name='M') + WLegendV('M',RightOf())
dev.off()

## plot DNMT and UHRF
geneselect <- c("TET1|80312", "TET2|54790", "TET3|200424", "DNMT1|1786", "DNMT3A|1788", "DNMT3B|1789", "DNMT3L|29947", "UHRF1|29128", "UHRF1BP1|54887", "UHRF1BP1L|23074", "UHRF2|115426")
for (gene in geneselect) {
  df01 <- apply(df0[,gene,drop=F], 2, function(x) {(x - mean(x,na.rm=T)) / sd(x, na.rm=T);})
  png(sprintf('~/gallery/2017_02_20_expression_by_pmd_%s.png', gene), width=3000, height=50)
  print(WHeatmap(t(df01), cmp=CMPar(stop.points=c('green', 'black', 'red'), dmin=-1.2, dmax=1.2)) + WCustomize(mar.left=0, mar.top=0, mar.right=0, mar.bottom=0))
  dev.off()
}

## too many cancer types to color properly
cancertypes <- df0$cancer
cancertypes[!(cancertypes %in% c('BLCA', 'BRCA', 'COAD', 'LGG', 'THCA', 'SKCM', 'LIHC'))] <- 'OTHERS'
png('~/gallery/2017_02_20_expression_by_pmd_cancertypes.png', width=3000, height=50)
WColorBarH(cancertypes, cmp=CMPar(label2color=setNames(c('#A6CEE3', '#1F78B4', '#B2DF8A', '#33A02C', '#B21E00', '#FFC887', '#8DD3C7', '#CCCCCC'), c('BLCA', 'BRCA', 'COAD', 'LGG', 'THCA', 'SKCM', 'LIHC', 'OTHERS')))) + WCustomize(mar.left=0, mar.top=0, mar.right=0, mar.bottom=0)
dev.off()
pdf('~/gallery/2017_02_20_expression_by_pmd_cancertypes_legend.pdf', width=3, height=3)
WColorBarH(t(cancertypes), cmp=CMPar(label2color=setNames(c('#A6CEE3', '#1F78B4', '#B2DF8A', '#33A02C', '#B21E00', '#FFC887', '#8DD3C7', '#CCCCCC'), c('BLCA', 'BRCA', 'COAD', 'LGG', 'THCA', 'SKCM', 'LIHC', 'OTHERS'))), name='a') + WLegendV('a',RightOf())
dev.off()

# line-plot
periodicmean <- rowMeans(df0[,allperiodic])
aa <- data.frame(delta=df0$commonHMD - df0$commonPMD, periodicmean=periodicmean)
aa$order <- seq_along(aa$delta)
bb <- apply(df0, 1, function(x)
  quantile(as.numeric(x[5:length(x)]), c(0.25,0.75,0.8,0.85), na.rm=T))
bb <- as.data.frame(t(bb))
colnames(bb) <- c('p25','p75','p80','p85')
aa <- cbind(aa, bb, df0[,1:4])
## aa$zeroexpression <- zeroexpressioncnt
## aa$zeroexpressionperiodic <- zeroexpressionperiodic
pdf('~/gallery/2017_02_15_periodicmean_HMDmPMD.pdf', width=5, height=2)
ggplot(aa) + geom_line(aes(order, delta), color='red')
dev.off()
pdf('~/gallery/2017_02_15_periodicmean_meanexp_rank.pdf', width=5, height=2)
ggplot(aa) + geom_smooth(aes(order, periodicmean), color='blue') + geom_smooth(aes(order, p75), color='orange') + geom_smooth(aes(order,p25), color='orange') + ylab('Rank')
dev.off()
pdf('~/gallery/2017_02_15_periodicmean_meanexp_rank_nopct.pdf', width=5, height=2)
ggplot(aa) + geom_smooth(aes(order, periodicmean), color='blue') + ylab('Rank')
dev.off()
pdf('~/gallery/2017_02_15_zeroexpressioncnt.pdf', width=5, height=2)
ggplot(aa) + geom_smooth(aes(order, zeroexpression), color='blue') + ylab('Zero Exp. Cnt')
dev.off()

## STAD
load('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/expression_correlation/allexp.rda')
ebvbarcode <- read.table('/Users/wandingzhou/projects/hs-tcga/data/2015_03_23_Hui_annotation/20161218_STAD_EBV.txt', header=T)$barcode
df <- df[df$cancer == 'STAD',]
df <- cbind(df[,1:4], t(apply(df[,5:ncol(df)],1,rank))) # rank expression within sample
aastad <- df[,1:4]
aastad$delta <- aastad$commonHMD - aastad$commonPMD
aastad$periodicmean <- rowMeans(df[,allperiodic])
aastad$EBV <- substr(rownames(aastad),1,12) %in% ebvbarcode
pdf('~/gallery/2017_02_16_PeriodicMean_vs_PMD_STAD.pdf', width=5.5, height=4)
ggplot(aastad) + geom_point(aes(delta, periodicmean, color=EBV), size=2) + xlab('HMD - PMD') + ylab('Mean Periodic Gene Expression')
dev.off()

## kidney
x <- read.table('/Users/wandingzhou/projects/hs-tcga/data/2015_03_23_Hui_annotation/nejmoa1505917_appendix_3_simple.txt', sep='\t', header=T)
kcimpbarcode <- x[x$dnamecluster=='CIMP','barcode']
aakidney <- aa[aa$cancer %in% c('KIRP', 'KICH', 'KIRC'),]
aakidney$cimp <- substr(rownames(aakidney),1,12) %in% kcimpbarcode
pdf('~/gallery/2017_02_16_PeriodicMean_vs_PMD_kidney.pdf', width=5.5, height=4)
ggplot(aakidney) + geom_point(aes(delta, periodicmean, shape=cancer, color=cimp), size=2) + xlab('HMD - PMD') + ylab('Mean Periodic Gene Expression')
dev.off()

## other cancer types
load('~/projects/hs-tcga/data/2015_03_23_Hui_annotation/subtypes.rda')
df0 <- cbind(df[,1:4], t(apply(df[,5:ncol(df)],1,rank))) # rank expression within sample
aa <- data.frame(
  delta = df0$commonHMD - df0$commonPMD,
  periodicmean = rowMeans(df0[,allperiodic]))
aa <- cbind(df0[,1:4], aa)
aa$subtype <- subtypes[match(substr(rownames(aa),1,15), substr(subtypes$id_TCGA,1,15)), 'subtype']
cancers <- unique(subtypes$cancer)
cancers <- cancers[cancers != 'OV']
for (cancer in cancers) {
  cat(cancer,'\n')
  aacancer <- aa[aa$cancer == cancer,]
  aacancer$subtype[is.na(aacancer$subtype)] <- 'NA'
  if (dim(aacancer)[1] > 5) {
    pdf(sprintf('~/gallery/2017_02_16_PeriodicMean_vs_PMD_auto_%s.pdf', cancer), width=5.5, height=4)
    print(ggplot(aacancer) + geom_point(aes(delta, periodicmean, color=subtype), size=2) + xlab('HMD - PMD') + ylab('Mean Periodic Gene Expression'))
    dev.off()
  }
}

## version 2 raw Log2 RPKM
# heatmap
## df0 <- df[df$purity >= 0.7,]
df0 <- df
df0 <- df0[, c(1:4, 4+which(colMeans(df0[,5:ncol(df0)]) >= 50))] # threshold RPKM
df0 <- df0[df0$cancer %in% c('CHOL','PRAD','STAD','ACC','SARC','THCA','KIRC','PCPG','KIRP','BLCA','UCEC','LUAD','UVM','LGG','MESO'),]
a <- WColorBarH(df0$commonHMD-df0$commonPMD)
for (phase in c('G1','G1/S','S','G2','G2/M','M')) {
  df1 <- df0[,intersect(colnames(df0), bypeak[[phase]])]
  df11 <- apply(df1, 2, function(x) {
    lx <- log(1+x)
    (lx - mean(lx,na.rm=T)) / sd(lx, na.rm=T)
  })
  a <- a + WHeatmap(row.cluster(t(df11))$mat, cmp=CMPar(dmin=-1.2, dmax=1.2), Beneath(pad=1))
}
png('~/gallery/2017_02_14_expression_by_pmd_highcancertypes.png', width=3000, height=1000)
print(a)
dev.off()

# line-plot
allperiodic <- intersect(unique(do.call(c, bypeak)), colnames(df0))
periodicmean <- log2(1+rowMeans(df0[,allperiodic]))
aa <- data.frame(delta=df0$commonHMD - df0$commonPMD, periodicmean=periodicmean)
aa$order <- seq_along(aa$delta)
bb <- apply(df0, 1, function(x)
  quantile(log(1+as.numeric(x[5:length(x)])), c(0.25,0.75,0.8,0.85), na.rm=T))
bb <- as.data.frame(t(bb))
colnames(bb) <- c('p25','p75','p80','p85')
aa <- cbind(aa, bb, df0[,1:4])
pdf('~/gallery/2017_02_15_periodicmean_HMDmPMD_highcancertypes.pdf', width=5, height=2)
ggplot(aa) + geom_line(aes(order, delta), color='red')
dev.off()
pdf('~/gallery/2017_02_15_periodicmean_meanexp_highcancertypes.pdf', width=5, height=2)
ggplot(aa[aa$purity > 0.7,]) + geom_smooth(aes(order, periodicmean), color='blue') + geom_smooth(aes(order, p75), color='orange') + ylab('Rank')
dev.off()

## by cancertyes
bycancer <- split(df, df$cancer)
for (cancer in names(bycancer)) {
  df0 <- bycancer[[cancer]]
  df0 <- df0[df0$purity >= 0.6,]
  a <- WColorBarH(df0$commonHMD-df0$commonPMD)
  for (phase in c('G1','G1/S','S','G2','G2/M','M')) {
    df1 <- df0[,bypeak[[phase]]]
    df11 <- apply(df1, 2, function(x) {
      lx <- log(x+1);
      (lx - mean(lx,na.rm=T)) / sd(lx, na.rm=T)
    })
    a <- a + WHeatmap(row.cluster(t(df11))$mat, cmp=CMPar(dmin=-1.2, dmax=1.2), Beneath(pad=1))
  }
  png(sprintf('~/gallery/2017_02_14_expression_by_pmd_%s.png', cancer), width=3000, height=1000)
  print(a)
  dev.off()
}

load('~/projects/hs-tcga/data/2015_03_23_Hui_annotation/subtypes.rda')
df0 <- bycancer[['LGG']]
df0$subtype <- subtypes$subtype[match(substr(rownames(df0),1,15), substr(subtypes$id_TCGA,1,15))]
df1 <- df0[, bypeak[['M']]]
df11 <- apply(df1, 2, function(x) {
  lx <- log(x+1);
  (lx - mean(lx,na.rm=T)) / sd(lx, na.rm=T)
})
mm <- both.cluster(t(df11))$mat
df01 <- df0[colnames(mm),]
WHeatmap(mm, cmp=CMPar(dmin=-1.2, dmax=1.2)) + WColorBarH(df01$purity, TopOf(), label='purity') + WColorBarH(df01$commonHMD-df01$commonPMD, TopOf(), label='HMD-PMD') + WColorBarH(df01$subtype, TopOf(), label='subtype', name='a') + WLegendV('a',TopRightOf())


## cor.allcancertypes.mean <- apply(cor.allcancertypes,1,function(x) mean(x, na.rm=T))
## cor.allcancertypes.mean <- apply(cor.allcancertypes[,-28],1,function(x) mean(x, na.rm=T))

## head(top.neg, n=30)
##     CELF2|10659        CD37|951    CNRIP1|25927    CYFIP2|26999     GAPT|202309
##               9               8               8               8               8
##       GNG7|2788     PTPRN2|5799  ARHGAP31|57514      BEX4|56271     CACNA1C|775
##               8               8               7               7               7
##       CCND2|894       CCR6|1235       FLI1|2313       GAS7|8522   GLT1D1|144423
##               7               7               7               7               7
##       GYPC|2995     IL10RA|3587  PPP1R16B|26051      PRKCQ|5588      RAI2|10742
##               7               7               7               7               7
##        SLA|6503     ZC4H2|55906  ARHGAP15|55843      BNC2|54796 C14orf64|388011
##               7               7               6               6               6
##  C8orf48|157773   C8orf79|57604     CBFA2T3|863        CD22|933       CD52|1043
##               6               6               6               6               6

cancers <- sort(colnames(cor.allcancertypes))
genes <- names(top.pos)
hits <- melt(cor.allcancertypes[rownames(cor.allcancertypes) %in% pos.genes,
                                colnames(cor.allcancertypes) %in% cancers])
hits.pval <- melt(cor.allcancertypes.pval[rownames(cor.allcancertypes.pval) %in% pos.genes,
                                          colnames(cor.allcancertypes.pval) %in% cancers])
df <- data.frame(
  cancerindex = match(hits$Var2, cancers),
  geneindex = match(hits$Var1, genes),
  estimate = hits$value)
df$direction <- df$estimate > 0
df$correlation <- abs(df$estimate)

t(apply(which(abs(cor.allcancertypes) > 0.9, arr.ind=T), 1, function(x) c(rownames(cor.allcancertypes)[x[1]], colnames(cor.allcancertypes)[x[2]], cor.allcancertypes[x[1],x[2]])))

corlist <- melt(cor.allcancertypes)
colnames(corlist) <- c('ID', 'cancer', 'r')
corlist$gene <- sapply(strsplit(as.character(corlist$ID),'\\|'), function(x) x[1])
corlist$absr <-  abs(corlist$r)
corlist <- corlist[order(corlist$absr, decreasing=T),]
bycancer <- split(corlist, corlist$cancer)
for (cancer in names(bycancer)) {
  write.table(bycancer[[cancer]], file=sprintf('/secondary/projects/laird/projects/2016_12_26_TCGA_WGBS/hm450/expression_correlation/allexpbycancerlist/allexpbycancerlist_%s.tsv', cancer), sep='\t', col.names=NA)
}
