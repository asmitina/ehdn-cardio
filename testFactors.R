library(GenomicRanges)
library(data.table)
library(plotrix)
library(scales)
library(readr)

ci <- function(vector, interval) {
  v.sd <- sd(vector) # Standard deviation 
  n <- length(vector) # Sample size
  v.mean <- mean(vector) # Sample mean
  error <- qt((interval + 1)/2, df = n - 1) * v.sd / sqrt(n) # Error
  result <- c('n' = as.numeric(n), 'mean' = v.mean, 'se' = error, "lower" = v.mean - error, "upper" = v.mean + error) # Confidence interval
  return(result)
}

#### Annotation and rare expansions files
ann <- fread("merged.trs.annotation.tsv", data.table = F)
ann <- ann[ann$chr %in% paste0("chr", 1:22), ] 
ann$typeseq_priority <- factor(ann$typeseq_priority, 
                               levels = c("exonic", "splicing", "exonic;splicing", "ncRNA_exonic", "ncRNA_splicing", "ncRNA_exonic;ncRNA_splicing", 
                                          "UTR5", "UTR3", "intronic", "ncRNA_intronic", "upstream", "downstream", "intergenic"))
ann <- ann[order(ann$typeseq_priority), ]
ann <- ann[!duplicated(ann$varid), ]
ann <- ann[ann$typeseq_priority %in% c("exonic", "splicing", "exonic;splicing",  
                                       "UTR5", "UTR3", "intronic", "upstream", "downstream", "intergenic"), ] # 32299 -> 41319 -> 32305

feat <- list("exonic" = "exonic",
             "splicing" = c("splicing", "exonic;splicing"),
             "UTR5" = "UTR5","UTR3" = "UTR3",
             "intronic" = "intronic","upstream" = "upstream", 
             "downstream" = "downstream", "intergenic" = "intergenic")

dt <- read.delim("merged.rare.expansions.tsv", stringsAsFactors = F)
dt$repeatID <- paste(dt$chr, dt$start, dt$end, dt$motif, sep="#")
dt$varid <- paste(dt$chr, dt$start, dt$end, sep="#")
ann <- ann[ann$varid %in% dt$varid, ]


#### Sampleinfo
sampleinfo <- read.delim("samples.with.TRcount.txt", stringsAsFactors = F)
outliers <- readLines("samples.with.exceed.trs.txt")
sampleinfo <- sampleinfo[(!sampleinfo$Var1 %in% outliers) & (sampleinfo$cohort != "1000G"), ]
names(sampleinfo) <- c("sample", "trs", "cohort")
sampleinfo <- sampleinfo[sampleinfo$sample %in% clean.sample,]
sampleinfo$status <- ifelse(sampleinfo$cohort == "CMP", 1, 0)
sampleinfo <- sampleinfo[!duplicated(sampleinfo$sample),]

# To analyse clinical factors in TSS-proximal rare TREs use this:
# tss.5k <- read.delim('5k.tss.cmp.child.tsv', stringsAsFactors = F)
# dt = dt[dt$varid %in% tss.5k$varid,]

#### Genotype
start <- Sys.time()
tmp <- do.call(rbind, lapply(1:nrow(dt), getOutlierMergeData, dt, sampleinfo))
end <- Sys.time()
difftime(end, start)
tmp <- tmp[, -2]
tmp$varid <- paste(tmp$chr, tmp$start, tmp$end, sep="#")

ehdn.genotype.data <- tmp$sample[tmp$varid %in% ann$varid]
ehdn.genotype.data <- data.frame(table(ehdn.genotype.data), stringsAsFactors = F)
names(ehdn.genotype.data) <- c("sample", "ncloci")

ehdn.genotype.data$cohort = ifelse(ehdn.genotype.data$sample %in% cmp.clean$Sample, 'CMP', 'CHILD')
ehdn.genotype.cmp <- ehdn.genotype.data[ehdn.genotype.data$cohort == 'CMP',]

### keep only CMP samples
ehdn.genotype.data <- ehdn.genotype.data[ehdn.genotype.data$sample %in% cmp.meta$Sample,]
cmp.meta$ehdn <- 0

for (i in 1:nrow(cmp.meta)){
  if (cmp.meta$Sample[i] %in% ehdn.genotype.data$sample){
    cmp.meta$ehdn[i] <- ehdn.genotype.data$ncloci[ehdn.genotype.data$sample == cmp.meta$Sample[i]]
  }
}

ehdn.genotype.cmp <- ehdn.genotype.cmp[!duplicated(ehdn.genotype.cmp$Sample),]
ehdn.genotype.cmp <- ehdn.genotype.cmp[ehdn.genotype.cmp$Relationship == 'PROBAND' 
                                      | ehdn.genotype.cmp$Relationship == 'INDEX',]
#### Test sex
wilcox.sex <- as.data.frame(rbind(ci(ehdn.genotype.cmp$ehdn[ehdn.genotype.cmp$Sex == 'F'], 0.95), 
                                 ci(ehdn.genotype.cmp$ehdn[ehdn.genotype.cmp$Sex == 'M'], 0.95)))
test <- wilcox.test(ehdn.genotype.cmp$ehdn[ehdn.genotype.cmp$Sex == 'F'],
                   ehdn.genotype.cmp$ehdn[ehdn.genotype.cmp$Sex == 'M'])
wilcox.sex$p.value <- rep(test$p.value, 2)
wilcox.sex$name <- rep('female', 2)
wilcox.sex$condition <- c(T,F)

#### Test ethnicity
# wilcox.eur <- as.data.frame(rbind(ci(ehdn.genotype.cmp$ehdn[ehdn.genotype.cmp$PredictedAncestry == 'EUR'], 0.95), 
#                                  ci(ehdn.genotype.cmp$ehdn[ehdn.genotype.cmp$PredictedAncestry != 'EUR'], 0.95)))
# test <- wilcox.test(ehdn.genotype.cmp$ehdn[ehdn.genotype.cmp$PredictedAncestry == 'EUR'],
#                    ehdn.genotype.cmp$ehdn[ehdn.genotype.cmp$PredictedAncestry != 'EUR'])[3]
# wilcox.eur$p.value <- rep(test$p.value, 2)
# wilcox.eur$name <- rep('eur', 2)
# wilcox.eur$condition <- c(T,F)
# wilcox.out <- rbind(wilcox.sex, wilcox.eur)

#### Test age
wilcox.onset <- as.data.frame(rbind(ci(ehdn.genotype.cmp$ehdn[ehdn.genotype.cmp$AgeAtPhenotypicOnset..years. <= 5], 0.95), 
                                   ci(ehdn.genotype.cmp$ehdn[ehdn.genotype.cmp$AgeAtPhenotypicOnset..years. > 5], 0.95)))
test <- wilcox.test(ehdn.genotype.cmp$ehdn[ehdn.genotype.cmp$AgeAtPhenotypicOnset..years. <= 5],
                   ehdn.genotype.cmp$ehdn[ehdn.genotype.cmp$AgeAtPhenotypicOnset..years. > 5])[3]
wilcox.onset$p.value <- rep(test$p.value, 2)
wilcox.onset$name <- rep('early onset', 2)
wilcox.onset$condition <- c(T,F)
wilcox.out <- rbind(wilcox.out, wilcox.onset)


#### Test confirmed variants
for (i in 1:nrow(ehdn.genotype.cmp)){
  ehdn.genotype.cmp$variant[i] <- cmp.meta$variant[cmp.meta$Sample == ehdn.genotype.cmp$Sample[i]]
}
ehdn.genotype.save <- ehdn.genotype.cmp

ehdn.genotype.cmp <- ehdn.genotype.cmp[!is.na(ehdn.genotype.cmp$variant),]
wilcox.variants <- as.data.frame(rbind(ci(ehdn.genotype.cmp$ehdn[ehdn.genotype.cmp$variant == T], 0.95), 
                                      ci(ehdn.genotype.cmp$ehdn[ehdn.genotype.cmp$variant == F], 0.95)))
test <- wilcox.test(ehdn.genotype.cmp$ehdn[ehdn.genotype.cmp$variant == T],
                   ehdn.genotype.cmp$ehdn[ehdn.genotype.cmp$variant == F])[3]
wilcox.variants$p.value <- rep(test$p.value, 2)
wilcox.variants$name <- rep('confirmed variant', 2)
wilcox.variants$condition <- c(T,F)
wilcox.out <- rbind(wilcox.out, wilcox.variants)


#### Test diagnosis
ehdn.genotype.cmp <- ehdn.genotype.save
table(ehdn.genotype.cmp$DX)
ehdn.genotype.cmp$dx.2 <- ehdn.genotype.cmp$DX
ehdn.genotype.cmp$dx.2[grep('DCM ', ehdn.genotype.cmp$DX)] <- 'DCM'
ehdn.genotype.cmp$dx.2[grep('HCM ', ehdn.genotype.cmp$DX)] <- 'HCM'
ehdn.genotype.cmp$dx.2[grep('LVNC ', ehdn.genotype.cmp$DX)] <- 'LVNC'
ehdn.genotype.cmp$dx.2[grep('LV dysfx ', ehdn.genotype.cmp$DX)] <- 'LVD'
ehdn.genotype.cmp$dx.2 <- as.character(ehdn.genotype.cmp$dx.2)

wilcox.dcm <- as.data.frame(rbind(ci(ehdn.genotype.cmp$ehdn[grep('DCM', ehdn.genotype.cmp$dx.2)], 0.95), 
                                 ci(ehdn.genotype.cmp$ehdn[-grep('DCM', ehdn.genotype.cmp$dx.2)], 0.95)))

test <- wilcox.test(ehdn.genotype.cmp$ehdn[grep('DCM', ehdn.genotype.cmp$dx.2)],
                   ehdn.genotype.cmp$ehdn[-grep('DCM', ehdn.genotype.cmp$dx.2)])[3]
wilcox.dcm$p.value <- rep(test$p.value, 2)
wilcox.dcm$name <- rep('have DCM', 2)
wilcox.dcm$condition <- c(T,F)
wilcox.out <- rbind(wilcox.out, wilcox.dcm)

#### Go back and regenerate ehdn (keep parents in metainfo)
ehdn.genotype.cmp <- cmp.meta
ehdn.genotype.other <- ehdn.genotype.cmp[ehdn.genotype.cmp$Relationship != 'PROBAND' 
                                        & ehdn.genotype.cmp$Relationship != 'INDEX',]
ehdn.genotype.proband <- ehdn.genotype.cmp[ehdn.genotype.cmp$Relationship == 'PROBAND' 
                                          | ehdn.genotype.cmp$Relationship == 'INDEX',]

families <- sort(intersect(ehdn.genotype.other$Family, ehdn.genotype.proband$Family))
ehdn.genotype.family <- rbind(ehdn.genotype.other[ehdn.genotype.other$Family %in% families,],
                             ehdn.genotype.proband[ehdn.genotype.proband$Family %in% families,])
ehdn.genotype.family <- ehdn.genotype.family[order(ehdn.genotype.family$Family),]
ehdn.genotype.family$short.rel <- sapply(ehdn.genotype.family$Relationship, switch, SISTER='sibling', BROTHER='sibling', 'HALF-BROTHER' = 'sibling',
                                        'MOTHER'='parent', 'FATHER'='parent', 'PROBAND - FATHER'='parent',
                                        'PROBAND - MOTHER'='parent', 'PROBAND-MOTHER'='parent', PROBAND='proband',INDEX='proband')
ehdn.genotype.family$short.rel <- as.character(ehdn.genotype.family$short.rel)
ehdn.genotype.burden <- ehdn.genotype.family[ehdn.genotype.family$short.rel != 'sibling',]

affected.parent <- data.frame()
for (family in unique(ehdn.genotype.burden$Family)){
  tmp <- ehdn.genotype.burden[ehdn.genotype.burden$Family == family,]
  affected.parent <- rbind(affected.parent, data.frame('sample' = tmp$Sample[tmp$short.rel == 'proband'],
                                                      'affected.parent' = 'affected' %in% tmp$short.dx[tmp$short.rel == 'parent'],
                                                      'ehdn' = tmp$ehdn[tmp$short.rel == 'proband']))
}

#### Test for an affected parent
wilcox.aff.parent <- as.data.frame(rbind(ci(affected.parent$ehdn[affected.parent$affected.parent == T], 0.95), 
                                        ci(affected.parent$ehdn[affected.parent$affected.parent == F], 0.95)))
test <- wilcox.test(affected.parent$ehdn[affected.parent$affected.parent == T],
                   affected.parent$ehdn[affected.parent$affected.parent == F])[3]
wilcox.aff.parent$p.value <- rep(test$p.value, 2)
wilcox.aff.parent$name <- rep('affected parent', 2)
wilcox.aff.parent$condition <- c(T,F)
wilcox.out <- rbind(wilcox.out, wilcox.aff.parent)

write.table(wilcox.out, 'test.factors.rare.all.tsv', quote = F, row.names = F, col.names = T, sep = ',')
# write.table(wilcox.out, 'test.factors.rare.TSS.tsv', quote = F, row.names = F, col.names = T, sep = ',')


wilcox.out <- read.delim('test.factors.rare.all.tsv', stringsAsFactors = F, sep = ',')
wilcox.out$col <- ifelse(wilcox.out$condition == T, 'orange', 'maroon')
wilcox.out$alpha <- ifelse(wilcox.out$p.value <= 0.05, 1, 0.5)
wilcox.out$name <- sapply(wilcox.out$name,switch,'eur'='European ancestry','confirmed variant'='Pathogenic variant',
                         'female'='Female','early onset'='Early onset','have DCM'='Have DCM','affected parent'='Affected parent')
# plotx = c(1:17)[-c(3*1:5)] 

wilcox.out <- wilcox.out[!wilcox.out$name == 'European ancestry',]
plotx <- c(1:15)[-c(3*1:5)]

#### Plot for all expansions
pdf('factors.burden.rare.all.pdf', w = 6, h = 4)
par(mar = c(8,4,3,6), mgp = c(3,1,0), mfrow=c(1,1)) 
plotCI(plotx,wilcox.out$mean, ui = wilcox.out$upper, li = wilcox.out$lower, ylab = 'Burden of rare TREs', xlab = '',
       col = alpha(wilcox.out$col, wilcox.out$alpha), 
       xaxt = 'n', yaxt = 'n', lwd = 1.2, pch = 15, cex = 1.5, sfrac = 0.01)
text(seq(1,length(plotx), by = 2)*1.5, par("usr")[3]-0.2, srt = 45, adj = 1, xpd = TRUE, labels = unique(wilcox.out$name))
axis(1, at = seq(1,length(plotx), by = 2)*1.5, labels = NA)
axis(2, at = seq(0,3.5,by = 0.5), las = 2)
legend('topright', pch=15, lwd = 1.2, lty = 1, 
       pt.cex=1.5, inset=c(-0.25,0), xpd = TRUE, bty = "n", 
       col = c('orange', 'maroon'), legend=c("TRUE","FALSE"), cex = 0.8)

# format(wilcox.out$p.value[1], scientific = T, digits = 2)# "3.5e-02"
text(x=1,y=2,expression(3.5~x~10^-02),srt=90, cex = 0.9)
text(x=1.2,y=1.3,paste('p = '),srt=90, cex = 0.9)

# format(wilcox.out$p.value[3], scientific = T, digits = 2)#  "2.9e-03"
text(x=4,y=2,expression(2.9~x~10^-03),srt=90, cex = 0.9)
text(x=4.2,y=1.3,paste('p = '),srt=90, cex = 0.9)# 

dev.off()

#### Plot for all expansions (Supplementary)
pdf('factors.burden.rare.all.eth.pdf', w = 6, h = 4)
par(mar = c(8,4,3,6), mgp = c(3,1,0), mfrow=c(1,1)) 
plotCI(plotx,wilcox.out$mean, ui = wilcox.out$upper, li = wilcox.out$lower, ylab = 'Burden of rare TREs', xlab = '',
       col = alpha(wilcox.out$col, wilcox.out$alpha), 
       xaxt = 'n', yaxt = 'n', lwd = 1.2, pch = 15, cex = 1.5, sfrac = 0.01)
text(seq(1,length(plotx), by = 2)*1.5, par("usr")[3]-0.2, srt = 45, adj = 1, xpd = TRUE, labels = unique(wilcox.out$name))
axis(1, at = seq(1,length(plotx), by = 2)*1.5, labels = NA)
axis(2, at = seq(0,3.5,by = 0.5), las = 2)
legend('topright', pch=15, lwd = 1.2, lty = 1, 
       pt.cex=1.5, inset=c(-0.25,0), xpd = TRUE, bty = "n", 
       col = c('orange', 'maroon'), legend=c("TRUE","FALSE"), cex = 0.8)

# format(wilcox.out$p.value[3], scientific = T, digits = 1)# "1.3e-05"
text(x=3.8,y=0.35,expression(1.3~x~10^-05),srt=90, cex = 0.9)
text(x=4.0,y=0.2,paste('p = '),srt=90, cex = 0.9)# 

# format(wilcox.out$p.value[13], scientific = T, digits = 2)#  "4.2e-02"
text(x=16,y=0.35,expression(4.2~x~10^-02),srt=90, cex = 0.9)
text(x=16.2,y=0.17,paste('p = '),srt=90, cex = 0.9)#
dev.off()
