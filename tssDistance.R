library(GenomicRanges)
library(readr)

refflat <- read.delim("hg38_refFlat.txt", stringsAsFactors = F, header = F)# 77227
ref.info = read.delim('refflat.tsv', stringsAsFactors = F, header = F, sep = '')
colnames(refflat) = gsub(';', '',ref.info$V2)
refflat = refflat[,1:8]
refflat$txSize = refflat$txEnd - refflat$txStart
refflat$cdsSize = refflat$cdsEnd - refflat$cdsStart
refflat$tis100kbstart <- ifelse(refflat$strand == "+", refflat$txStart-10000, refflat$txEnd-10000)
refflat$tis100kbend <- ifelse(refflat$strand == "+", refflat$txStart+10000, refflat$txEnd+10000)
refflat.g <- GRanges(refflat$chrom, IRanges(refflat$tis100kbstart, refflat$tis100kbend))


getdistance <- function(refflat.g, rare.expansion.g, refflat, rare.expansion){
  rare.olap <- data.frame(findOverlaps(rare.expansion.g, refflat.g))
  rare.olap$strain <- refflat$strand[rare.olap$subjectHits]
  rare.olap$tss <- ifelse(rare.olap$strain == "+", 
                          refflat$tis100kbstart[rare.olap$subjectHits] + 10000, 
                          refflat$tis100kbend[rare.olap$subjectHits] - 10000)
  rare.olap$expansion.start <- rare.expansion$start[rare.olap$queryHits]
  rare.olap$expansion.end <- rare.expansion$end[rare.olap$queryHits]
  rare.olap$mid.point <- (rare.olap$expansion.start + rare.olap$expansion.end) / 2
  
  rare.olap$distance <- rare.olap$tss - rare.olap$mid.point
  rare.olap$distance <- ifelse(rare.olap$strain == "+", rare.olap$distance, -rare.olap$distance)
  rare.olap <- rare.olap[order(abs(rare.olap$distance)), ]
  rare.olap <- rare.olap[!duplicated(rare.olap$queryHits), ]
  return(rare.olap)
}

rare.expansion <- read.delim("merged.rare.expansions.tsv", stringsAsFactors = F)
rare.expansion <- rare.expansion[rare.expansion$chr %in% paste0("chr", c(1:22)), ]
rare.expansion$size <- rare.expansion$end - rare.expansion$start

for (i in 1:nrow(rare.expansion)){ 
  samples = strsplit(rare.expansion$outliers[i], ';')[[1]]
  rare.expansion$cmp_clean[i] = paste0(samples[samples %in% sampleinfo$sample[sampleinfo$cohort == 'CMP']], collapse = ';')
  rare.expansion$control_clean[i] = paste0(samples[samples %in% sampleinfo$sample[sampleinfo$cohort == 'CHILD']], collapse = ';')
  
  rare.expansion$cmp_cl.count[i] = length(samples[samples %in% sampleinfo$sample[sampleinfo$cohort == 'CMP']])
  rare.expansion$control_cl.count[i] = length(samples[samples %in% sampleinfo$sample[sampleinfo$cohort == 'CHILD']])
}

rare.cmp = rare.expansion[rare.expansion$cmp_cl.count > 0 ,]
# median.size.cmp <- median(rare.cmp$size)# 1233 -> 1234.5
# iqr.size.cmp <- IQR(rare.cmp$size)# 713 -> 685.25
rare.cmp <- data.frame(reduce(GRanges(rare.cmp$chr.x, IRanges(as.numeric(rare.cmp$start.x), as.numeric(rare.cmp$end.x)), "*")))
rare.cmp.g <- GRanges(rare.cmp$seqnames, IRanges(rare.cmp$start, rare.cmp$end), "*")
rare.distance.cmp <- getdistance(refflat.g, rare.cmp.g, refflat, rare.cmp) # 59
rare.distance.cmp$type <- "rare cmp"

rare.child = rare.expansion[rare.expansion$control_clean > 0 ,]
# median.size.child <- median(rare.child$size)# 1085 -> 1072.5
# iqr.size.child <- IQR(rare.child$size)# 584 -> 580
rare.child <- data.frame(reduce(GRanges(rare.child$chr.x, IRanges(as.numeric(rare.child$start.x), as.numeric(rare.child$end.x)), "*")))
rare.child.g <- GRanges(rare.child$seqnames, IRanges(rare.child$start, rare.child$end), "*")
rare.distance.child <- getdistance(refflat.g, rare.child.g, refflat, rare.child) # 64
rare.distance.child$type <- "rare child"

detected.expansion <- read.delim("merged.ehdn.tsv", stringsAsFactors = F)
detected.expansion <- detected.expansion[detected.expansion$chr %in% paste0("chr", c(1:22)), ]
detected.expansion$size <- detected.expansion$end - detected.expansion$start
detected.expansion <- data.frame(reduce(GRanges(detected.expansion$chr, IRanges(detected.expansion$start, detected.expansion$end), "*")))
detected.expansion.g <- GRanges(detected.expansion$seqnames, IRanges(detected.expansion$start, detected.expansion$end), "*")
detected.distance <- getdistance(refflat.g, detected.expansion.g, refflat, detected.expansion)
detected.distance$type <- "all expansions"

# Plot the distance to TSS from rare TREs in CA control
test = wilcox.test(abs(rare.distance.child$distance)[abs(rare.distance.child$distance) <= 5000], 
                   abs(detected.distance$distance)[abs(detected.distance$distance) <= 5000], 
                   alternative = "less")$p.value 
pdf('density.plot.ca.control.pdf', w = 4, h = 3.5)
par(mar = c(5,5,3,2), mgp = c(2,1,0), mfrow=c(1,1)) 
plot(density(rare.distance.child$distance, bw = 100), col = 'lightskyblue3', lwd = 2, xlim = c(-5000,5000), 
     ylab = '', xaxt="n", yaxt = "n", main = '', xlab = 'Distance from TSS', ylim = c(0,0.0005)) # plots the results
polygon(density(rare.distance.child$distance, bw = 100), col = alpha('lightskyblue3', 0.3))
lines(density(rare.distance.child$distance, bw = 100), col = "lightskyblue3", lwd = 2, xlim = c(-5000,5000))
polygon(density(detected.distance$distance, bw = 100), col = alpha('grey28', 0.3))
lines(density(detected.distance$distance, bw = 100), col = "grey28", lwd = 2, xlim = c(-5000,5000))
axis(2, at=0, labels = 0, las = 2)
axis(2, at=0.0001, labels = expression(1~x~10^-04), las = 2, font = 5)
axis(2, at=0.0003, labels = expression(3~x~10^-04), las = 2, font = 5)
axis(2, at=0.0005, labels = expression(5~x~10^-04), las = 2, font = 5)

axis(1, at=seq(-5000, 5000, by=5000), labels = TRUE, las = 1)
abline(v = 0, lty = 2)
legend('topright', lwd = 2, lty = 1, pt.cex=1, inset=c(-0,-0.3), xpd = TRUE, bty = "n",
       col = c('lightskyblue3'), legend =c('rare TREs in CA control'))
legend('topright', lwd = 2, lty = 1, pt.cex=1, inset=c(-0,-0.2), xpd = TRUE, bty = "n",
       col = c('grey28' ), legend =c('all detected TREs'))
legend('topright', legend = paste0('p = ', signif(test, digit = 2)), bty = "n", inset = c(0.1, 0), xpd = TRUE)
dev.off()

# Plot the distance to TSS from rare TREs in CA CMP
test = wilcox.test(abs(rare.distance.cmp$distance)[abs(rare.distance.cmp$distance) <= 5000], 
                   abs(detected.distance$distance)[abs(detected.distance$distance) <= 5000], 
                   alternative = "less")$p.value
pdf('density.plot.ca.control.pdf', w = 4, h = 3.5)
par(mar = c(5,5,3,2), mgp = c(2,1,0), mfrow=c(1,1)) 
plot(density(rare.distance.cmp$distance, bw = 100), col = 'dodgerblue3', lwd = 2, xlim = c(-5000,5000), 
     ylab = '', xaxt="n", yaxt = "n", main = '', xlab = 'Distance from TSS')
polygon(density(rare.distance.cmp$distance, bw = 100), col = alpha('dodgerblue3', 0.3))
lines(density(rare.distance.cmp$distance, bw = 100), col = "dodgerblue3", lwd = 2, xlim = c(-5000,5000))
polygon(density(detected.distance$distance, bw = 100), col = alpha('grey28', 0.3))
lines(density(detected.distance$distance, bw = 100), col = "grey28", lwd = 2, xlim = c(-5000,5000))
axis(2, at=0, labels = 0, las = 2)
axis(2, at=0.0001, labels = expression(1~x~10^-04), las = 2, font = 5)
axis(2, at=0.0003, labels = expression(3~x~10^-04), las = 2, font = 5)
axis(2, at=0.0005, labels = expression(5~x~10^-04), las = 2, font = 5)

axis(1, at=seq(-5000, 5000, by=5000), labels = TRUE, las = 1)
abline(v = 0, lty = 2)
legend('topright', lwd = 2, lty = 1, pt.cex=1, inset=c(-0,-0.3), xpd = TRUE, bty = "n",
       col = c('dodgerblue3'), legend =c('rare TREs in CMP'))
legend('topright', lwd = 2, lty = 1, pt.cex=1, inset=c(-0,-0.2), xpd = TRUE, bty = "n",
       col = c('grey28' ), legend =c('all detected TREs'))
legend('topright', legend = paste0('p = '), bty = "n", inset = c(0.3, 0), xpd = TRUE)
# format(test, scientific = T, digits = 2)# "1.2e-02"
legend('topright', legend = expression(1.2~x~10^-02), bty = "n", inset = c(-0.0, -0.008), xpd = TRUE)
dev.off()


#### get regions within 5K from TSS
cmp.dist = rare.distance.cmp[abs(rare.distance.cmp$distance) <= 5000,]
cmp.select = rare.cmp[rownames(rare.cmp) %in% cmp.dist$queryHits,]
cmp.select$varid = paste(cmp.select$seqnames, cmp.select$start, cmp.select$end, sep="#")
cmp.tss.5k = all[all$varid %in% cmp.select$varid,]

child.dist = rare.distance.child[abs(rare.distance.child$distance) <= 5000,]
child.select = rare.child[rownames(rare.child) %in% child.dist$queryHits,]
child.select$varid = paste(child.select$seqnames, child.select$start, child.select$end, sep="#")# 30
child.tss.5k = all[all$varid %in% child.select$varid,]# 24;;; 30
child.tss = sort(unique(child.tss.5k$gene_symbol))
child.tss.old = readLines('gene_set/core.child.genes.txt')
setdiff(child.tss.old, child.tss)
write.table(child.tss, "gene_set/core.child.genes.updated.22.txt", sep="\t", row.names = F, col.names=F, quote=F)
write.table(child.tss, "Documents/_Ryan/_Seema_Myocardium/_INOVA_CHILD_EHdn/gene_set/core.child.genes.updated.22.txt", sep="\t", row.names = F, col.names=F, quote=F)

tss.5k <- rbind(cmp.tss.5k, child.tss.5k)
tss.5k = tss.5k[order(tss.5k$varid),]
tss.5k = tss.5k[!duplicated(tss.5k$varid),]
write.table(tss.3.5k, "5k.tss.cmp.child.tsv", sep="\t", row.names = F, col.names=T, quote=F)
