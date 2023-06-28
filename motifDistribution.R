library(scales)
library(readr)

getCGContent <- function(seq){
  cg.len <- length(grep("C|G", strsplit(seq, "")[[1]]))
  return(round(cg.len/nchar(seq), digits = 4))
}

# CA CMP data
all.sk <- read.delim('all.table.sk.txt', stringsAsFactors = F)
for (i in 1:nrow(all.sk)){
  motifs <- strsplit(all.sk$motif[i], ';')[[1]]
  all.sk$d.motif[i] <- motifs[nchar(motifs) <= min(nchar(motifs))][1]
  all.sk$gc.cont[i] <- getCGContent(all.sk$d.motif[i])
}
all.sk$motif.length <- nchar(all.sk$d.motif)

all.sk.cmp <- all.sk[all.sk$cmp_cl.count > 0,]# 305
all.sk.child <- all.sk[all.sk$control_cl.count > 0,]# 320


# UK CMP data
all.ge <- read.delim('all.table.ge.txt', stringsAsFactors = F)
for (i in 1:nrow(all.ge)){
  motifs <- strsplit(all.ge$motif[i], ';')[[1]]
  all.ge$d.motif[i] <- motifs[nchar(motifs) <= min(nchar(motifs))][1]
  all.ge$gc.cont[i] <- getCGContent(all.ge$d.motif[i])
}
all.ge$motif.length <- nchar(all.ge$d.motif)

all.ge.cmp <- all.ge[all.ge$cmp_count >0,]# 280
all.ge.child <- all.ge[all.ge$cntrl_count >0,]# 232


h.sk.cmp.gc <- hist(all.sk.cmp$gc.cont)
h.ge.cmp.gc <- hist(all.ge.cmp$gc.cont)

h.sk.child.gc <- hist(all.sk.child$gc.cont)
h.ge.child.gc <- hist(all.ge.child$gc.cont)

label.pts <- h.sk.cmp.gc$mids*46
label.val <- h.sk.cmp.gc$breaks[-1]

test.sk <- wilcox.test(all.sk.child$gc.cont, all.sk.cmp$gc.cont)
format(test.sk$p.value, scientific = T, digits = 2) # "1.3e-06"
## GC hist one by one plus boxplots CA
pdf('figure_1_gc-composition_sk_.pdf', w = 4, h = 4)
par(mar = c(1,5,3,3), mgp = c(3,1,0), mfrow=c(2,1)) 
b <- barplot(h.sk.cmp.gc$counts/sum(h.sk.cmp.gc$counts),col = alpha('dodgerblue3', 0.8), 
        xlab = '', ylab = 'Percent of regions', ylim = c(0, 0.3),yaxt = 'n')
axis(1, at = b[c(1,5,10,15,20)], labels = c('0.0', '0.25', '0.5', '0.75','1.0'), las = 1, padj = -1, tick = F)
axis(2, at = seq(0,0.3, 0.1), labels = T, las = 2)

par(mar = c(6,5,1,3), mgp = c(1,1,0)) 
boxplot(all.sk.child$gc.cont, all.sk.cmp$gc.cont,
        col = c(alpha('grey28', 0.6), alpha('dodgerblue3', 0.8)), 
        yaxt = 'n', xaxt = 'n', ylab ='', xlab = 'GC-composition (%)', horizontal = T)

axis(2, at = seq(1,2), las = 2, labels = c('CA control', 'CA CMP'))
text(0.8, 1.5, expression(1.3~x~10^-6), cex = 0.8)
text(0.6, 1.35, paste('p = '), cex = 0.8)
dev.off()


test.ge <- wilcox.test(all.ge.child$gc.cont, all.ge.cmp$gc.cont)
format(test.ge$p.value, scientific = T, digits = 2) #  "4.3e-02"
## GC hist one by one plus boxplots UK
pdf('figure_1_gc-composition_ge_.pdf', w = 4, h = 4)
par(mar = c(1,5,3,3), mgp = c(3,1,0), mfrow=c(2,1)) 
b <- barplot(h.ge.cmp.gc$counts/sum(h.ge.cmp.gc$counts),col = alpha('lightskyblue1', 0.8), 
            xlab = '', ylab = 'Percent of regions', ylim = c(0, 0.3),yaxt = 'n')
axis(1, at = b[c(1,5,10,15,20)], labels = c('0.0', '0.25', '0.5', '0.75','1.0'), las = 1, padj = -1, tick = F)
axis(2, at = seq(0,0.3, 0.1), labels = T, las = 2)

par(mar = c(6,5,1,3), mgp = c(1,1,0)) 
boxplot(all.ge.child$gc.cont, all.ge.cmp$gc.cont,
        col = c(alpha('grey58', 0.6), alpha('lightskyblue1', 0.8)), 
        yaxt = 'n', xaxt = 'n', ylab ='', xlab = 'GC-composition (%)', horizontal = T)

axis(2, at = seq(1,2), las = 2, labels = c('UK control', 'UK CMP'))
text(0.85, 1.5, expression(4.3~x~10^-2), cex = 0.8)
text(0.65, 1.35, paste('p = '), cex = 0.8)
dev.off()


h.sk.cmp.gc <- hist(all.sk.cmp$motif.length)
h.ge.cmp.gc <- hist(all.ge.cmp$motif.length)

h.sk.child.gc <- hist(all.sk.child$motif.length)
h.ge.child.gc <- hist(all.ge.child$motif.length)

test.sk <- wilcox.test(all.sk.child$motif.length, all.sk.cmp$motif.length)
format(test.sk$p.value, scientific = T, digits = 2) 
## Motif hist one by one plus boxplots CA
pdf('figure_1_motif_length_sk_.pdf', w = 4, h = 4)
par(mar = c(1,5,3,3), mgp = c(3,1,0), mfrow=c(2,1)) 
b <- barplot(h.sk.cmp.gc$counts/sum(h.sk.cmp.gc$counts),col = alpha('dodgerblue3', 0.8), 
            xlab = '', ylab = 'Percent of regions', ylim = c(0, 0.3),yaxt = 'n')
axis(1, at = b[c(1,3,5,9,12,15,18)], labels = c(h.sk.cmp.gc$breaks[c(1,3,5,9,12,15,18)]), las = 1, padj = -1, tick = F)
axis(2, at = seq(0,0.3, 0.1), labels = T, las = 2)

par(mar = c(6,5,1,3), mgp = c(1,1,0)) 
boxplot(all.sk.child$motif.length, all.sk.cmp$motif.length,
        col = c(alpha('grey28', 0.6), alpha('dodgerblue3', 0.8)), 
        yaxt = 'n', xaxt = 'n', ylab ='', xlab = 'Motif size (bp)', horizontal = T)

axis(2, at = seq(1,2), las = 2, labels = c('CA control', 'CA CMP'))
text(15, 1.35, paste('p = ', format(test.sk$p.value, digits = 2)), cex = 0.8)
dev.off()

test.ge <- wilcox.test(all.ge.child$motif.length, all.ge.cmp$motif.length)
format(test.ge$p.value, scientific = T, digits = 2) 
## Motif hist one by one plus boxplots UK
pdf('figure_1_motif_length_ge_.pdf', w = 4, h = 4)#, h = 7
par(mar = c(1,5,3,3), mgp = c(3,1,0), mfrow=c(2,1)) 
b <- barplot(h.ge.cmp.gc$counts/sum(h.ge.cmp.gc$counts),col = alpha('lightskyblue1', 0.8), 
            xlab = '', ylab = 'Percent of regions', ylim = c(0, 0.3),yaxt = 'n')
axis(1, at = b[c(1,3,5,9,12,15,18)], labels = c(h.sk.cmp.gc$breaks[c(1,3,5,9,12,15,18)]), las = 1, padj = -1, tick = F)
axis(2, at = seq(0,0.3, 0.1), labels = T, las = 2)

par(mar = c(6,5,1,3), mgp = c(1,1,0)) 
boxplot(all.sk.child$motif.length, all.sk.cmp$motif.length,
        col = c(alpha('grey58', 0.6), alpha('lightskyblue1', 0.8)), 
        yaxt = 'n', xaxt = 'n', ylab ='', xlab = 'Motif size (bp)', horizontal = T)

axis(2, at = seq(1,2), las = 2, labels = c('UK control', 'UK CMP'))
text(15, 1.35, paste('p = ', format(test.ge$p.value, digits = 2)), cex = 0.8)
dev.off()
