library(BSgenome.Hsapiens.UCSC.hg38)
library(phastCons100way.UCSC.hg38)
library(GenomicRanges)
library(Biostrings)
library(data.table)
library(RSQLite)
library(ordinal)
library(scales)

bin.dt <- read.delim("subset_phylop_1kb_bins.tsv", stringsAsFactors = F)# 3050056
bin.dt <- bin.dt[bin.dt$seqnames %in% paste0("chr", c(1:22, "X", "Y")), ]# 2937992
bin.g <- GRanges(bin.dt$seqnames, IRanges(bin.dt$start, bin.dt$end), "*")

expansion <- read.delim('merged.ehdn.tsv', stringsAsFactors = F)
expansion$repeatID <- paste(expansion$chr, expansion$start, expansion$end, expansion$motif, sep="#")

dt <- read.delim("merged.rare.expansions.tsv", stringsAsFactors = F)
dt$repeatID <- paste(dt$chr, dt$start, dt$end, dt$motif, sep="#")
dt$varid <- paste(dt$chr, dt$start, dt$end, sep="#")

for (i in 1:nrow(dt)){ 
  samples <- strsplit(dt$outliers[i], ';')[[1]]
  dt$cmp_clean[i] <- paste0(samples[samples %in% sampleinfo$sample[sampleinfo$cohort == 'CMP']], collapse = ';')
  dt$control_clean[i] <- paste0(samples[samples %in% sampleinfo$sample[sampleinfo$cohort == 'CHILD']], collapse = ';')
  
  dt$cmp_cl.count[i] <- length(samples[samples %in% sampleinfo$sample[sampleinfo$cohort == 'CMP']])
  dt$control_cl.count[i] <- length(samples[samples %in% sampleinfo$sample[sampleinfo$cohort == 'CHILD']])
}

dt.cmp <- dt[dt$cmp_cl.count > 0,]
dt.child <- dt[dt$control_cl.count > 0,]

expansion <- expansion[expansion$repeatID %in% dt.cmp$repeatID,]
# expansion <- expansion[expansion$repeatID %in% dt.child$repeatID,]

expansion.g <- GRanges(expansion$chr, IRanges(expansion$start, expansion$end), "*")
olap <- data.frame(findOverlaps(bin.g, expansion.g))
olap <- aggregate(subjectHits ~ queryHits, olap, length)
bin.dt$expan <- 0
bin.dt$expan[olap$queryHits] <- olap$subjectHits

fragile <- read.delim("fragileSites_hg38.tsv", stringsAsFactors = F)# 148
fragile <- fragile[!(is.na(fragile$start.hg38.)), ]# 146
fragile.g <- GRanges(fragile$chr.hg38., IRanges(as.numeric(fragile$start.hg38.), as.numeric(fragile$end.hg38.)), "*")
olap <- data.frame(findOverlaps(bin.g, fragile.g))
bin.dt$fragile <- 0
bin.dt$fragile[unique(olap$queryHits)] <- 1

bin.dt <- bin.dt[bin.dt$seqnames %in% paste0("chr", c(1:22, "X", "Y")), ]
features <- c("fragile", "GC", "phylop.mean")

bin.dt[is.na(bin.dt)] <- 0
bin.dt$EHdn <- factor(bin.dt$expan > 0)
p <- list()

dt.out <- data.frame()
type <- "EHdn"
for(f in features){
  bin.dt$feat <- bin.dt[, f]#scale(bin.dt[, f])#
  
  lm <- glm(sprintf("%s ~ feat", type), data = bin.dt, family = binomial(link = "logit"))
  pvalue <- summary(lm)$coefficients["feat", "Pr(>|z|)"]
  
  conf <- confint.default(lm)
  dt.out <- rbind(dt.out, data.frame("feature" = f, "Odds ratio" = exp(lm$coefficients["feat"]),
                                     "OR.lower" = exp(conf["feat", 1]),
                                     "OR.upper" = exp(conf["feat", 2]),
                                     "type" = type,
                                     "pvalue" = pvalue, stringsAsFactors = F))
}

dt.out.cmp <- dt.out
dt.out.cmp$type <- 'CA CMP'

dt.out.child <- dt.out
dt.out.child$type <- 'child'

# write.table(dt.out.cmp, 'bins.out.cmp.tsv', quote = F, sep = ',', row.names = F, col.names = T)
# write.table(dt.out.child, 'bins.out.child.upd.tsv', quote = F, sep = ',', row.names = F, col.names = T)

dt.out.cmp <- read.delim('bins.out.cmp.upd.tsv', stringsAsFactors = F, sep = ',')
dt.out.child <- read.delim('bins.out.child.upd.tsv', stringsAsFactors = F, sep = ',')

dt.out <- rbind(dt.out.cmp, dt.out.child)
dt.out$col <- c(rep(alpha('dodgerblue3', 0.6), 3), rep(alpha('grey28', 0.6), 3))

dt.plot <- dt.out[c(1,4,2,5,3,6),] 
plotx <- c(1:8)[-c(3*1:5)]

dt.plot$feature <- sapply(dt.plot$feature,switch,fragile='Fragile sites',GC="GC content",phylop="PhyloP")

pdf('fragile.cmp.child.pdf', w = 4, h = 4)
par(mar = c(8,5,3,6), mgp = c(2,1,0), mfrow=c(1,1)) 
b <- barplot(dt.plot$Odds.ratio, col = dt.plot$col, space = c(rep(c(1,0.2),3)),
            yaxt = 'n', ylim = c(0,1.6), ylab = 'Odds Ratio')
arrows(b, dt.plot$OR.lower, b, dt.plot$OR.upper, col='black', lwd = 1, code = 3, angle = 90, length = 0.05)# code=1
abline(h = 1, lty = 2)
axis(2, at = c(0,1,5,12), las = 2)
text(seq(1.4,length(plotx), by = 2)*1.6, par("usr")[3]-0.1, srt = 45, adj = 1, xpd = TRUE, labels = unique(dt.plot$feature))
legend(x = "topright", legend = c('CA CMP','CA control'), fill = c(alpha('dodgerblue3', 0.6), alpha('grey28', 0.6)), 
       inset=c(-0,-0.2), bty = "n", xpd = TRUE)
dev.off()
