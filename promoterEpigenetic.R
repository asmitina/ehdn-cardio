library(GenomicRanges)
library(data.table)
library(plotrix)
library(readr)

#### Get annotation and expansions tables
ann <- fread("merged.expansions.annotation.tsv", data.table = F)
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
             "UTR5" = "UTR5",
             "UTR3" = "UTR3",
             "intronic" = "intronic",
             "upstream" = "upstream",
             "downstream" = "downstream",
             "intergenic" = "intergenic")

dt <- read.delim("merged.rare.expansions.tsv", stringsAsFactors = F)
dt$repeatID <- paste(dt$chr, dt$start, dt$end, dt$motif, sep="#")
dt$varid <- paste(dt$chr, dt$start, dt$end, sep="#")

tss.5k <- read.delim('5k.tss.cmp.child.tsv', stringsAsFactors = F)
dt <- dt[dt$varid %in% tss.5k$varid,]
ann <- ann[ann$varid %in% dt$varid, ]


#### PROMOTER ####
promoter <- read.delim('promoter.coordinates.bed', stringsAsFactors = F, sep = ' ')
promoter.g <- GRanges(promoter$chr, IRanges(as.numeric(promoter$start), as.numeric(promoter$end)), "*")
expansion.g <- GRanges(ann$chr, IRanges(ann$start, ann$end), "*")

olap <- data.frame(findOverlaps(promoter.g, expansion.g))
olap$queryHits[duplicated(olap$subjectHits)]
olap <- aggregate(queryHits ~ subjectHits, olap, length)

ann$promoter <- 0
ann$promoter[olap$subjectHits] <- olap$queryHits

#### H3K4me1 ####
h3k4me1.locus <- read.delim('LV.h3k4me1.locus.bed', stringsAsFactors = F, sep = ' ')
# h3k4me1.locus <- read.delim('RV.h3k4me1.locus.bed', stringsAsFactors = F, sep = ' ')
# h3k4me1.locus <- read.delim('RA.h3k4me1.locus.bed', stringsAsFactors = F, sep = ' ')

h3k4me1.g <- GRanges(h3k4me1.locus$chr, IRanges(as.numeric(h3k4me1.locus$start), as.numeric(h3k4me1.locus$end)), "*")
expansion.g <- GRanges(all$chr, IRanges(all$start, all$end), "*")
olap <- data.frame(findOverlaps(h3k4me1.g, expansion.g))# 
olap <- aggregate(queryHits ~ subjectHits, olap, length)

ann$h3k4me1 <- 0
ann$h3k4me1[olap$subjectHits] <- olap$queryHits

#### H3K27ac ####
h3k27ac.locus <- read.delim('LV.h3k27ac.locus.bed', stringsAsFactors = F, sep = ' ')
# h3k27ac.locus <- read.delim('RV.h3k27ac.locus.bed', stringsAsFactors = F, sep = ' ')
# h3k27ac.locus <- read.delim('RA.h3k27ac.locus.bed', stringsAsFactors = F, sep = ' ')

h3k27ac.g <- GRanges(h3k27ac.locus$chr, IRanges(as.numeric(h3k27ac.locus$start), as.numeric(h3k27ac.locus$end)), "*")
olap <- data.frame(findOverlaps(h3k27ac.g, expansion.g))
olap <- aggregate(queryHits ~ subjectHits, olap, length)

ann$h3k27ac <- 0
ann$h3k27ac[olap$subjectHits] <- olap$queryHits

#### H3K27me3 ####
h3k27me3.locus <- read.delim('LV.h3k27me3.locus.bed', stringsAsFactors = F, sep = ' ')
# h3k27me3.locus <- read.delim('RV.h3k27me3.locus.bed', stringsAsFactors = F, sep = ' ')
# h3k27me3.locus <- read.delim('RA.h3k27me3.locus.bed', stringsAsFactors = F, sep = ' ')

h3k27me3.g <- GRanges(h3k27me3.locus$chr, IRanges(as.numeric(h3k27me3.locus$start), as.numeric(h3k27me3.locus$end)), "*")
olap <- data.frame(findOverlaps(h3k27me3.g, expansion.g))# 
olap <- aggregate(queryHits ~ subjectHits, olap, length)

ann$h3k27me3 <- 0
ann$h3k27me3[olap$subjectHits] <- olap$queryHits

#### Prepare sample info table
sampleinfo <- read.delim("samples.with.TRcount.txt", stringsAsFactors = F)
outliers <- readLines("samples.with.exceed.trs.txt")
sampleinfo <- sampleinfo[(!sampleinfo$Var1 %in% outliers) & (sampleinfo$cohort != "1000G"), ]
names(sampleinfo) <- c("sample", "trs", "cohort")
sampleinfo = sampleinfo[sampleinfo$sample %in% clean.sample,]
sampleinfo$status = ifelse(sampleinfo$cohort == 'CMP', 1,0)

start <- Sys.time()
tmp <- do.call(rbind, lapply(1:nrow(dt), getOutlierMergeData, dt, sampleinfo))
end <- Sys.time()
difftime(end, start)
tmp <- tmp[, -2]
tmp$varid <- paste(tmp$chr, tmp$start, tmp$end, sep="#")

#ehdn.genotype.data <- tmp$sample[tmp$varid %in% ann$varid[ann$typeseq_priority == "intergenic"]]
# If intergenic result in global burden is significant, change to this

ehdn.genotype.data <- tmp$sample[tmp$varid %in% ann$varid]
ehdn.genotype.data <- data.frame(table(ehdn.genotype.data), stringsAsFactors = F)
names(ehdn.genotype.data) <- c("sample", "ncloci")
ehdn.genotype.data$dummy <- 1
cov <- c("ncloci")

genic.count <- tmp$sample[tmp$varid %in% ann$varid[ann$typeseq_priority != "intergenic"]]
genic.count <- data.frame(table(genic.count), stringsAsFactors = F) 
names(genic.count) <- c("sample", "ncloci") 

feat <- c("promoter","h3k27ac","h3k4me1","h3k27me3")
promoter.burden <- data.frame("mark" = "promoter", 
                              testByAggregate(dt$repeatID[dt$varid %in% ann$varid[ann[,'promoter'] > 0]], 
                                              tmp, sampleinfo, ehdn.genotype.data, cov))
for(f in feat[2:4]){
  tmp.burden <- data.frame("mark" = f, testByAggregate(dt$repeatID[dt$varid %in% ann$varid[ann[, f] > 0]],
                                                       tmp, sampleinfo, ehdn.genotype.data, cov))
  promoter.burden <- rbind(promoter.burden, tmp.burden)
}

promoter.plot = promoter.burden

write.table(promoter.plot, 'promoter.burden.all.txt', quote = F, sep = '\t')
# write.table(promoter.plot, 'promoter.burden.eur.txt', quote = F) 
# write.table(promoter.plot, 'promoter.burden.non-eur.txt', quote = F)


#### Plot CA CMP and UK CMP promoter and epigenetic marks burden
prom.burden.all.ca = read.delim('promoter.burden.all.txt', stringsAsFactors = F, sep = '\t')
prom.burden.all.ca = prom.burden.all.sk[,c(1:10)]
prom.burden.all.ca$col = 'dodgerblue3'
prom.burden.all.ca$mark[2:4] = c('H3K27ac', 'H3K4me1', 'H3K27me3')
rownames(prom.burden.all.ca) = paste0('ca', c(1:length(rownames(prom.burden.all.ca))))

prom.burden.all.ge = read.delim('promoter.burden.all.ge.tsv', stringsAsFactors = F, sep = '\t')
prom.burden.all.ge$col = 'lightskyblue1'
prom.burden.all.ge$mark[2:4] = c('H3K27ac', 'H3K4me1', 'H3K27me3')
rownames(prom.burden.all.ge) = paste0('ge', c(1:length(rownames(prom.burden.all.ge))))

p.numbers = c(rbind(1:4, 5:8))
promoter.plot = rbind(prom.burden.all.ca, prom.burden.all.ge)[p.numbers,]

pdf('promoter.burden.all.SK.GE.LV.pdf', w = 4, h = 2.5)
par(mar = c(4,4,1,5), mgp = c(3,1,0), mfrow=c(1,1)) 
b = barplot(height = promoter.plot$OR, space = c(rep(c(1,0.2),4)) ,
            col = alpha(promoter.plot$col, alpha = 0.8),
            yaxt = 'n', ylim = c(0,12), ylab = 'Odds ratio')
arrows(b, promoter.plot$lower, b, promoter.plot$upper, col='black', lwd = 1, code = 3, angle = 90, length = 0.05)
abline(h = 1, lty = 2)
axis(2, at = c(0,1,5,12), las = 2)
text(c(0.7:4)*3.2, par("usr")[3]-0.5, srt = 45, adj = 1, xpd = TRUE, labels = paste(promoter.plot$mark)[seq(1,length(promoter.plot$mark), 2)])
# format(promoter.plot$perm_p[3], scientific = T, digits = 3)# "3.7e-02"
text(x=4.25,y=7.8,expression(3.7~x~10^-02),srt=90, cex = 0.8)
text(x=4.5,y=4.5,paste('p = '),srt=90, cex = 0.8)

# format(promoter.plot$perm_p[4], scientific = T, digits = 3)# "1.1e-02"
text(x=5.45,y=8.2,expression(1.1~x~10^-02),srt=90, cex = 0.8)
text(x=5.7,y=5,paste('p = '),srt=90, cex = 0.8)# "1.1e-02"
legend(x = "topright", legend = c('CA CMP','UK CMP'), fill = c(alpha('dodgerblue3', 0.8), alpha('lightskyblue1', 0.8)), 
       inset=c(-0,-0.2), bty = "n", xpd = TRUE)
dev.off()
