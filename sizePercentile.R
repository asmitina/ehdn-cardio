library(GenomicRanges)
library(data.table)
library(stringr)
library(plotrix)
library(readr)

kgp <- read.delim('1000G_metadata.tsv', stringsAsFactors = F)

#### Rare expansions 
dt <- read.delim("merged.rare.expansions.tsv", stringsAsFactors = F)
dt$repeatID <- paste(dt$chr, dt$start, dt$end, dt$motif, sep="#")
dt$varid <- paste(dt$chr, dt$start, dt$end, sep="#")

sampleinfo <- read.delim("samples.with.TRcount.txt", stringsAsFactors = F)
outliers <- readLines("samples.with.exceed.trs.txt")
sampleinfo <- sampleinfo[(!sampleinfo$Var1 %in% outliers) & (sampleinfo$cohort != "1000G"), ]
names(sampleinfo) <- c("sample", "trs", "cohort")
sampleinfo = sampleinfo[sampleinfo$sample %in% clean.sample,]
sampleinfo = sampleinfo[!duplicated(sampleinfo$sample),]
sampleinfo$status <- ifelse(sampleinfo$cohort == "CMP", 1, 0)

for (i in 1:nrow(dt)){ 
  samples = strsplit(dt$outliers[i], ';')[[1]]
  dt$cmp_clean[i] = paste0(samples[samples %in% sampleinfo$sample[sampleinfo$cohort == 'CMP']], collapse = ';')# 247
  dt$control_clean[i] = paste0(samples[samples %in% sampleinfo$sample[sampleinfo$cohort != 'CMP']], collapse = ':')# 386
  
  dt$cmp_cl.count[i] = length(samples[samples %in% sampleinfo$sample[sampleinfo$cohort == 'CMP']])
  dt$control_cl.count[i] = length(samples[samples %in% sampleinfo$sample[sampleinfo$cohort != 'CMP']])
}

dt.cmp <- dt[dt$cmp_cl.count >0,]
dt.cmp.g <- GRanges(dt.cmp$chr, IRanges(as.numeric(dt.cmp$start), as.numeric(dt.cmp$end)), "*")

#### All EHdn expansions
exp.cmp  <- read.delim("output_regions.min2.txt", stringsAsFactors = F, header = F)# 46642
colnames(exp.cmp)[c(1:4,7)] <- c('chr','start','end','motif','outliers')
exp.cmp <- exp.cmp[exp.cmp$chr %in% paste0("chr", c(1:22)), ]# 41477 restrict to autosomes
for(i in 1:nrow(exp.cmp)){
  outliers <- strsplit(gsub(",", ":", exp.cmp$outliers[i]), ":")[[1]]
  samples <- outliers[seq(1, length(outliers), 2)]
  size  <- outliers[seq(2, length(outliers), 2)]
  exp.cmp$samples[i] <- paste(samples[grep('CMP',samples)], collapse = ';')
  exp.cmp$size[i] <- paste(size[grep('CMP',samples)], collapse = ';')
  exp.cmp$freq[i] <- sum(samples %in% kgp$Sample.ID)/nrow(kgp)*100
}

exp.cmp = exp.cmp[exp.cmp$freq <= 0.2,]
ehdn.g <- GRanges(exp.cmp$chr, IRanges(exp.cmp$start, exp.cmp$end), "*")
olap <- data.frame(findOverlaps(all.cmp.g, ehdn.g))

for (i in unique(olap$queryHits)){
  ehdn.tmp <- exp.cmp[olap$subjectHits[olap$queryHits == i],]
  ehdn.table <- data.frame('samples' = c(), 'size.ehdn' = c())
  for (line in 1:nrow(ehdn.tmp)){
    table.tmp <- data.frame('samples' = strsplit(ehdn.tmp$samples[line], ';')[[1]], 
                           'size.ehdn' = strsplit(ehdn.tmp$size[line], ';')[[1]])
    ehdn.table <- rbind(ehdn.table, table.tmp)
    ehdn.table$samples <- as.character(ehdn.table$samples)
    ehdn.table$size.ehdn <- as.numeric(as.character(ehdn.table$size.ehdn))
  }
  ehdn.table <- ehdn.table[order(ehdn.table$size.ehdn, decreasing = T),]
  ehdn.table <- ehdn.table[!duplicated(ehdn.table$samples),]
  
  for (one in 1:nrow(ehdn.table)){
    if (ehdn.table$samples[one] %in% cmp.meta$Sample){
      ehdn.table$relation[one] <- cmp.meta$Relationship[cmp.meta$Sample == ehdn.table$samples[one]] 
      ehdn.table$family[one] <- cmp.meta$Family[cmp.meta$Sample == ehdn.table$samples[one]]
      ehdn.table$outlier[one] <- ifelse(ehdn.table$samples[one] %in% strsplit(all.cmp$cmp_clean[i], ';')[[1]], T, F)
    }
  }
  if (length(ehdn.table$samples[ehdn.table$family %in% ehdn.table$family[ehdn.table$outlier == T]]) > 1){
    print(paste(count, all.cmp$varid[i], collapse = ' '))
    write.table(ehdn.table, paste0('_ehdn_percentile/', count, '_', all.cmp$varid[i], '.txt'),
                quote = F, col.names = T, row.names = F, sep = '\t')
    count <- count + 1
  }
}

outpath <- '_final_results/_ehdn_percentile/'
families <- data.frame()

#### Calculate and plot percentile vs EHdn size 
pdf('ehdn.percentile_0.2_cutoff.pdf', h = 7, w = 10)
par(mar = c(4,5,3,1), mgp = c(3,1,0), mfrow=c(3,4))
for (file in sort(list.files(outpath))){
  name <- strsplit(gsub('.txt', '', file), '_')[[1]][2]
  size.tmp <- read.delim(sprintf("%s%s", outpath, file), header = T, stringsAsFactors = F)
  no.geno <- nrow(cmp.meta) - nrow(size.tmp)
  size.tmp$perc <- (sapply(lapply(size.tmp$size.ehdn, ">", size.tmp$size.ehdn), sum) + no.geno)/(nrow(cmp.meta) - 1)
  size.tmp$short.rel <- sapply(size.tmp$relation, switch, SISTER='sibling', BROTHER='sibling', 'HALF-BROTHER' = 'sibling',
                              'MOTHER'='parent', 'FATHER'='parent', 'PROBAND - FATHER'='parent',
                              'PROBAND - MOTHER'='parent', 'PROBAND-MOTHER'='parent', PROBAND='proband',INDEX='proband')
  size.tmp <- size.tmp[size.tmp$short.rel != 'NULL',]
  size.tmp$pch <- lapply(size.tmp$short.rel, switch, proband=21,sibling=22,parent=24,'NULL'='snow')
  size.tmp$bg <- ifelse(size.tmp$outlier, 'orangered1', 'grey58')
  size.tmp$bg <- unlist(size.tmp$bg)
  size.tmp$pch <- unlist(size.tmp$pch)
  size.tmp$short.rel <- unlist(size.tmp$short.rel)
  size.tmp$region <- name
  
  fam.tmp <- size.tmp[size.tmp$family %in% size.tmp$family[size.tmp$outlier == T],]
  
  plot(size.tmp$perc, size.tmp$size.ehdn, pch = size.tmp$pch, bg = size.tmp$bg, main = name, 
       ylab = 'EHdn size', xlab = 'percentile', col = 'black', cex = size.tmp$cex)
  
  families <- rbind(families, fam.tmp)
}
dev.off()


families.selected <- data.frame()
for (one in unique(families$region)){
  fam.tmp <- families[families$region == one,]
  for (fam in unique(fam.tmp$family)){
    if ( nrow(fam.tmp[fam.tmp$family == fam,]) >= 2 ){
      families.selected <- rbind(families.selected, fam.tmp[fam.tmp$family == fam,])
    }
  }
}

families.selected$percentile <- families.selected$perc*100
write.table(families.selected, 'families.selected.tsv', quote = F, col.names = T, row.names = F, sep = '\t')


families.selected <- read.delim('families.selected.tsv', stringsAsFactors = F, sep = '\t')
families.selected <- families.selected[families.selected$short.rel != 'sibling',]
families.selected$short.rel <- str_to_title(families.selected$short.re)

#### Boxplot for percentile in proband and parent 
pdf('boxplot.perc.ehdn.pdf', h = 3.2, w = 5.5)
par(mar = c(3,5,3,1), mgp = c(3,1,0), mfrow=c(1,2))

boxplot(families.selected$perc ~ families.selected$short.rel, yaxt = 'n',
        col = c('honeydew', 'orangered1'),main = '', xlab = '', ylab = 'Percentile')
axis(2, at = seq(88,100,4), las = 2)

boxplot(families.selected$perc ~ families.selected$short.rel, yaxt = 'n', 
        col = c('honeydew', 'orangered1'), main = '', xlab = '', 
        ylab = 'Percentile', ylim = c(99,100))
axis(2, at = seq(0,100,1), las = 2)
dev.off()
