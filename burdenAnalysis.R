library(GenomicRanges)
library(data.table)
library(plotrix)
library(readr)

ann <- fread("merged.expansions.forannotation.tsv", data.table = F)
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
ann <- ann[ann$varid %in% dt$varid, ] 

sampleinfo <- read.delim("samples.with.TRcount.txt", stringsAsFactors = F)
outliers <- readLines("samples.with.exceed.trs.txt")
sampleinfo <- sampleinfo[(!sampleinfo$Var1 %in% outliers) & (sampleinfo$cohort != "1000G"), ]
names(sampleinfo) <- c("sample", "trs", "cohort")
sampleinfo = sampleinfo[sampleinfo$sample %in% clean.sample,]
sampleinfo = sampleinfo[!duplicated(sampleinfo$sample),]
sampleinfo$status <- ifelse(sampleinfo$cohort == "CMP", 1, 0)

start <- Sys.time()
tmp <- do.call(rbind, lapply(1:nrow(dt), getOutlierMergeData, dt, sampleinfo))
end <- Sys.time()
difftime(end, start)
tmp <- tmp[, -2]

tmp$varid <- paste(tmp$chr, tmp$start, tmp$end, sep="#")
ehdn.genotype.data <- tmp$sample[tmp$varid %in% ann$varid]
ehdn.genotype.data <- data.frame(table(ehdn.genotype.data), stringsAsFactors = F)
names(ehdn.genotype.data) <- c("sample", "ncloci")

cov <- c("ncloci")

genic.burden <- data.frame("gs" = "genic", testByAggregate(dt$repeatID[dt$varid %in% ann$varid[ann$typeseq_priority != "intergenic"]], tmp, sampleinfo, ehdn.genotype.data, cov))

for(f in 1:length(feat)){
  if(sum(tmp$varid %in% ann$varid[ann$typeseq_priority %in% feat[f]]) > 0){
    tmp.burden <- data.frame("gs" = names(feat)[f], 
                             testByAggregate(dt$repeatID[dt$varid %in% ann$varid[ann$typeseq_priority %in% feat[[f]]]], 
                                             tmp, sampleinfo, ehdn.genotype.data, cov))
    # for significant intergenic situation, remove the else condition
    
    #else
    #tmp.burden <- data.frame("gs" = names(feat)[f], 
    # testByAggregate(dt$repeatID[dt$varid %in% ann$varid[ann$typeseq_priority %in% feat[[f]]]], 
    #  tmp, sampleinfo, ehdn.genotype.data, "dummy"),
    #fisher(dt$outliers[dt$varid %in% ann$varid[ann$typeseq_priority %in% feat[[f]]]], 
    #genic.count, sampleinfo))
    genic.burden <- rbind(genic.burden, tmp.burden)
  }
}

# splicing
for(f in c(2)){
  tmp.burden <- data.frame("gs" = names(feat)[f],
                           testByAggregate(dt$repeatID[dt$varid %in% ann$varid[ann$typeseq_priority %in% feat[[f]][1]]], 
                                           tmp, sampleinfo, ehdn.genotype.data, cov))
}
genic.burden <- rbind(genic.burden, tmp.burden)

