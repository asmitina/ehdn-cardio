library(scales)
library(readr)

### Constraint scores
constraint <- read.delim('gnomAD_pLI_v2.1_Feb8.txt', stringsAsFactors = F)# 80953
constraint <- constraint[constraint$canonical == 'true',]# 19705

rare.ca.cmp <- readLines('all.genic.ca.cmp.127.txt')
const.cmp.exp <- constraint$oe_lof_upper[constraint$gene %in% rare.ca.cmp]# 127
const.all.genes <- constraint$oe_lof_upper[!constraint$gene %in% rare.ca.cmp]# 19588

rare.ge.cmp <- readLines('all.genic.ge.cmp.142.txt')
const.ge.exp <- constraint$oe_lof_upper[constraint$gene %in% rare.ge.cmp]# 134
const.all.genes.ge <- constraint$oe_lof_upper[!constraint$gene %in% const.ge.exp]# 19705

test.sk <- wilcox.test(const.cmp.exp, const.all.genes, alternative = 'less')
test.ge <- wilcox.test(const.ge.exp, const.all.genes.ge, alternative = 'less')

pdf('boxplot.gnomAD.pdf', w = 3.2, h = 4)
par(mar = c(8,4,3,1), mgp = c(3,1,0), mfrow=c(1,1)) 

boxplot(at = c(1,2,4,5), const.cmp.exp, const.all.genes, 
        const.ge.exp, const.all.genes.ge, xaxt = 'n', yaxt = 'n', 
        ylab = 'Constraint score', col = c(alpha('dodgerblue3', 0.8), alpha('snow',0.8), 
                                           alpha('lightskyblue1', 0.8), alpha('snow',0.8)))
axis(1, at = c(1,2,4,5), labels = NA)
text(c(1,2,4,5), par("usr")[3]-0.2, srt = 45, adj = 1, xpd = TRUE, 
     labels = c('with TRE in CA', 'without TRE','with TRE in UK', 'without TRE'))
axis(2, at = seq(0,2,by = 0.5), las = 2)

# format(test.sk$p.value, scientific = T, digits = 2)# "2.4e-07"
text(x=2.0,y=2.22,expression(2.4~x~10^-7),srt=0, cex = 1, xpd = T)
text(x=0.75,y=2.15,paste('p = '),srt=0, cex = 1, xpd = T)

# format(test.ge$p.value, scientific = T, digits = 2)# 1e-14
text(x=4.7,y=2.22,expression(1~x~10^-14),srt=0, cex = 1, xpd = T)
text(x=3.6,y=2.15,paste('p = '),srt=0, cex = 1, xpd = T)
dev.off()
