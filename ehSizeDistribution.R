library(scales)
library(readr)

###### DIP2B genotype and pedigrees
dip2b_eh <- read.delim('dip2b.eh.cmp.tsv', sep = ' ', stringsAsFactors = F)   
allele_freq <- dip2b_eh[which (dip2b_eh$chr12.50505001.50505022.GGC != './.' ),]
all_freq <- as.data.frame(table(allele_freq$chr12.50505001.50505022.GGC))

for (i in 1:length(all_freq$Var1)){
  if(all_freq$Freq[i] == 1 ){
    all_freq$name[i] <- allele_freq$Sample[which (allele_freq$chr12.50505001.50505022.GGC == all_freq$Var1[i])]
  }else{all_freq$name[i] = ''}
}

for (i in 1:length(all_freq$Var1)){
  all_freq$allele1[i] <- as.numeric(strsplit(all_freq$Var1[i], '/',fixed = TRUE)[[1]][1])
  all_freq$allele2[i] <- as.numeric(strsplit(all_freq$Var1[i], '/',fixed = TRUE)[[1]][2])
}

all_freq$size <- 1.75
all_freq$size[all_freq$Freq <= 50] <- 1.5
all_freq$size[all_freq$Freq <= 20] <- 1.25
all_freq$size[all_freq$Freq <= 10] <- 1.1
all_freq$size[all_freq$Freq == 1] <- 1

all_freq$color <- 'grey'
all_freq$relationship <- ''
all_freq$col <- 'black'

all_freq$color[all_freq$name == 'Family-1_PROBAND'] <- alpha('firebrick', 0.4)
all_freq$color[all_freq$name == 'Family-1_MOTHER'] <- alpha('firebrick', 0.4)
all_freq$color[all_freq$name == 'Family-1_HALF_BROTHER'] <- alpha('firebrick', 0.4)

all_freq$relationship[all_freq$name == 'Family-1_PROBAND'] <- 'PROBAND'
all_freq$relationship[all_freq$name == 'Family-1_MOTHER'] <- 'MOTHER'
all_freq$relationship[all_freq$name == 'Family-1_HALF_BROTHER'] <- 'HALF_BROTHER'

all_freq$color[all_freq$name == 'Family-2_PROBAND'] <- alpha('darkorange', 0.4)
all_freq$color[all_freq$name == 'Family-2_FATHER'] <- alpha('darkorange', 0.4)

all_freq$relationship[all_freq$name == 'Family-2_PROBAND'] <- 'PROBAND'
all_freq$relationship[all_freq$name == 'Family-2_FATHER'] <- 'FATHER'

all_freq$color[all_freq$name == 'Family-3_PROBAND'] <- 'yellow'
all_freq$relationship[all_freq$name == 'Family-3_PROBAND'] <- 'PROBAND'

all_freq$x <- as.numeric(all_freq$allele1)
all_freq$y <- as.numeric(all_freq$allele2)

### PLOT 
pdf('dip2b_CCG_label__.pdf', w = 5.5, h = 4.5)
par(mar = c(8,6,3,7), mgp = c(2.5,1,0))
plot(all_freq$allele1, all_freq$allele2, col = 'black', bg = all_freq$color, pch = all_freq$pch, cex = all_freq$size,
     xlab = 'Smaller allele', ylab = 'Larger allele', yaxt = 'n', xaxt = 'n', xlim = c(5,25), ylim = c(0,145),
     main = 'DIP2B CGG')
axis(2, at=seq(0, 140, by=20), las = 2)
axis(1, at=c(7,10,15,20,25), las = 1)
abline(h = 7, lty = 2, lwd = 0.5)
abline(v = 7, lty = 2, lwd = 0.5)

text(all_freq$x, all_freq$y+5,paste0(all_freq$relationship), cex = 0.7, col = 'black')
legend('topright', legend = paste0(c('1','10','20','50','> 50')),col='black',pt.bg = 'grey',pch=21,pt.cex=c(unique(all_freq$size)),
       cex = 0.8,title = 'Frequency', xpd = TRUE, inset = c(-0.3,-0), bty = 'n')
legend('topright', legend = paste0(c('PROBAND with large expansion')),col='black',
       pt.bg = c('grey'),pch=c(23),pt.cex=1,title = '',
       cex = 0.8,xpd = TRUE, inset = c(-0,-0.2), bty = 'n')

legend('bottomright', legend = paste0(c('Family-1', 'Family-2','Family-3')),col='black',
       pt.bg = c('firebrick', 'darkorange', 'yellow'),pch=21,pt.cex=1,
       cex = 0.8,title = '', xpd = TRUE, inset = c(-0.3,0.3), bty = 'n')
dev.off()



size <- read.delim('Table.3-Table.3: Individuals with DIP2B expansions.csv', 
                  stringsAsFactors = F, sep = ',')[-c(7,9),-c(5,6)]
size$Participant <- toupper(size$Participant)

size$smaller.eh <- as.numeric(sapply(strsplit(size$EH...repeat.size, '/'), '[', 1))
size$larger.eh <- as.numeric(sapply(strsplit(size$EH...repeat.size, '/'), '[', 2))

size$smaller.np <- as.numeric(sapply(strsplit(size$Nanopore.repeat.size, '/'), '[', 1))
size$larger.np <- as.numeric(sapply(strsplit(size$Nanopore.repeat.size, '/'), '[', 2))

size$color <- unlist(sapply(size$Family, switch,'Family-1'='firebrick', 'Family-2'='darkorange', 'Family-3'='yellow', 'Family-4'='grey'))
size$pch <- ifelse(size$Participant == 'PROBAND', 23,21)
size$alpha <- 1
size$col <- 'black'

pdf('eh.nanopore.allele.pdf', w = 6, h = 3)#, h = 7
par(mar = c(5,4,3,3), mgp = c(3,1,0), mfrow=c(1,2)) 
cor.test.sm <- cor.test(size$smaller.np, size$smaller.eh)
plot(size$smaller.np, size$smaller.eh, xlab = 'LRS allele size', ylab = 'SRS allele size', col = size$col,
     bg = alpha(size$color, size$alpha), xaxt = 'n', yaxt = 'n', pch = size$pch, xlim = c(5,25), ylim = c(5,25), cex = 1.5)
axis(1, at = c(7,15, 25))
axis(2, at = c(7,15, 25), las=2)
abline(lm(size$smaller.eh ~ size$smaller.np), col = 'black', lty = 2)

legend('topright', legend = paste0('r = ', round(cor.test.sm$estimate, 2) ),
       bty = "n", inset = c(0, -0.4), xpd = TRUE)
# format(cor.test.sm$p.value, scientific = T, digits = 2)  # "2.2e-02"
text(x=22, y=28.5, expression(2.2~x~10^-02), bty = "n", xpd = TRUE)
text(x=15, y=27.55, paste('p = '), bty = "n", xpd = TRUE) 

cor.test.lg <- cor.test(size$larger.np, size$larger.eh)
plot(size$larger.np, size$larger.eh, xlab = 'LRS allele size', ylab = 'SRS allele size', cex = 1.5,
     bg = alpha(size$color, size$alpha), xaxt = 'n', yaxt = 'n', pch = size$pch, xlim = c(0,750), ylim = c(0,700))
axis(1, at = c(0,300,700))
axis(2, at = c(0,300,700), las = 2)
abline(lm(size$larger.eh ~ size$larger.np), col = 'black', lty = 2)

legend('topright', legend = paste0('r = ', round(cor.test.lg$estimate, 2)),
       bty = "n", inset = c(0, -0.4), xpd = TRUE)
legend('topright', legend = paste0('p = ', round(cor.test.lg$p.value, digits = 2)),
       bty = "n", inset = c(0, -0.25), xpd = TRUE)
dev.off()
