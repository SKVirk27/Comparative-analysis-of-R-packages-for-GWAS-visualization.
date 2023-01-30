# Project 
#Comparative analysis of R packages for GWAS Visualization

library(qqman)
library(data.table)
library(readr)
library(fastman)
library(tictoc)
setwd("~/Desktop")
#load the sumstats
data= fread("daner_adhd_meta_filtered_NA_iPSYCH23_PGC11_sigPCs_woSEX_2ell6sd_EUR_Neff_70 (1).meta")
head(data)
tic();
png("fastman.png", width=10, height=6, units="in", res=300)
#Manhattan plot
fastman(data,chr = 'CHR',ps='BP',p="P",color=c("darkslateblue","gold2","#ff5733"))
dev.off()
toc();# 25.843 sec elapsed
#QQ plot
tic();
png("QQman.png", width=10, height=6, units="in", res=300)
manhattan(data, main = "Manhattan Plot", ylim = c(0, 15), cex = 0.6, cex.axis = 0.9, 
          col = c("blue4", "orange3","darkslateblue") )
dev.off()
toc();  #124.929 sec

#QQplot by Fastmaan
tic();
png("Fastmanplot.png", width=10, height=6, units="in", res=300)
par(pty="s")
fastqq(data)
dev.off()
toc();#35.054 sec elapsed
#QQplot by qqman
tic();
png("QQman_plot.png", width=10, height=6, units="in", res=300)
qq(data$P, main = "Q-Q plot of GWAS p-values", xlim = c(0, 7), 
   ylim = c(0, 12), pch = 18, col = "blue4", cex = 1.5, las = 1)
dev.off()
toc();#187.906 sec
#QQplot by Qqman

# Chromosome3 of the data.
png("QQman_ch3.png", width=10, height=6, units="in", res=300)
manhattan(subset(data,CHR==3), main = "Manhattan Plot", ylim = c(0, 10), cex = 0.6, cex.axis = 0.9, 
          col = c("blue4", "orange3","darkslateblue") )
dev.off()

#top SNP per chromosome that exceeds the annotatePval threshold.
png("QQman_topsnp.png", width=10, height=6, units="in", res=300)
manhattan(data, main = "Manhattan Plot",annotatePval = 0.01, ylim = c(0, 15),col = c("blue4", "orange3","darkslateblue") )
dev.off()



