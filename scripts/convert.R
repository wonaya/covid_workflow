library(ape)
library(seqinr)
f<-read.fasta("orf_6_aligned_rmstop.fasta")
write.dna(f,"orf_6_aligned_rmstop.phy", nbcol=1,colsep="", colw=1000000)
