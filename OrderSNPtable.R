#!/usr/bin/Rscript --vanilla

#get command line arguments
args <- commandArgs(trailingOnly = TRUE)

snptable <- read.table(args[1])
snptable_ordered <- unique(snptable[order(snptable[,3]),])
filename <- paste('ordered_',args[1], sep = "")
write.table(snptable_ordered, filename, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

