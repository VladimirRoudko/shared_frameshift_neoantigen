#!/usr/bin/env Rscript

library("biomaRt")

args = commandArgs(trailingOnly=TRUE)

mart <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="dec2016.archive.ensembl.org")
data <- read.table(args[1], sep=" ")
ENSTID <- t(data$V6)
ENSTID <- unique(ENSTID)
seq = getSequence(id=ENSTID, type="ensembl_transcript_id", seqType=args[2], mart = mart, verbose=FALSE)
exportFASTA(seq, args[3])
coord <- getBM(attributes = c("ensembl_transcript_id", "strand"), filters="ensembl_transcript_id", values=ENSTID, mart=mart)
write.table(coord, "frameshift.strand.txt", sep="\t")


# library("biomaRt")
# mart <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="dec2016.archive.ensembl.org")
# seq = getSequence(id="ENST00000444749", type="ensembl_transcript_id", seqType="coding", mart = mart, verbose=FALSE)
# exportFASTA(seq, "filename")
