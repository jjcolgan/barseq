library(ensembldb)
library(biomartr)

meta <- meta.retrieval(kingdom = "bacteria", db = "refseq")
bifido <- meta[grep("Bifidobacterium breve", meta$organism), ]
View(bifido)

bBreveDb <- ensDbFromGtf(gtf = "/Users/johnjamescolgan/Downloads/bBreve/ncbi_dataset/data/GCF_000220135.1/genomic.gtf",
                         organism = 'Bifidobacterium breve UCC2003')
?ensDbFromGtf

txdb <- makeTxDbFromGFF("/Users/johnjamescolgan/Downloads/bBreve/ncbi_dataset/data/GCF_000220135.1/genomic.gtf", format = "gtf")
genes(txdb)