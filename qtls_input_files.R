message("Loading input files...")

geneinfo   = read.table("/publicdata/gencode_v19_20151104/gene_info.tsv"      , header = TRUE, sep = "\t")
chromsizes = read.table("/publicdata/gatk_bundle_2.8/hg19/ucsc.hg19.fasta.fai", header = FALSE, col.names = c("chrom", "size", "start", "na1" , "na2"))[,1:3]


