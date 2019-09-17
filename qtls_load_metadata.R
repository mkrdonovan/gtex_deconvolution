message("Loading metadata...")

metadata           = read.table("input/metadata/metadata.txt"            , header = TRUE , check.names = FALSE)
covariates_subject = read.table("input/covariates/covariates.subject.txt", header = TRUE , check.names = FALSE)
covariates_rna     = read.table("input/covariates/covariates.rna.txt"    , header = TRUE , check.names = FALSE)
rna_list           = read.table("input/phenotypes/rna_list.txt"          , header = FALSE, check.names = FALSE)$V1
unrelated          = read.table("input/metadata/unrelated.txt"           , header = FALSE, check.names = FALSE)$V1
