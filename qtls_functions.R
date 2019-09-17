### Version history:
# V05: add conditional QTLs
# V04: initialize QTL analysis: create folders and links for all analyses, so that we can run multiple QTL analyses in different folders
# V03: from Alasoo, NatGenet 2018 (PMID: 29379200), run LMM, get p-value for each variant/gene pair, then run all other LMMs with interactions and compare with ANOVA (model is better if AIC is smaller)
# V02: changed the linear model for QTLs: first pass, just run a linear model without interactions

message("Loading functions...")

### Functions to read database data

getTableFromDb = function(query)
{    
    con = dbConnect(MySQL(), user="cardips", password="G00dC@kes", dbname="cardips", host="10.0.16.10")# works from compute nodes

    rs <- dbSendQuery(con, query)
    db_table <- fetch(rs, n = -1)

    dbClearResult(rs)
    dbDisconnect(con)
    
    return(db_table)
}

addDashUUID = function(uuid)
{
    splitted = unlist(strsplit(uuid, ""))

    out = paste(paste(splitted[ 1: 8], collapse = ""),
                paste(splitted[ 9:12], collapse = ""),
                paste(splitted[13:16], collapse = ""),
                paste(splitted[17:20], collapse = ""),
                paste(splitted[21:32], collapse = ""),
                sep = "-"
               )
    return(out)
}

transformUuid = function(uuid)
{
	uuid = gsub("-", ".", uuid)
	uuid[grepl("^[[:digit:]]", uuid, perl = TRUE)] = paste("X", uuid[grepl("^[[:digit:]]", uuid, perl = TRUE)], sep = "")
	
	return(uuid)
}

### Functions to handle VCF files:
### Filter VCF files and return a genotype matrix

filter_vcf = function(wgs_list, bcftools, snp_list_file, vcf_input, name)
{
    gttable_vcf = paste("private/gttable", name, "vcf", "gz", sep = ".")

    command1 = paste(bcftools, "view", "-R", snp_list_file, "-s", paste(wgs_list, collapse = ","), vcf_input, "-O z", "-o", gttable_vcf)
    command2 = paste(bcftools, "index", "--tbi", "-f", gttable_vcf)

    system(command1)
    system(command2)
    
    return(gttable_vcf)
}

find_genotype_matrix = function(wgs_list, bcftools, vcf_input)
{
    gttable_out = "private/gttable.txt" 
    command     = paste(bcftools, "query", "-H", "-s", paste(wgs_list, collapse = ","), vcf_input, "-f", '"%ID[\\t%GT]\\n"' , "-o", gttable_out)

    system(command)

    gttable           = fread(gttable_out, header = TRUE, sep = "\t", data.table = FALSE)
    colnames(gttable) = c("rsid", wgs_list)
    gttable           = gttable[gttable$rsid != ".",]
    gttable$rsid      = NULL
    gtmatrix          = as.matrix(gttable)
    gtmatrix          = gsub("\\.", "0", gtmatrix)
    gtmatrix          = gsub("2"  , "1", gtmatrix)
    gtmatrix          = gsub("3"  , "1", gtmatrix)
    gtmatrix          = gsub("4"  , "1", gtmatrix)
    gtmatrix          = gsub("5"  , "1", gtmatrix)
    gtmatrix          = suppressMessages(mapvalues(gtmatrix, from=c("0/0", "0/1", "1/0", "1/1"), to=c(0, 0.5, 0.5, 1)))
    class(gtmatrix)   = "numeric" 
    gtmatrix          = gtmatrix[unlist(apply(gtmatrix, 1, FUN = sd)) > 0,]
    
    return(as.data.frame(gtmatrix))
}


### Calculate kinship matrix
normalize_kinmat <- function(kinmat)
{
    #normalize kinship so that Kij \in [0,1]
    tmp=kinmat - min(kinmat)
    tmp=tmp/max(tmp)
    tmp[1:9,1:9]
    #fix eigenvalues to positive
    diag(tmp)=diag(tmp)-min(eigen(tmp)$values)
    tmp[1:9,1:9]  
    return(tmp)
}

calculate_kinship = function(vcf_file, write_output = FALSE)
{
    gds_file = sub("vcf", "gds", vcf_file)

    snpgdsVCF2GDS(vcf_file, gds_file, method = "biallelic.only")

    genofile          = snpgdsOpen  (gds_file)
    ibd               = snpgdsIBDMoM(genofile  , kinship = TRUE)
    kinship           = as.data.frame(ibd$kinship)
    colnames(kinship) = ibd$sample.id
    rownames(kinship) = ibd$sample.id

    snpgdsClose(genofile)
    
    if (write_output == TRUE)
    {
        write.table(kinship, "input/genotypes/kinship.txt", quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA)
    }

    # https://www.r-bloggers.com/fixing-non-positive-definite-correlation-matrices-using-r-2/

    origMat = kinship

    cholStatus <- try(u <- chol(origMat), silent = FALSE)
    cholError  <- ifelse(class(cholStatus) == "try-error", TRUE, FALSE)
    newMat     <- origMat

    iter <- 0
    while (cholError) 
    {
        iter <- iter + 1
        cat("iteration ", iter, "\n")

        # replace -ve eigen values with small +ve number
        newEig <- eigen(newMat)
        newEig2 <- ifelse(newEig$values < 0, 0, newEig$values)

        # create modified matrix eqn 5 from Brissette et al 2007, inv = transp for
        # eig vectors
        newMat <- newEig$vectors %*% diag(newEig2) %*% t(newEig$vectors)

        # normalize modified matrix eqn 6 from Brissette et al 2007
        newMat <- newMat/sqrt(diag(newMat) %*% t(diag(newMat)))

        # try chol again
        cholStatus <- try(u <- chol(newMat), silent = TRUE)
        cholError <- ifelse(class(cholStatus) == "try-error", TRUE, FALSE)
    }
    
    rownames(newMat) = rownames(kinship)
    colnames(newMat) = colnames(kinship)
    
    newMat = normalize_kinmat(newMat)
    
    return(newMat)
}

### Merge RNA-seq data and create one single TPM table (rows = genes; columns = RNA samples)

read_rna = function(ii, rnadata)
{
    id     = rnadata[ii, "rna_id"]
    infile = rnadata[ii, "file"  ]
    
    if (file.exists(infile))
    {   
        indata = read.table(infile, header = TRUE)[,c("gene_id", "TPM")]
        colnames(indata) = c("gene_id", id)
        return(indata)
    }
}

### Divide genotypes by gene
filter_vcf_by_region = function(wgs_list, bcftools, region, maf_threshold, vcf_input, outfile)
{
    txt     = sub("vcf.gz", "txt", outfile)
    command1 = paste(bcftools, "view", 
                     "-r", region, 
                     "-q", paste(maf_threshold, "minor", sep = ":"), 
                     "-s", paste(wgs_list, collapse = ","), 
                     vcf_input, 
                     "-O z", 
                     "-o", outfile)
    
    command2    = paste(bcftools, "index", "--tbi", "-f", outfile)
    vcf_to_text = paste(bcftools, "query", "-H", "-s", paste(wgs_list, collapse = ","), outfile, "-f", '"%CHROM\\t%POS\\t%REF\\t%ALT{0}\\t%ID[\\t%GT]\\n"', "-o", txt)
    
    system(command1)
    system(command2)
    system(vcf_to_text)
    return(txt)
}

read_gttable = function(infile, wgs_list)
{
    gttable           = read.table(infile, header = TRUE, sep = "\t")
    colnames(gttable) = c("chrom", "pos", "ref", "alt", "rsid",  wgs_list)
    gttable$id        = paste(gttable$chrom, gttable$pos, gttable$ref, gttable$alt, sep = "_") 
    rownames(gttable) = gttable$id
    return(gttable)
}

write_genotypes_rna = function(gene_id, outfolder, geneinfo_file)
{
    # to work with rparallel
    suppressMessages(source("analysis/cardiac_qtls_packages.R"      ))
    suppressMessages(source("analysis/cardiac_qtls_input_files.R"   ))
    suppressMessages(source("analysis/cardiac_qtls_functions.R"     ))
    suppressMessages(source("analysis/cardiac_qtls_input_data.R"    ))
    suppressMessages(source("analysis/cardiac_qtls_load_metadata.R" ))
    
    geneinfo = read.table(geneinfo_file, header = TRUE)
    
    # Get WGS sample lists
    wgs_list_ipscore = sort(unique(metadata[metadata$subject_id %in% covariates_subject[covariates_subject$study == "iPSCORE", "subject_id"], "wgs_id"]))
    wgs_list_gtex    = sort(unique(metadata[metadata$subject_id %in% covariates_subject[covariates_subject$study == "GTEx"   , "subject_id"], "wgs_id"]))

    # Get region coordinates
    chrom     = gsub("chr", "", geneinfo[geneinfo$gene_id == gene_id, "chrom"])
    gene_from =                 geneinfo[geneinfo$gene_id == gene_id, "start"]
    gene_to   =                 geneinfo[geneinfo$gene_id == gene_id, "end"  ]

    region_from = max(0                                            , gene_from - qtl_distance)
    region_to   = min(chromsizes[chromsizes$chrom == chrom, "size"], gene_to   + qtl_distance)

    region      = paste(chrom, ":", region_from, "-", region_to, sep = "")

    # Merge data from iPSCORE and GTEx and return the genotype tables
    tmp_vcf_ipscore     = paste("input//genotypes/rna/tmp", gene_id, "ipscore", "vcf", "gz", sep = ".")
    tmp_vcf_gtex        = paste("input//genotypes/rna/tmp", gene_id, "gtex"   , "vcf", "gz", sep = ".")
    tmp_vcf_merged      = paste("input//genotypes/rna/tmp", gene_id, "merged" , "vcf"      , sep = ".")
    tmp_gttable_ipscore = filter_vcf_by_region(wgs_list_ipscore, bcftools, region, maf_threshold, ipscore_vcf_input, tmp_vcf_ipscore)
    tmp_gttable_gtex    = filter_vcf_by_region(wgs_list_gtex   , bcftools, region, maf_threshold, gtex_vcf_input   , tmp_vcf_gtex   )

    gttable_ipscore = read_gttable(tmp_gttable_ipscore, wgs_list_ipscore)
    gttable_gtex    = read_gttable(tmp_gttable_gtex   , wgs_list_gtex   )
    common_snps     = intersect(gttable_ipscore$id, gttable_gtex$id)
    gt_info         = gttable_ipscore[common_snps, c("chrom", "pos", "ref", "alt", "rsid", "id")]
    gt_info         = gt_info[order(gt_info$chrom, gt_info$pos, gt_info$ref, gt_info$alt),]
    gttable         = cbind(gttable_ipscore[common_snps, wgs_list_ipscore], gttable_gtex[common_snps, wgs_list_gtex])
    gtmatrix        = as.matrix(gttable)
    gtmatrix        = gsub("\\.", "0", gtmatrix)
    gtmatrix        = suppressMessages(mapvalues(gtmatrix, from=c("0/0", "0/1", "1/0", "1/1"), to=c(0, 0.5, 0.5, 1)))
    
    if (length(gtmatrix) > 0)
    {
        class(gtmatrix) = "numeric" 
        gtmatrix        = gtmatrix[unlist(apply(gtmatrix, 1, FUN = sd)) > 0,]
        gttable         = as.data.frame(gtmatrix)

        write.table(gt_info , paste("input/genotypes/rna/gt_info", gene_id, "txt", sep = "."), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
        write.table(gtmatrix, paste("input/genotypes/rna/gt_data", gene_id, "txt", sep = "."), quote = FALSE, sep = "\t", row.names = TRUE , col.names = NA  )
    }

    file.remove(list.files(path = "input/genotypes/rna/", pattern = paste("tmp", gene_id, sep = "."), full.names = TRUE))
}

write_genotypes = function(ii, outfolder, geneinfo_file)
{
    source("analysis/cardiac_qtls_packages.R"      )
    source("analysis/cardiac_qtls_input_files.R"   )
    source("analysis/cardiac_qtls_functions.R"     )
    source("analysis/cardiac_qtls_input_data.R"    )
    source("analysis/cardiac_qtls_load_metadata.R" )
    
    geneinfo = read.table(geneinfo_file, header = TRUE)
    gene_id  = geneinfo[ii, "gene_id"]
    
    if (grepl("ENSG", gene_id) == TRUE ){qtl_distance = 1e6}
    if (grepl("ENSG", gene_id) == FALSE){qtl_distance = 1e5}
    
    # Get WGS sample lists
    wgs_list = sort(unique(metadata[metadata$subject_id %in% covariates_subject$subject_id, "wgs_id"]))

    # Get region coordinates
    chrom     = gsub("chr", "", geneinfo[geneinfo$gene_id == gene_id, "chrom"])
    gene_from =                 geneinfo[geneinfo$gene_id == gene_id, "start"]
    gene_to   =                 geneinfo[geneinfo$gene_id == gene_id, "end"  ]

    region_from = max(0                                            , gene_from - qtl_distance)
    region_to   = min(chromsizes[chromsizes$chrom == chrom, "size"], gene_to   + qtl_distance)
    region      = paste(chrom, ":", region_from, "-", region_to, sep = "")

    # Merge data from iPSCORE and GTEx and return the genotype tables
    tmp_vcf     = paste("input//genotypes/rna/tmp", gene_id, "gtex"   , "vcf", "gz", sep = ".")  
    tmp_gttable = filter_vcf_by_region(wgs_list, bcftools, region, maf_threshold, gtex_vcf_input, tmp_vcf)
    gttable     = read_gttable(tmp_gttable, wgs_list)
    common_snps = gttable$id
    gt_info     = gttable[common_snps, c("chrom", "pos", "ref", "alt", "rsid", "id")]
    gt_info     = gt_info[order(gt_info$chrom, gt_info$pos, gt_info$ref, gt_info$alt),]
    gtmatrix    = as.matrix(gttable[, wgs_list])
    gtmatrix    = gsub("\\.", "0", gtmatrix)
    gtmatrix    = suppressMessages(mapvalues(gtmatrix, from=c("0/0", "0/1", "1/0", "1/1"), to=c(0, 0.5, 0.5, 1)))
    
    if (length(gtmatrix) > 0)
    {
        class(gtmatrix) = "numeric" 
        gtmatrix        = gtmatrix[unlist(apply(gtmatrix, 1, FUN = sd)) > 0,]
        gttable         = as.data.frame(gtmatrix)

        write.table(gt_info , paste(outfolder, paste("gt_info", gene_id, "txt", sep = "."), sep = "/"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
        write.table(gtmatrix, paste(outfolder, paste("gt_data", gene_id, "txt", sep = "."), sep = "/"), quote = FALSE, sep = "\t", row.names = TRUE , col.names = NA  )
    }

    suppressMessages(file.remove(list.files(path = outfolder, pattern = paste("tmp", gene_id, sep = "."), full.names = TRUE)))
    closeAllConnections()
}

# Run qsub to divide genotypes
create_sh_divide_genotypes = function(rscript, folder, geneinfo_file)
{
    command = c("#!/usr/bin/sh", "", 
                paste(rscript, 
                      paste(getwd(), "analysis", "cardiac_qtls_run_divide_genotypes.R", sep = "/"),
                      "--taskid $SGE_TASK_ID", 
                      "--folder", folder, 
                      "--geneinfo_file", geneinfo_file
                     ),
                sep = "\n")
    
    sh_file = paste(folder, "sh", sep = ".")
    writeLines(text = command, con = sh_file, sep = "\n")
    
    return(paste(getwd(), sh_file, sep = "/"))
}

run_qsub_genotypes = function(qtl_folder, rscript, geneinfo_file, n_genes, run_qsub = FALSE, queue = "week")
{
    sh_file  = create_sh_divide_genotypes(rscript, qtl_folder, geneinfo_file)
    
    err_file = gsub("\\.sh", ".err", sh_file)
    out_file = gsub("\\.sh", ".out", sh_file)
    
    suppressWarnings(file.remove(err_file))
    suppressWarnings(file.remove(out_file))

    qsub_options = paste("-l", queue, 
                         "-pe", "smp 1", 
                         "-t", paste(1, n_genes, sep = "-"), 
                         "-tc", 300,
                         "-o", err_file,
                         "-e", out_file
                        )

    command = paste("qsub", qsub_options, sh_file)
    
    message(command)

    if (run_qsub == TRUE ){system(command)}
    if (run_qsub == FALSE){return(command)}
}

### Gene expression normalization
my.invnorm = function(x)
{
    res = rank(x)
    res = qnorm(res/(length(res)+0.5))
    return(res)
}

transform_standard_normal = function(df)
{
    data_valid_expressed_full_qn = normalize.quantiles(as.matrix(df), copy=FALSE)

    input_mat = as.data.frame(t(apply(t(data_valid_expressed_full_qn), 2, my.invnorm)))
    
    return(input_mat)
}

normalize_tpm = function(tpm, min_tpm = 2, min_samples = 0.1) # add PEER factors?
{
    message("Normalizing expression...")
    message(paste("Total genes/peaks", nrow(tpm), sep = " = "))
    message(paste("Total samples"    , ncol(tpm), sep = " = "))
    
    expressed = as.matrix(tpm)

    expressed[as.matrix(tpm) <  min_tpm] = 0
    expressed[as.matrix(tpm) >= min_tpm] = 1

    tpm_f          = tpm[rowSums(expressed) >= (min_samples * ncol(tpm)),]
    message(paste("Total expressed genes/peaks"    , nrow(tpm_f), sep = " = "))
    tpm_f_std_norm = transform_standard_normal(tpm_f)
    
    return(list(tpm_f = tpm_f, tpm_f_std_norm = tpm_f_std_norm, gene_ids = rownames(tpm_f), sample_ids = colnames(tpm_f)))
}

divide_phenotypes_by_gene = function(expdata, outfolder)
{
    message("Dividing phenotype data by gene/peak...")
    gene_ids   = expdata$gene_ids
    sample_ids = expdata$sample_ids
    rawexp     = expdata$tpm_f
    normexp    = expdata$tpm_f_std_norm
    
    for(gene in gene_ids)
    {
        outdata = data.frame(sample_id = sample_ids, raw = as.numeric(rawexp[gene,sample_ids]), norm = as.numeric(normexp[gene,sample_ids]))
        write.table(outdata, paste(outfolder, paste(gene, "txt", sep = "."), sep = "/"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
    }
}

divide_phenotypes_by_study = function(ids, prefix)
{
    message(prefix)
    suppressMessages(source("analysis/cardiac_qtls_packages.R"      ))
    suppressMessages(source("analysis/cardiac_qtls_input_files.R"   ))
    suppressMessages(source("analysis/cardiac_qtls_functions.R"     ))
    suppressMessages(source("analysis/cardiac_qtls_input_data.R"    ))
    suppressMessages(source("analysis/cardiac_qtls_load_metadata.R" ))
    
    tpm           = fread("input/phenotypes/tpm.txt", sep = "\t", header = TRUE, data.table = FALSE)
    rownames(tpm) = tpm$V1
    tpm$V1        = NULL
    tpm           = tpm[, ids] 
    expdata       = normalize_tpm(tpm, min_tpm = phenotype_min_value, min_samples = phenotype_min_samples)
    outfolder     = paste("input/phenotypes/rna", prefix, sep = "_")

    dir.create(outfolder, showWarnings = FALSE)
    divide_phenotypes_by_gene(expdata, outfolder)

    write.table(expdata$gene_ids      , paste(outfolder, "list.txt", sep = "_"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
    write.table(expdata$tpm_f         , paste(outfolder, "raw.txt" , sep = "_"), quote = FALSE, sep = "\t", row.names = TRUE , col.names = NA   )
    write.table(expdata$tpm_f_std_norm, paste(outfolder, "norm.txt", sep = "_"), quote = FALSE, sep = "\t", row.names = TRUE , col.names = NA   )
}


# Initialize QTL analysis
write_read_table = function(folder, x, header = TRUE, one_column = FALSE)
{
    out = paste(x, "=", 'read.table(', paste(paste(paste('"', paste(folder, paste(x,"txt", sep = "."), sep = "/"), '"', sep = ""), 
                                                   paste("header"     , header , sep = " = "), 
                                                   paste("check.names", "FALSE", sep = " = "), 
                                                   sep = ", "), ")", sep = ""), sep = "")
    
    if (one_column == TRUE){out = paste(out, "$V1", sep = "")}
    
    return(out)
}


initialize_qtl_analysis = function(analysis_name, phenotype, phenotype_folder, sample_list, gene_list, vars0_assay, vars1_assay, 
                                   qtl_distance          = 1e6   , # Distance threshold between each gene/peak and variants
                                   maf_threshold         =   0.01, # MAF threshold for variants used for QTL analysis
                                   phenotype_min_value   =   0.5 , # Threshold to consider a gene/peak expressed
                                   phenotype_min_samples =   0.1 , # Fraction of samples that have expression above phenotype_min_value
                                   geneinfo_file         = "/publicdata/gencode_v19_20151104/gene_info.tsv", # for ATAC and ChIP, use the appropriate info files
                                   n_perm                = 10,
                                   primary               = TRUE
                                  )
{
    suppressMessages(source("analysis/cardiac_qtls_packages.R"      ))
    suppressMessages(source("analysis/cardiac_qtls_input_files.R"   ))
    suppressMessages(source("analysis/cardiac_qtls_functions.R"     ))
    suppressMessages(source("analysis/cardiac_qtls_input_data.R"    ))
    suppressMessages(source("analysis/cardiac_qtls_load_metadata.R" ))
    
    geneinfo = read.table(geneinfo_file, header = TRUE)

    # Create output folders
    dir.create(paste("qtls", analysis_name, sep = "/"), showWarnings = FALSE)
    invisible(lapply(c("input", "qtls", "analysis", "log", "sh"), function(x){dir.create(paste("qtls", analysis_name, x, sep = "/"), showWarnings = FALSE)}))
    
    phenotype2 = phenotype
    if (grepl("atac", phenotype) == TRUE){phenotype2 = "atac"}
    
    # Rewrite input files and metadata script
    metadata           = metadata[metadata[,paste(phenotype2, "id", sep = "_")] %in% sample_list,]
    metadata           = metadata[,paste(c("subject", "wgs", phenotype2), "id", sep = "_")]
    colnames(metadata) = paste(c("subject", "wgs", "assay"), "id", sep = "_")
    covariates_subject = covariates_subject[covariates_subject$subject_id %in% metadata$subject_id,]
    covariates_assay   = read.table(paste("input/covariates/covariates", phenotype2, "txt", sep = "."), header = TRUE , check.names = FALSE)
    covariates_assay   = covariates_assay[covariates_assay$assay_id %in% metadata$assay_id,]
    unrelated          = unrelated[unrelated %in% metadata$subject_id]
    
    #gene_list = read.table(paste(phenotype_folder, "list.txt", sep = "_"), header = FALSE)$V1
    
    outfolder = paste(getwd(), "qtls", analysis_name, sep = "/")

    write.table(metadata          , paste(outfolder, "input", "metadata.txt"          , sep = "/"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE )
    write.table(covariates_subject, paste(outfolder, "input", "covariates_subject.txt", sep = "/"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE )
    write.table(covariates_assay  , paste(outfolder, "input", "covariates_assay.txt"  , sep = "/"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE )
    write.table(unrelated         , paste(outfolder, "input", "unrelated.txt"         , sep = "/"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
    write.table(gene_list         , paste(outfolder, "input", "gene_list.txt"         , sep = "/"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
    
    # link input data
    invisible(suppressWarnings(file.symlink(paste(getwd(),phenotype_folder                     , sep = "/"), paste(outfolder, "input", "phenotypes"   , sep = "/"))))
    invisible(suppressWarnings(file.symlink(paste(getwd(),"input", "genotypes" , phenotype     , sep = "/"), paste(outfolder, "input", "genotypes"    , sep = "/"))))
    invisible(suppressWarnings(file.symlink(paste(getwd(),"analysis/cardiac_qtls_packages.R"   , sep = "/"), paste(outfolder, "input", "packages.R"   , sep = "/"))))
    invisible(suppressWarnings(file.symlink(paste(getwd(),"analysis/cardiac_qtls_input_files.R", sep = "/"), paste(outfolder, "input", "input_files.R", sep = "/"))))
    invisible(suppressWarnings(file.symlink(paste(getwd(),"analysis/cardiac_qtls_functions.R"  , sep = "/"), paste(outfolder, "input", "functions.R"  , sep = "/"))))
    invisible(suppressWarnings(file.symlink(paste(getwd(),"analysis/cardiac_qtls_input_data.R" , sep = "/"), paste(outfolder, "input", "input_data.R" , sep = "/"))))
    invisible(suppressWarnings(file.symlink(              geneinfo_file                                    , paste(outfolder, "input", "geneinfo.txt" , sep = "/"))))

    # write metadata.R file
    write_read_table_list = c(write_read_table(paste(outfolder, "input", sep = "/"), "metadata"          , header = TRUE , one_column = FALSE),
                              write_read_table(paste(outfolder, "input", sep = "/"), "covariates_subject", header = TRUE , one_column = FALSE),
                              write_read_table(paste(outfolder, "input", sep = "/"), "covariates_assay"  , header = TRUE , one_column = FALSE),
                              write_read_table(paste(outfolder, "input", sep = "/"), "unrelated"         , header = FALSE, one_column = TRUE ),
                              write_read_table(paste(outfolder, "input", sep = "/"), "gene_list"         , header = FALSE, one_column = TRUE ),
                              write_read_table(paste(outfolder, "input", sep = "/"), "geneinfo"          , header = TRUE , one_column = FALSE),
                              "",
                              paste("vars0_assay"          , paste('c("', paste(vars0_assay, collapse = '", "'), '")', sep = ""), sep = " = "),
                              paste("vars1_assay"          , paste('c("', paste(vars1_assay, collapse = '", "'), '")', sep = ""), sep = " = "),
                              paste("qtl_distance"         , qtl_distance         , sep = " = "),
                              paste("maf_threshold"        , maf_threshold        , sep = " = "),
                              paste("phenotype_min_value"  , phenotype_min_value  , sep = " = "),
                              paste("phenotype_min_samples", phenotype_min_samples, sep = " = "),
                              paste("n_perm"               , n_perm               , sep = " = "),
                              paste("primary"              , primary              , sep = " = "),
                              "",
                              paste('setwd("', outfolder, '")', sep = "")
                             )

    writeLines(text = write_read_table_list, con = paste("qtls", analysis_name, "input", "load_metadata.R", sep = "/"), sep = "\n")
    
    # write SH file
    to_sh = c('#!/usr/bin/sh',
              "",
              paste(rscript, "/frazer01/projects/CARDIPS/analysis/cardiac_qtls/analysis/cardiac_qtls_run_eqtls.R", "--taskid $SGE_TASK_ID", "--folder", outfolder)
             )
    
    writeLines(text = to_sh, con = paste(outfolder, "sh", "cardiac_qtls_run_eqtls.sh", sep = "/"), sep = "\n")
    
    
    return(list(folder = outfolder, n_genes = length(gene_list)))
}

# Run QTLs
get_lmm_pval = function(lmm, n, gene_id, snp)
{
    anova0   = anova(lmm)
    fval     = anova0["gt", "F value"]
    pval     = pf(fval, df1 = 1, df2 = n - 1, lower.tail = FALSE)
    beta     = summary(lmm)$coefficients[2,1]
    se       = summary(lmm)$coefficients[2,2]
    
    return(data.frame(gene_id = gene_id, id = snp, beta = beta, se = se, pval = pval))
}

run_lmm_by_variant = function(gene_id, snp, expdata, gtdata, meta, covariates_assay, covariates_subject, vars0, vars1, compare = FALSE, type = "", n_perm = 0)
{
    meta           = meta   [ meta$assay_id %in% intersect(rownames(expdata), colnames(gtdata)),]
    expdata        = expdata[meta$assay_id,]
    gtdata         = gtdata [             , meta$assay_id]
    input          = expdata
    input$gt       = as.numeric(gtdata[snp, meta$assay_id])
    input$assay_id = rownames(input)
    input          = merge(input, meta)
    input          = merge(input, covariates_assay)
    input          = merge(input, covariates_subject)
    input          = input[order(input$wgs_id),]
    #input$study    = as.numeric(factor(input$study ))
    input$sex      = as.numeric(factor(input$sex   ))
    
    lmm0   = suppressMessages(suppressWarnings(lmer(paste("norm", paste(vars0, collapse = "+"), sep = "~"), data = input, REML = FALSE)))
    
    if (compare == FALSE)
    {
        out_lm = get_lmm_pval(lmm0, nrow(input), gene_id, snp)
        return(out_lm)
    }
    if (compare == TRUE)
    {
        compare_lmm         = suppressMessages(suppressWarnings(find_best_model(snp, do.call("cbind", lapply(vars1, function(var1){compare_lmms(lmm0, vars0, var1, input, n_perm)})), vars1)))
        compare_lmm$gene_id = gene_id
        compare_lmm$type    = type
        return(compare_lmm)
    }
}

calculate_r2 = function(told, x1, x2)
{
    x1  = as.numeric(told[x1,])
    x2  = as.numeric(told[x2,])
    xx  = data.frame(x1 = x1, x2 = x2)
    p1  = sum(xx$x1) / nrow(xx)
    q1  = sum(xx$x2) / nrow(xx)
    x11 = sum(unlist(apply(xx, 1, min))) / nrow(xx)
    d   = x11 - p1 * q1
    r2  = d^2 / (p1 * q1 * (1 - p1) * (1 - q1))
    
    return(r2)
}

find_independent_variants = function(lmm_data, gtdata, unrelated, meta)
{
    lmm_signif     = lmm_data[lmm_data$bonferroni < 0.1,]
    told           = gtdata[lmm_signif$id, meta[meta$subject_id %in% unrelated, "assay_id"]]
    lmm_signif$top = ""

    start_time = Sys.time()

    for(ii in 1:nrow(lmm_signif))
    {
        snp1   = lmm_signif[ii, "id"]
        if (lmm_signif[ii, "top"] == "")
        {
            snps2  = lmm_signif[(ii + 1): nrow(lmm_signif), "id"]
            snp2r2 = data.frame(snp = snps2, r2 = unlist(lapply(snps2, function(snp2){calculate_r2(told, snp1, snp2)})))

            lmm_signif[lmm_signif$id %in% snp2r2[snp2r2$r2 > 0.6, "snp"], "top"] = snp1
        }
    }
    return(lmm_signif)
    
}

compare_lmms = function(lmm0, vars0, var1, input, n_perm)
{
    lmm1        = suppressWarnings(lmer(paste("norm", paste(c(vars0, var1), collapse = "+"), sep = "~"), data = input, REML = FALSE))
    anova_data  = as.data.frame(anova(lmm0, lmm1))
    out_compare = data.frame(aic0 = anova_data[1,"AIC"], aic1 = anova_data[2,"AIC"], delta_aic = anova_data[2,"AIC"] - anova_data[1,"AIC"], pval = anova_data[2, "Pr(>Chisq)"])
    out_compare$pval_perm = PBmodcomp(lmm1, lmm0, nsim = n_perm)$test$p.value[[2]]
    
    colnames(out_compare) = paste(var1, colnames(out_compare), sep = ":")
    
    return(out_compare)
}

find_best_model = function(snp, compare_lmm, vars1)
{
    min_pos       = which.min(c(compare_lmm[,paste(vars1, "delta_aic", sep = ":")]))
    best          = vars1[min_pos]
    min_pval      = compare_lmm[,paste(best, "pval"     , sep = ":")]
    min_pval_perm = compare_lmm[,paste(best, "pval_perm", sep = ":")]
    min_aic       = compare_lmm[,paste(best, "delta_aic", sep = ":")]
    
    if (min_aic > 0)
    {
        best     = "gt"
        min_pval = 1
    }
    
    compare_lmm$best          = best
    compare_lmm$min_aic       = min_aic
    compare_lmm$min_pval      = min_pval
    compare_lmm$min_pval_perm = min_pval_perm
    compare_lmm$id            = snp
    
    return(compare_lmm)
}

run_qtl_analysis = function(gene_id, wd, outfolder, sample_list, covariates_assay, vars0, vars1, var_list = NULL)
{
    setwd(wd)
    
    suppressMessages(source("input/packages.R"      ))
    suppressMessages(source("input/input_files.R"   ))
    suppressMessages(source("input/functions.R"     ))
    suppressMessages(source("input/input_data.R"    ))
    suppressMessages(source("input/load_metadata.R" ))
    
    # Read gene-specific data
    gtdata  = read.table(paste(wd, "input/genotypes" , paste("gt_data", gene_id, "txt", sep = "."), sep = "/"), header = TRUE, check.names = FALSE, row.names = 1)
    gtinfo  = read.table(paste(wd, "input/genotypes" , paste("gt_info", gene_id, "txt", sep = "."), sep = "/"), header = TRUE, check.names = FALSE)
    expdata = read.table(paste(wd, "input/phenotypes", paste(           gene_id, "txt", sep = "."), sep = "/"), header = TRUE, check.names = FALSE, row.names = 1)
    
    if (is.null(var_list) == FALSE)
    {
        gtdata = gtdata[var_list,]
    }
    
    gtdata   = gtdata[is.na(rowSums(gtdata)) == FALSE, ]
    gtinfo   = gtinfo[gtinfo$id %in% rownames(gtdata),]
    var_list = rownames(gtdata)
    
    # Transform to have the same sample IDs
    meta              = metadata[metadata$wgs_id %in% colnames(gtdata) & metadata$assay_id %in% rownames(expdata),]
    meta              = meta[meta$assay_id %in% sample_list,]
    meta              = meta[order(meta$wgs_id),]
    gtdata            = gtdata [             , meta$wgs_id]
    expdata           = expdata[meta$assay_id,            ]
    colnames(gtdata ) = meta$assay_id
    rownames(expdata) = meta$assay_id
    
    # First pass: run QTLs for all variants
    min_p      = 0
    cond       = 0
    vars0_cond = vars0

    lmm_data_list   = list()
    top_hit_list    = list()
    top_hits        = c()
    this_covariates = covariates_assay

    while(min_p < 0.1)
    {
        if (cond == 0)
        {
            type            = "primary"
            new_exp         = expdata
            this_covariates = covariates_assay
        }else
        {
            type            = paste("conditional", cond, sep = "")
            this_hit        = top_hit_list[[cond]]
            this_hit        = this_hit[1, "id"]
            this_covariates = merge(this_covariates, data.frame(assay_id = colnames(gtdata), cond = as.numeric(gtdata[this_hit,])))

            this_covariates[, paste("var", this_hit, sep = "")] = this_covariates$cond

            vars0_cond        = c(vars0_cond, paste("var", this_hit, sep = ""))
            new_exp           = merge(expdata, this_covariates[,c("assay_id", paste("var", top_hits, sep = ""))], by.x = "row.names", by.y = "assay_id")
            rownames(new_exp) = new_exp$Row.names
            res               = residuals(lm(paste("norm", paste(paste("var", top_hits, sep = ""), collapse = "+"), sep = "~"), data = new_exp))
            new_exp$norm      = mean(new_exp$norm) + res
        }

        lmm_data            = do.call("rbind", lapply(var_list, function(snp){run_lmm_by_variant(gene_id, snp, new_exp, gtdata, meta, this_covariates, covariates_subject, vars0_cond, vars1, compare = FALSE, type = "")}))
        lmm_data            = lmm_data[order(lmm_data$pval),]
        lmm_data$bonferroni = p.adjust(lmm_data$pval)
        lmm_data$type       = type
        top_hit             = lmm_data[1,]

        if ((top_hit[1, "bonferroni"] < 0.1)|(cond == 0))
        {
            lmm_data_list[[cond + 1]] = lmm_data
            top_hit_list [[cond + 1]] = top_hit
            top_hits                  = c(top_hits, top_hit[1, "id"])
        }
        
        if (primary == TRUE){break}
        if (cond    == 5   ){break}

        cond  = cond + 1
        min_p = top_hit[1, "bonferroni"]
    }

    lmm_data = do.call("rbind", lmm_data_list)
    top_hit  = do.call("rbind", top_hit_list )
    top_hit  = top_hit[order(top_hit$pval),]
    
    # Write output files
    write.table(merge(gtinfo, lmm_data   ), file = paste(outfolder, paste("all"    , gene_id, "txt", sep = "."), sep = "/"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
    write.table(merge(gtinfo, top_hit    ), file = paste(outfolder, paste("top_hit", gene_id, "txt", sep = "."), sep = "/"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
    
    top_hit = top_hit[top_hit$type == "primary",]
    
    #if ((nrow(lmm_data[lmm_data$bonferroni < 0.1,]) > 0)&(file.exists(paste(outfolder, paste("signif", gene_id, "txt", sep = "."), sep = "/")) == FALSE))
    if ((nrow(top_hit[top_hit$bonferroni < 0.1,]) > 0)) # rerun lmm_signif
    {
        #lmm_compare = do.call("rbind", lapply(1:nrow(top_hit), function(ii){run_lmm_by_variant(gene_id, top_hit[ii,"id"], expdata, gtdata, meta, covariates_assay, covariates_subject, vars0, vars1, compare = TRUE, type = top_hit[ii,"type"])})) # too slow
        lmm_compare = run_lmm_by_variant(gene_id, top_hit[1,"id"], expdata, gtdata, meta, covariates_assay, covariates_subject, vars0, vars1, compare = TRUE, type = top_hit[1,"type"], n_perm)

        # Write output files
        write.table(merge(gtinfo, lmm_compare), file = paste(outfolder, paste("lmm_compare", gene_id, "txt", sep = "."), sep = "/"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
    }
    closeAllConnections()
}

# Run qsub for all genes/peaks
run_qsub = function(qtl_folder, n_genes, run_qsub = FALSE, queue = "week", tc = 300)
{
    err_file = paste(qtl_folder, "log", "cardiac_qtls_run_eqtls.out", sep = "/")
    out_file = paste(qtl_folder, "log", "cardiac_qtls_run_eqtls.err", sep = "/")
    
    suppressWarnings(file.remove(err_file))
    suppressWarnings(file.remove(out_file))


    qsub_options = paste("-l", queue, 
                         "-pe", "smp 1", 
                         "-t", paste(1, n_genes, sep = "-"), 
                         "-tc", tc,
                         "-o", err_file,
                         "-e", out_file
                        )


    command = paste("qsub", qsub_options, paste(qtl_folder, "sh", "cardiac_qtls_run_eqtls.sh", sep = "/"))

    if (run_qsub == TRUE ){message("Running qsub")}
    if (run_qsub == TRUE ){system(command)}
    if (run_qsub == FALSE){return(command)}
}

# Monitor QTL progression
monitor_qtls = function(analysis_name)
{
    genes     = length(readLines(paste("qtls", analysis_name, "input", "gene_list.txt", sep = "/")))
    genes_run = length(gsub("all\\.", "", gsub("\\.txt", "", list.files(paste("qtls", analysis_name, "qtls", sep = "/"), pattern = "^all"        ))))
    lmms_run  = length(gsub("all\\.", "", gsub("\\.txt", "", list.files(paste("qtls", analysis_name, "qtls", sep = "/"), pattern = "^lmm_compare"))))

    message(paste(Sys.time(),
                  analysis_name,
                  paste("Genes analyzed"  ,         genes_run, sep = " = "), 
                  paste("Genes to analyze", genes - genes_run, sep = " = "), 
                  paste("LMMs compared"   , lmms_run         , sep = " = "), 
                  sep = "\n"))
}

# Merge eGenes to a single file
merge_files = function(infolder, pattern, geneinfo)
{
    qtls = do.call("rbind", lapply(list.files(infolder, pattern = pattern, full.names = TRUE), function(x){read.table(x, header = TRUE, check.names = FALSE)}))
    if (sum(colnames(qtls) %in% c("chrom")) > 0){geneinfo = geneinfo[,colnames(geneinfo) != "chrom"]}
    
    qtls = merge(geneinfo, qtls, by.x = "gene_id", by.y = "gene_id")
    
    return(qtls)
}

fdr_by_type = function(egenes, column, type)
{
    this       = egenes[egenes$type == type,]
    this$fdr   = p.adjust(this[,column], method = "BH", n = length(unique(egenes$gene_id)))
    this$egene = FALSE
    this[this$fdr < 0.1, "egene"] = TRUE
    return(this)
}

merge_qtls = function(wd, phenotype)
{
    old_wd = getwd()
    setwd(wd)
    
    suppressMessages(source("input/packages.R"      ))
    suppressMessages(source("input/input_files.R"   ))
    suppressMessages(source("input/functions.R"     ))
    suppressMessages(source("input/input_data.R"    ))
    suppressMessages(source("input/load_metadata.R" ))
    
    egenes             = merge_files("qtls", "top_hit"    , geneinfo)
    lmm_compare        = merge_files("qtls", "lmm_compare", geneinfo)
    n_egenes           = length(unique(egenes$gene_id))
    egenes             = do.call("rbind", lapply(sort(unique(egenes$type)), function(type){fdr_by_type(egenes, "bonferroni", type)}))
    egenes_tested      = unique(paste(egenes[egenes$egene == TRUE, "gene_id"], egenes[egenes$egene == TRUE, "type"]))
    lmm_compare$tested = paste(lmm_compare$gene_id, lmm_compare$type)
    lmm_compare        = lmm_compare[lmm_compare$gene_id %in% egenes[egenes$egene == TRUE, "gene_id"],]
    min_pval           = "min_pval_perm"
    
    if (grepl("rna", phenotype) == FALSE){min_pval = "min_pval"}
    
    lmm_compare        = do.call("rbind", lapply(sort(unique(lmm_compare$type)), function(type){fdr_by_type(lmm_compare, min_pval, type)}))
    lmm_compare[lmm_compare$fdr > 0.1, "best"] = "gt"
    
    lmm_compare        = lmm_compare[lmm_compare$tested %in% egenes_tested, ]
    lmm_compare$egene  = NULL
    lmm_compare$tested = NULL

    write.table(egenes     , paste("analysis/egenes"     , phenotype, "txt", sep = "."), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
    write.table(lmm_compare, paste("analysis/lmm_compare", phenotype, "txt", sep = "."), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

    message(paste(phenotype, paste("Tested genes", n_egenes, sep = " = "), paste("eGenes", nrow(egenes[egenes$egene == TRUE, ]), sep = " = "), "", sep = "\n"))
    
    setwd(old_wd)
    #return(list(egenes = egenes, lmm_compare = lmm_compare))
}

# COLOC
read_gwas_table = function(infile)
{
    indata = read.table(infile, header = TRUE, sep = "\t")[,1]
    y      = unique(unlist(lapply(indata, function(z){unlist(strsplit(z, " "))[[1]]})))
    
    return(y)
}

get_ukbb_data = function(gene_id, main, chrom, from, to, traits, type)
{
    outfile = paste("private", "coloc", main, paste(gene_id, type, "txt", sep = "."), sep = "/")
    command = paste(bcftools, "query", 
                    "-r", paste(gsub("chr", "", chrom), ":", from, "-", to, sep = ""),
                    "-s", paste(traits, collapse = ","),
                    "-f", paste('"%POS\t%ID\t%AF[\t%', type, ']\n"', sep = ""),
                    paste("/frazer01/publicdata/ukbiobank_20180801/ukbb2vcf/vcf", paste(chrom, "vcf", "gz", sep = "."), sep = "/"),
                    ">", 
                    outfile
                   )
    
    system(command)
    
    #return(command)
    
    indata           = fread(outfile, sep = "\t", header = FALSE, data.table = FALSE)
    colnames(indata) = c("pos", "id", "af", traits)
    
    #write.table(indata, outfile, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
    file.remove(outfile)
    
    return(indata)
}


qtl2df = function(gene_id, tissue)
{
    infile = paste("qtls", tissue, "qtls", paste("all", gene_id, "txt", sep = "."), sep = "/")

    indata = fread(infile, sep = "\t", header = TRUE, data.table = FALSE)
    outdata = indata[,c("id", "beta", "se", "pval")]

    return(outdata)
}

trait2df = function(trait, pv)
{
    this         = data.frame(id = pv$id, af = pv$af, pval = pv[,trait])
    
    return(this)
}

run_coloc = function(to_coloc, id1, id2)
{
    b1 = to_coloc[[id1]]
    b2 = to_coloc[[id2]]
    
    rownames(b1) = b1$id
    rownames(b2) = b2$id

    ids = intersect(b1$id, b2$id)

    b1 = b1[ids,]
    b2 = b2[ids,]

    coloc_mapped = coloc.abf(dataset1=list(snp = b1$id, pvalues = b1$pval, N = length(ids), type = "quant"),
                             dataset2=list(snp = b2$id, pvalues = b2$pval, N = length(ids), type = "quant"),
                             MAF     =b2$af
                            )   
    
    probs           = as.data.frame(t(coloc_mapped$summary))
    myres           = coloc_mapped$results
    myres           = coloc_mapped$results
    myres           = myres[, c(which(colnames(myres) == "snp"), ncol(myres))]
    colnames(myres) = c("snp_id", "pp_snp")
    myres           = myres[order(myres$pp_snp, decreasing = TRUE),]
    
    return(cbind(data.frame(id1 = id1, id2 = id2), probs, myres[1,]))
}

run_coloc_all_interactions = function(gene_id, to_coloc, to_test, main)
{
    coloc         = do.call("rbind", lapply(1:nrow(to_test), function(ii){run_coloc(to_coloc, to_test[ii, "qtl"], to_test[ii, "trait"])}))
    coloc$gene_id = gene_id
    
    write.table(coloc, paste("private", "coloc", main, paste(gene_id, "txt", sep = "."), sep = "/"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}

coloc_by_main = function(main, taskid)
{
    suppressMessages(source("analysis/cardiac_qtls_packages.R"      ))
    suppressMessages(source("analysis/cardiac_qtls_input_files.R"   ))
    suppressMessages(source("analysis/cardiac_qtls_functions.R"     ))
    suppressMessages(source("analysis/cardiac_qtls_input_data.R"    ))
    suppressMessages(source("analysis/cardiac_qtls_load_metadata.R" ))

    gene2main = read.table(paste("private", "coloc", "egenes.txt", sep = "/"), header = TRUE)
    gene_ids  = sort(unique(gene2main$gene_id))
    gene_id   = gene_ids[[taskid]]
    chrom     = geneinfo[geneinfo$gene_id == gene_id, "chrom"]
    from      = geneinfo[geneinfo$gene_id == gene_id, "start"] - qtl_distance
    to        = geneinfo[geneinfo$gene_id == gene_id, "start"] + qtl_distance
    
    if (file.exists(paste("/frazer01/projects/GTEx_v7/analysis/eqtls_deconvolution/qtls/", main, "_original/input/phenotypes/", gene_id, ".txt", sep = "")) == TRUE)
    {
        traits             = read_gwas_table(paste("/home/mdonovan/gtex_deconvolution/tables/", main, "_gwas.txt", sep = ""))
        pv                 = get_ukbb_data(gene_id, main, chrom, from, to, traits, "PV")
        trait2coloc        = lapply(traits, function(trait){trait2df(trait, pv)})
        qtl2coloc          = lapply(tissue2name[tissue2name$main == main, "tissue"], function(tissue){qtl2df(gene_id, tissue)})
        names(trait2coloc) = traits
        names(qtl2coloc)   = tissue2name[tissue2name$main == main, "tissue"]
        to_coloc           = c(qtl2coloc, trait2coloc)
        to_test            = expand.grid(qtl = names(qtl2coloc), trait = names(trait2coloc), stringsAsFactors = FALSE)

        run_coloc_all_interactions(gene_id, to_coloc, to_test, main)
    }
}

run_qsub_coloc = function(qtl_folder, n_genes, run_qsub = FALSE, queue = "week", tc = 500)
{
    err_file = paste(qtl_folder, "log", "coloc.out", sep = "/")
    out_file = paste(qtl_folder, "log", "coloc.err", sep = "/")
    
    suppressWarnings(file.remove(err_file))
    suppressWarnings(file.remove(out_file))


    qsub_options = paste("-l", queue, 
                         "-pe", "smp 1", 
                         "-t", paste(1, n_genes, sep = "-"), 
                         "-tc", tc,
                         "-o", err_file,
                         "-e", out_file
                        )

    command = paste("qsub", qsub_options, paste(qtl_folder, "analysis/cardiac_qtls_run_coloc.sh", sep = "/"))

    if (run_qsub == TRUE ){message("Running qsub")}
    if (run_qsub == TRUE ){system(command)}
    if (run_qsub == FALSE){return(command)}
}


