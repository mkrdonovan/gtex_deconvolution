
setwd("/frazer01/projects/GTEx_v7/analysis/eqtls_deconvolution")

source("/frazer01/projects/GTEx_v7/analysis/eqtls_deconvolution/qtls/liver_cells/input/packages.R"      )
source("/frazer01/projects/GTEx_v7/analysis/eqtls_deconvolution/qtls/liver_cells/input/input_files.R"   )
source("/frazer01/projects/GTEx_v7/analysis/eqtls_deconvolution/qtls/liver_cells/input/functions.R"     )
source("/frazer01/projects/GTEx_v7/analysis/eqtls_deconvolution/qtls/liver_cells/input/input_data.R"    )
source("/frazer01/projects/GTEx_v7/analysis/eqtls_deconvolution/qtls/liver_cells/input/load_metadata.R" )

option_list = list(make_option("--taskid", type="integer"  , default=0, help="SGE task ID"                                   , metavar="character")) 

opt_parser  = OptionParser(option_list=option_list)
opt         = parse_args(opt_parser)
taskid      = opt$taskid

qtls = fread("analysis//egenes.liver_cells.txt", header = TRUE, sep = "\t", data.table = FALSE)

subject_ids = as.character(fread(cmd = paste("head", "-n 1", paste("input//genotypes", paste("gt_data", qtls[1,"gene_id"], "txt", sep = "."), sep = "/")), header = FALSE, data.table = FALSE))
subject_ids[[1]] = "id"

run_lmm_by_variant_perm = function(gene_id, gtdata, meta, covariates_assay, covariates_subject, vars0, vars1, compare = FALSE, type = "", n_perm = 0)
{
    expdata           = fread(paste("input/phenotypes", paste(gene_id, "txt", sep = "."), sep = "/"), header = TRUE, sep = "\t", data.table = FALSE)
    rownames(expdata) = expdata$sample_id
    meta              = meta   [ meta$assay_id %in% intersect(rownames(expdata), colnames(gtdata)),]
    gtdata            = gtdata [             , meta$assay_id]
    input             = expdata
    input$gt          = as.numeric(gtdata[gene_id, meta$assay_id])
    input$assay_id    = rownames(input)
    input             = merge(input, meta)
    input             = merge(input, covariates_assay)
    input             = merge(input, covariates_subject)
    input             = input[order(input$wgs_id),]
    input$sex         = as.numeric(factor(input$sex   ))
    
    lmm0   = suppressMessages(suppressWarnings(lmer(paste("norm", paste(vars0, collapse = "+"), sep = "~"), data = input, REML = FALSE)))
    out_lm = get_lmm_pval(lmm0, nrow(input), gene_id, gene_id)
    
    out_lm$bonferroni = min(c(1,out_lm$pval * as.numeric(unlist(strsplit(system(paste("wc -l", paste("input/genotypes", paste("gt_data", gene_id, "txt", sep = "."), sep = "/")), intern = TRUE), split = " "))[[1]])))
    
    return(out_lm)
}

run_permute_qtls = function(perm = 1)
{
    covariates_assay  = fread(paste("perm/input", paste("covariates_assay", perm, "txt", sep = "."), sep = "/"), sep = "\t", header = TRUE, data.table = FALSE)
    gtdata            = fread("perm/input/gtdata.txt"                                                          , sep = "\t", header = TRUE, data.table = FALSE)
    sample_list       = covariates_assay$assay_id
    rownames(gtdata)  = gtdata$V1
    
    meta              = metadata[metadata$wgs_id %in% colnames(gtdata),]
    meta              = meta[meta$assay_id %in% sample_list,]
    meta              = meta[order(meta$wgs_id),]
    gtdata            = gtdata [             , meta$wgs_id]
    colnames(gtdata ) = meta$assay_id
    
    lmm_data          = as.data.frame(rbindlist(lapply(rownames(gtdata), function(gene_id){run_lmm_by_variant_perm(gene_id, gtdata, meta, covariates_assay, covariates_subject, vars0_assay, vars1_assay, compare = FALSE, type = "")})), stringsAsFactors = FALSE)
    lmm_data$fdr   = p.adjust(lmm_data$bonferroni, method = "BH")
    lmm_data$perm     = perm
    
    fwrite(lmm_data, file = paste("perm/qtls", paste("qtls", perm, "txt", sep = "."), sep = "/"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}

run_permute_qtls(taskid)
