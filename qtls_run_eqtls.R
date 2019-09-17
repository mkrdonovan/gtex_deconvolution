### Version history:
# V02: changed to pair with functions_V03

setwd("/frazer01/projects/GTEx_v7/analysis/eqtls_deconvolution")

suppressMessages(source("analysis/cardiac_qtls_packages.R"      ))
suppressMessages(source("analysis/cardiac_qtls_input_files.R"   ))
suppressMessages(source("analysis/cardiac_qtls_functions.R"     ))
suppressMessages(source("analysis/cardiac_qtls_input_data.R"    ))
#suppressMessages(source("analysis/cardiac_qtls_load_metadata.R" ))

option_list = list(make_option("--taskid", type="integer"  , default=0, help="SGE task ID"                                   , metavar="character"),
                   make_option("--folder", type="character", default=0, help="Analysis folder (from initialize_qtl_analysis)", metavar="character")) 

opt_parser  = OptionParser(option_list=option_list)
opt         = parse_args(opt_parser)
taskid      = opt$taskid
infolder    = opt$folder

suppressMessages(source(paste(infolder, "input", "load_metadata.R", sep = "/")))

suppressMessages(source("input/packages.R"      ))
suppressMessages(source("input/input_files.R"   ))
suppressMessages(source("input/functions.R"     ))
suppressMessages(source("input/input_data.R"    ))

gene_id     = gene_list[[taskid]]

run_qtl_analysis(gene_id          = gene_id, 
                 wd               = infolder, 
                 outfolder        = paste(infolder, "qtls", sep = "/"), 
                 sample_list      = metadata$assay_id, 
                 covariates_assay = covariates_assay, 
                 vars             = vars0_assay, 
                 vars1            = vars1_assay, 
                 var_list         = NULL
                )
