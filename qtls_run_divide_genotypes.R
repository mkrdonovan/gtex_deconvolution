### Version history:
# V02: changed to pair with functions_V03

setwd("/frazer01/projects/GTEx_v7/analysis/eqtls_deconvolution")

suppressMessages(source("analysis/cardiac_qtls_packages.R"      ))
suppressMessages(source("analysis/cardiac_qtls_input_files.R"   ))
suppressMessages(source("analysis/cardiac_qtls_functions.R"     ))
suppressMessages(source("analysis/cardiac_qtls_input_data.R"    ))
suppressMessages(source("analysis/cardiac_qtls_load_metadata.R" ))

option_list = list(make_option("--taskid"       , type="integer"  , default=0, help="SGE task ID"                                   , metavar="character"),
                   make_option("--folder"       , type="character", default=0, help="Analysis folder (from initialize_qtl_analysis)", metavar="character"),
                   make_option("--geneinfo_file", type="character", default=0, help="geneinfo file"                                 , metavar="character")) 

opt_parser    = OptionParser(option_list=option_list)
opt           = parse_args(opt_parser)
taskid        = opt$taskid
infolder      = opt$folder
geneinfo_file = opt$geneinfo_file

write_genotypes(taskid, infolder, geneinfo_file)
