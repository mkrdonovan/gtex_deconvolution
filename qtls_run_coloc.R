### Version history:
# V02: changed to pair with functions_V03

setwd("/frazer01/projects/CARDIPS/analysis/cardiac_qtls")

suppressMessages(source("/frazer01/projects/CARDIPS/analysis/cardiac_qtls/analysis/cardiac_qtls_packages.R"      ))
suppressMessages(source("/frazer01/projects/CARDIPS/analysis/cardiac_qtls/analysis/cardiac_qtls_input_files.R"   ))
suppressMessages(source("/frazer01/projects/CARDIPS/analysis/cardiac_qtls/analysis/cardiac_qtls_functions.R"     ))
suppressMessages(source("/frazer01/projects/CARDIPS/analysis/cardiac_qtls/analysis/cardiac_qtls_input_data.R"    ))
suppressMessages(source("/frazer01/projects/CARDIPS/analysis/cardiac_qtls/analysis/cardiac_qtls_load_metadata.R" ))

option_list = list(make_option("--taskid", type="integer", default=0, help="SGE task ID", metavar="character")) 

opt_parser  = OptionParser(option_list=option_list)
opt         = parse_args(opt_parser)
taskid      = opt$taskid

nodes       = read.table("private/qtl_network/intersect/nodes.txt", header = TRUE, sep = "\t")
intersected = read.table("private/qtl_network/intersect/edges.txt", header = TRUE, sep = "\t")
#gene_ids    = nodes[nodes$type == "gene", "id"]
gene_ids    = sort(unique(intersected$gene_id1))
gene_id     = gene_ids[[taskid]]

run_coloc_all_interactions(gene_id, nodes, intersected)
