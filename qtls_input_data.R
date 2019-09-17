message("Loading input data...")

bcftools = "/frazer01/software/bcftools-1.9/bcftools"
vcftools = "/frazer01/software/vcftools-0.1.14/bin/vcftools"
rscript  = "/frazer01/home/matteo/software/R-3.5.1/bin/Rscript"

ipscore_vcf_input = "/frazer01/projects/reference_files/cellType_Invariant/IPSCORE_WGS.biallelic.b37.vcf.gz"
gtex_vcf_input    = "/frazer01/projects/GTEx_v7/decrypted/GenotypeFiles/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_AllVar_QC_metrics.vcf.gz"

# Constants
qtl_distance          = 1e6    # Distance threshold between each gene/peak and variants
maf_threshold         =   0.01 # MAF threshold for variants used for QTL analysis
phenotype_min_value   =   0.5  # Threshold to consider a gene/peak expressed
phenotype_min_samples =   0.1  # Fraction of samples that have expression above phenotype_min_value

# For QTL analysis, divided by assay:
#vars0_rna = c("gt" , "population1", "age_sample", "sex", paste("PC", 1:10, sep = ""), "(1|wgs_id)", "(1|family_id)") # list of variants for LMM formula
#vars1_rna = c("gt:population1", "gt:age_sample", "gt:population1:age_sample")                                        # list of interaction terms to test vs LMM without interactions
