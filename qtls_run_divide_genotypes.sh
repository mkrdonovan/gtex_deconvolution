#!/usr/bin/sh

/frazer01/home/matteo/software/R-3.5.1/bin/Rscript /frazer01/projects/CARDIPS/analysis/eqtls_deconvolution_gtex/analysis/cardiac_qtls_run_divide_genotypes.R -e /frazer01/projects/CARDIPS/analysis/eqtls_deconvolution_gtex/log/divide_genotypes.err -o /frazer01/projects/CARDIPS/analysis/eqtls_deconvolution_gtex/log/divide_genotypes.out --taskid $SGE_TASK_ID
