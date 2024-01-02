# conda environment is /zfs/omics/projects/metatools/TOOLS/dadasnake/conda/38b38ff488abf3036df0ac8c2c2c369e 
condap <- Sys.getenv("CONDA_PREFIX")
.libPaths(paste0(condap,"/lib/R/library"))
require(Tax4Fun2)
#only required the first time around:
#system(paste0("mkdir -p /zfs/omics/projects/metatools/DB/amplicon/Functions/tax4fun2/Tax4Fun2_ReferenceData_v2/TOOLS/blast_bin/bin/"))
#system(paste0("ln -s ",condap,"/bin/makeblastdb /zfs/omics/projects/metatools/DB/amplicon/Functions/tax4fun2/Tax4Fun2_ReferenceData_v2/TOOLS/blast_bin/bin/makeblastdb"))

#set every time:
inputSeqs <- "TIFN_Tax4Fun2_rarefied_10k_OTUs.fa"
inputTab <- "TIFN_Tax4Fun2_rarefied_10k.RDS"
outputFolder <- "TIFN_Tax4Fun2_rarefied_10k"
database <- "Ref99NR"
reference_database_folder <- "/zfs/omics/projects/metatools/DB/amplicon/Functions/tax4fun2/Tax4Fun2_ReferenceData_v2"

#this can stay as it is:
threads <- 1
include_user_data <- T
path_to_user_data <- "/zfs/omics/projects/metatools/DB/amplicon/Functions"
name_of_user_data <- "HOMD_cut"
runRefBlast(path_to_otus = inputSeqs, 
            path_to_reference_data = reference_database_folder,
            path_to_temp_folder = outputFolder, 
            database_mode = database,
            use_force = T,
            num_threads = threads,
            include_user_data = include_user_data,
            path_to_user_data = path_to_user_data,
            name_of_user_data = name_of_user_data)

source("/zfs/omics/projects/metatools/TOOLS/dadasnake/workflow/scripts/functionalPredictionCustom.R")

makeFunctionalPredictionCustom(inputTab,
                               path_to_reference_data = reference_database_folder,
                               path_to_temp_folder = outputFolder,
                               database_mode = database,
                               normalize_by_copy_number =T, 
                               min_identity_to_reference = 0.97,
                               normalize_pathways = FALSE,
                               include_user_data = include_user_data,
                               path_to_user_data = path_to_user_data,
                               name_of_user_data = name_of_user_data)

calculateFunctionalRedundancy(path_to_otu_table = inputSeqs,
                              path_to_reference_data = reference_database_folder,
                              path_to_temp_folder = outputFolder,
                              database_mode = database)
