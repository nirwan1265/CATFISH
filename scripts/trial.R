# # Packages
# #BiocManager::install(version = "3.14")
# #BiocManager::install(c("GenomicRanges", "GenomeInfoDb", "TailRank", "IRanges", "Gviz"))
# #install.packages("JuliaCall")
# library(JuliaCall)
# library(dplyr)
# #install.packages("RTIGER")
# library(RTIGER)

# # 2. Code Organization                                                                                                                          
                                                                                                                                                
# #   - Split the large OMNIBUS_test.R (~1400 lines) into smaller logical files                                                                     
# #   - Standardize naming conventions (some use _, some use camelCase)                                                                             
# #   - Add consistent input validation across wrapper functions                                                                                    
                                                                                                                                                
  
                                                                                                                                                
#   # 4. Vignettes                                                                                                                                  
                                                                                                                                                
#   # - Create a "Getting Started" vignette                                                                                                         
#   # - Create a methods/statistical background vignette                                                                                            
#   # - Add worked examples with sample data                                                                                                        
                                                                                                                                                
#   # 5. CRAN/Bioconductor Readiness                                                                                                                
                                                                                                                                                
#   # - Run R CMD check and fix any NOTEs/WARNINGs                                                                                                  
#   # - Ensure all dependencies are properly declared                                                                                               
#   # - Add @importFrom statements to reduce namespace pollution                                                                                    
                                                                                                                                                
#   # 6. Error Handling                                                                                                                             
                                                                                                                                                
#   # - Add informative error messages with cli or rlang                                                                                            
#   # - Validate inputs at function entry points                                                                                                    
                                                                    

# ### SETUP
# # Done once
# #setupJulia(JULIA_HOME="/Applications/Julia-1.0.app/Contents/Resources/julia/bin")
# setupJulia(JULIA_HOME="/Applications/Julia-1.10.app/Contents/Resources/julia/bin")
# # Needs to be run everytime we load RTIGER
# sourceJulia()


# library(doParallel)
# library(foreach)
# library(doParallel)
# library(foreach)
# library(tools)
# library(RTIGER)


# # ----- RTIGER constants -----
# chr_len <- c(308452471,243675191,238017767,250330460,226353449,
#              181357234,185808916,182411202,163004744,152435371)
# names(chr_len) <- paste0("chr", 1:10)


# # chr_len <- c(243675191)
# # names(chr_len) <- paste0("chr", c(2))

# post_post.processing <- TRUE


# in_dir  <- "/Users/nirwantandukar/Documents/Research/data/BZea/angsd_genotyping/allele_counts"
# files   <- list.files(in_dir, pattern="\\.tsv$", full.names=TRUE)
# files <- files[c(51:100)]

# expDesign <- data.frame(
#   files = files,
#   name  = sub("\\.rtiger.tsv$", "", basename(files)),
#   stringsAsFactors = FALSE
# )

# outdir <- "/Users/nirwantandukar/Documents/Research/data/BZea/rtiger_results/rtiger_results_angsd"

# myres <- RTIGER(
#   expDesign = expDesign,
#   outputdir = outdir,
#   seqlengths = chr_len,
  
#   # START CONSERVATIVE for your pruned marker density
#   rigidity = 100,
#   # Best rigidity = 512
  
#   # Let RTIGER tune (but give realistic depth)
#   autotune = F,
#   average_coverage = 0.8,
#   crossovers_per_megabase = 0.05,
  
#   nstates = 3,
#   post.processing = post_post.processing,
#   save.results = TRUE,
#   verbose = TRUE
# )