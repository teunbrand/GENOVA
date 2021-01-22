# Declare files -----------------------------------------------------------

# Data from the Haarhuis *et al*. (2017) paper.
# Raw data is available at the Gene Expression Omnibus (GEO) accession GSE95014.
# Specifically, the file paths below are merged data from the following samples 
# noted by their name on the GEO accession:
# Hap1-A Hi-C, 
# Hap1-B Hi-C, 
# Hap1-C Hi-C.
# Raw sequence data was processed and mapped against the hg19 human genome build
# using HiC-Pro v2.7.7

base_dir <- "/shared/dewit/oidBackUp/WAPL_Project/Hi-C/matrices/combined/WT/combined/output/hic_results/matrix/WT/"

signal_150k  = paste0(base_dir, "iced/150000/WT_150000_iced.matrix")
indices_150k = paste0(base_dir, "raw/150000/WT_150000_abs.bed")

signal_40k   = paste0(base_dir, "iced/40000/WT_40000_iced.matrix")
indices_40k  = paste0(base_dir, "raw/40000/WT_40000_abs.bed")

# Data import -------------------------------------------------------------

library(GENOVA)

exp_150k <- load_contacts(
  signal_path  = signal_150k,
  indices_path = indices_150k,
  sample_name  = "Hap1_WT_150k",
  colour = "dodgerblue"
)

exp_40k <- load_contacts(
  signal_path  = signal_40k,
  indices_path = indices_40k,
  sample_name  = "Hap1_WT_40k",
  colour = "tomato"
)

# Subsetting data ---------------------------------------------------------

library(data.table)

subset_exp <- function(exp, chromosomes = NULL) {
  match <- match(chromosomes, unique(exp$IDX$V1))
  chromosomes <- chromosomes[!is.na(match)]
  
  if (is.null(chromosomes) || length(chromosomes) == 0) {
    return(exp)
  }
  
  # Subset relevant components
  idx <- exp$IDX[chromosomes]
  
  mat <- exp$MAT[V1 %in% idx$V4 & V2 %in% idx$V4]
  
  centro <- exp$CENTROMERES[chrom %in% chromosomes]
  
  chrs <- exp$CHRS[exp$CHRS %in% chromosomes]
  
  # We need to unset data.table because I imagine `.internal.selfref` attribute 
  # can cause problems.
  setDF(idx)
  setDF(mat)
  setDF(centro)
  
  out <- structure(list(
    MAT = mat,
    IDX = idx,
    CHRS = chrs,
    CENTROMERES = centro
  ))
  attributes(out) <- attributes(exp)
  out
}

exp_150k <- subset_exp(exp_150k, c("chr21", "chr22"))
exp_40k  <- subset_exp(exp_40k,  c("chr21", "chr22"))

# Saving data -------------------------------------------------------------

saveRDS(object = exp_150k, here::here("data", "test_150k.rds"))
saveRDS(object = exp_40k,  here::here("data", "test_40k.rds"))

# Checksums ---------------------------------------------------------------

tools::md5sum(here::here("data", "test_150k.rds"))
# > "8b89fa2344a544ec01f7218d8dcbde86"

tools::md5sum(here::here("data", "test_40k.rds"))
# > "6d011f0fe5632f24bc97497d26df791b"