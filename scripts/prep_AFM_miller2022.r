# Create X_bin from Allele frequency matrix as in Miller et al., 2022

library(tidyverse)
library(data.table)


##


## Load function from Peter https://github.com/petervangalen/MAESTER-2021/blob/main/Auxiliary_files/210215_FunctionsGeneral.R
# Function that computes all heteroplasmic variants from MAEGATK output (from Caleb Lareau).

computeAFMutMatrix <- function(SE){
  
  cov <- assays(SE)[["coverage"]]+ 0.000001
  ref_allele <- as.character(rowRanges(SE)$refAllele)
  
  getMutMatrix <- function(letter){
    mat <- (assays(SE)[[paste0(letter, "_counts_fw")]] + assays(SE)[[paste0(letter, "_counts_rev")]]) / cov
    rownames(mat) <- paste0(as.character(1:dim(mat)[1]), "_", toupper(ref_allele), ">", letter)
    return(mat[toupper(ref_allele) != letter,])
  }
  
  rbind(getMutMatrix("A"), getMutMatrix("C"), getMutMatrix("G"), getMutMatrix("T"))
}


##


cutf <- function(x, f=1, d="/") sapply(strsplit(x, d), function(i) paste(i[f], collapse=d))


##


# Read data
path <-'/Users/IEO5505/Desktop/MI_TO/MI_TO_analysis_repro/data/MI_TO_bench/AFMs/maegatk/MDA_PT/tables/'
maegatk <- readRDS(paste0(path, 'AFM.rds'))
af.dm <- suppressWarnings(data.matrix(computeAFMutMatrix(maegatk))*100)

af.dm %>% dim
af.dm[1:10,1:10]

# Compute stats, extract variants tables, and binarize (af>=t, or >0)

# Qualities
assays.ls <- suppressWarnings(lapply(maegatk@assays@data, function(x) as.matrix(x)))
qual.num <- sapply(rownames(af.dm), function(x) {
  #x <- "2141_T>C"
  pos <- as.numeric( cutf(x, d = "_") )
  mut <- cutf(x, d = ">", f = 2)
  # Get the mean quality of reads for this call (only use cells in which the base was sequenced) - forward
  covered_fw <- assays.ls[[str_c(mut, "_counts_fw")]][pos,] > 0
  qual_fw <- assays.ls[[str_c(mut, "_qual_fw")]][pos, covered_fw]
  # Same for reverse
  covered_rev <- assays.ls[[str_c(mut, "_counts_rev")]][pos,] > 0
  qual_rev <- assays.ls[[str_c(mut, "_qual_rev")]][pos, covered_rev]
  qual <- mean(c(qual_fw, qual_rev))
  return(qual)
})

# Make tibble
vars.tib <- tibble(var = rownames(af.dm),
                   mean_af = rowMeans(af.dm),
                   mean_cov = rowMeans(assays(maegatk)[["coverage"]])[as.numeric(cutf(rownames(af.dm), d = "_"))],
                   quality = qual.num)

# Calculate the number of cells that exceed VAF thresholds 0, 1, 5, 10, 50 (3 minutes)
vars.tib <- vars.tib %>%
  mutate(n0 = apply(af.dm, 1, function(x) sum(x == 0))) %>%
  mutate(n1 = apply(af.dm, 1, function(x) sum(x > 1))) %>%
  mutate(n5 = apply(af.dm, 1, function(x) sum(x > 5))) %>%
  mutate(n10 = apply(af.dm, 1, function(x) sum(x > 10))) %>%
  mutate(n50 = apply(af.dm, 1, function(x) sum(x > 50)))

Variant_CellN<-apply(af.dm,1,function(x){length(which(x>0))})
vars.tib<-cbind(vars.tib,Variant_CellN)

# Extract X
n1.2 <- vars.tib %>% filter(mean_cov >= 10, quality >= 30, n0 > 0.2*ncol(af.dm), n1>2) %>% .$var
n1.5 <- vars.tib %>% filter(mean_cov >= 10, quality >= 30, n0 > 0.2*ncol(af.dm), n1>5) %>% .$var
n5.2 <- vars.tib %>% filter(mean_cov >= 10, quality >= 30, n0 > 0.2*ncol(af.dm), n5>2) %>% .$var
n5.5 <- vars.tib %>% filter(mean_cov >= 10, quality >= 30, n0 > 0.2*ncol(af.dm), n5>5) %>% .$var
n10.2 <- vars.tib %>% filter(mean_cov >= 10, quality >= 30, n0 > 0.2*ncol(af.dm), n10>2) %>% .$var
n10.5 <- vars.tib %>% filter(mean_cov >= 10, quality >= 30, n0 > 0.2*ncol(af.dm), n10>=5) %>% .$var
var.list <- list(n1.2, n1.5, n5.2, n5.5, n10.2, n10.5)
names(var.list) <- c("n1.2", "n1.5", "n5.2", "n5.5", "n10.2", "n10.5")

for (name in names(var.list)) {
  write.csv(af.dm[var.list[[name]],] %>% t, paste0(path, name, ".csv"))
}


##

