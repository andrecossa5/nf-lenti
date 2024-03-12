# DecontX 

# Code
library(readr)
library(celda)
library(SingleCellExperiment)
library(Matrix)
library(umap)
library(ggplot2)


##


# path
path_main <- '/Users/ieo6943/Documents/data/8_clones/count_matrix_&_DecontX'


##



# Read the CSV file into a data frame for different coverage threshold
for (coverage_threshold in seq(0, 100, by=20)) {
  
  count_matrix <- sprintf("M_%s.csv", coverage_threshold)
  path_in <- file.path(path_main, count_matrix)
  data <- read.csv(path_in, header=TRUE, row.names=1)
  data<-as.matrix(data)
  
  #Decontamination
  sce_dec<-decontX(t(data))

  #Decontaminated count matrix and saving
  m_dec<-round(sce_dec$decontXcounts)
  dense_matrix <- as.matrix(m_dec)
  dense_matrix <- as.matrix(m_dec)
  df <- as.data.frame(t(dense_matrix))
  count_matrix_dec <- paste("M_", as.character(coverage_threshold), "_dec.csv", sep = "")
  output_path <- file.path(path_main, count_matrix_dec)
  write.csv(df, file =output_path, row.names = FALSE)
  
  #Contamination
  contamination=sce_dec$estimates$all_cells$contamination
  contamination_df<-as.data.frame(contamination)
  contamination_file <- sprintf("contamination_%s.csv", coverage_threshold)
  path_contamination <-file.path(path_main, contamination_file)
  write.csv(contamination_df, file =path_contamination, row.names = FALSE)
  
  
  ### UMAP and CLUSTERING 
  m_umap <- sce_dec$estimates$all_cells$UMAP
  umap_file <- sprintf("umap_%s.csv", coverage_threshold)
  path_umap <- file.path(path_main, umap_file)
  write.csv(m_umap, file=path_umap, row.names = FALSE)
  
  z <- sce_dec$estimates$all_cells$z
  z_file <- sprintf("z%s.csv", coverage_threshold)
  path_z <- file.path(path_main, z_file)
  write.csv(z, file=path_z , row.names = FALSE)
  
}


##