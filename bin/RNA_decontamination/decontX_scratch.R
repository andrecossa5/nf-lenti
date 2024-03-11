# DecontX 

# Code
library(readr)
library(celda)
library(SingleCellExperiment)
library(Matrix)
library(umap)
library(ggplot2)


##


# Args
path_results <- '/Users/ieo6943/Documents/Guido/mito_preprocessing/bin/RNA_decontamination/results/M_'




# Read the CSV file into a data frame
for (coverage_threshold in seq(0, 100, by=20)) {
  
  path_in <- paste("/Users/ieo6943/Documents/Guido/mito_preprocessing/bin/RNA_decontamination/results/M_", as.character(coverage_threshold), ".csv", sep = "")
  #path<-paste("/Users/ieo6943/Documents/Guido/NTS/M_", as.character(coverage_threshold), ".csv", sep = "")
  #MM<- list()
  #for (i in seq(1,10)){
  data <- read.csv(path_in, header=TRUE, row.names=1)
  data<-as.matrix(data)
  

  
  ############## USANDO SEMPLICEMENTE LA MATRICE]
  sce_dec<-decontX(t(data))

  
  #output matrix decontaminated
  m_dec<-round(sce_dec$decontXcounts)
  dense_matrix <- as.matrix(m_dec)
  #MM[[i]]=dense_matrix
#}
# Convert dgCMatrix to regular matrix
  dense_matrix <- as.matrix(m_dec)
  
  # Convert the dense matrix to a data frame 
  df <- as.data.frame(t(dense_matrix))
  output_path<-paste("/Users/ieo6943/Documents/Guido/mito_preprocessing/bin/RNA_decontamination/results/M_", as.character(coverage_threshold), "_dec.csv", sep = "")
  #output_path<-paste("/Users/ieo6943/Documents/Guido/NTS/M_", as.character(coverage_threshold), "_dec.csv", sep = "")
  #write.csv(df, file =output_path, row.names = FALSE)
  
  contamination=sce_dec$estimates$all_cells$contamination
  contamination_df<-as.data.frame(contamination)
  path_contamination<-paste("/Users/ieo6943/Documents/Guido/mito_preprocessing/bin/RNA_decontamination/results/contamination_", as.character(coverage_threshold), ".csv", sep = "")
  #write.csv(contamination_df, file =path_contamination, row.names = FALSE)
  
  
  ### UMAP and CLUSTERING 
  
  m_umap <- sce_dec$estimates$all_cells$UMAP
  path_out<-"/Users/ieo6943/Documents/Guido/mito_preprocessing/bin/RNA_decontamination/results/"
  #path="/Users/ieo6943/Documents/Guido/NTS/"
  
  umap_path <- paste(path_out,"umap_",as.character(coverage_threshold), ".csv", sep = "")
  #write.csv(m_umap, file = umap_path , row.names = FALSE)
  
  z_path <- paste(path_out,"z",as.character(coverage_threshold), ".csv", sep = "")
  #write.csv(sce_dec$estimates$all_cells$z, file = z_path , row.names = FALSE)
  plotDimReduceCluster(x = sce_dec$estimates$all_cells$z,
                       dim1 = m_umap[, 1], dim2 = m_umap[, 2], )
}


##