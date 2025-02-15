### --------------------------------------------------------------------------------------------------
#   HEATMAP - PRESENCE & ABSENCE
#
#   Author : Eduarda Valério de Jesus 
#   Last updated 2024/02/14
#
#   Script to create a heatmap based in binary data (presence/absence)
### -------------------------------------------------------------------------------------------------- 


# Install pheatmap if not already installed
if (!require(pheatmap)) install.packages("pheatmap")
library(pheatmap)

# Load data frame 
data <- read.csv("draft - heatmap_porifera.csv", row.names=1)

# Create a matrix 
data_matrix <- data.matrix(data)
mybinarymap <- heatmap(data_matrix, Rowv = NA, Colv = NA, col = c("white", "black"), scale = "none")

# Create the heatmap
pheatmap(
  data_matrix,
  color = c("gray", "black"),      # Binary color scheme
  cluster_rows = FALSE,            # Disable row clustering
  cluster_cols = FALSE,            # Disable column clustering             
  fontsize_row = 4,                # Row label size
  fontsize_col = 10,               # Column label size 
  labels_col = c("COI NCBI", "COI BOLD", "18S SILVA", "18S NCBI", "28S NCBI"),  # Add custom column labels
  angle_col = 315,
  legend = FALSE                   # Disable legend 
)

# Save graph in TIFF file 
dev.copy(tiff, filename = "heatmap_porifera.tiff", width = 8, height = 15, units = "in", res = 600)
dev.off() 



