setwd("C:/Users/gabri/Documents/Thom_Lab_Project/R_confusion_matrix")
library(readxl)
library(reshape)
library(ggplot2)
# 1. Creating list of confusion matrices
file_list <- list.files(
  pattern = "confusion_matrix\\.csv$", 
  full.names = TRUE
)

confusion_matrices <- list()

for (file in file_list) {
  clean_name <- tools::file_path_sans_ext(basename(file))
  
  confusion_matrices[[clean_name]] <- read.csv(file)
}

#Cleaning Up
confusion_matrices <- lapply(confusion_matrices, \(x) x[, -1, drop = FALSE])

confusion_matrices <- lapply(confusion_matrices, \(x) {
  colnames(x) <- c("Model 1", "Model 2", "Model 3", "Model 4")
  x
})
confusion_matrices <- lapply(confusion_matrices, \(x) {
  rownames(x) <- c("Model 1", "Model 2", "Model 3", "Model 4")
  x
})

confusion_matrices <- lapply(confusion_matrices, \(x) {
 as.matrix(x)
})

#Creating Heat Maps
confusion_matrices_melted <- lapply(confusion_matrices, \(x) {
  melt(x)
})

confusion_matrices_melted <- lapply(confusion_matrices_melted, \(x) {
  colnames(x) <- c("Pseudo-Observed", "Estimated", "value")
  x
})

heat_maps <- lapply(confusion_matrices_melted, \(x) {
  ggplot(x, aes(x = `Pseudo-Observed`,
                        y = Estimated,
                        fill = value)) +                          
    geom_tile() +
    scale_fill_gradient(low = "yellow", high = "red") +
    coord_fixed()
})

for (nm in names(heat_maps)) {
  ggsave(
    filename = paste0(nm, ".pdf"),
    plot = heat_maps[[nm]],
    width = 8,
    height = 6
  )
}
#Combining matrices
combined_confusion_matrix <- Reduce("+", confusion_matrices)

#Creating Combined Heat Map
combined_melt <- melt(combined_confusion_matrix)                                           # Reorder data
colnames(combined_melt) <- c("Pseudo-Observed", "Estimated", "value")
 
combined_heat_map <- ggplot(combined_melt, aes(x = `Pseudo-Observed`,
                            y = Estimated,
                            fill = value)) +                           # Create heatmap with ggplot2
  geom_tile() +
  scale_fill_gradient(low = "yellow", high = "red") +
  coord_fixed()

ggsave(filename = "combined_confusion_matrix.pdf",
       plot = combined_heat_map,
    width = 8,
    height = 6
  )


