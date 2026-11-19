library(LEA)
library(devtools)
library(tess3r)

setwd("c:/Users/gabri/Documents/R_VCF/3_vcfs/filtered/")
files <- list.files(path="./", pattern=".recode.vcf", full.names=F, recursive=FALSE)

files <- files[grep(".vcf.", files, invert = T)]

for (file in files) {
  result <- snmf(paste(file, ".geno", sep = ""), 
                 K = 1:3, 
                 ploidy = 2, 
                 entropy = TRUE, 
                 CPU = 5, 
                 project = "new", 
                 alpha = 1, 
                 repetitions = 250)
  
  assign(paste(file, ".snmf.a1", sep = ""), result)
}
for (file in files) {
  result <- snmf(paste(file, ".geno", sep = ""), 
                 K = 1:3, 
                 ploidy = 2, 
                 entropy = TRUE, 
                 CPU = 5, 
                 project = "new", 
                 alpha = 10, 
                 repetitions = 250)
  
  assign(paste(file, ".snmf.a10", sep = ""), result)
}

for (file in files) {
  result <- snmf(paste(file, ".geno", sep = ""), 
                 K = 1:3, 
                 ploidy = 2, 
                 entropy = TRUE, 
                 CPU = 5, 
                 project = "new", 
                 alpha = 100, 
                 repetitions = 250)
  
  assign(paste(file, ".snmf.a100", sep = ""), result)
}

for (file in files) {
  result <- summary(file, ".snmf.a1", sep = "")
  assign(paste(file, ".snmf.a1.summary", sep = ""), result)
}
# Loop over your files to summarize each snmf result and assign it to a new object
for (file in files) {
  
  # Construct the object name from the file name (this should match the names created in your first loop)
  object_name <- paste(file, ".geno.snmf.a1", sep = "")
  
  # Use get() to access the dynamically created object
  snmf_obj <- get(object_name)
  
  # Compute the summary of the snmf object
  summary_result <- summary(snmf_obj)
  
  # Assign the summary to a new object, using the file name
  assign(paste(file, "_a1_summary", sep = ""), summary_result)
  
  # Optionally, store it in a list as well if you want to collect all summaries in one object
  # summaries[[file]] <- summary_result
}

# Define alpha values for your objects
alpha_values <- c(1, 10, 100)

# Loop over each file and alpha value
for (file in files) {
  for (alpha_val in alpha_values) {
    # Construct the original object name
    original_name <- paste(file, ".snmf.a", alpha_val, sep = "")
    
    # Generate summary and save it as a new object
    summary_result <- summary(get(original_name))
    
    # Assign the summary result to a new object name
    summary_name <- paste(original_name, ".summary", sep = "")
    assign(summary_name, summary_result)
  }
}

# Loop over each file to combine and save crossEntropy values
for (file in files) {
  # Initialize an empty list to store crossEntropy data for each alpha
  crossEntropy_list <- list()
  
  # Loop through each alpha value
  for (alpha_val in alpha_values) {
    # Construct the summary object name
    summary_name <- paste(file, ".snmf.a", alpha_val, ".summary", sep = "")
    
    # Extract crossEntropy values
    crossEntropy_data <- get(summary_name)$crossEntropy
    
    # Add a column for alpha value
    crossEntropy_data <- as.data.frame(crossEntropy_data)
    crossEntropy_data$alpha <- alpha_val
    
    # Store in list
    crossEntropy_list[[as.character(alpha_val)]] <- crossEntropy_data
  }
  
  # Combine all crossEntropy data into one data frame
  combined_crossEntropy <- do.call(rbind, crossEntropy_list)
  
  # Write to CSV, with a unique filename for each file
  output_filename <- paste(file, "_crossEntropy.csv", sep = "")
  write.csv(combined_crossEntropy, file = output_filename, row.names = FALSE)
}

###"Species" is used as a default name in this code
write.csv(rbind(species_0.7_DP4.recode.vcf.snmf.a1.summary$crossEntropy,species_0.7_DP4.recode.vcf.snmf.a10.summary$crossEntropy,species_0.7_DP4.recode.vcf.snmf.a100.summary$crossEntropy), paste(file, "species_cross_entropy.csv"))

species_cross_entropy <- read.csv("C:/Users/gabri/Documents/R_VCF/species_cross_entropy.csv")

#Selecting runs w/ lowest cross-entropy

species_ce = cross.entropy(species_0.7_DP4.recode.vcf.snmf.a1, K = 1)
best_species = which.min(species_ce) #species: K1, A1

## Indivuals for each species

ind_species <- read.csv("C:/Users/gabri/Documents/R_VCF/3_vcfs/filtered/species_FINAL.recode.vcf.csv")


##Plot
barplot(t(Q(species_0.7_DP4.recode.vcf.snmf.a1, K = 1, run = best_species)))
Q.matrix = Q(species_0.7_DP4.recode.vcf.snmf.a1, K = 1, run = best_species)
x <- cbind(ind_species, Q.matrix)
x <- x[order(x$V1),]
y <- cbind(x[,4])

barplot(t(y), col=c("coral","aquamarine3"),border = NA, space = .05,
        ylab = "Admixture coefficients", main = "Ancestry matrix", horiz = FALSE, 
        names.arg = x[,1] ,cex.names=1.0, las = 2)

pdf("species_k1.pdf", width=18, height=9)
par(mar=c(9.1,5.1,4.1,3.1))
barplot(t(y), col=c("coral","aquamarine3"),border = NA, space = .05,
        ylab = "Admixture coefficients", main = "Ancestry matrix", horiz = FALSE, 
        names.arg = x[,1] ,cex.names=1.0, las = 2)
dev.off()

# k2 runs with lowest cross entropy

species_ce_k2 = cross.entropy(species_0.7_DP4.recode.vcf.snmf.a1, K = 2)
best_species_k2 = which.min(species_ce_k2) #species: K2, A1

#Plotting K2
##Plot
barplot(t(Q(species_0.7_DP4.recode.vcf.snmf.a1, K = 2, run = best_species_k2)))
Q.matrix = Q(species_0.7_DP4.recode.vcf.snmf.a1, K = 2, run = best_species)
x <- cbind(ind_species, Q.matrix)
x <- x[order(x$V1),]
y <- cbind(x[,4], x[,5])

barplot(t(y), col=c("coral","aquamarine3"),border = NA, space = .05,
        ylab = "Admixture coefficients", main = "Ancestry matrix", horiz = FALSE, 
        names.arg = x[,1] ,cex.names=1.0, las = 2)

pdf("species_k2.pdf", width=18, height=9)
par(mar=c(9.1,5.1,4.1,3.1))
barplot(t(y), col=c("coral","aquamarine3"),border = NA, space = .05,
        ylab = "Admixture coefficients", main = "Ancestry matrix", horiz = FALSE, 
        names.arg = x[,1] ,cex.names=1.0, las = 2)
dev.off()
