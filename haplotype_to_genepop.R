### Stacks haplotype conversion facility ###
# Mark Ravinet, University of Gothenburg - October 2013

rm(list = ls())
library(stringr)
source("~/Dropbox/R/mark_source/haplotype_conversion_functions.R")

### NB - for now this function works on haplotype file output from POPULATIONS module of the stacks
### pipeline

# Specify command line arguments
args <- commandArgs(TRUE)
infile <- read.delim(args[1])
outfile <- args[2]
header <- args[3]

# Assign RAD_loci names
loci <- infile$Catalog.ID
cat(sprintf("Converting %s\n", header))
cat(sprintf("Haplotype data for %s loci read in\n", length(loci)))

# Create a haplotype matrix
haplo <- file_convert(infile)

# First sweep through haplo for any individuals with > 2 alleles
haplo <- apply(haplo, 2, function(x){
  sapply(x, function(y){
    if(str_count(y, "/") > 1){
      y <- "-"
    } else{
      y <- y
    }
  })
})

# remove all loci where an allele is consensus
# and there are only two alleles
consensus <- loci[apply(haplo, 2, function(x){
  "consensus" %in% x & length(unique(x)) == 2
})]

# then remove all loci where an allele is "-"
# and there is only a single allele otherwise
nonpoly_null <- loci[apply(haplo, 2, function(x){
  "-" %in% x && length(unique(x)) == 2
})]

# perform a second scan for nonpoly loci that are not
# covered by either of the previous two filters
nonpoly_other <- loci[apply(haplo, 2, function(x){
  length(unique(x)) == 1
})]

# a vector of loci to remove
remove <- sort(unique(c(consensus, nonpoly_null, nonpoly_other)))

# remove incompatible loci
# or leave alone if necessary
if(length(remove > 0)){
  haplo <- haplo[, -which(colnames(haplo) %in% remove)]
  # update loci vector
  loci <- loci[-which(loci %in% remove)]
} else {
  haplo <- haplo
  loci <- loci
}

# Print to screen how many removed
cat(sprintf("A total of %s loci removed\n", length(remove)))
cat(sprintf("%s loci to be written to file\n", length(loci)))

# Split data into alleles
alleles <- apply(haplo, 2, strsplit, "/")

# Produce a list of alleles - copy them twice too
alleles_list <- lapply(alleles, function(x){
  lapply(x, function(y){
    if(length(y) < 2){
      y <- c(y, y)
      } else{
        y <- y
      }
    })
  })

# Scan alleles list for individuals with > 3 alleles
# remove them from the dataset
# alleles_list <- lapply(alleles_list, function(x){
#   lapply(x, function(y){
#     if(length(y) != 2){
#       y <- c(y[1] <- "-", y[2] <- "-") 
#     } else {
#       y <- c(y[1], y[2])
#     }
#   })
# })

# Convert to a dataframe
allele_df<- lapply(alleles_list, function(x){
  x <- do.call(rbind, x)
  x <- data.frame(x)
  data.frame(lapply(x, function(y){
     levels(y)[levels(y) == "-"] <- NA
     y
  }))
})

# Convert to genepop style
new_all <- lapply(allele_df, function(x){
  y <- as.numeric(unlist(x))
  y[is.na(y)] <- 0
  y <- formatC(y, width = 3, flag = "0")
  matrix(y, ncol = 2)
})

# Combine "alleles"
comb_all <- lapply(new_all, function(x){
  paste(x[, 1], x[, 2], sep = "")
})

# sort out ind
ind <- paste(rownames(haplo), ",", sep = "")
# Set pop by using first two characters of ind name
pop <- factor(substr(ind, 1, 2))

# produce a dataframe - all ind all loci
hap_mat <- data.frame(ind, do.call(cbind, comb_all))
names(hap_mat) <- NULL

# split that dataset by population
pop_data <- split(hap_mat, pop)

# Now write the output 
newfile <- file(outfile, "w")
cat(header, "\n", file = newfile, append = TRUE)
cat(loci, sep = "\n", file = newfile, append = TRUE)
for(i in 1:length(pop_data)){
  cat("pop\n", file = newfile, append = TRUE)
  write.table(pop_data[[i]], file = newfile, append = TRUE, col.names = FALSE,
              row.names = FALSE, sep = " ", quote = FALSE)
}
close(newfile)

cat("Done\n\n")






