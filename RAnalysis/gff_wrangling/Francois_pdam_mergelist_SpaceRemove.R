# reading data in - list of paths and sampleIDs gtf
my_data <- read.delim("~/Desktop/pdam_mergelist.txt", header = FALSE)

# Subbing the space in each line for a '/'
my_data$V1 <- sub(" ", "/", my_data$V1)
head(my_data)

# Write out new file 
write_delim(my_data, "~/Desktop/pdam_mergelist_no_spaces.txt")
