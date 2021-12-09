# import data
otu_table <- read.delim("otu_table_no_hem.txt")

#drop columns

otu <- otu_table[ , ! names(otu_table) %in% c("label", "numOtus")]

#Set coulumn as label
#install.packages("tidyverse")
library(tidyverse)
otu_data <- otu %>% remove_rownames %>% column_to_rownames(var="Group")

#rarefy data
#https://search.r-project.org/CRAN/refmans/GUniFrac/html/Rarefy.html
#install.packages("GUniFrac")
#library(GUniFrac)

#rare_otu <- Rarefy(otu_data,374)

# relative abundance
#https://rdrr.io/cran/funrar/man/make_relative.html

install.packages("funrar")
library(funrar)

otu_matrix <- as.matrix(otu_data, rownames.force = NA)

relabun_otu <- make_relative(otu_matrix)

write.csv(relabun_otu,"rel_abund_otu.csv", row.names = TRUE)