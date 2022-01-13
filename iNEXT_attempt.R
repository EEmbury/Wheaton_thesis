
###### iNEXT attempt ######
#install.packages("iNEXT")

library(iNEXT)

#set otu table in iNEXT format

otu_inext_form <- as.abucount(otu_data)

#calculate coverage

coverage <- iNEXT(otu_data, q = 0, datatype = "abundance")

#calculate estimate d
estimate_d <- estimateD(otu_data, datatype = "abundance", base = "coverage")

##### phyloseq ######

library(phyloseq)

#read data
otudata <- read.csv("otu_refined.csv", header = TRUE)
otutax <- read.csv("taxa_metadata.csv")
samdata <- read.csv("Metadata_updated.csv")


#change row labels
#install.packages("tidyverse")
library(tidyverse)

otudata$OUT.ID <- sub(otudata$OUT.ID, 
                        pattern = "OTU_", replacement = "")

otu_data <- otudata %>% remove_rownames %>% column_to_rownames(var="OUT.ID")


otutax$OTU_ID <- sub(otutax$OTU_ID, 
                  pattern = "OTU", replacement = "")

otu_taxa <- otutax %>% remove_rownames %>% column_to_rownames(var="ID")

rownames(samdata) <-samdata$InputFileName

#coerce matrix (taxa needs to be in a matrix format)

test <- as.matrix(otu_taxa, rownames.force = NA)

# setting data to phyloseq format

library(phyloseq)
OTU <- otu_table(otu_data, taxa_are_rows = TRUE)
TAX <-  tax_table(test)
SAM <- sample_data(samdata)
OTU
TAX
SAM

# combine physeq tables into 1 unit
table <- phyloseq(OTU, TAX)

physeq1 <- merge_phyloseq(table, SAM)
