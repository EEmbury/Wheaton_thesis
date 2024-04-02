rm(list = ls())
###### iNEXT attempt ######
#install.packages("iNEXT")

library(iNEXT)

#set otu table in iNEXT format

otu_inext_form <- as.abucount(otu_data)

#calculate coverage

coverage <- iNEXT(otu_data, q = 0, datatype = "abundance")

#calculate estimate d
estimate_d <- estimateD(otu_data, datatype = "abundance", base = "coverage")


###### rarefy #######

estimate_d_values <- read.csv("estimate_d.csv", header = TRUE)

t_otu <- t(otu_data)

##cl20r -- 307

cl20r <- otu_data[c("CL20r")]
typeof(cl20r)
final_cl20r <- as.data.frame(t(cl20r))
typeof(final_cl20r)
CL20r_rare <- rrarefy(final_cl20r, sample = 307)

# for loop

vectooor <- c(1)
for (val in 2:314){
  vectooor <- append(vectooor,val)
}

for (q in 1:314) {
  table_val <- estimate_d_values[q,1]
  otu <- otu_data[, table_val]
  typeof(table_val)
  t_val <- as.data.frame(t(otu))
  typeof(estimate_d_values[q, 2])
  rare <- rrarefy(t_val, sample = (estimate_d_values[q, 2]))
  t_rare <- t(rare)
  art <- as.data.frame(t_rare)
  vectooor[q] <- art
  }

vectar <- as.data.frame(vectooor)

write.csv(vectar, "test.csv")

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


