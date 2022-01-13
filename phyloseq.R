# https://joey711.github.io/phyloseq/import-data.html#_microbio_me_qiime_(defunct)

# load packages 
library(phyloseq)
library(ggplot2)
library(vegan)

### all data ####
#read data
otudata <- read.csv("rel_abund_otu_transpose.csv.csv")
otutax <- read.csv("new_taxa.csv")
samdata <- read.csv("Metadata2.csv")


library(tidyverse)

otudata$X <- sub(otudata$X, 
                   pattern = "OTU", replacement = "")

otu_data <- otudata %>% remove_rownames %>% column_to_rownames(var="X")




otutax$OTU_ID <- sub(otutax$OTU_ID, 
                  pattern = "OTU", replacement = "")

otu_taxa <- otutax %>% remove_rownames %>% column_to_rownames(var="OTU_ID")

#coerce matrix

taxa_matrix <- as.matrix(otu_taxa, rownames.force = NA)

# setting data to phyloseq format

library(phyloseq)
OTU <- otu_table(otu_data, taxa_are_rows = TRUE)
TAX <-  tax_table(taxa_matrix)
SAM <- sample_data(samdata)
OTU
TAX
SAM

AllData <- phyloseq(OTU, TAX)

#mystery answer to error in validobject

rownames(samdata) <-samdata$InputFileName


physeq1 <- merge_phyloseq(AllData, SAM)


#setting up ggplot themes

theme_set(theme_bw())
pal = "Set1"
scale_colour_discrete <-  function(palname=pal, ...){
  scale_colour_brewer(palette=palname, ...)
}
scale_fill_discrete <-  function(palname=pal, ...){
  scale_fill_brewer(palette=palname, ...)
}

#fix the random column

df = subset(otu_data, select = -c(X) )


#Plotting

p <- plot_richness(physeq1, x="Treatment", color="Treatment", measures=c("Shannon"))
p + geom_point(size=1, alpha=0.05)



#### Taxa Reorganization ####

taxa <- read.csv("taxa.csv")

library(tidyr)
y <- taxa %>% separate(Taxonomy, c("kingdom", "phylum", "class", "order", "family", "genus"), "[,][a-z ]:")

y$kingdom <- sub(y$kingdom,  pattern = "k:", replacement = "")

y$phylum <- sub(y$phylum, pattern = ",p:", replacement = "")

y$class <- sub(y$class, pattern = ",c:", replacement = "")

y$order <- sub(y$order, pattern = ",o:", replacement = "")

y$family <- sub(y$family, pattern = ",f:", replacement = "")

y$genus <- sub(y$genus,  pattern = ",g:", replacement = "")


write.csv(y,"new_taxa.csv", row.names = FALSE)

