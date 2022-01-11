
#####1st attempt -- older otu table formats#######
#read otu table
otu_table <- read.table("otu_table_no_hem.txt", header = TRUE)


##Remove extra columns
df_bop <- otu_table[ -c(1,3) ]

##rename rows
rownames(df_bop) <-df_bop$Group

df_bop_bop <- df_bop[ -c(1)]



###### iNEXT attempt ######
install.packages("iNEXT")

library(iNEXT)


OTU_t <- read.csv("t_otu_updated.csv", header = TRUE)

rownames(OTU_t) <-OTU_t$X
OTU <- OTU_t[ -c(1)]


otu_inext_form <- as.abucount(OTU)

coverage <- iNEXT(otu_inext_form, q = 0, datatype = "abundance", size = NULL, 
      endpoint = NULL, knots = 40, se = TRUE, conf = 0.95, nboot = 50)

plot(coverage, type = 1, se = TRUE, show.legend = TRUE,
     show.main = TRUE, col = NULL)

estimate_d <- DataInfo(otu_inext_form, datatype = "abundance")



##### phyloseq ######


library(phyloseq)

#read data
otudata <- read.csv("t_otu.csv", header = TRUE)
otutax <- read.csv("taxa_divided.csv")
samdata <- read.csv("Metadata.csv")

#coerce matrix

test <- as.matrix(otu_taxa, rownames.force = NA)

#change row labels

otudata$X <- sub(otudata$X, 
                        pattern = "OTU", replacement = "")

install.packages("tidyverse")
library(tidyverse)
otu_data <- otudata %>% remove_rownames %>% column_to_rownames(var="X")


otutax$OTU_ID <- sub(otutax$OTU_ID, 
                  pattern = "OTU", replacement = "")

otu_taxa <- otutax %>% remove_rownames %>% column_to_rownames(var="OTU_ID")

rownames(samdata) <-samdata$InputFileName

# setting data to phyloseq format

library(phyloseq)
OTU <- otu_table(otu_data, taxa_are_rows = TRUE)
TAX <-  tax_table(test)
SAM <- sample_data(samdata)
OTU
TAX
SAM

table <- phyloseq(OTU, TAX)


physeq1 <- merge_phyloseq(table, SAM)
