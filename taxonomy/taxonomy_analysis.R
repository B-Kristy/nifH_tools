setwd("C:/KBS_Code/taxonomy")
library(readxl)
library(phyloseq)
library(tidyverse)
library(microbiome)

metadata <- read_excel("nifH_soil_root_metadata.xlsx")
metadata <- metadata %>% remove_rownames %>% column_to_rownames(var="Sample_ID")

metadata <- sample_data(metadata)

otutable <- read.delim("otu_table.txt", header = TRUE, sep="\t")

colnames(otutable) <-gsub("_R1_001\\.fastq$", "", colnames(otutable))
colnames(otutable)[1] <- "Feature.ID"

otutable <- otutable %>% remove_rownames %>% column_to_rownames(var="Feature.ID")

otutable <- as.matrix(otutable)
otutable <- otu_table(otutable, taxa_are_rows = TRUE)

##NFIXDB##
nfix_tax <- read.delim("nfixdb_taxonomy.tsv", header = TRUE, sep="\t")
#remove first column
nfix_tax <- nfix_tax [-1,]
#get rid of sample label
nfix_tax$Feature.ID <- gsub(";[^;]*;", "", nfix_tax$Feature.ID)
nfix_tax <- nfix_tax %>% remove_rownames %>% column_to_rownames(var="Feature.ID")
nfix_tax_sep <- nfix_tax %>% 
  mutate(Taxon = str_split(Taxon, ";")) %>%
  unnest_wider(Taxon, names_sep = "_") %>%
  mutate(across(everything(), ~ ifelse(is.na(.), "", .))) %>%
  transmute(
    Domain = str_remove(Taxon_1, "d__"),
    Phylum = str_remove(Taxon_2, "p__"),
    Class = str_remove(Taxon_3, "c__"),
    Order = str_remove(Taxon_4, "o__"),
    Family = str_remove(Taxon_5, "f__"),
    Genus = str_remove(Taxon_6, "g__"),
    Species = str_remove(Taxon_7, "s__")
  )

fill_unclassified <- function(Domain, Phylum, Class, Order, Family, Genus, Species) {
  tax_vec <- c(Domain, Phylum, Class, Order, Family, Genus, Species)
  levels <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  last_known <- NA
  for (i in seq_along(tax_vec)) {
    if (!is.na(tax_vec[i]) && tax_vec[i] != "") {
      last_known <- tax_vec[i]
    } else {
      tax_vec[i] <- paste("Unclassified", last_known)
    }
  }
  # Return as a named list to ensure it's treated as a row
  set_names(as.list(tax_vec), levels)
}
nfix_tax_filled <- nfix_tax_sep %>% 
  mutate(across(everything(), as.character)) %>%
  pmap_dfr(fill_unclassified)

nfixdb_tax_final <- cbind(nfix_tax, nfix_tax_filled)

nfixdb_tax_final <- as.matrix(nfixdb_tax_final)
nfixdb_tax_final <- tax_table(nfixdb_tax_final)

nfixdb_all_phyloseq <- merge_phyloseq(otutable, metadata, nfixdb_tax_final)

###ZEHR###
zehr_tax <- read.delim("zehr_taxonomy.tsv", header = TRUE, sep="\t")
#remove first column
zehr_tax <- zehr_tax [-1,]
#get rid of sample label
zehr_tax$Feature.ID <- gsub(";[^;]*;", "", zehr_tax$Feature.ID)
zehr_tax <- zehr_tax %>% remove_rownames %>% column_to_rownames(var="Feature.ID")
zehr_tax_sep <- zehr_tax %>% 
  mutate(Taxon = str_split(Taxon, ";")) %>%
  unnest_wider(Taxon, names_sep = "_") %>%
  mutate(across(everything(), ~ ifelse(is.na(.), "", .))) %>%
  transmute(
    Domain = str_remove(Taxon_1, "d__"),
    Phylum = str_remove(Taxon_2, "p__"),
    Class = str_remove(Taxon_3, "c__"),
    Order = str_remove(Taxon_4, "o__"),
    Family = str_remove(Taxon_5, "f__"),
    Genus = str_remove(Taxon_6, "g__"),
    Species = str_remove(Taxon_7, "s__")
  )

zehr_tax_filled <- zehr_tax_sep %>% 
  mutate(across(everything(), as.character)) %>%
  pmap_dfr(fill_unclassified)

zehr_tax_final <- cbind(zehr_tax, zehr_tax_filled)

zehr_tax_final <- as.matrix(zehr_tax_final)
zehr_tax_final <- tax_table(zehr_tax_final)

zehr_all_phyloseq <- merge_phyloseq(otutable, metadata, zehr_tax_final)

##BUCKLEY##
buckley_tax <- read.delim("buckley_taxonomy.tsv", header = TRUE, sep="\t")
buckley_tax <- buckley_tax [-1,]
buckley_tax$Feature.ID <- gsub(";[^;]*;", "", buckley_tax$Feature.ID)
buckley_tax <- buckley_tax %>% remove_rownames %>% column_to_rownames(var="Feature.ID")
buckley_tax_sep <- buckley_tax %>% 
  mutate(Taxon = str_split(Taxon, ";")) %>%
  unnest_wider(Taxon, names_sep = "_") %>%
  mutate(across(everything(), ~ ifelse(is.na(.), "", .))) %>%
  transmute(
    Domain = str_remove(Taxon_1, "d__"),
    Kingdom = str_remove(Taxon_2, "k__"),
    Phylum = str_remove(Taxon_3, "p__"),
    Class = str_remove(Taxon_4, "c__"),
    Order = str_remove(Taxon_5, "o__"),
    Family = str_remove(Taxon_6, "f__"),
    Genus = str_remove(Taxon_7, "g__"),
    Species = str_remove(Taxon_8, "s__")
  )
fill_unclassified_buck <- function(Domain, Kingdom, Phylum, Class, Order, Family, Genus, Species) {
  tax_vec <- c(Domain, Kingdom, Phylum, Class, Order, Family, Genus, Species)
  levels <- c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  last_known <- NA
  for (i in seq_along(tax_vec)) {
    if (!is.na(tax_vec[i]) && tax_vec[i] != "") {
      last_known <- tax_vec[i]
    } else {
      tax_vec[i] <- paste("Unclassified", last_known)
    }
  }
  # Return as a named list to ensure it's treated as a row
  set_names(as.list(tax_vec), levels)
}

buckley_tax_filled <- buckley_tax_sep %>% 
  mutate(across(everything(), as.character)) %>%
  pmap_dfr(fill_unclassified_buck)

buckley_tax_filled_no_king <- buckley_tax_filled %>% dplyr::select(-c(Kingdom))

buckley_tax_final <- cbind(buckley_tax, buckley_tax_filled_no_king)

buckley_tax_final <- as.matrix(buckley_tax_final)
buckley_tax_final <- tax_table(buckley_tax_final)

buckley_all_phyloseq <- merge_phyloseq(otutable, metadata, buckley_tax_final)

#MAKE MERGED#
nfixdb_tax_final <- as.data.frame(nfixdb_tax_final)
nfixdb_tax_final$nfixdbtaxon <- paste(nfixdb_tax_final$Domain,
                                      nfixdb_tax_final$Phylum,
                                      nfixdb_tax_final$Class,
                                      nfixdb_tax_final$Order,
                                      nfixdb_tax_final$Family,
                                      nfixdb_tax_final$Genus,
                                      nfixdb_tax_final$Species,
                                      sep=";")
nfixdb_tax_final <- nfixdb_tax_final %>%
  rownames_to_column(var="OTU_ID")
##zehr##
zehr_tax_final <- as.data.frame(zehr_tax_final)
zehr_tax_final$zehrtaxon <- paste(zehr_tax_final$Domain,
                                  zehr_tax_final$Phylum,
                                  zehr_tax_final$Class,
                                  zehr_tax_final$Order,
                                  zehr_tax_final$Family,
                                  zehr_tax_final$Genus,
                                  zehr_tax_final$Species,
                                      sep=";")
zehr_tax_final <- zehr_tax_final %>%
  rownames_to_column(var="OTU_ID")
#buckley#
buckley_tax_final <- as.data.frame(buckley_tax_final)
buckley_tax_final$buckleytaxon <- paste(buckley_tax_final$Domain,
                                        buckley_tax_final$Phylum,
                                        buckley_tax_final$Class,
                                        buckley_tax_final$Order,
                                        buckley_tax_final$Family,
                                        buckley_tax_final$Genus,
                                        buckley_tax_final$Species,
                                      sep=";")
buckley_tax_final <- buckley_tax_final %>%
  rownames_to_column(var="OTU_ID")

mini_merge <- merge(nfixdb_tax_final,zehr_tax_final, by="OTU_ID")
merged_tax <- merge(mini_merge, buckley_tax_final, by="OTU_ID")
merged_tax <- merged_tax %>% dplyr::select(c(OTU_ID, nfixdbtaxon,zehrtaxon,buckleytaxon))
