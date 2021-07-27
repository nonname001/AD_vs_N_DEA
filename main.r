# install packages
# install.packages("tidyverse")
library(tidyverse)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("limma")
library(limma)
source("GSEAenricher.R")

matrixx <- readRDS("pgdata_matrix_2.rds") # AD vs N dataset only contains AD data??
metadata <- readRDS("pgdata_metadata_2.rds")
doe <- read_tsv("manifest2.tsv")

Gender <- factor(make.names(doe$Gender))
Age <- doe$AGE
Braak <- factor(make.names(doe$Braak.stage))
Batch <- factor(make.names(doe$Gel.Batch))

design <- model.matrix(~0+Gender+Age+Braak+Batch)
contrast <- makeContrasts(GenderComparison = GenderFemale - GenderMale, levels = design)

fit <- lmFit(matrixx, design) # what do these do?
# fit <- eBayes(fit)
fit2 <- contrasts.fit(fit, contrast) # what is the difference between fit and fit2?
fit2 <- eBayes (fit2)

# top <- topTable(fit, coef = "GenderComparison", adjust.method = "BH",n=Inf, sort.by = "P")
top2 <- topTable(fit2, coef = "GenderComparison", adjust.method = "BH",n=Inf, sort.by = "P")
# top <- merge(top, metadata, by=0, all=TRUE)
top2 <- merge(top2, metadata, by=0, all=TRUE)

top2_sig <- top2 %>%
  filter(P.Value <= 0.05)

print (summary(decideTests(fit,adjust.method = "BH", p.value = 0.05)))
print (summary(decideTests(fit,adjust.method = "none", p.value = 0.05)))
print (summary(decideTests(fit2,adjust.method = "BH", p.value = 0.05)))
print (summary(decideTests(fit2,adjust.method = "none", p.value = 0.05)))

write.table(top2_sig, file='sig.tsv', quote=FALSE, sep='\t', col.names = NA)

load("gmt-go.Rdata")

g1 <- do_GSEA(top2, 0.05, gmt.go, group_name = "F_vs_M", Protein.ID.Column ="Protein.IDs")
