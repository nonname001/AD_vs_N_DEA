# install packages
# install.packages("tidyverse")
library(tidyverse)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("limma")
library(limma)

matrixx <- readRDS("pgdata_matrix_2.rds") # AD vs N dataset only contains AD data??
metadata <- readRDS("pgdata_metadata.rds")
doe <- read_tsv("manifest2.tsv")

Gender <- factor(make.names(doe$Gender)) # dummy variable

design <- model.matrix(~0+Gender)
contrast <- makeContrasts(GenderComparison = GenderFemale - GenderMale, levels = design)

fit <- lmFit(matrixx, design) # what do these do?
fit <- eBayes(fit)
fit2 <- contrasts.fit(fit, contrast) # what is the difference between fit and fit2?
fit2 <- eBayes (fit2)

top <- topTable(fit, coef = 1, adjust.method = "BH",n=Inf, sort.by = "P")
top2 <- topTable(fit2, coef = 1, adjust.method = "BH",n=Inf, sort.by = "P")
top <- merge(top, metadata, by=0, all=TRUE)
top2 <- merge(top2, metadata, by=0, all=TRUE)

top2_sig <- top2 %>%
  filter(P.Value <= 0.05)

print (summary(decideTests(fit,adjust.method = "BH", p.value = 0.05)))
print (summary(decideTests(fit,adjust.method = "none", p.value = 0.05)))
print (summary(decideTests(fit2,adjust.method = "BH", p.value = 0.05)))
print (summary(decideTests(fit2,adjust.method = "none", p.value = 0.05)))