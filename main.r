# install packages
# install.packages("tidyverse")
library(tidyverse)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("limma")
library(limma)
library(plotly)
library(htmlwidgets)
source("GSEAenricher.R")

top_tables <- list()
pea_tables <- list()
bubble_plots <- list()

run_dea_pea <- function(matrix_file, metadata_file, doe_file, index, save_path_dea, save_path_pea) {
  matrixx <- readRDS(matrix_file) 
  metadata <- readRDS(metadata_file)
  doe <- read_tsv(doe_file)
  
  Gender <- factor(make.names(doe$Gender))
  Age <- doe$AGE
  Braak <- factor(make.names(doe$Braak.stage))
  Batch <- factor(make.names(doe$Gel.Batch))
  
  design <- model.matrix(~0+Gender+Age+Braak+Batch)
  contrast <- makeContrasts(GenderComparison = GenderFemale - GenderMale, levels = design)
  
  fit <- lmFit(matrixx, design) 
  # fit <- eBayes(fit)
  fit2 <- contrasts.fit(fit, contrast) 
  fit2 <- eBayes (fit2)
  
  # top <- topTable(fit, coef = "GenderComparison", adjust.method = "BH",n=Inf, sort.by = "P")
  top2 <- topTable(fit2, coef = "GenderComparison", adjust.method = "BH",n=Inf, sort.by = "P")
  # top <- merge(top, metadata, by=0, all=TRUE)
  top2 <- merge(top2, metadata, by=0, all=TRUE)
  top2$Shown.ID <- top2$hgnc_symbol
  top2$Shown.ID[top2$Shown.ID==""] <- top2$Protein.IDs[top2$Shown.ID==""]
  
  top_tables[[index]] <<- top2
  
  top2_sig <- top2 %>%
     filter(P.Value <= 0.05)
  
  print (summary(decideTests(fit,adjust.method = "BH", p.value = 0.05)))
  print (summary(decideTests(fit,adjust.method = "none", p.value = 0.05)))
  print (summary(decideTests(fit2,adjust.method = "BH", p.value = 0.05)))
  print (summary(decideTests(fit2,adjust.method = "none", p.value = 0.05)))
  
  # write.table(top2_sig, file='sig.tsv', quote=FALSE, sep='\t', col.names = NA)
  
  load("gmt-go.Rdata")
  
  g1 <- do_GSEA(top2, 1, gmt.go, group_name = "F_vs_M", Protein.ID.Column ="Protein.IDs")
  
  g1_filtered <- g1$df %>%
    filter(pvalue < 0.05)
  
  pea_tables[[index]] <<- g1
  
  hcal_dea <- max ((length (unique (top2$Protein.IDs)) * 16.5 + 25 + 10 + 100), 500)
  hcal_pea <- max ((length (unique (g1$df$ID)) * 16.5 + 25 + 10 + 100), 500)
  #gp.pt <- ggplotly(p, tooltip = "text", height = hcal, source="DEAPlotSource") %>% 
  #  layout (yaxis=(list(automargin = F)), margin=list (l=200))

  # ggplot(top2_sig, aes(y=Protein.IDs, x=logFC)) + 
  #   geom_point(aes(size=P.Value, color=AveExpr)) + 
  #   labs(y="Protein IDs", x="logFC", color="Average expression", size="P-value") +
  #   theme(axis.text.x=element_text(angle=45, hjust=1), axis.text.y=element_text(size = 8, angle = 0, hjust = 1, face = "plain"))
  
  p1 <- ggplot (data=top2_sig, aes (x=logFC, y=Protein.IDs)) +
    geom_point(aes (size=AveExpr, fill=logFC, color=-log10(P.Value), stroke=0.5, text=paste('<b>Log fold change:</b>', logFC, '<br>',
                                                                                         '<b>Average expression:</b>', AveExpr,'<br>',
                                                                                         '<b>t-statistic:</b>', t,'<br>',
                                                                                         '<b>P-value:</b>',P.Value,'<br>',
                                                                                         '<b>Negative log10-P-value:</b>',-log10(P.Value),'<br>',
                                                                                         '<b>B-statistic:</b>',B,'<br>'))) +
    scale_fill_gradient2(limits=c(-max(abs(top2_sig$logFC)), max(abs(top2_sig$logFC))), low = "#0199CC", mid = "white", high = "#FF3705", midpoint = 0) +
    scale_color_gradient2(limits=c(0, max(abs(log10(top2_sig$P.Value)))), low = "white", mid = "white", high = "black", midpoint = 1) +
    theme_minimal () + xlab("Comparison") + ylab("") + 
    theme(axis.text.x=element_text(angle=45, hjust=1), axis.text.y=element_text(size = 8, angle = 0, hjust = 1, face = "plain"))
  
  gp.pt.1 <- ggplotly(p1, tooltip = "text", height = hcal_dea) %>% 
    layout (yaxis=(list(automargin = F)), margin=list (l=200))
  
  saveWidget(gp.pt.1, save_path_dea)
  
  # ggsave(save_path_dea)
  
  # ggplot(g1_filtered, aes(y=ID, x=rank)) + 
  #   geom_point(aes(size=pvalue, color=enrichmentScore)) + 
  #   labs(y="Protein IDs", x="logFC", color="Enrichment score", size="P-value") +
  #   theme(axis.text.x=element_text(angle=45, hjust=1), axis.text.y=element_text(size = 8, angle = 0, hjust = 1, face = "plain"))
  
  p2 <- ggplot (data=g1_filtered, aes (x=enrichmentScore, y=ID)) +
    geom_point(aes (size=setSize, fill=NES, color=-log10(pvalue), stroke=0.5, text=paste('<b>Core_enrichment:</b>', core_enrichment, '<br>',
                                                                                         '<b>Comparison:</b>', Group,'<br>',
                                                                                         '<b>Setsize:</b>', setSize,'<br>',
                                                                                         '<b>NES:</b>',NES,'<br>',
                                                                                         '<b>pvalue:</b>',pvalue,'<br>',
                                                                                         '<b>rank:</b>',rank,'<br>',
                                                                                         '<b>leading_edge:</b>',leading_edge,'<br>'))) +
    scale_fill_gradient2(limits=c(-max(abs(g1_filtered$NES)), max(abs(g1_filtered$NES))), low = "#0199CC", mid = "white", high = "#FF3705", midpoint = 0) +
    scale_color_gradient2(limits=c(0, max(abs(log10(g1_filtered$pvalue)))), low = "white", mid = "white", high = "black", midpoint = 1) +
    theme_minimal () + xlab("Comparison") + ylab("") + 
    theme(axis.text.x=element_text(angle=45, hjust=1), axis.text.y=element_text(size = 8, angle = 0, hjust = 1, face = "plain"))
  
  gp.pt.2 <- ggplotly(p2, tooltip = "text", height = hcal_pea) %>% 
    layout (yaxis=(list(automargin = F)), margin=list (l=200))
  
  saveWidget(gp.pt.2, save_path_pea)
  
  # ggsave(save_path_pea)
}

run_dea_pea("pgdata_matrix_1.rds", "pgdata_metadata_1.rds", "manifest2.tsv", "filter_mindet", "dea_1.html", "pea_1.html")
run_dea_pea("pgdata_matrix_2.rds", "pgdata_metadata_2.rds", "manifest2.tsv", "filter_bpca", "dea_2.html", "pea_2.html")
run_dea_pea("pgdata_matrix_3.rds", "pgdata_metadata_3.rds", "manifest2.tsv", "nofilter_mindet", "dea_3.html", "pea_3.html")
run_dea_pea("pgdata_matrix_4.rds", "pgdata_metadata_4.rds", "manifest2.tsv", "nofilter_bpca", "dea_4.html", "pea_4.html")