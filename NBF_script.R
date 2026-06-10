#-------------------------------
# EXPRESSÃO DIFERENCIAL COM VIA METILAÇÃO MUTADOS (NBF)
# Lorenzo Santos Furtado
#-------------------------------

#-------------------------------
# Passo 0: Dependências
#-------------------------------

library(BiocManager)
library(tidyverse)
library(TCGAbiolinks)
library(maftools)
library(SummarizedExperiment)
library(sesameData)
library(sesame)
library(pheatmap)
library(ggpubr)
library(DESeq2)
library(ggrepel)
library(biomaRt)
library(apeglm)
library(org.Hs.eg.db)
library(clusterProfiler)
library(pathview)
library(gage)
library(gageData)


#-------------------------------
# Passo 1: Baixando dados
#-------------------------------


## 1.2: Começamos obtendo dados clínicos dos pacientes do projeto TCGA - LAML
query_clinic <- GDCquery(
  project = "TCGA-LAML", 
  data.category = "Clinical",
  data.type = "Clinical Supplement", 
  data.format = "BCR Biotab"
)

GDCdownload(query_clinic)

laml <- GDCprepare(query_clinic)
laml.clinico <- laml[["clinical_patient_laml"]]

names(laml.clinico)

laml.clinico <- laml.clinico[-c(1:2),]  
laml.clinico <- laml.clinico %>% mutate("Tumor_Sample_Barcode" = bcr_patient_barcode) # Entrada do maf deve ser Tumor_Sample_Barcode

## 1.3: Obtendo dados de mutação dos pacientes do projeto 
laml.maf = system.file('extdata', 'tcga_laml.maf.gz', package = 'maftools')

laml.maf = read.maf(maf = laml.maf,
                clinicalData = laml.clinico,
                verbose = FALSE)

gistic_res_folder <- system.file("extdata", package = "maftools")
laml.gistic = readGistic(gisticDir = gistic_res_folder, isTCGA = TRUE)

## 1.4: Obtendo dados de metilação dos pacientes do projeto 
query_met <- GDCquery(
  project= "TCGA-LAML", 
  data.category = "DNA Methylation", 
  data.type = "Methylation Beta Value",
  platform = "Illumina Human Methylation 450"
)

GDCdownload(query_met,
            files.per.chunk = 20)

laml.met <- GDCprepare(query = query_met, summarizedExperiment = TRUE)

## 1.5: Obtendo dados de expressão dos pacientes do projeto 
query_exp <- GDCquery(
  project = "TCGA-LAML", 
  data.category = "Transcriptome Profiling", 
  data.type = "Gene Expression Quantification", 
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  sample.type = c("Blood Derived Normal","Primary Blood Derived Cancer - Peripheral Blood")
)

GDCdownload(query_exp,
            files.per.chunk = 20,
            method = "api")

laml.exp <- GDCprepare(query_exp, summarizedExperiment = T)

#-------------------------------
# Passo 2: Separando mutação dos genes alvo
#-------------------------------

somaticInteractions(maf = laml.maf, top = 6, pvalue = c(0.05, 0.1)) # verificando co-ocorrencia 
                                                                    #(caso queira escolher vias independentes (co-ocorrentes))
genes <- ("DNMT3A")

maf.data <- laml.maf@data
maf.data <- maf.data[maf.data$Hugo_Symbol %in% genes,] # mantem apenas linhas com os genes escolhidos
maf.data <- maf.data[!duplicated(maf.data$Tumor_Sample_Barcode),] # mantem linhas quem nao (!) tem TRUE 
duplicated(maf.data) # deve retornar FALSE para todos
maf.data <- maf.data %>% pivot_wider(names_from = Hugo_Symbol, values_from = Variant_Type) # transforma valores de hugo para colunas e preenche com valores de variante


levels(maf.data$DNMT3A)

maf.data<- maf.data %>% 
  mutate(DNMT3A = case_when(
    DNMT3A == "SNP" ~ 1,
    DNMT3A == "INS" ~ 1,
    DNMT3A == "DEL" ~ 1,
    is.na(DNMT3A) ~ 0
))

maf.data <- maf.data[,-c(1:10)]
maf.data <- maf.data[,-c(1,4,5)]
names(maf.data)

# Passo 2.1: Juntando a informação de mutação aos dados clinicos
laml.dnmt <- merge(maf.data, laml.clinico, by.x = "Tumor_Sample_Barcode", by.y = "bcr_patient_barcode", all.y = T) #juntando a mutação aos dados clinicos
laml.dnmt<- laml.dnmt %>% # NAs viram zero
  mutate(DNMT3A = case_when(
    DNMT3A == 1 ~ 1,
    is.na(DNMT3A) ~ 0
))

# Passo 2.2: Limpando data frame clinico
colunas_interesse <- c("Tumor_Sample_Barcode", "Protein_Change","DNMT3A","gender",
"race","ethnicity","vital_status","last_contact_days_to",
"death_days_to","fab_category","cyto_risk_group") 

laml.dnmt <- laml.dnmt[,colunas_interesse]

names(laml.dnmt)

laml.dnmt<- laml.dnmt %>%
  mutate(count_days = case_when(
    vital_status == "Dead" ~ death_days_to,
    vital_status == "Alive" ~ last_contact_days_to
))

laml.dnmt<- laml.dnmt %>%
  mutate(vital_status = case_when(
    vital_status == "Dead" ~ 1,
    vital_status == "Alive" ~ 0
))

laml.dnmt$DNMT3A <- as.factor(laml.dnmt$DNMT3A)
laml.dnmt$Tumor_Sample_Barcode <- as.character(laml.dnmt$Tumor_Sample_Barcode)


#-------------------------------
# Passo 3: Metilação Diferencial
#-------------------------------

met.s_metadata <- as.data.frame(colData(laml.met)) # Extraindo os metadados do SummarizedExp
names(met.s_metadata) 

sampless <- met.s_metadata$samples # formato com mais caractere que o tumor_sample_barcode
met.s_metadata <- met.s_metadata %>% mutate(samples = substr(sampless, 1, 12)) # diminuindo os caracteres
laml.dnmt <- laml.dnmt[match(met.s_metadata$samples, laml.dnmt$Tumor_Sample_Barcode), ] # ordenando laml.dnmt$tumorsample pela ordem da matriz de contagem

merged_data <- merge(met.s_metadata, laml.dnmt[, c("Tumor_Sample_Barcode","DNMT3A", "vital_status")], by.x = "samples", by.y = "Tumor_Sample_Barcode", all.x = TRUE)

colData(laml.met) <- DataFrame(merged_data)

laml.met <- laml.met[rowSums(is.na(assay(laml.met))) == 0, ] # de 485,577 para 378,705 elementos

data <- TCGAanalyze_DMC(
  data = laml.met, 
  groupCol = "DNMT3A",
  group1 = "1",
  group2 = "0",
  p.cut = 10^-3,
  diffmean.cut = 0.2,
  legend = "State",
  plot.filename = "METVOLCANODNMT.png"
)

TCGAvisualize_meanMethylation(
  laml.met,
  groupCol = "DNMT3A",
  print.pvalue = T,
  plot.jitter = TRUE,
  jitter.size = 3,
  filename = "groupMeanMet.pdf",
  ylab = expression(paste("Mean DNA methylation (", beta, "-values)")),
  xlab = NULL,
  title = "Mean DNA methylation"
)

met.mean <- data.frame(
  "Sample.mean" = colMeans(assay(laml.met), na.rm = TRUE),
  "groups" = laml.met$DNMT3A
)


met.boxplot <- ggboxplot(
  data = met.mean,
  y = "Sample.mean",
  x = "groups",
  color = "groups",
  add = "jitter",
  ylab = "Mean DNA methylation (beta-values)",
  xlab = ""
) + stat_compare_means() 


# 4. DNA Methylation heatmap -------------------------
status.col <- "status"
probes <- rownames(dmc_s1)[grep("hypo|hyper",dmc_s1$status,ignore.case = TRUE)]
sig.met <- met_s[probes,]


# top annotation, which samples are LGG and GBM
# We will add clinical data as annotation of the samples
# we will sort the clinical data to have the same order of the DNA methylation matrix
clinical.ordered <- leuk.clinic[match(substr(colnames(sig.met),1,12),leuk.clinic$Patient.ID),]

ta <- HeatmapAnnotation(
  df = clinical.ordered[, c("gender", "vital_status", "race")],
  col = list(
    gender = c("male" = "blue", "female" = "pink")
  )
)

# row annotation: add the status for LGG in relation to GBM
# For exmaple: status.gbm.lgg Hypomethyated means that the
# mean DNA methylation of probes for lgg are hypomethylated
# compared to GBM ones.
ra = rowAnnotation(
  df = dmc_s1[probes, status.col],
  col = list(
    "status" =
      c(
        "Hypomethylated in Dead" = "orange",
        "Hypermethylated in Dead" = "darkgreen"
      )
  ),
  width = unit(1, "cm")
)

heatmap  <- Heatmap(
  matrix = assay(sig.met),
  name = "DNA methylation",
  col = matlab::jet.colors(200),
  show_row_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  bottom_annotation = ta,
  column_title = "DNA Methylation"
) 
# Save to pdf
png("heatmap.png",width = 600, height = 400)
draw(heatmap, annotation_legend_side =  "bottom")
dev.off()
#-------------------------------
# Passo 4: Expressão Diferencial
#-------------------------------

laml.exp_matrix <- assay(laml.exp,"unstranded") 
laml.exp_df <- as.data.frame(laml.exp_matrix)
listazerodeseq <- rowSums(laml.exp_df) == 0
laml.exp_df <- laml.exp_df[!listazerodeseq, ]

exp.ordem_metadado <- laml.exp@colData@listData[["cases"]]
exp.ordem_metadado <- substr(exp.ordem_metadado, 1, 12) # diminuindo os caracteres
colnames(laml.exp_df) <- exp.ordem_metadado


laml.dnmt <- laml.dnmt[match(exp.ordem_metadado, laml.dnmt$Tumor_Sample_Barcode), ] # ordenando laml.dnmt$tumorsample pela ordem da matriz de contagem
identical(laml.dnmt$Tumor_Sample_Barcode, exp.ordem_metadado) # checando a ordem
laml.dnmt <- laml.dnmt %>% column_to_rownames("Tumor_Sample_Barcode")
colData(laml.exp) <- DataFrame(laml.dnmt)


all(colnames(laml.exp_df) %in% rownames(colData(laml.exp))) # esperado TRUE (coluna de counts alinhadas com linha de metadados)


barplot(colSums(laml.exp_df)/1e7, las =3)

# Criando objeto DESeq2
dds_laml <- DESeqDataSetFromMatrix(countData =  laml.exp_df,    
                                   colData =  colData(laml.exp),
                                   design =~ DNMT3A
) 

dds_laml <- DESeq(dds_laml) # rodando a expressão diferencial 

rld = vst(dds_laml)
plotPCA(rld, intgroup = "DNMT3A")

resultsNames(dds_laml) # extraindo resultados (string) da analise para fazer o contraste

res <- results(dds_laml, contrast=c("DNMT3A","1","0")) # ordena pelo controle (com mutação)

laml_res_all <- as.data.frame(res)
laml_res_all <- laml_res_all %>% arrange(laml_res_all$padj) # ordenando os dados pelo menor pajustado
sum(laml_res_all$padj < 0.01, na.rm=TRUE) # quantos genes estão dif. exp para pajus <0.05?

with(laml_res_all, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", xlim=c(-3,3)))
with(subset(laml_res_all, padj<.01 ), points(log2FoldChange, -log10(padj), pch=20, col="blue"))
with(subset(laml_res_all, padj<.01 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(padj), pch=20, col="red"))


# Estringir mais 
laml_res_all <- laml_res_all %>%
  filter(padj < 0.05) # out = 1186 genes

rownames(laml_res_all) <- sub("\\..*", "", rownames(laml_res_all)) 
  
enslista <- rownames(laml_res_all)
enslista <- sub("\\..*", "", enslista) 


ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl") # adicionando 

resultssENS <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                     filters = "ensembl_gene_id",
                     values = enslista,
                     mart = ensembl)

sum(resultssENS$hgnc_symbol == "") # 188

faltantes <- resultssENS$ensembl_gene_id[
  resultssENS$hgnc_symbol == ""]

faltantes

laml_res_all$genelabels <- resultssENS$hgnc_symbol[match(rownames(laml_res_all), resultssENS$ensembl_gene_id)]

vetor <- laml_res_all$genelabels[!is.na(laml_res_all$genelabels) & laml_res_all$genelabels != ""]
vetor

ggplot(laml_res_all, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(alpha = 0.4) +
  geom_point(data = subset(laml_res_all, padj < 0.05 & abs(log2FoldChange) > 2),
             color = "red") +
  geom_text_repel(
    data = subset(laml_res_all, padj < 0.01 & abs(log2FoldChange) > 1),
    aes(label = genelabels),
    size = 3
  ) +
  xlim(c(-3, 3)) +
  theme_minimal() +
  labs(title = "Volcano plot")

#-------------------------------
# Passo 5: Enriquecimento Funcional
#-------------------------------

laml_res_all <- laml_res_all %>%
  drop_na()

laml_res_all <- tibble::rownames_to_column(laml_res_all, "ENSEMBL")


gene_conversion <- bitr(laml_res_all$ENSEMBL, fromType = "ENSEMBL",
                        toType = "ENTREZID", OrgDb = org.Hs.eg.db) # out: 18.6% of input gene IDs are fail to map...

gene_list <- unique(na.omit(gene_conversion$ENTREZID))

#-------------------------------
# GO Analysis
#-------------------------------

## Criando a lista para o argumento universo 

# Obtenha todos os termos GO e seus genes anotados (Entrez IDs)

# Aqui, "BP" significa Biological Process, que é onde GO:0030097 está.
go_map_genes <- AnnotationDbi::as.list(org.Hs.eg.db::org.Hs.egGO2EG)

# Definindo o Termo GO que você quer buscar (0030097: Hematopoietic cell lineage)
termo_go_alvo <- "GO:0030097"

# Vetor de Entrez IDs para o termo alvo
entrez_ids_hemopoiesis <- go_map_genes[[termo_go_alvo]]

# Remove possíveis NAs ou IDs duplicados 
entrez_ids_hemopoiesis <- unique(na.omit(entrez_ids_hemopoiesis))

input.hemato <- bitr(entrez_ids_hemopoiesis, fromType = "ENTREZID",
                                 toType = "ENSEMBL", OrgDb = org.Hs.eg.db) # 1.64% of input gene IDs are fail to map...


# Molecular Function (MF), Cellular Component (CC), and Biological Process (BP)
go_BP <- enrichGO(gene = gene_conversion$ENSEMBL,
                OrgDb = "org.Hs.eg.db",
                ont = "BP",
                universe = input.hemato$ENSEMBL,
                pAdjustMethod = "BH",
                keyType = 'ENSEMBL')

barplot(go_BP, showCategory = 20, font.size = 5)  
go_BP.df <- as.data.frame(go_BP@result)


go_MF <- enrichGO(gene = gene_conversion$ENSEMBL,
                  OrgDb = "org.Hs.eg.db",
                  ont = "MF",
                  universe = input.hemato$ENSEMBL,
                  pAdjustMethod = "BH",
                  keyType = 'ENSEMBL')

barplot(go_MF, showCategory = 20, font.size = 5)  
go_MF.df <- as.data.frame(go_MF@result)


## Analise sem definir argumento universo 

# Molecular Function (MF), Cellular Component (CC), and Biological Process (BP)
go_BP <- enrichGO(gene = gene_conversion$ENSEMBL,
                  OrgDb = "org.Hs.eg.db",
                  ont = "BP",
                  pAdjustMethod = "BH",
                  keyType = 'ENSEMBL')

go_BP <- setReadable(go_BP, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL") # transformando ENS p HUGO
go_BP.df <- as.data.frame(go_BP@result)
go_BP.df <- go_BP.df[go_BP.df$p.adjust <= 0.01,]

cnetplot(go_BP)
heatplot(go_BP)
go_bp2 <- pairwise_termsim(go_BP)


barplot(go_BP2, showCategory = 5, font.size = 5)  


go_MF <- enrichGO(gene = gene_conversion$ENSEMBL,
                  OrgDb = "org.Hs.eg.db",
                  ont = "MF",
                  pAdjustMethod = "BH",
                  keyType = 'ENSEMBL')

go_MF <- setReadable(go_MF, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL") # transformando ENS p HUGO

barplot(go_MF, showCategory = 20, font.size = 5)  
go_MF.df <- as.data.frame(go_MF@result)
go_MF.df <- go_MF.df[go_MF.df$p.adjust <= 0.01,]
write_csv(go_MF.df, file ='go_MF.df.csv')


#-------------------------------
# KEGG Analysis
#-------------------------------

foldchanges <- res$log2FoldChange
names(foldchanges) <- res$ENTREZID

kegg_results <- enrichKEGG(gene = gene_list, 
                           organism = "hsa",  
                           pvalueCutoff = 0.05,
)

# KEGG hematopoiesis (hsa04640: Hematopoietic cell lineage - Homo sapiens (human))

# Feito manuealmente a partir de: https://www.kegg.jp/entry/hsa04640 > Gene)

ids_genes_hsa04640 <- c(
  "100133941", "102723407", "1378", "1379", "1380", "1435", "1436", "1437", 
  "1438", "1440", "1441", "1604", "1791", "2056", "2057", "2208", "2209", 
  "2322", "2323", "2811", "2812", "2814", "2815", "290", "2993", "3108", 
  "3109", "3111", "3112", "3113", "3115", "3117", "3118", "3119", "3120", 
  "3122", "3123", "3125", "3126", "3127", "3552", "3553", "3554", "3559", 
  "3562", "3563", "3565", "3566", "3567", "3568", "3569", "3570", "3574", 
  "3575", "3581", "3589", "3590", "3655", "3672", "3673", "3674", "3675", 
  "3676", "3678", "3684", "3690", "3815", "4254", "4311", "7037", "7066", 
  "7124", "7850", "909", "910", "911", "912", "913", "914", "915", "916", 
  "917", "920", "921", "924", "925", "926", "927", "928", "929", "930", 
  "931", "933", "945", "947", "948", "951", "952", "960", "966"
)

# transformar lista com os IDS e adicionar 
gene_list_named <- laml_res_all %>%
  inner_join(gene_conversion, by = c("ENSEMBL" = "ENSEMBL")) %>%  
  dplyr::select(ENTREZID, log2FoldChange) %>% 
  distinct(ENTREZID, .keep_all = TRUE) %>% 
  tibble::deframe()

# A partir do head plotei os genes listados

#cell growth and death (24)
hsa04640 <- pathview(gene.data  = gene_list_named,
                     pathway.id = "hsa04640",
                     species    = "hsa",
                     limit      = list(gene=max(abs(gene_list_named)), cpd=1))








# visualisando dps de fazer o encolhimento por ln
# resLFC <- lfcShrink(dds_laml, coef="DNMT3A_1_vs_0", type="apeglm")
# resLFC
# resLFC <- as_data_frame(resLFC)
# resLFC <- resLFC %>% arrange(resLFC$padj) # ordenando os dados pelo menor pajustado
# 
# with(resLFC, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", xlim=c(-3,3)))
# with(subset(resLFC, padj<.01 ), points(log2FoldChange, -log10(padj), pch=20, col="blue"))
# with(subset(resLFC, padj<.01 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(padj), pch=20, col="red"))
# ficou bizarro

# Referências

#https://bioconductor.org/packages/release/bioc/html/TCGAbiolinks.html pipelines
#https://bioconductor.org/packages/devel/bioc/vignettes/maftools/inst/doc/maftools.html pipeline maftools
#https://www.bioconductor.org/packages/devel/bioc/manuals/TCGAbiolinks/man/TCGAbiolinks.pdf documentacao tcga biolinks
#https://bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/download_prepare.html download para diferentes dados
#https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html