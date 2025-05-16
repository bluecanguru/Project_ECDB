# ╔════════════════════════════════════════════════════════════════╗
# ║                        TRABALHO DE EXTRAÇÃO                    ║
# ║                      Brain Lower Grade Glioma                  ║
# ╠════════════════════════════════════════════════════════════════╣
# ║ Autores:                                                       ║
# ║   • Cátia Rosário (pg57791)                                    ║
# ║   • Vanessa Rodriguez (pg49131)                                ║
# ║   • André Dias (55127)                                         ║
# ║                                                                ║
# ║ Data: 10/04/2025                                               ║
# ╚════════════════════════════════════════════════════════════════╝

### 0. Importação das packages e dos dados ----

# Grouped Packages by Purpose

cran_packages <- c(
  "ggplot2",      # Visualization
  "gplots",       # Enhanced plots
  "gridExtra",    # Combine multiple plots
  "DT",           # Interactive datatables
  "dplyr",        # Data manipulation
  "plyr",         # Data manipulation (older)
  "stringr",      # String operations
  "knitr",        # RMarkdown rendering
  "data.table",   # Efficient data handling
  "viridis",      # Color palettes
  "openxlsx",     # Excel export
  "corrplot",     # Correlation plots
  "mgcv",         # GAM models
  "ggpubr",       # Publication-ready plots
  "GGally",       # Pair plots (extension of ggplot2)
  "htmltools",    # To print text in specific format
  "caret",        # Modelation
  "randomForest", # Random Forest
  "smotefamily",  # SMOTE
  "xgboost",      # xgboost
  "e1071",        # SVM
  "pROC",         # ROC CURVE
  "MLmetrics",    # F1_score
  "pander",       # Tables
  "tibble",       # Interactive datatables
  "Rtsne",        # UMAP
  "uwot",         # T-SNE
  "factoextra"    # T-SNE
)

bioc_packages <- c(
  "edgeR",             # RNA-seq analysis
  "limma",             # Linear models for microarray/RNA-seq
  "clusterProfiler",   # Functional enrichment
  "org.Hs.eg.db",      # Human gene annotation
  "AnnotationDbi",     # Annotation infrastructure
  "biomaRt",           # Interface to BioMart databases
  "enrichplot"         # Visualization of enrichment results
)

other_packages <- c(
  "gprofiler2"         # g:Profiler interface (CRAN but bio-focused)
)


# Installation Block

# Install BiocManager if missing
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install missing CRAN packages
cran_missing <- cran_packages[!cran_packages %in% rownames(installed.packages())]
if (length(cran_missing) > 0) install.packages(cran_missing)

# Install missing Bioconductor packages
bioc_missing <- bioc_packages[!bioc_packages %in% rownames(installed.packages())]
if (length(bioc_missing) > 0) BiocManager::install(bioc_missing)

# Install other packages (gprofiler2 is from CRAN)
other_missing <- other_packages[!other_packages %in% rownames(installed.packages())]
if (length(other_missing) > 0) install.packages(other_missing)

# Load All Packages
all_packages <- c(cran_packages, bioc_packages, other_packages)
invisible(lapply(all_packages, function(pkg) {
  suppressMessages(library(pkg, character.only = TRUE))
}))

# Import data
clini_data <- read.delim("lgg_tcga_pan_can_atlas_2018_clinical_data.tsv", stringsAsFactors=TRUE)
gene_data <- read.delim("all_genes_mrna.txt")

palette_3 <- c("#fd7f6f", "#7eb0d5", "#b2e061", "#bd7ebe", "#ffb55a", "#ffee65", "#beb9db", "#fdcce5", "#8bd3c7")

# 1. Tratamento e Vizualização de dados ----
# 1.a. Preparação e pré-processamento dos dados ----
# 1.1. Dados Clínicos ----

## 2.1. Estrutura dos dados ----

### Preview dos dados
datatable(head(clini_data), options = list(scrollX = TRUE))

### Estrutura e formato dos dados
str(clini_data)

## 2.2. Manipulação de dados ----

### Change column names
nomes <- names(clini_data)
nomes <- gsub("[\\.]", " ", nomes)
nomes <- gsub("  ", " ", nomes)
nomes <- gsub(" $", "", nomes, perl=T)
names(clini_data) <- nomes

### Filter by columns of interest
novas_colunas<-c("Sample ID","Diagnosis Age","Cancer Type Detailed","Months of disease specific survival","Disease specific Survival status","Fraction Genome Altered","Genetic Ancestry Label","Sex","Tumor Break Load")
clini_data <- clini_data[,novas_colunas]

### Change variable names
colnames(clini_data) <- c("sample_ID","diagnosis_age","cancer_type","surv_months","surv_status","genome_alt","ancestry","sex","tbl")

### Convert sample ID to same format as expression data
sample_id <- clini_data$sample_ID
sample_id <- gsub("-", "\\.", sample_id)
clini_data$sample_ID <- sample_id

### Convert survival months to years
clini_data$surv_months <- clini_data$surv_months/12
colnames(clini_data) <- c("sample_ID","diagnosis_age","cancer_type","surv_years","surv_status","genome_alt","ancestry","sex","tbl")

## 2.3. Limpeza de Dados - Remoção de NAs ----

### Remoção de NAs
clini_data_no_na <- na.omit(clini_data)

### Proporção de NAs no dataset
(nrow(clini_data) - nrow(clini_data_no_na))/nrow(clini_data)

## 2.4. Sumarização e Visualização de Dados ----

### Medidas de localização e dispersão
summary(clini_data_no_na)

# 1.a. Vizualização de dados ----

# 1.1 Visualização univariada ----

# 1.1.1 Variáveis numéricas ----
a <- ggplot(alpha=0.8) + geom_histogram(data=clini_data_no_na, aes(x=diagnosis_age), fill="#7eb0d5") + theme_bw()
b <- ggplot(alpha=0.8) + geom_histogram(data=clini_data_no_na, aes(x=surv_years), fill="#7eb0d5") + theme_bw()
ggarrange(a, b)

a <- ggplot(alpha=0.8) + geom_histogram(data=clini_data_no_na, aes(x=genome_alt), fill="#7eb0d5") + theme_bw()
b <- ggplot(alpha=0.8) + geom_histogram(data=clini_data_no_na, aes(x=tbl), fill="#7eb0d5") + theme_bw()
grid.arrange(a, b, ncol=2)

a <- ggplot(alpha=0.8) + geom_bar(data=clini_data_no_na, aes(x=sex, fill=sex)) + scale_fill_manual(values = palette_3) + theme_bw()
b <- ggplot(alpha=0.8) + geom_bar(data=clini_data_no_na, aes(x=cancer_type, fill=cancer_type)) + scale_fill_manual(values = palette_3) + theme_bw()
grid.arrange(a, b, ncol=2)

a <- ggplot(alpha=0.8) + geom_bar(data=clini_data_no_na, aes(x=ancestry, fill=ancestry)) + scale_fill_manual(values = palette_3) + theme_bw()
b <- ggplot(alpha=0.8) + geom_bar(data=clini_data_no_na, aes(x=surv_status, fill=surv_status)) + scale_fill_manual(values = palette_3) + theme_bw()
grid.arrange(a, b, ncol=2)

# 1.2 Visualização multivariada ----
## Variável numérica vs. numérica
ggpairs(select_if(clini_data_no_na, is.numeric), lower = list(continuous = "smooth"))

## 1.2.1 Variável numérica vs. categórica ----

### Remover outliers
### get outliers
out_age <- which(clini_data_no_na$diagnosis_age>(quantile(clini_data_no_na$diagnosis_age, .75) + 1.5*IQR(clini_data_no_na$diagnosis_age)))
out_surv <- which(clini_data_no_na$surv_years>(quantile(clini_data_no_na$surv_years, .75) + 1.5*IQR(clini_data_no_na$surv_years)))
out_gen <- which(clini_data_no_na$genome_alt>(quantile(clini_data_no_na$genome_alt, .75) + 1.5*IQR(clini_data_no_na$genome_alt)))
out_tbl <- which(clini_data_no_na$tbl>(quantile(clini_data_no_na$tbl, .75) + 1.5*IQR(clini_data_no_na$tbl)))

### remove outliers
no_out_data <- clini_data_no_na[-c(out_age, out_surv, out_gen, out_tbl),]

## 1.2.1.a. Variável sex ----
### sex vs. diagnosis age
a <- ggplot(alpha=0.8, data=no_out_data, aes(x=sex, y=diagnosis_age, fill=sex)) + geom_violin() + geom_boxplot(width=.1) + scale_fill_manual(values = palette_3) + theme_bw()

### sex vs. genome altered
b <- ggplot(alpha=0.8, data=no_out_data, aes(x=sex, y=genome_alt, fill=sex))+ geom_violin()  + geom_boxplot(width=.1)  + scale_fill_manual(values = palette_3) + theme_bw()
ggarrange(a, b, common.legend=T)

### sex vs. tbl
a <- ggplot(alpha=0.8, data=no_out_data, aes(x=sex, y=tbl, fill=sex)) + geom_violin() + geom_boxplot(width=.1) + scale_fill_manual(values = palette_3) + theme_bw()

### sex vs. survival
b <- ggplot(alpha=0.8, data=no_out_data, aes(x=sex, y=surv_years, fill=sex)) + geom_violin() + geom_boxplot(width=.1) + scale_fill_manual(values = palette_3) + theme_bw()
ggarrange(a, b, common.legend=T)

## 1.2.1.b. Variável Cancer Type ----
### cancer type vs. diagnosis age
ggplot(alpha=0.8, data=no_out_data, aes(x=cancer_type, y=diagnosis_age, fill=sex)) + geom_violin() + geom_boxplot(width=.1,position = position_dodge(width =0.9)) +  scale_fill_manual(values = palette_3) + theme_bw()

### cancer type vs. genome altered
ggplot(alpha=0.8, data=no_out_data, aes(x=cancer_type, y=genome_alt, fill=sex)) + geom_violin() + geom_boxplot(width=.1,position = position_dodge(width =0.9)) +  scale_fill_manual(values = palette_3) + theme_bw()

### cancer type vs. tbl
a <- ggplot(alpha=0.8, data=no_out_data, aes(x=cancer_type, y=tbl, fill=cancer_type)) + geom_violin() + geom_boxplot(width=.1,position = position_dodge(width =0.9)) +  scale_fill_manual(values = palette_3) + theme_bw()

### cancer type vs. years survival
b <- ggplot(alpha=0.8, data=no_out_data, aes(x=cancer_type, y=surv_years, fill=cancer_type)) + geom_violin() + geom_boxplot(width=.2,position = position_dodge(width =0.9)) +  scale_fill_manual(values = palette_3) + theme_bw()
ggarrange(a, b, common.legend=T)


## 1.2.1.c. Variável Ancestry ----
### ancestry vs. genome altered
a <- ggplot(no_out_data, aes(x = ancestry, y = genome_alt)) +
  geom_boxplot(fill = palette_3[1:length(unique(no_out_data$ancestry))], alpha = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.8, color = "grey10")   +
  xlab("") + 
  theme_minimal()

### ancestry vs. age
b <- ggplot(no_out_data, aes(x = ancestry, y = diagnosis_age)) +
  geom_boxplot(fill = palette_3[1:length(unique(no_out_data$ancestry))], alpha = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.8, color = "grey10") +
  xlab("") +
  theme_minimal()

ggarrange(a, b, common.legend=T)

### ancestry vs. tbl
a <- ggplot(no_out_data, aes(x = ancestry, y = tbl)) +
  geom_boxplot(fill = palette_3[1:length(unique(no_out_data$ancestry))], alpha = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.8, color = "grey10")  + 
  xlab("") +
  theme_minimal()

### ancestry vs. survival
b <- ggplot(no_out_data, aes(x = ancestry, y = surv_years)) +
  geom_boxplot(fill = palette_3[1:length(unique(no_out_data$ancestry))], alpha = 0.5) + 
  geom_jitter(width = 0.2, alpha = 0.8, color = "grey10")  + 
  xlab("") +
  theme_minimal()

ggarrange(a, b, common.legend=T)


## 1.2.1.d. Variável surv_status ----
### sex vs. diagnosis age
a <- ggplot(alpha=0.8, data=no_out_data, aes(x=surv_status, y=diagnosis_age, fill=surv_status)) + geom_violin() + geom_boxplot(width=.1) + scale_fill_manual(values = palette_3) + theme_bw()

### sex vs. genome altered
b <- ggplot(alpha=0.8, data=no_out_data, aes(x=surv_status, y=genome_alt, fill=surv_status))+ geom_violin()  + geom_boxplot(width=.1)  + scale_fill_manual(values = palette_3) + theme_bw()

ggarrange(a, b, common.legend=T)

### sex vs. tbl
a <- ggplot(alpha=0.8, data=no_out_data, aes(x=surv_status, y=tbl, fill=surv_status)) + geom_violin() + geom_boxplot(width=.1) + scale_fill_manual(values = palette_3) + theme_bw()

### sex vs. survival
b <- ggplot(alpha=0.8, data=no_out_data, aes(x=surv_status, y=surv_years, fill=surv_status)) + geom_violin() + geom_boxplot(width=.1) + scale_fill_manual(values = palette_3) + theme_bw()

ggarrange(a, b, common.legend=T)

# 1.b. Preparação e pré-processamento dos dados ----

# 1.1. Dados de Expressão Génica ----

## 2.1. Estrutura dos dados ----

# Preview dos dados
head(gene_data[,1:5])

# Estrutura dos dados
str(gene_data[,1:10])

## 2.2. Remoção de duplicados e NAs ----

# Gene descontinuado DGCR9, row nº 4886 - será removido
gene_data_nd <- gene_data[-c(4886),]

# Remove duplicates
gene_data_nd <- gene_data_nd[-c(which(duplicated(gene_data_nd$Entrez_Gene_Id))),]

# Remoção de NAs
gene_data_nona <- na.omit(gene_data_nd)

# Proporção de NAs no dataset
(nrow(gene_data_nd) - nrow(gene_data_nona))/nrow(gene_data_nd)

# Converter Gene ID em rownames
rownames(gene_data_nona) <- gene_data_nona$Entrez_Gene_Id

# Exclusão de pacientes sem dados clínicos
complete_ids <- clini_data_no_na$sample_ID
gene_data_nona <- gene_data_nona[,complete_ids]

## 2.3. Tratamento de dados ----

# Converter em formato DGEList (package edgeR)
## sem agrupamento
d0 <- DGEList(gene_data_nona) # no grouping

# criação de grupos
f_cancer <- revalue(clini_data_no_na$cancer_type, c("Astrocytoma" = "A", "Oligoastrocytoma" = "B", "Oligodendroglioma" = "C", "Low-Grade Glioma (NOS)" = "D"))
f_survival <- revalue(clini_data_no_na$surv_status, c("0:ALIVE OR DEAD TUMOR FREE" = "A", "1:DEAD WITH TUMOR" = "B"))
f_ancestry <- revalue(clini_data_no_na$ancestry, c("EUR" = "A", "AFR" = "B", "AFR_ADMIX" = "C", "EAS" = "D", "EUR_ADMIX" = "E", "AMR" = "F", "SAS" = "G", "SAS_ADMIX" = "H"))
f_sex <- revalue(clini_data_no_na$sex, c("Male" = "A", "Female" = "B"))

# agrupamento por fatores
# d0_g <- DGEList(gene_data_nona, group=f_cancer) # or f_sex, etc.

# Visualizar estrutura dos dados
head(d0$samples)

## 2.4. Remoção de genes pouco expressos ----

# Excluir genes pouco expressos
keep <- filterByExpr(d0) 
d0_f <- d0[keep, , keep.lib.sizes=FALSE] # passámos de uma lista de 20504 elementos para 14357

## 2.5. Filtragem por níveis de expressão ----

# Filtros flat patterns - 2*mediana

gene_exp_sd <- apply(gene_data_nona, 1, sd)
d_med <- 2*median(gene_exp_sd)
high_exp <- which(gene_exp_sd > d_med)

mean_exp <- data.frame(gene_id=high_exp,exp=apply(gene_data_nona[high_exp,], 1, mean), row.names = c())
mean_exp$exp_max <- mean_exp$exp/max(mean_exp$exp)
mean_exp$exp_log <- log(mean_exp$exp)

# Filtros flat patterns - max/min > 2
gene_min <- apply(gene_data_nona, 1, min)
gene_max <- apply(gene_data_nona, 1, max)

rat <- which(gene_max/gene_min > 2)
mean_exp_f2 <- data.frame(id=rat, exp=apply(gene_data_nona[rat,-1], 1, mean))
mean_exp_f2$exp_max <- mean_exp_f2$exp/max(mean_exp_f2$exp)
mean_exp_f2$exp_log <- log(mean_exp_f2$exp)

# 1.b. Vizualização de dados ----
## 2.1. Sumarização e visualização de dados ----
# Filtro flat pattern 1
a <- ggplot(mean_exp) + geom_histogram(aes(x=exp), bins=400, fill="#7eb0d5") +ggtitle("Flat pattern 1")
b <- ggplot(mean_exp) + geom_histogram(aes(x=exp_log), bins=400, fill="#7eb0d5") + ggtitle("Flat pattern 1 - log transformation")

grid.arrange(a, b, ncol=2)

# Filtro flat pattern 2
a <- ggplot(mean_exp_f2) + geom_histogram(aes(x=exp), bins=200, fill="#7eb0d5") + ggtitle("Flat pattern 2")
b <- ggplot(mean_exp_f2) + geom_histogram(aes(x=exp_log), bins=200, fill="#7eb0d5") + ggtitle("Flat pattern 2 - log transformation")

grid.arrange(a, b, ncol=2)


## 2.2. Ontologia de genes ----
# 20 genes mais expressos
all_genes_means <- data.frame(gene_id=c(rownames(gene_data_nona)),exp=apply(gene_data_nona, 1, mean))
all_genes_means <- all_genes_means[order(all_genes_means$exp, decreasing=T),]
high_20 <- head(all_genes_means, n=20)

# Visualização das funções
go_prof <- gost(high_20$gene_id, organism = "hsapiens")
gostplot(go_prof, capped=F, interactive = T)


# 2. Análise Univariada ---- 

# Preparar dados clínicos e fatores
clean_data_for_2part <- clini_data_no_na

f_sex <- as.factor(clean_data_for_2part$sex)
names(f_sex) <- clean_data_for_2part$sample_ID

f_cancer <- as.factor(clean_data_for_2part$cancer_type)
names(f_cancer) <- clean_data_for_2part$sample_ID

f_ancestralidade <- as.factor(clean_data_for_2part$ancestry)
names(f_ancestralidade) <- clean_data_for_2part$sample_ID

f_survival <- as.factor(clean_data_for_2part$surv_status)
names(f_survival) <- clean_data_for_2part$sample_ID

fatores <- list(
  Sexo = f_sex,
  TipoTumor = f_cancer,
  Ancestralidade = f_ancestralidade,
  Sobrevivencia = f_survival
)

# Expressão log2-CPM
log_cpm <- cpm(d0_f, log = TRUE)

# Inicializar variáveis
resultados_df <- data.frame()
plot_list <- list()
i <- 1

# Loop por genes e variáveis clínicas
for (gene_id in high_20$gene_id) {
  for (nome_var in names(fatores)) {
    grupo <- fatores[[nome_var]]
    amostras_comuns <- intersect(colnames(log_cpm), names(grupo))
    if (length(amostras_comuns) < 2) next
    
    grupo_valid <- grupo[amostras_comuns]
    expr <- log_cpm[as.character(gene_id), amostras_comuns]
    
    if (length(unique(grupo_valid)) < 2) next
    
    df_plot <- data.frame(expr = expr, grupo = grupo_valid)
    
    # Boxplot
    p <- ggplot(df_plot, aes(x = grupo, y = expr, fill=grupo)) +
      geom_boxplot() +
      scale_fill_manual(values = palette_3) + theme_bw() +
      labs(title = paste(gene_id, "-", nome_var), x = nome_var, y = "Log2 CPM")
    
    plot_list[[i]] <- p
    i <- i + 1
    
    # Teste estatístico
    p_val <- if (length(unique(grupo_valid)) == 2) {
      t.test(expr ~ grupo_valid)$p.value
    } else {
      summary(aov(expr ~ grupo_valid))[[1]][["Pr(>F)"]][1]
    }
    
    resultados_df <- rbind(resultados_df, data.frame(
      gene = gene_id,
      variavel_clinica = nome_var,
      p_value = p_val
    ))
  }
}

# Mostrar os boxplots em grupos de 4 com espaços para comentários
text_chunks <- c(
  "Análise de tendências observadas nesta secção.",
  "Exploração da variabilidade entre diferentes grupos de amostras.",
  "Investigação de correlações e possíveis padrões.",
  "Revisão de resultados e implicações."
)

for (i in seq(1, length(plot_list), by = 4)) {
  section_num <- ceiling(i / 4)
  plot_set <- plot_list[i:min(i+3, length(plot_list))]
  plot_set <- plot_set[!sapply(plot_set, is.null)]
  
  if (length(plot_set) > 0) {
    grid.arrange(grobs = plot_set, ncol = 2, nrow = 2)
  }
  if (section_num <= length(text_chunks)) {
    html_chunk <- HTML(paste0(
      '<div style="background-color:#e6f2ff; padding:10px; border-left:5px solid #337ab7; margin-top:10px; margin-bottom:10px;">',
      ">", text_chunks[section_num],
      '</div>'
    ))
    print(html_chunk)
  }
}

# Mostrar p-values em tabela
datatable(resultados_df, caption = "p-values da análise univariada", options = list(pageLength = 10))

# 3. Análise de Expressão Diferencial ----

for (nome_var in names(fatores)) {
cat("\n\n### Analisando expressão diferencial para:", nome_var, "###\n")

fator <- fatores[[nome_var]]
amostras_comuns <- intersect(colnames(d0_f), names(fator))
grupo <- fator[amostras_comuns]
y <- d0_f[, amostras_comuns]
grupo <- droplevels(grupo[!is.na(grupo)])
y <- y[, names(grupo)]

if (length(unique(grupo)) < 2) next

design <- model.matrix(~ grupo)
colnames(design) <- make.names(colnames(design))
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = 2)
top_genes <- topTags(qlf, n = Inf)
topGenesTable <- top_genes$table
topGenesTable$threshold <- as.factor(abs(topGenesTable$logFC) > 1 & topGenesTable$FDR < 0.05)

print(ggplot(topGenesTable, aes(x = logFC, y = -log10(FDR), color = threshold)) +
        geom_point(alpha = 0.6) +
        scale_color_manual(values = c("grey", "red")) +
        theme_bw() +
        ggtitle(paste("Volcano plot -", nome_var)))

cat("Total de DEGs com FDR < 0.05 e |logFC| > 1:", sum(topGenesTable$threshold == TRUE), "\n")
datatable(head(topGenesTable[topGenesTable$threshold == TRUE, ], 20), caption = paste("Top DEGs -", nome_var))
}

# Loop por cada variável categórica
for (nome_var in names(fatores)) {
  cat(paste0("\n\n### Enriquecimento funcional para: ", nome_var, " ###\n"))
  
  # Subconjunto de amostras e fatores
  fator <- fatores[[nome_var]]
  amostras_comuns <- intersect(colnames(d0_f), names(fator))
  grupo <- droplevels(fator[amostras_comuns])
  
  if (length(unique(grupo)) < 2) {
    cat("Fator sem pelo menos 2 grupos distintos. Pulando.\n")
    next
  }

  y <- d0_f[, names(grupo)]
  design <- model.matrix(~ grupo)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef = 2)
  
  top_genes <- topTags(qlf, n = Inf)
  deg_table <- top_genes$table
  deg_genes <- deg_table[deg_table$FDR < 0.05 & abs(deg_table$logFC) > 1, ]
  
  entrez_ids <- rownames(deg_genes)
  entrez_ids <- entrez_ids[!is.na(entrez_ids)]

  cat("Total de genes diferenciais para enriquecimento:", length(entrez_ids), "\n\n")

  if (length(entrez_ids) > 10) {
    ## GO: Biological Processes
    ego <- enrichGO(
      gene = entrez_ids,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      qvalueCutoff = 0.05,
      readable = TRUE
    )

    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
      cat("**Dotplot GO - Processos Biológicos:**\n")
      print(dotplot(ego, showCategory = 20) +
              theme(axis.text.y = element_text(size = 10)) +
              scale_y_discrete(labels = function(x) str_wrap(x, width = 40)))
      
      datatable(as.data.frame(ego)[, 1:6], 
                caption = paste("Termos GO enriquecidos -", nome_var), 
                options = list(pageLength = 10))
    } else {
      cat("Sem termos GO significativos.\n")
    }

    ## KEGG Pathways
    ekegg <- enrichKEGG(
      gene = entrez_ids,
      organism = 'hsa',
      pAdjustMethod = "BH",
      qvalueCutoff = 0.05
    )

    if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
      cat("**Dotplot KEGG - Vias metabólicas:**\n")
      print(dotplot(ekegg, showCategory = 20) +
              theme(axis.text.y = element_text(size = 10)) +
              scale_y_discrete(labels = function(x) str_wrap(x, width = 40)))
      
      datatable(as.data.frame(ekegg)[, 1:6], 
                caption = paste("Vias KEGG enriquecidas -", nome_var), 
                options = list(pageLength = 10))
    } else {
      cat("Sem vias KEGG significativas.\n")
    }

  } else {
    cat("Genes insuficientes para enriquecimento funcional.\n")
  }

  cat("\n\n---\n\n")
}

# 4. Enriquecimento Funcional ----

# Loop por cada variável categórica
for (nome_var in names(fatores)) {
  cat(paste0("\n\n### Enriquecimento funcional para: ", nome_var, " ###\n"))
  
  # Subconjunto de amostras e fatores
  fator <- fatores[[nome_var]]
  amostras_comuns <- intersect(colnames(d0_f), names(fator))
  grupo <- droplevels(fator[amostras_comuns])
  
  if (length(unique(grupo)) < 2) {
    cat("Fator sem pelo menos 2 grupos distintos. Pulando.\n")
    next
  }

  y <- d0_f[, names(grupo)]
  design <- model.matrix(~ grupo)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef = 2)
  
  top_genes <- topTags(qlf, n = Inf)
  deg_table <- top_genes$table
  deg_genes <- deg_table[deg_table$FDR < 0.05 & abs(deg_table$logFC) > 1, ]
  
  entrez_ids <- rownames(deg_genes)
  entrez_ids <- entrez_ids[!is.na(entrez_ids)]

  cat("Total de genes diferenciais para enriquecimento:", length(entrez_ids), "\n\n")

  if (length(entrez_ids) > 10) {
    ## GO: Biological Processes
    ego <- enrichGO(
      gene = entrez_ids,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      qvalueCutoff = 0.05,
      readable = TRUE
    )

    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
      cat("**Dotplot GO - Processos Biológicos:**\n")
      print(dotplot(ego, showCategory = 20) +
              theme(axis.text.y = element_text(size = 10)) +
              scale_y_discrete(labels = function(x) str_wrap(x, width = 40)))
      
      datatable(as.data.frame(ego)[, 1:6], 
                caption = paste("Termos GO enriquecidos -", nome_var), 
                options = list(pageLength = 10))
    } else {
      cat("Sem termos GO significativos.\n")
    }

    ## KEGG Pathways
    ekegg <- enrichKEGG(
      gene = entrez_ids,
      organism = 'hsa',
      pAdjustMethod = "BH",
      qvalueCutoff = 0.05
    )

    if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
      cat("**Dotplot KEGG - Vias metabólicas:**\n")
      print(dotplot(ekegg, showCategory = 20) +
              theme(axis.text.y = element_text(size = 10)) +
              scale_y_discrete(labels = function(x) str_wrap(x, width = 40)))
      
      datatable(as.data.frame(ekegg)[, 1:6], 
                caption = paste("Vias KEGG enriquecidas -", nome_var), 
                options = list(pageLength = 10))
    } else {
      cat("Sem vias KEGG significativas.\n")
    }

  } else {
    cat("Genes insuficientes para enriquecimento funcional.\n")
  }

  cat("\n\n---\n\n")
}

# 5. Análise PCA ----

# Usar todos os genes da matriz de expressão
expr_all <- log_cpm  # log2 CPM

# Transpor: linhas = amostras, colunas = genes
expr_all_t <- t(expr_all)

# PCA
pca_res <- prcomp(expr_all_t, scale. = TRUE)

# Criar diretório para salvar resultados se necessário
dir.create("resultados_PCA", showWarnings = FALSE)

# Loop pelas variáveis clínicas
for (nome_var in names(fatores)) {
  cat("### PCA para:", nome_var, "\n\n")
  
  grupo_var <- fatores[[nome_var]]
  
  # Alinhar amostras
  amostras_comuns <- intersect(rownames(expr_all_t), names(grupo_var))
  
  if (length(amostras_comuns) < 2) {
    cat("Variável", nome_var, "não tem amostras suficientes. Pulando.\n\n")
    next
  }
  
  grupo_valid <- grupo_var[amostras_comuns]
  grupo_valid <- as.factor(grupo_valid)
  
  if (length(unique(grupo_valid)) < 2) {
    cat("Variável", nome_var, "não tem pelo menos 2 grupos distintos. Pulando.\n\n")
    next
  }
  
  # Preparar dados PCA para o gráfico
  pca_data <- as.data.frame(pca_res$x[amostras_comuns, 1:2])
  pca_data$grupo <- grupo_valid
  
  # Gráfico PCA
  p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = grupo)) +
    geom_point(size = 3, alpha = 0.8) +scale_fill_manual(values = palette_3) + theme_bw() +
    labs(title = paste("PCA -", nome_var),
         x = "PC1", y = "PC2")
  
  
  print(p)  # Mostrar no HTML
  
  cat("\n\n---\n\n")  # Separador entre gráficos no relatório
}

#7. Análise UMAP com PC´s**
                               
# Reduzir para 30 componentes antes do UMAP
pca_matrix <- pca_res$x[, 1:30]

# UMAP diretamente sobre os PCs
set.seed(123)  # para reprodutibilidade
umap_res <- umap(pca_matrix, n_neighbors = 15, min_dist = 0.1, metric = "euclidean")

# Converter resultado UMAP para data frame com nomes de amostras
rownames(umap_res) <- rownames(expr_all_t)
umap_df <- as.data.frame(umap_res)
colnames(umap_df) <- c("UMAP1", "UMAP2")


# Loop pelas variáveis clínicas 
for (nome_var in names(fatores)) {
  cat("### UMAP para:", nome_var, "\n\n")
  
  grupo_var <- fatores[[nome_var]]
  amostras_comuns <- intersect(rownames(umap_df), names(grupo_var))
  
  grupo_valid <- as.factor(grupo_var[amostras_comuns])
  umap_plot_df <- umap_df[amostras_comuns, ]
  umap_plot_df$Grupo <- grupo_valid

  # Gráfico
  p <- ggplot(umap_plot_df, aes(x = UMAP1, y = UMAP2, color = Grupo)) +
    geom_point(size = 3, alpha = 0.8) +
    theme_bw() +
    labs(title = paste("UMAP -", nome_var))
  
  print(p)
  
  cat("\n\n---\n\n")
}
## **8. Análise T-SNE com PC´s**

# t-SNE
set.seed(123)
tsne_res <- Rtsne(pca_matrix, dims = 2, perplexity = 30, verbose = FALSE)
tsne_df <- as.data.frame(tsne_res$Y)
rownames(tsne_df) <- rownames(pca_matrix)
colnames(tsne_df) <- c("tSNE1", "tSNE2")

# Loop pelas variáveis clínicas
for (nome_var in names(fatores)) {
  cat("### t-SNE para:", nome_var, "\n\n")
  
  grupo_var <- fatores[[nome_var]]
  amostras_comuns <- intersect(rownames(tsne_df), names(grupo_var))
  
  grupo_valid <- as.factor(grupo_var[amostras_comuns])
  tsne_plot_df <- tsne_df[amostras_comuns, ]
  tsne_plot_df$Grupo <- grupo_valid

  # Gráfico
  p <- ggplot(tsne_plot_df, aes(x = tSNE1, y = tSNE2, color = Grupo)) +
    geom_point(size = 3, alpha = 0.8) +
    theme_bw() + labs(title = paste("t-SNE -", nome_var))

  # Mostrar e salvar
  print(p)
  
  cat("\n\n---\n\n")
}


## **9. UMAP vs T-SNE**
## **10. Clustering**
# Selecionar os 30 genes com menor p-value
top_30 <- head(deg_table[order(deg_table$PValue), ], 30)

# Usar rownames como identificador dos genes (caso do PC da colega)
gene_ids <- rownames(top_30)

# Obter dados de expressão para os genes válidos
genes_matrix <- expr_all[gene_ids, ]

# Normalizar por gene (z-score)
genes_matrix_scaled <- t(scale(t(genes_matrix)))

# Clustering hierárquico
dist_matrix <- dist(genes_matrix_scaled)
hc <- hclust(dist_matrix, method = "complete")

# Plotar dendrograma
plot(hc,
     main = paste("Clustering hierárquico de", length(gene_ids), "genes mais significativos"),
     ylab = "Distância Euclidiana",
     xlab = "",
     sub = "")



library(factoextra)

# Selecionar os 30 genes com menor p-value
top_30 <- head(deg_table[order(deg_table$PValue), ], 30)

# Usar rownames como identificador dos genes (caso do PC da colega)
gene_ids <- rownames(top_30)

# Filtrar dados de expressão
genes_matrix <- expr_all[rownames(expr_all) %in% gene_ids, ]

# Normalizar por gene (z-score por linha)
genes_matrix_scaled <- t(scale(t(genes_matrix)))

# Método do cotovelo para escolher o número ideal de clusters
fviz_nbclust(genes_matrix_scaled, kmeans, method = "wss") +
  labs(title = "Número ótimo de clusters\nMétodo do Cotovelo")


set.seed(123)
k <- 3
kmeans_result <- kmeans(genes_matrix_scaled, centers = k, nstart = 25)

# Visualizar o resultado do clustering
fviz_cluster(kmeans_result, data = genes_matrix_scaled,
             main = "Clustering K-means dos 30 genes mais significativos",
             geom = "point", ellipse.type = "convex", 
             palette = palette_3,repel=TRUE)

## **11. Análise supervisionada**

# 1. Modelos possíveis 
# A. Random Forest
# B. Support Vector Machine
# c. XGBoost

# 2. Treino/teste
# A. train_test_split
# B. cross-validation

# 3. Avaliação (métricas)
# A. Accuracy
# B. F1-score
# C. AUT-ROC


# Preparar dados clínicos e fatores
clean_data_for_2part <- clini_data_no_na

f_sex <- as.factor(clean_data_for_2part$sex)
names(f_sex) <- clean_data_for_2part$sample_ID

f_cancer <- as.factor(clean_data_for_2part$cancer_type)
names(f_cancer) <- clean_data_for_2part$sample_ID

f_ancestralidade <- as.factor(clean_data_for_2part$ancestry)
names(f_ancestralidade) <- clean_data_for_2part$sample_ID

f_survival <- as.factor(clean_data_for_2part$surv_status)
names(f_survival) <- clean_data_for_2part$sample_ID

fatores <- list(
  Sexo = f_sex,
  TipoTumor = f_cancer,
  Ancestralidade = f_ancestralidade,
  Sobrevivencia = f_survival
)

# Base trainControl parameters
base_ctrl <- list(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  savePredictions = "all"
)

ctrl <- do.call(trainControl, base_ctrl)
ctrl1 <- do.call(trainControl, c(base_ctrl, list(summaryFunction = multiClassSummary)))

### **11.1 Random Forest**
## A. Random Forest com smote ---- SURVIVAL

# Verificar consistência dos dados
common_samples <- intersect(colnames(expr_all), names(f_survival))
logCPM_sub <- expr_all[, common_samples]
labels <- f_survival[common_samples]

# Criar dataframe para modelar
data_model <- as.data.frame(t(logCPM_sub))  # Genes como colunas, amostras como linhas
data_model$SurvivalStatus <- factor(labels)

# Garantir que os levels são válidos para SMOTE
levels(data_model$SurvivalStatus) <- make.names(levels(data_model$SurvivalStatus))

# Dividir treino e teste
set.seed(123)
trainIndex <- createDataPartition(data_model$SurvivalStatus, p = 0.5, list = FALSE)
trainData <- data_model[trainIndex, ]
testData  <- data_model[-trainIndex, ]

# Aplicar SMOTE
# Remover colunas não numéricas
trainData_numeric <- trainData[, sapply(trainData, is.numeric)]

# SMOTE para balancear classe
smote_result_rf <- SMOTE(trainData_numeric, trainData$SurvivalStatus, K = 3, dup_size = 2)
smote_data_rf <- smote_result_rf$data

# Ajustar nomes das colunas
colnames(smote_data_rf)[ncol(smote_data_rf)] <- "SurvivalStatus"
smote_data_rf$SurvivalStatus <- factor(smote_data_rf$SurvivalStatus, 
                                       levels = levels(trainData$SurvivalStatus))

# Treinar modelo Random Forest
set.seed(123)
fit_rf1 <- train(SurvivalStatus ~ ., data = smote_data_rf, method = "rf", 
                trControl = ctrl,
                tuneGrid = expand.grid(mtry = 2),  # Ajuste do número de variáveis
                ntree = 100)  # Número de árvores ajustado

# Avaliação no teste
pred_rf <- predict(fit_rf1, testData)
conf_matrix <- confusionMatrix(pred_rf, testData$SurvivalStatus)
pander::pander(conf_matrix, caption = "Confusion Matrix Statistics (per class)")

# Precision, Recall, F1
positive_class <- "X0.ALIVE.OR.DEAD.TUMOR.FREE"

precision <- posPredValue(pred_rf, testData$SurvivalStatus, positive = positive_class)
recall <- sensitivity(pred_rf, testData$SurvivalStatus, positive = positive_class)
f1 <- 2 * (precision * recall) / (precision + recall)

# Print precision, recall,F1
cat("Precision:", precision, "\n")
cat("Recall:", recall, "\n")
cat("F1 Score:", round(f1, 3), "\n")

# Calcular AUC e Curva ROC
roc_curve <- roc(testData$SurvivalStatus, as.numeric(pred_rf))
auc_score <- auc(roc_curve)
print(paste("AUC Score:", auc_score))

plot(roc_curve, col = palette_3, main = "ROC Curve for Random Forest (com SMOTE)")

# Ver importância dos genes
importance <- varImp(fit_rf1)

# Extract variable importance as data frame
importance_df <- as.data.frame(importance$importance)
importance_df$Gene <- rownames(importance_df)
importance_df <- importance_df[order(-importance_df$Overall), ]  # Sort by importance

# Select top 20
top_20 <- head(importance_df, 20)

# Tabela
knitr::kable(top_20[, c("Gene", "Overall")], caption = "Top 20 Most Important Genes")

plot(importance, top = 20, main = "Top 20 Important Features")
# Visualizar distribuição de probabilidades
prob_rf <- predict(fit_rf1, testData, type = "prob")
boxplot(prob_rf[,1] ~ testData$SurvivalStatus, col = palette_3, 
        main = "Class Probability Distribution", ylab = "Predicted Probability")

## A. Random Forest com SMOTE ---- SEX

# Verificar consistência dos dados
common_samples <- intersect(colnames(expr_all), names(f_sex))
logCPM_sub <- expr_all[, common_samples]
labels <- f_sex[common_samples]

# Criar dataframe para modelagem
data_model <- as.data.frame(t(logCPM_sub))
data_model$Sex <- factor(labels)

# Garantir que os niveis são válidos para SMOTE
levels(data_model$Sex) <- make.names(levels(data_model$Sex))

# Dividir treino e teste
set.seed(123)
trainIndex <- createDataPartition(data_model$Sex, p = 0.5, list = FALSE)
trainData <- data_model[trainIndex, ]
testData  <- data_model[-trainIndex, ]

# Aplicar SMOTE
# Remover colunas não numéricas
trainData_numeric <- trainData[, sapply(trainData, is.numeric)]

# SMOTE para balancear classe
smote_result_rf <- SMOTE(trainData_numeric, trainData$Sex, K = 3, dup_size = 2)
smote_data_rf <- smote_result_rf$data

# Ajustar nomes das colunas
colnames(smote_data_rf)[ncol(smote_data_rf)] <- "Sex"
smote_data_rf$Sex <- factor(smote_data_rf$Sex, levels = levels(trainData$Sex))

# Treinar modelo Random Forest
set.seed(123)
fit_rf2 <- train(Sex ~ ., data = smote_data_rf, method = "rf", 
                trControl = ctrl, 
                tuneGrid = expand.grid(mtry = 2),
                ntree = 100)

# Avaliação no teste
pred_rf <- predict(fit_rf2, testData)
conf_matrix <- confusionMatrix(pred_rf, testData$Sex)
pander::pander(conf_matrix, caption = "Confusion Matrix Statistics (per class)")
# F1, Recall, precision
testData$Sex <- factor(testData$Sex)                       
pred_rf <- factor(pred_rf, levels = levels(testData$Sex))  

# Define positive class
positive_class <- "Female"

# Check precision e recall para classe positiva
precision <- posPredValue(pred_rf, testData$Sex, positive = positive_class)
recall <- sensitivity(pred_rf, testData$Sex, positive = positive_class)

# Print precision e recall
cat("Precision:", precision, "\n")
cat("Recall:", recall, "\n")
# F1 score
if (!is.na(precision) && !is.na(recall) && (precision + recall) > 0) {
  f1 <- 2 * (precision * recall) / (precision + recall)
  cat("F1 Score:", round(f1, 3), "\n")
} else {
  cat("F1 Score cannot be computed due to missing precision or recall.\n")
}

# Calcular AUC e Curva ROC
roc_curve <- roc(testData$Sex, as.numeric(pred_rf))
auc_score <- auc(roc_curve)
print(paste("AUC Score:", auc_score))

plot(roc_curve, col = palette_3, main = "ROC Curve for Random Forest (com SMOTE)")
# Ver importância dos genes
importance <- varImp(fit_rf2)

# Extrair variavel importance como data frame
importance_df <- as.data.frame(importance$importance)
importance_df$Gene <- rownames(importance_df)
importance_df <- importance_df[order(-importance_df$Overall), ]  # Sort by importance

# Seleção top 20
top_20 <- head(importance_df, 20)

# Tabela
knitr::kable(top_20[, c("Gene", "Overall")], caption = "Top 20 Most Important Genes")

plot(importance, top = 20, main = "Top 20 Important Features")
# Visualizar distribuição de probabilidades
prob_rf <- predict(fit_rf2, testData, type = "prob")
boxplot(prob_rf[,1] ~ testData$Sex, 
        col = palette_3, 
        main = "Class Probability Distribution", ylab = "Predicted Probability")



## A. Random Forest com SMOTE ---- Ancestralidade

# Verificar consistência dos dados
common_samples <- intersect(colnames(expr_all), names(f_ancestralidade))
logCPM_sub <- expr_all[, common_samples]
labels <- f_ancestralidade[common_samples]

# Criar dataframe para modelagem
data_model <- as.data.frame(t(logCPM_sub))  
data_model$Ancestralidade <- factor(labels)

# Garantir que os levels são válidos para SMOTE
levels(data_model$Ancestralidade) <- make.names(levels(data_model$Ancestralidade))

# Dividir treino e teste
set.seed(123)
trainIndex <- createDataPartition(data_model$Ancestralidade, p = 0.5, list = FALSE)
trainData <- data_model[trainIndex, ]
testData  <- data_model[-trainIndex, ]

# Conta ocorrencias por classe
class_counts <- table(trainData$Ancestralidade)

# Filtrar classes <= 5 samples
valid_classes <- names(class_counts[class_counts >= 5])
trainData <- trainData[trainData$Ancestralidade %in% valid_classes, ]

# Niveis fatoriais atualizados
trainData$Ancestralidade <- factor(trainData$Ancestralidade)

# Smote
trainData_numeric <- trainData[, sapply(trainData, is.numeric)]
min_class_size <- min(table(trainData$Ancestralidade))
K_value <- max(1, min_class_size - 1)

smote_result_rf <- SMOTE(trainData_numeric, trainData$Ancestralidade, K = K_value, dup_size = 2)
smote_data_rf <- smote_result_rf$data

colnames(smote_data_rf)[ncol(smote_data_rf)] <- "Ancestralidade"
smote_data_rf$Ancestralidade <- factor(smote_data_rf$Ancestralidade, levels = levels(trainData$Ancestralidade))
smote_data_rf <- na.omit(smote_data_rf)


# Treinar modelo Random Forest
set.seed(123)
fit_rf3 <- train(Ancestralidade ~ ., data = smote_data_rf, method = "rf", 
                trControl = ctrl1,
                tuneGrid = expand.grid(mtry = 2),
                ntree = 100)

# Avaliação no teste
pred_rf <- predict(fit_rf3, testData)
true_labels <- testData$Ancestralidade

# Confusion Matrix
all_levels <- union(levels(factor(pred_rf)), levels(factor(true_labels)))
pred_rf <- factor(pred_rf, levels = all_levels)
true_labels <- factor(true_labels, levels = all_levels)
conf_matrix <- confusionMatrix(pred_rf, true_labels)
pander::pander(conf_matrix, caption = "Confusion Matrix Statistics (per class)")

# Calcular métricas por classe
precision_vec <- diag(conf_matrix$table) / colSums(conf_matrix$table)
recall_vec <- diag(conf_matrix$table) / rowSums(conf_matrix$table)
f1_vec <- ifelse((precision_vec + recall_vec) > 0,
                 2 * (precision_vec * recall_vec) / (precision_vec + recall_vec),
                 NA)

f1_scores <- data.frame(
  Class = rownames(conf_matrix$table),
  Precision = round(precision_vec, 3),
  Recall = round(recall_vec, 3),
  F1 = round(f1_vec, 3)
)

# Mostrar tabela de métricas
knitr::kable(f1_scores, caption = "Precision, Recall, and F1 Score per Class")
roc_curve <- roc(testData$Ancestralidade, as.numeric(pred_rf))
auc_score <- auc(roc_curve)
print(paste("AUC Score:", auc_score))
plot(roc_curve, col = palette_3, main = "ROC Curve for Random Forest (com SMOTE)")
# Importância dos genes
importance <- varImp(fit_rf3)

# Extração variavel importance como data frame
importance_df <- as.data.frame(importance$importance)
importance_df$Gene <- rownames(importance_df)
importance_df <- importance_df[order(-importance_df$Overall), ]  # Sort by importance

# Seleção top 20
top_20 <- head(importance_df, 20)

# Tabela
knitr::kable(top_20[, c("Gene", "Overall")], caption = "Top 20 Most Important Genes")

plot(importance, top = 20, main = "Top 20 Important Features")
# Visualização da distribuição de probabilidades 
prob_rf <- predict(fit_rf3, testData, type = "prob")
boxplot(prob_rf[,1] ~ testData$Ancestralidade, col = palette_3, 
        main = "Class Probability Distribution", ylab = "Predicted Probability")
## A. Random Forest com SMOTE ---- Cancer

# Verificar consistência dos dados
common_samples <- intersect(colnames(expr_all), names(f_cancer))
logCPM_sub <- expr_all[, common_samples]
labels <- f_cancer[common_samples]  # Ajustando para usar f_cancer corretamente

# Criar dataframe para modelagem
data_model <- as.data.frame(t(logCPM_sub))  # Genes como colunas, amostras como linhas
data_model$Cancer <- factor(labels)  # Alterando nome para refletir f_cancer

# Garantir que os levels são válidos para SMOTE
levels(data_model$Cancer) <- make.names(levels(data_model$Cancer))

# Dividir treino e teste
set.seed(123)
trainIndex <- createDataPartition(data_model$Cancer, p = 0.5, list = FALSE)
trainData <- data_model[trainIndex, ]
testData  <- data_model[-trainIndex, ]

# Aplicar SMOTE
# Remover colunas não numéricas
trainData_numeric <- trainData[, sapply(trainData, is.numeric)]

# Encontrar o tamanho mínimo da classe
min_class_size <- min(table(trainData$Cancer))

# Ajustar K para evitar erro de SMOTE
K_value <- max(1, min_class_size - 1)  # Garante que K seja pelo menos 1

# Aplicar SMOTE com K ajustado
smote_result_rf <- SMOTE(trainData_numeric, trainData$Cancer, K = K_value, dup_size = 2)
smote_data_rf <- smote_result_rf$data

# Ajustar nomes das colunas
colnames(smote_data_rf)[ncol(smote_data_rf)] <- "Cancer"
smote_data_rf$Cancer <- factor(smote_data_rf$Cancer, levels = levels(trainData$Cancer))

# Remover valores ausentes
smote_data_rf <- na.omit(smote_data_rf)

# Treinar modelo Random Forest
set.seed(123)
fit_rf4 <- train(Cancer ~ ., data = smote_data_rf, method = "rf", 
                trControl = ctrl1,
                tuneGrid = expand.grid(mtry = 2),  # Ajuste do número de variáveis
                metric = "ROC",
                ntree = 100)  # Número de árvores ajustado

# Avaliação no teste
pred_rf <- predict(fit_rf4, testData)
conf_matrix <- confusionMatrix(pred_rf, testData$Cancer)
pander::pander(conf_matrix, caption = "Confusion Matrix Statistics (per class)")

# Tabela confusion matrix 
conf_table <- conf_matrix$table

# Inicializa data frame para armazenar metricas
class_metrics <- data.frame(Class = levels(testData$Cancer),
                            Precision = numeric(length(levels(testData$Cancer))),
                            Recall = numeric(length(levels(testData$Cancer))),
                            F1_Score = numeric(length(levels(testData$Cancer))))

# Ciclo sobre cada classe para calcular precisão, sensibilidade e F1 score
for (class in levels(testData$Cancer)) {
  # Verdadeiros Positivos (VP): Previsões corretas da classe
  TP <- conf_table[class, class]
  
  # Falsos Positivos (FP): Previsões da classe, mas pertencem a outra
  FP <- sum(conf_table[, class]) - TP
  
  # Falsos Negativos (FN): Pertencem à classe, mas foram previstos como outra
  FN <- sum(conf_table[class, ]) - TP
  
  # Precisão: VP / (VP + FP)
  precision <- ifelse((TP + FP) == 0, NA, TP / (TP + FP))
  
  # Sensibilidade (Recall): VP / (VP + FN)
  recall <- ifelse((TP + FN) == 0, NA, TP / (TP + FN))
  
  # F1 Score: 2 * (Precisão * Sensibilidade) / (Precisão + Sensibilidade)
  f1 <- ifelse(is.na(precision) | is.na(recall) | (precision + recall) == 0, NA, 2 * (precision * recall) / (precision + recall))
  
  # Guardar os resultados
  class_metrics[class_metrics$Class == class, ] <- c(class, precision, recall, f1)
}


# Converter colunas numéricas
class_metrics[, 2:4] <- lapply(class_metrics[, 2:4], as.numeric)

# Round
class_metrics_rounded <- class_metrics
class_metrics_rounded[, 2:4] <- round(class_metrics_rounded[, 2:4], 3)

# Tabela
knitr::kable(class_metrics_rounded, caption = "Per-Class Performance Metrics (Precision, Recall, F1 Score)")

# Calcular AUC e Curva ROC
roc_curve <- roc(testData$Cancer, as.numeric(pred_rf))
auc_score <- auc(roc_curve)
print(paste("AUC Score:", auc_score))

plot(roc_curve, col = palette_3, main = "ROC Curve for Random Forest (com SMOTE)")

# Ver importância dos genes
importance <- varImp(fit_rf4)

# Extração variavel importance como data frame
importance_df <- as.data.frame(importance$importance)
importance_df$Gene <- rownames(importance_df)
importance_df <- importance_df[order(-importance_df$Overall), ]  # Sort by importance

# Seleção top 20
top_20 <- head(importance_df, 20)

# Tabela
knitr::kable(top_20[, c("Gene", "Overall")], caption = "Top 20 Most Important Genes")

plot(importance, top = 20, main = "Top 20 Important Features")
# Visualizar distribuição de probabilidades
prob_rf <- predict(fit_rf4, testData, type = "prob")
boxplot(prob_rf[,1] ~ testData$Cancer, col = palette_3, 
        main = "Class Probability Distribution", ylab = "Predicted Probability")
### **11.2 Support Vector Machine**

## B. Support Vector Machine com SMOTE ---- SURVIVAL

# Preparar dados
common_samples <- intersect(colnames(expr_all), names(f_survival))
logCPM_sub <- expr_all[, common_samples]
labels <- f_survival[common_samples]

# Criar data frame com amostras como linhas e genes como colunas
data_model <- as.data.frame(t(logCPM_sub))
data_model$SurvivalStatus <- factor(labels)
levels(data_model$SurvivalStatus) <- make.names(levels(data_model$SurvivalStatus))  # tornar nomes válidos

# Dividir treino e teste
set.seed(123)
trainIndex <- createDataPartition(data_model$SurvivalStatus, p = 0.5, list = FALSE)
trainData <- data_model[trainIndex, ]
testData  <- data_model[-trainIndex, ]

# Aplicar SMOTE
trainData_numeric <- trainData[, sapply(trainData, is.numeric)]
smote_result <- SMOTE(trainData_numeric, target = trainData$SurvivalStatus, K = 3, dup_size = 2)
smote_data <- smote_result$data
colnames(smote_data)[ncol(smote_data)] <- "SurvivalStatus"
smote_data$SurvivalStatus <- factor(smote_data$SurvivalStatus,
                                    levels = levels(trainData$SurvivalStatus))

# Treinar modelo SVM
set.seed(123)
fit_svm1 <- train(SurvivalStatus ~ ., data = smote_data, method = "svmRadial",
                  trControl = ctrl,
                  preProcess = c("center", "scale"),
                  tuneLength=5,
                  metric = "ROC")

# Previsões no conjunto de teste
# Garantir que os níveis são iguais
testData$SurvivalStatus <- factor(testData$SurvivalStatus,
                                  levels = levels(fit_svm1$trainingData$.outcome))
pred_class <- predict(fit_svm1, testData)
pred_class <- factor(pred_class, levels = levels(testData$SurvivalStatus))

# Matriz de confusão
confusion <- confusionMatrix(pred_class, testData$SurvivalStatus)
pander::pander(confusion, caption = "Confusion Matrix Statistics (per class)")

# Precision, Recall, F1
positive_class <- "X0.ALIVE.OR.DEAD.TUMOR.FREE"
precision <- posPredValue(pred_class, testData$SurvivalStatus, positive = positive_class)
recall <- sensitivity(pred_class, testData$SurvivalStatus, positive = positive_class)
f1 <- 2 * (precision * recall) / (precision + recall)

cat("Precision:", precision, "\n")
cat("Recall:", recall, "\n")
cat("F1 Score:", round(f1, 3), "\n")

# Prever probabilidades
pred_prob <- predict(fit_svm1, testData, type = "prob")

# AUC
correct_col <- grep("ALIVE.OR.DEAD.TUMOR.FREE", colnames(pred_prob), value = TRUE)
roc_curve <- roc(response = testData$SurvivalStatus, predictor = pred_prob[[correct_col]])
auc_score <- auc(roc_curve)
cat("AUC Score:", auc_score, "\n")

# Curva ROC
plot(roc_curve, col = palette_3, main = "ROC Curve for SVM (com SMOTE)")

# Variavel importance
importance <- varImp(fit_svm1, scale = FALSE)

# importance em data frame
importance_df <- as.data.frame(importance$importance)
importance_df$Gene <- rownames(importance_df)

importance_df$Overall <- rowMeans(importance_df[, c("X0.ALIVE.OR.DEAD.TUMOR.FREE", "X1.DEAD.WITH.TUMOR")])

# Sort por importance
importance_df <- importance_df[order(-importance_df$Overall), ]

# Seleciona top 20 
top_20 <- head(importance_df, 20)

# tabela
knitr::kable(top_20[, c("Gene", "Overall")], caption = "Top 20 Most Important Genes")

# plot
plot(importance, top = 20, main = "Top 20 Important Features")
# Boxplot da distribuição de probabilidades
df_plot <- data.frame(Prob = pred_prob[[correct_col]], Status = testData$SurvivalStatus)
boxplot(Prob ~ Status, data = df_plot,
        col = palette_3,
        main = "Class Probability Distribution", ylab = "Predicted Probability")



## B. Support Vector Machine com SMOTE ---- SEX

# Preparar dados
common_samples <- intersect(colnames(expr_all), names(f_sex))
logCPM_sub <- expr_all[, common_samples]
labels <- f_sex[common_samples]

# Criar data frame com amostras como linhas e genes como colunas
data_model <- as.data.frame(t(logCPM_sub))
data_model$Sex <- factor(labels)

# Garantir que os levels são válidos para uso com class probabilities
levels(data_model$Sex) <- make.names(levels(data_model$Sex))

# Dividir treino e teste
set.seed(123)
trainIndex <- createDataPartition(data_model$Sex, p = 0.5, list = FALSE)
trainData <- data_model[trainIndex, ]
testData  <- data_model[-trainIndex, ]

# Corrigir aplicação do SMOTE
# Remover colunas não numéricas
trainData_numeric <- trainData[, sapply(trainData, is.numeric)]  

# Encontrar o tamanho mínimo da classe
min_class_size <- min(table(trainData$Sex))

# Ajustar K para evitar erro de SMOTE
K_value <- max(1, min_class_size - 1)  

# Aplicar SMOTE corretamente (usando factor numérico)
smote_result <- SMOTE(trainData_numeric, target = trainData$Sex, K = K_value, dup_size = 2)
smote_data <- smote_result$data

# Ajustar nomes das colunas
colnames(smote_data)[ncol(smote_data)] <- "Sex"
smote_data$Sex <- factor(smote_data$Sex, levels = levels(trainData$Sex))

# Remover valores ausentes
smote_data <- na.omit(smote_data)

# Treinar modelo SVM com validação cruzada e métrica ROC
set.seed(123)
fit_svm2 <- train(Sex ~ ., data = smote_data, method = "svmRadial",
                 trControl = ctrl,
                 preProcess = c("center", "scale"),
                 tuneLength=5,
                 metric = "ROC")

# Definir a classe positiva para avaliação
positive_class <- "Female"

# Previsão do modelo no conjunto de teste
pred_class <- predict(fit_svm2, testData)

# Criar matriz de confusão para avaliar o desempenho do modelo
confusion <- confusionMatrix(pred_class, testData$Sex)
pander::pander(confusion, caption = "Estatísticas da Matriz de Confusão (por classe)")

# Garantir que testData$Sex é um fator
testData$Sex <- factor(testData$Sex)                       

# Garantir que as previsões têm os mesmos níveis que os dados reais
pred_svm <- factor(pred_class, levels = levels(testData$Sex)) 


# Calcular Precisão, Recall e F1 Score
precision <- posPredValue(pred_svm, testData$Sex, positive = positive_class)
recall <- sensitivity(pred_svm, testData$Sex, positive = positive_class)

# Exibir os valores calculados
cat("Precisão:", precision, "\n")
cat("Recall:", recall, "\n")

# Calcular o F1 Score apenas se houver valores válidos de precisão e recall
if (!is.na(precision) && !is.na(recall) && (precision + recall) > 0) {
  f1 <- 2 * (precision * recall) / (precision + recall)
  cat("F1 Score:", round(f1, 3), "\n")
} else {
  cat("O F1 Score não pode ser calculado devido à ausência de precisão ou recall.\n")
}

# Previsão das probabilidades para cada classe
pred_prob <- predict(fit_svm2, testData, type = "prob")

# Garantir que testData$Sex tem os mesmos níveis que smote_data$Sex
testData$Sex <- factor(testData$Sex, levels = levels(smote_data$Sex))

# Encontrar dinamicamente a coluna correta de probabilidades
sex_levels <- levels(testData$Sex)
correct_col <- grep(paste(sex_levels, collapse = "|"), colnames(pred_prob), value = TRUE)

# Verificar se a coluna correta foi encontrada
if (length(correct_col) == 0) {
  stop(paste("Erro: Coluna de probabilidade para a classificação de Sex não encontrada. Nomes disponíveis:", 
             paste(colnames(pred_prob), collapse = ", ")))
}

# Garantir que apenas uma coluna é selecionada
if (length(correct_col) > 1) {
  correct_col <- correct_col[1]  # Selecionar a primeira correspondência
}

# Garantir que o preditor e a resposta têm o mesmo comprimento
if (length(testData$Sex) != nrow(pred_prob)) {
  stop("Erro: Incompatibilidade de tamanhos entre testData$Sex e as probabilidades previstas.")
}

# Converter preditor para numérico
numeric_pred <- as.numeric(pred_prob[[correct_col]])  

# Calcular a curva ROC
roc_curve <- roc(response = testData$Sex, predictor = numeric_pred)

# Calcular a AUC (Área sob a curva ROC)
auc_score <- auc(roc_curve)
cat("AUC Score:", auc_score, "\n")

# Plotar a curva ROC
plot(roc_curve, col = palette_3, main = "Curva ROC para SVM (com SMOTE)")
# Calcular a importância das variáveis (genes)
importance <- varImp(fit_svm2, scale = FALSE)

# Converter importância para um data frame
importance_df <- as.data.frame(importance$importance)
importance_df$Gene <- rownames(importance_df)

# Criar uma métrica geral de importância combinando as classes "Female" e "Male"
importance_df$Overall <- rowMeans(importance_df[, c("Female", "Male")])

# Ordenar por importância
importance_df <- importance_df[order(-importance_df$Overall), ]

# Selecionar as 20 características mais importantes
top_20 <- head(importance_df, 20)

# Exibir tabela com os genes mais importantes
knitr::kable(top_20[, c("Gene", "Overall")], caption = "Top 20 Genes Mais Importantes")

# Plot da importância das características
plot(importance, top = 20, main = "Top 20 Características Mais Importantes")

# Criar um data frame para o boxplot
df_plot <- data.frame(Prob = pred_prob[, correct_col], Status = testData$Sex)

# Gerar o boxplot para visualizar a distribuição das probabilidades previstas
boxplot(pred_prob[, correct_col] ~ testData$Sex, 
        col = palette_3, 
        main = "Distribuição das Probabilidades por Classe", ylab = "Probabilidade Prevista")



# --- B. Support Vector Machine com SMOTE ---- ANCESTRALIDADE

# Preparar dados
clean_data_for_2part <- clini_data_no_na

# Definir variável resposta
f_ancestralidade <- as.factor(clini_data_no_na$ancestry)
names(f_ancestralidade) <- clini_data_no_na$sample_ID

# Alinhar dados de expressão com amostras com dados clínicos
common_samples <- intersect(colnames(expr_all), names(f_ancestralidade))
logCPM_sub <- expr_all[, common_samples]
labels <- f_ancestralidade[common_samples]

# Criar data frame com amostras como linhas e genes como colunas
data_model <- as.data.frame(t(logCPM_sub))
data_model$Ancestralidade <- factor(labels)

# Garantir que os levels são válidos para uso com class probabilities
levels(data_model$Ancestralidade) <- make.names(levels(data_model$Ancestralidade))

# Dividir treino e teste
set.seed(123)
trainIndex <- createDataPartition(data_model$Ancestralidade, p = 0.5, list = FALSE)
trainData <- data_model[trainIndex, ]
testData  <- data_model[-trainIndex, ]

# --- FILTRAR CLASSES RARAS (< 5 amostras) ---
min_samples <- 5
table_counts <- table(trainData$Ancestralidade)
keep_levels <- names(table_counts[table_counts >= min_samples])

trainData_filtered <- trainData[trainData$Ancestralidade %in% keep_levels, ]
trainData_filtered$Ancestralidade <- droplevels(trainData_filtered$Ancestralidade)

# Atualizar levels no conjunto de teste para corresponder
testData <- testData[testData$Ancestralidade %in% keep_levels, ]
testData$Ancestralidade <- droplevels(testData$Ancestralidade)

# Preparar dados numéricos para SMOTE
trainData_numeric <- trainData_filtered[, sapply(trainData_filtered, is.numeric)]  

# Encontrar novo tamanho mínimo da classe
min_class_size <- min(table(trainData_filtered$Ancestralidade))
K_value <- max(1, min_class_size - 1)

# Aplicar SMOTE
smote_result <- SMOTE(trainData_numeric, target = trainData_filtered$Ancestralidade, K = K_value, dup_size = 2)
smote_data <- smote_result$data
colnames(smote_data)[ncol(smote_data)] <- "Ancestralidade"
smote_data$Ancestralidade <- factor(smote_data$Ancestralidade, levels = levels(trainData_filtered$Ancestralidade))

# Remover NAs
smote_data <- na.omit(smote_data)

# Treinar modelo SVM
set.seed(123)
fit_svm3 <- train(Ancestralidade ~ ., data = smote_data, method = "svmRadial",
                 trControl = ctrl1,
                 preProcess = c("center", "scale"),
                 tuneLength=5,
                 metric = "Accuracy")

# Previsão do modelo no conjunto de teste
pred_class <- predict(fit_svm3, testData)

# Criar matriz de confusão
conf_matrix <- confusionMatrix(pred_class, testData$Ancestralidade)
pander::pander(conf_matrix, caption = "Estatísticas da Matriz de Confusão (por classe)")
# Calcular Precisão, Recall e F1 Score por classe
precision_vec <- diag(conf_matrix$table) / colSums(conf_matrix$table)
recall_vec <- diag(conf_matrix$table) / rowSums(conf_matrix$table)

# Calcular F1-score para cada classe
f1_vec <- ifelse((precision_vec + recall_vec) > 0,
                 2 * (precision_vec * recall_vec) / (precision_vec + recall_vec),
                 NA)

# Criar um data frame com os valores de F1-score
f1_scores <- data.frame(
  Classe = rownames(conf_matrix$table),
  Precisão = round(precision_vec, 3),
  Recall = round(recall_vec, 3),
  F1 = round(f1_vec, 3)
)

cat("F1, Recall e Precisão por classe:\n")
knitr::kable(f1_scores, caption = "Precisão, Recall e F1 Score por Classe")

# Previsão de probabilidades (para AUC e boxplot)
pred_prob <- predict(fit_svm3, testData, type = "prob")

# Garantir que testData$Ancestralidade tem os mesmos níveis que pred_prob
testData$Ancestralidade <- factor(testData$Ancestralidade, levels = colnames(pred_prob))


# Criar listas para armazenar curvas ROC e valores de AUC
roc_list <- list()
auc_scores <- c()
classes <- colnames(pred_prob)

# Calcular curvas ROC para cada classe (um-contra-todos)
for (class_name in classes) {
  # Criar rótulos binários: 1 para a classe atual, 0 para as restantes
  binary_response <- factor(ifelse(testData$Ancestralidade == class_name, class_name, paste0("not_", class_name)))
  
  # Calcular curva ROC
  roc_curve <- roc(response = binary_response,
                   predictor = pred_prob[, class_name],
                   levels = c(paste0("not_", class_name), class_name))
  
  # Guardar valores de AUC e curva ROC
  auc_score <- auc(roc_curve)
  auc_scores <- c(auc_scores, auc_score)
  roc_list[[class_name]] <- roc_curve
  
  cat(paste("AUC para", class_name, "vs todos:", round(auc_score, 3), "\n"))
}

# Plot a primeira curva ROC
plot(roc_list[[1]], col = palette_3[1], main = "Curvas ROC Multi-classe para SVM (com SMOTE)")
# Calcular a importância das variáveis (genes)
importance <- varImp(fit_svm3, scale = FALSE)

# Converter importância para um data frame
importance_df <- as.data.frame(importance$importance)
importance_df$Gene <- rownames(importance_df)

# Criar uma métrica geral de importância
importance_df$Overall <- rowMeans(importance_df[, c("AFR", "EUR")])

# Ordenar por importância
importance_df <- importance_df[order(-importance_df$Overall), ]

# Selecionar as 20 características mais importantes
top_20 <- head(importance_df, 20)

# Exibir tabela com os genes mais importantes
knitr::kable(top_20[, c("Gene", "Overall")], caption = "Top 20 Genes Mais Importantes")

# Plotar a importância das características
plot(importance, top = 20, main = "Top 20 Características Mais Importantes")
### Correção do Boxplot da Distribuição de Probabilidades
# Encontrar dinamicamente a coluna correta de probabilidades
correct_col1 <- grep(paste(levels(testData$Ancestralidade), collapse = "|"), colnames(pred_prob), value = TRUE)

# Garantir que apenas uma coluna é selecionada
if (length(correct_col1) > 1) {
    correct_col1 <- correct_col1[1]  # Selecionar a primeira correspondência
}

# Verificar se a coluna correta foi encontrada
if (!(correct_col1 %in% colnames(pred_prob))) {
    stop("Erro: Coluna de probabilidade para a classificação de Ancestralidade não encontrada.")
}

# Converter preditor para numérico
numeric_pred <- as.numeric(pred_prob[[correct_col1]])  

# Criar data frame para o boxplot
df_plot <- data.frame(Prob = numeric_pred, Status = testData$Ancestralidade)

# Gerar o boxplot
boxplot(df_plot$Prob ~ df_plot$Status, 
        col = palette_3, 
        main = "Distribuição das Probabilidades por Classe", ylab = "Probabilidade Prevista")
## B. Support Vector Machine com SMOTE ---- CANCER

# Preparar dados
common_samples <- intersect(colnames(expr_all), names(f_cancer))
logCPM_sub <- expr_all[, common_samples]
labels <- f_cancer[common_samples]

# Criar data frame com amostras como linhas e genes como colunas
data_model <- as.data.frame(t(logCPM_sub))
data_model$Cancer <- factor(labels)

# Garantir que os levels são válidos para uso com class probabilities
levels(data_model$Cancer) <- make.names(levels(data_model$Cancer))

# Dividir treino e teste
set.seed(123)
trainIndex <- createDataPartition(data_model$Cancer, p = 0.5, list = FALSE)
trainData <- data_model[trainIndex, ]
testData  <- data_model[-trainIndex, ]

# Remover colunas não numéricas
trainData_numeric <- trainData[, sapply(trainData, is.numeric)]  

# Encontrar o tamanho mínimo da classe
min_class_size <- min(table(trainData$Cancer))

# Ajustar K para evitar erro de SMOTE
K_value <- max(1, min_class_size - 1)  # Garante que K seja pelo menos 1

# Aplicar SMOTE corretamente (usando factor numérico)
smote_result <- SMOTE(trainData_numeric, target = trainData$Cancer, K = K_value, dup_size = 2)
smote_data <- smote_result$data

# Ajustar nomes das colunas
colnames(smote_data)[ncol(smote_data)] <- "Cancer"
smote_data$Cancer <- factor(smote_data$Cancer, levels = levels(trainData$Cancer))

# Remover valores ausentes
smote_data <- na.omit(smote_data)

# Treinar modelo SVM com validação cruzada e métrica ROC
set.seed(123)
fit_svm4 <- train(Cancer ~ ., data = smote_data, method = "svmRadial",
                 trControl = ctrl1,  # Use multiClassSummary
                 preProcess = c("center", "scale"),
                 tuneLength=5,
                 metric = "Accuracy")  # Use Accuracy for multi-class problems


# Avaliação no conjunto de teste
pred_class <- predict(fit_svm4, testData)
confusion <- confusionMatrix(pred_class, testData$Cancer)
pander::pander(confusion, caption = "Confusion Matrix Statistics (per class)")
# Extrair tabela confusion matrix 
conf_table <- confusion$table

# Initializa data frame para armazenar métricas
class_metrics <- data.frame(Class = levels(testData$Cancer),
                            Precision = numeric(length(levels(testData$Cancer))),
                            Recall = numeric(length(levels(testData$Cancer))),
                            F1_Score = numeric(length(levels(testData$Cancer))))

# Ciclo sobre cada classe para calcular precisão, sensibilidade e F1 score
for (class in levels(testData$Cancer)) {
  # Verdadeiros Positivos (VP): Previsões corretas da classe
  TP <- conf_table[class, class]
  
  # Falsos Positivos (FP): Previsões da classe, mas pertencem a outra
  FP <- sum(conf_table[, class]) - TP
  
  # Falsos Negativos (FN): Pertencem à classe, mas foram previstos como outra
  FN <- sum(conf_table[class, ]) - TP
  
  # Precisão: VP / (VP + FP)
  precision <- ifelse((TP + FP) == 0, NA, TP / (TP + FP))
  
  # Sensibilidade (Recall): VP / (VP + FN)
  recall <- ifelse((TP + FN) == 0, NA, TP / (TP + FN))
  
  # F1 Score: 2 * (Precisão * Sensibilidade) / (Precisão + Sensibilidade)
  f1 <- ifelse(is.na(precision) | is.na(recall) | (precision + recall) == 0, NA, 2 * (precision * recall) / (precision + recall))
  
  # Guardar os resultados
  class_metrics[class_metrics$Class == class, ] <- c(class, precision, recall, f1)
}

# Converter colunas numéricas
class_metrics[, 2:4] <- lapply(class_metrics[, 2:4], as.numeric)

# Round
class_metrics_rounded <- class_metrics
class_metrics_rounded[, 2:4] <- round(class_metrics_rounded[, 2:4], 3)

# Tabela
knitr::kable(class_metrics_rounded, caption = "Per-Class Performance Metrics (Precision, Recall, F1 Score)")

# Prever probabilidades
pred_prob <- predict(fit_svm4, testData, type = "prob")

# Garantir que testData$Cancer tem os mesmos níveis que smote_data$Cancer
testData$Cancer <- factor(testData$Cancer, levels = levels(smote_data$Cancer))

# Encontrar a coluna correta para probabilidades
cancer_levels <- levels(testData$Cancer)
correct_col <- grep(paste(cancer_levels, collapse = "|"), colnames(pred_prob), value = TRUE)

# Verificar se a coluna correta foi encontrada
if (length(correct_col) == 0) {
  stop(paste("Erro: Coluna de probabilidade para classificação de Cancer não encontrada. Nomes disponíveis:", 
             paste(colnames(pred_prob), collapse = ", ")))
}

# Garantir que apenas uma coluna seja selecionada
if (length(correct_col) > 1) {
  correct_col <- correct_col[1] 
}

# Verificar se os tamanhos das variáveis são compatíveis
if (length(testData$Cancer) != nrow(pred_prob)) {
  stop("Erro: Diferente número de elementos entre testData$Cancer e pred_prob.")
}
# Converter preditor para numérico
numeric_pred <- as.numeric(pred_prob[[correct_col]])  # Usa [[ ]] para extrair como vetor

# Calcular AUC e Curva ROC
roc_curve <- roc(response = testData$Cancer, predictor = numeric_pred)

auc_score <- auc(roc_curve)
cat("AUC Score:", auc_score, "\n")

# Plot Curva ROC
plot(roc_curve, col = palette_3, main = "ROC Curve for SVM (com SMOTE)")
# Calcular a importância das variáveis no modelo
importance <- varImp(fit_svm4, scale = FALSE)

# Converter a importância para um data frame
importance_df <- as.data.frame(importance$importance)

# Adicionar os nomes dos genes como uma coluna separada
importance_df$Gene <- rownames(importance_df)

# Criar uma métrica geral de importância combinando as classes "Atrocytoma", "Oligoastrocytoma" e "Oligodendroglioma"
importance_df$Overall <- rowMeans(importance_df[, c("Astrocytoma", "Oligoastrocytoma", "Oligodendroglioma")])

# Ordenar os genes por importância (do maior para o menor)
importance_df <- importance_df[order(-importance_df$Overall), ]

# Selecionar os 20 genes mais importantes
top_20 <- head(importance_df, 20)

# Exibir tabela com os genes mais importantes
knitr::kable(top_20[, c("Gene", "Overall")], caption = "Top 20 Genes Mais Importantes")

# Plotar a importância das características
plot(importance, top = 20, main = "Top 20 Características Mais Importantes")

# Boxplot da distribuição de probabilidades
df_plot <- data.frame(Prob = pred_prob[[correct_col]], Status = testData$Cancer)

boxplot(pred_prob[[correct_col]] ~ testData$Cancer, 
        col = palette_3, 
        main = "Class Probability Distribution", ylab = "Predicted Probability")

### **11.3 XGBoost**
## C. XGBoost com SMOTE ---- SURVIVAL

# Preparar dados
clean_data_for_2part <- clini_data_no_na

# Definir variável resposta
f_survival <- as.factor(clini_data_no_na$surv_status)
names(f_survival) <- clini_data_no_na$sample_ID

# Alinhar dados de expressão com amostras com dados clínicos
common_samples <- intersect(colnames(expr_all), names(f_survival))
logCPM_sub <- expr_all[, common_samples]
labels <- f_survival[common_samples]

# Construir dataframe para modelagem
data_model <- as.data.frame(t(logCPM_sub))  # genes como colunas
data_model$SurvivalStatus <- factor(labels)

# Ajustar levels para uso com SMOTE e xgboost
levels(data_model$SurvivalStatus) <- make.names(levels(data_model$SurvivalStatus))

# Dividir treino e teste
set.seed(123)
trainIndex <- createDataPartition(data_model$SurvivalStatus, p = 0.5, list = FALSE)
trainData <- data_model[trainIndex, ]
testData  <- data_model[-trainIndex, ]

# SMOTE
trainData_numeric <- trainData[, sapply(trainData, is.numeric)]
smote_result <- SMOTE(trainData_numeric, trainData$SurvivalStatus, K = 3, dup_size = 2)
smote_data <- smote_result$data
colnames(smote_data)[ncol(smote_data)] <- "SurvivalStatus"
smote_data$SurvivalStatus <- factor(smote_data$SurvivalStatus, 
                                    levels = levels(trainData$SurvivalStatus))

# Treinar modelo XGBoost
set.seed(123)
capture.output({
  fit_xgb1 <- train(SurvivalStatus ~ ., data = smote_data, method = "xgbTree",
                   trControl = ctrl, 
                   metric ="ROC", 
                   tuneLength=3)
}, file = "NUL")

# Avaliação no conjunto de teste
pred_class <- predict(fit_xgb1, testData)
 
conf_matrix <- confusionMatrix(pred_class, testData$SurvivalStatus)
pander::pander(conf_matrix, caption = "Confusion Matrix Statistics (per class)")
# Previsão de probabilidades
pred_prob <- suppressWarnings(predict(fit_xgb1, testData, type = "prob"))
correct_col <- grep("ALIVE.OR.DEAD.TUMOR.FREE", colnames(pred_prob), value = TRUE)

f1_smote <- F1_Score(y_pred = pred_class, y_true = testData$SurvivalStatus, positive = "ALIVE.OR.DEAD.TUMOR.FREE")

# Precision, Recall
positive_class <- "X0.ALIVE.OR.DEAD.TUMOR.FREE"
precision <- posPredValue(pred_class, testData$SurvivalStatus, positive = positive_class)
recall <- sensitivity(pred_class, testData$SurvivalStatus, positive = positive_class)
f1 <- 2 * (precision * recall) / (precision + recall)

cat("Recall:", recall, "\n")
cat("Precision:", precision, "\n")
cat("F1 Score:", round(f1, 3), "\n")
# ROC
roc_curve <- roc(response = testData$SurvivalStatus, predictor = pred_prob[, correct_col])
auc_score <- auc(roc_curve)
cat("AUC Score:", auc_score, "\n")

# Curva ROC
plot(roc_curve, col = palette_3, main = "ROC Curve for XGBoost (com SMOTE)")
# Importância das features
importance <- varImp(fit_xgb1, scale = FALSE)

# Extrair variavel importance como data frame
importance_df <- as.data.frame(importance$importance)
importance_df$Gene <- rownames(importance_df)
importance_df <- importance_df[order(-importance_df$Overall), ]  # Sort by importance

# Seleção top 20
top_20 <- head(importance_df, 20)

# Tabela
knitr::kable(top_20[, c("Gene", "Overall")], caption = "Top 20 Most Important Genes")

plot(importance, top = 20, main = "Top 20 Important Features")
# Boxplot das probabilidades
boxplot(pred_prob[, correct_col] ~ testData$SurvivalStatus,
        col = palette_3,
        main = "Class Probability Distribution", ylab = "Predicted Probability")
## C. XGBoost com SMOTE - Sex

# Preparar dados
clean_data_for_2part <- clini_data_no_na

# Definir variável resposta
f_sex <- as.factor(clean_data_for_2part$sex)
names(f_sex) <- clean_data_for_2part$sample_ID

# Alinhar dados de expressão com amostras com dados clínicos
common_samples <- intersect(colnames(expr_all), names(f_sex))
logCPM_sub <- expr_all[, common_samples]
labels <- f_sex[common_samples]

# Construir dataframe para modelagem
data_model <- as.data.frame(t(logCPM_sub))  # genes como colunas
data_model$Sex <- factor(labels)

# Ajustar levels para uso com SMOTE e xgboost
levels(data_model$Sex) <- make.names(levels(data_model$Sex))

# Dividir treino e teste
set.seed(123)
trainIndex <- createDataPartition(data_model$Sex, p = 0.5, list = FALSE)
trainData <- data_model[trainIndex, ]
testData  <- data_model[-trainIndex, ]

# SMOTE
trainData_numeric <- trainData[, sapply(trainData, is.numeric)]
smote_result <- SMOTE(trainData_numeric, trainData$Sex, K = 3, dup_size = 2)
smote_data <- smote_result$data
colnames(smote_data)[ncol(smote_data)] <- "Sex"
smote_data$Sex <- factor(smote_data$Sex, levels = levels(trainData$Sex))

# Treinar modelo XGBoost
set.seed(123)
capture.output({
  fit_xgb2 <- train(Sex ~ ., data = smote_data, method = "xgbTree",
                    trControl = ctrl, 
                    metric = "ROC", tuneLength=3)
}, file = "NUL")

# Avaliação no conjunto de teste
pred_class <- predict(fit_xgb2, testData)
conf_matrix <- confusionMatrix(pred_class, testData$Sex)
pander::pander(conf_matrix, caption = "Confusion Matrix Statistics (per class)")
# Previsão de probabilidades
pred_prob <- suppressWarnings(predict(fit_xgb2, testData, type = "prob"))
correct_col <- grep("Female", colnames(pred_prob), value = TRUE)  

# Precision, Recall, F1
positive_class <- "Female"

precision <- posPredValue(pred_class, testData$Sex, positive = positive_class)
recall <- sensitivity(pred_class, testData$Sex, positive = positive_class)
f1 <- 2 * (precision * recall) / (precision + recall)

# Print precision, recall,F1
cat("Precision:", precision, "\n")
cat("Recall:", recall, "\n")
cat("F1 Score:", round(f1, 3), "\n")
# Curva ROC e AUC
roc_curve <- roc(response = testData$Sex, predictor = pred_prob[, correct_col])
auc_score <- auc(roc_curve)
cat("AUC Score:", auc_score, "\n")

# Curva ROC
plot(roc_curve, col = palette_3, main = "ROC Curve for XGBoost (com SMOTE) - Sexo")
# Importância das features
importance <- varImp(fit_xgb2, scale = FALSE)

# Extrair variavel importance como data frame
importance_df <- as.data.frame(importance$importance)
importance_df$Gene <- rownames(importance_df)
importance_df <- importance_df[order(-importance_df$Overall), ]  # Sort by importance

# Seleção top 20
top_20 <- head(importance_df, 20)

# Tabela
knitr::kable(top_20[, c("Gene", "Overall")], caption = "Top 20 Most Important Genes")

plot(importance, top = 20, main = "Top 20 Important Features")

# Boxplot das probabilidades
boxplot(pred_prob[, correct_col] ~ testData$Sex,
        col = palette_3,
        main = "Distribuição das Probabilidades - Sexo", ylab = "Probabilidade Prevista")
## C. XGBoost com SMOTE - Ancestralidade
# Preparar dados
clean_data_for_2part <- clini_data_no_na

# Preparar f_ancestralidade sem NA
f_ancestralidade <- as.factor(clean_data_for_2part$ancestry)
names(f_ancestralidade) <- clean_data_for_2part$sample_ID
f_ancestralidade <- f_ancestralidade[!is.na(f_ancestralidade)]

# Alinhar com dados de expressão
common_samples <- intersect(colnames(expr_all), names(f_ancestralidade))
logCPM_sub <- expr_all[, common_samples]
labels <- f_ancestralidade[common_samples]

# Criar dataframe para modelagem
data_model <- as.data.frame(t(logCPM_sub))  # genes como colunas
data_model$Ancestralidade <- factor(labels)
data_model <- na.omit(data_model)  # Garantir que não há NAs
levels(data_model$Ancestralidade) <- make.names(levels(data_model$Ancestralidade))

# Dividir treino e teste
set.seed(123)
trainIndex <- createDataPartition(data_model$Ancestralidade, p = 0.5, list = FALSE)
trainData <- data_model[trainIndex, ]
testData  <- data_model[-trainIndex, ]

# Remover classes com poucas amostras no treino
min_class_size <- 5
keep_classes <- names(which(table(trainData$Ancestralidade) >= min_class_size))
trainData <- trainData[trainData$Ancestralidade %in% keep_classes, ]
trainData$Ancestralidade <- droplevels(trainData$Ancestralidade)

# Aplicar SMOTE
trainData_numeric <- trainData[, sapply(trainData, is.numeric)]
smote_result <- SMOTE(trainData_numeric, trainData$Ancestralidade, K = 3, dup_size = 2)
smote_data <- smote_result$data
colnames(smote_data)[ncol(smote_data)] <- "Ancestralidade"
smote_data$Ancestralidade <- factor(smote_data$Ancestralidade, levels = levels(trainData$Ancestralidade))

# Treinar modelo XGBoost
set.seed(123)
capture.output({
  fit_xgb3 <- train(Ancestralidade ~ ., data = smote_data, method = "xgbTree",
                    trControl = ctrl1,
                    metric = "ROC",
                    tuneLength=3)
}, file = "NUL")

# Avaliação no conjunto de teste
pred_class <- predict(fit_xgb3, testData)
conf_matrix <- confusionMatrix(pred_class, testData$Ancestralidade)
pander::pander(conf_matrix, caption = "Confusion Matrix Statistics (per class)")
# Previsão de probabilidades (para AUC e boxplot)
pred_prob <- suppressWarnings(predict(fit_xgb3, testData, type = "prob"))

# Filtrar apenas as classes "EUR" e "AFR"
testData <- testData[testData$Ancestralidade %in% c("EUR", "AFR"), ]
pred_prob <- predict(fit_xgb3, testData, type = "prob")

# Escolher uma classe como "positiva" arbitrária (por exemplo, a primeira)
positive_class <- colnames(pred_prob)[1]

# Calcular Precision e Recall
precision_vec <- diag(conf_matrix$table) / colSums(conf_matrix$table)
recall_vec <- diag(conf_matrix$table) / rowSums(conf_matrix$table)

# Calcular Macro F1-score manualmente (sem yardstick)
f1_per_class <- function(conf_mat) {
  by_class <- conf_mat$byClass
  if (is.matrix(by_class)) {
    f1 <- 2 * by_class[, "Sensitivity"] * by_class[, "Precision"] /
            (by_class[, "Sensitivity"] + by_class[, "Precision"])
    f1[is.na(f1)] <- 0
    return(mean(f1))
  } else {
    # Only two classes
    f1 <- 2 * by_class["Sensitivity"] * by_class["Precision"] /
            (by_class["Sensitivity"] + by_class["Precision"])
    return(ifelse(is.na(f1), 0, f1))
  }
}

f1_macro <- f1_per_class(conf_matrix)
cat("F1 Score:", round(f1_macro, 3), "\n")
cat("Recall:", round(recall_vec, 3), "\n")
cat("Precision:", round(precision_vec, 3), "\n")
roc_curve <- roc(response = testData$Ancestralidade, 
                 predictor = pred_prob[, positive_class], 
                 levels = c("EUR", "AFR"))

auc_score <- auc(roc_curve)
cat("AUC Score:", auc_score, "\n")

# Plot ROC curve
plot(roc_curve, col = palette_3, lwd = 2, main = "ROC Curve - EUR vs AFR")

# cálculo de AUC para multiclasse
multiclass_auc <- multiclass.roc(response = testData$Ancestralidade,
                                 predictor = as.matrix(pred_prob))
auc_score <- auc(multiclass_auc)
cat("Multiclass AUC Score:", round(auc_score, 3), "\n")


# Importância das features
importance <- varImp(fit_xgb3, scale = FALSE)

# Extração variavel importance como data frame
importance_df <- as.data.frame(importance$importance)
importance_df$Gene <- rownames(importance_df)
importance_df <- importance_df[order(-importance_df$Overall), ]  # Sort by importance

# Seleção top 20
top_20 <- head(importance_df, 20)

# Tabela
knitr::kable(top_20[, c("Gene", "Overall")], caption = "Top 20 Most Important Genes")

plot(importance, top = 20, main = "Top 20 Important Features")
# Boxplot das probabilidades para cada classe
positive_class <- colnames(pred_prob)[1]
boxplot(pred_prob[, positive_class] ~ testData$Ancestralidade,
        col = palette_3,
        main = paste("Probabilidade predita - Classe:", positive_class),
        ylab = "Predicted Probability")


## C. XGBoost com SMOTE - Cancer

# Preparar dados
clean_data_for_2part <- clini_data_no_na

# Definir variável resposta
f_cancer <- as.factor(clean_data_for_2part$cancer_type)
names(f_cancer) <- clean_data_for_2part$sample_ID

# Alinhar dados de expressão com amostras com dados clínicos
common_samples <- intersect(colnames(expr_all), names(f_cancer))
logCPM_sub <- expr_all[, common_samples]
labels <- f_cancer[common_samples]

# Construir dataframe para modelagem
data_model <- as.data.frame(t(logCPM_sub))  # genes como colunas
data_model$TipoTumor <- factor(labels)

# Ajustar levels para uso com SMOTE e xgboost
levels(data_model$TipoTumor) <- make.names(levels(data_model$TipoTumor))

# Dividir treino e teste
set.seed(123)
trainIndex <- createDataPartition(data_model$TipoTumor, p = 0.5, list = FALSE)
trainData <- data_model[trainIndex, ]
testData  <- data_model[-trainIndex, ]

# SMOTE
trainData_numeric <- trainData[, sapply(trainData, is.numeric)]
smote_result <- SMOTE(trainData_numeric, trainData$TipoTumor, K = 3, dup_size = 2)
smote_data <- smote_result$data
colnames(smote_data)[ncol(smote_data)] <- "TipoTumor"
smote_data$TipoTumor <- factor(smote_data$TipoTumor, levels = levels(trainData$TipoTumor))

# Treinar modelo XGBoost
set.seed(123)
capture.output({
  fit_xgb4 <- train(TipoTumor ~ ., data = smote_data, method = "xgbTree",
                    trControl = ctrl1, 
                    metric = "ROC", tuneLength=3)
}, file = "NUL")

# Avaliação no conjunto de teste
pred_class <- predict(fit_xgb4, testData)
conf_matrix <- confusionMatrix(pred_class, testData$TipoTumor)
pander::pander(conf_matrix, caption = "Confusion Matrix Statistics (per class)")

# Previsão de probabilidades
pred_prob <- suppressWarnings(predict(fit_xgb4, testData, type = "prob"))

# Escolher uma classe como "positiva" arbitrária (por exemplo, a primeira)
positive_class <- levels(testData$TipoTumor)[1]
roc_curve <-roc(response = testData$TipoTumor, predictor = pred_prob[, positive_class])
auc_score <- auc(roc_curve)
cat("AUC Score:", auc_score, "\n")

# Curva ROC
plot(roc_curve, col = palette_3, main = "ROC Curve for XGBoost (com SMOTE) - TipoTumor")
# Calcular F1-score, Recall e Precision (para a mesma classe positiva)
f1_smote <- F1_Score(y_pred = pred_class, y_true = testData$TipoTumor, positive = positive_class)

precision_vec <- diag(conf_matrix$table) / colSums(conf_matrix$table)
recall_vec <- diag(conf_matrix$table) / rowSums(conf_matrix$table)

cat("F1 Score:", f1_smote, "\n")
cat("Recall:", round(recall_vec, 3), "\n")
cat("Precision:", round(precision_vec, 3), "\n")



# Importância das features
importance <- varImp(fit_xgb4, scale = FALSE)

# Extrair variável importance como data frame
importance_df <- as.data.frame(importance$importance)
importance_df$Gene <- rownames(importance_df)
importance_df <- importance_df[order(-importance_df$Overall), ]  # Sort by importance

# Seleção top 20
top_20 <- head(importance_df, 20)

# Tabela
knitr::kable(top_20[, c("Gene", "Overall")], caption = "Top 20 Most Important Genes")

plot(importance, top = 20, main = "Top 20 Important Features")

# Boxplot das probabilidades
boxplot(pred_prob[, positive_class] ~ testData$TipoTumor,
        col = palette_3,
        main = "Distribuição das Probabilidades - TipoTumor", ylab = "Probabilidade Prevista")
## **12. Comparação dos modelos por variável**
# Variável Survival

min_resamples <- min(sapply(list(fit_rf1, fit_svm1, fit_xgb1), function(x) nrow(x$resample)))

fit_rf1$resample <- fit_rf1$resample[order(fit_rf1$resample$Resample), ]
fit_svm1$resample <- fit_svm1$resample[order(fit_svm1$resample$Resample), ]
fit_xgb1$resample <- fit_xgb1$resample[order(fit_xgb1$resample$Resample), ]

# Comparar resultados
results_survival <- resamples(list(RF = fit_rf1, SVM = fit_svm1, XGB = fit_xgb1))

summary_df <- as.data.frame(summary(results_survival)$statistics)
summary_df <- rownames_to_column(summary_df, "Metric_Model")

datatable(summary_df, 
          options = list(pageLength = 10, autoWidth = TRUE), 
          caption = "Summary of Resampling Metrics")

# Comparação Visual
bwplot(results_survival,
       col = palette_3,
       main = "Medidas ", ylab = "Model")

# Variável Sex

min_resamples <- min(sapply(list(fit_rf2, fit_svm2, fit_xgb2), function(x) nrow(x$resample)))

fit_rf2$resample <- fit_rf2$resample[1:min_resamples, ]
fit_svm2$resample <- fit_svm2$resample[1:min_resamples, ]
fit_xgb2$resample <- fit_xgb2$resample[1:min_resamples, ]

# Comparar resultados
results_sex <- resamples(list(RF = fit_rf2, SVM = fit_svm2, XGB = fit_xgb2))

summary_df_sex <- as.data.frame(summary(results_sex)$statistics)
summary_df_sex <- rownames_to_column(summary_df_sex, "Metric_Model")

datatable(summary_df_sex, 
          options = list(pageLength = 10, autoWidth = TRUE), 
          caption = "Summary of Resampling Metrics (Sex)")

# Comparação Visual
bwplot(results_sex,
       col = palette_3,
       main = "Medidas ", ylab = "Model")

# Variável Ancestralidade

min_resamples <- min(sapply(list(fit_rf3, fit_svm3, fit_xgb3), function(x) nrow(x$resample)))

fit_rf3$resample <- fit_rf3$resample[1:min_resamples, ]
fit_svm3$resample <- fit_svm3$resample[1:min_resamples, ]
fit_xgb3$resample <- fit_xgb3$resample[1:min_resamples, ]

# Comparar resultados
results_ancestralidade <- resamples(list(RF = fit_rf3, SVM = fit_svm3, XGB = fit_xgb3))

summary_df_ancestralidade <- as.data.frame(summary(results_ancestralidade)$statistics)
summary_df_ancestralidade <- rownames_to_column(summary_df_ancestralidade, "Metric_Model")

datatable(summary_df_ancestralidade, 
          options = list(pageLength = 10, autoWidth = TRUE), 
          caption = "Summary of Resampling Metrics (Ancestralidade)")

# Comparação Visual
bwplot(results_ancestralidade,
       col = palette_3,
       main = "Medidas ", ylab = "Model")


# Variável Cancer

min_resamples <- min(sapply(list(fit_rf4, fit_svm4, fit_xgb4), function(x) nrow(x$resample)))

fit_rf4$resample <- fit_rf4$resample[1:min_resamples, ]
fit_svm4$resample <- fit_svm4$resample[1:min_resamples, ]
fit_xgb4$resample <- fit_xgb4$resample[1:min_resamples, ]

# Comparar resultados
results_cancer <- resamples(list(RF = fit_rf4, SVM = fit_svm4, XGB = fit_xgb4))

summary_df_cancer <- as.data.frame(summary(results_cancer)$statistics)
summary_df_cancer <- rownames_to_column(summary_df_cancer, "Metric_Model")

datatable(summary_df_cancer, 
          options = list(pageLength = 10, autoWidth = TRUE), 
          caption = "Summary of Resampling Metrics (Cancer)")

# Comparação Visual
bwplot(results_cancer,
       col = palette_3,
       main = "Medidas ", ylab = "Model")

                               
