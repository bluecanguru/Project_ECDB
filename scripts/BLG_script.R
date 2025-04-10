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
  "htmltools"     # To print text in specific format
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
